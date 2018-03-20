/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2017 Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "utils/core/logging.hpp"

#include "sequence/counts.hpp"
#include "sequence/formats/fasta_input_iterator.hpp"
#include "sequence/formats/fasta_reader.hpp"
#include "sequence/formats/fasta_writer.hpp"
#include "sequence/functions/consensus.hpp"
#include "sequence/functions/entropy.hpp"
#include "sequence/functions/functions.hpp"
#include "sequence/sequence_set.hpp"
#include "sequence/sequence.hpp"

#include "taxonomy/formats/taxonomy_reader.hpp"
#include "taxonomy/formats/taxopath_generator.hpp"
#include "taxonomy/formats/taxopath_parser.hpp"
#include "taxonomy/functions/split.hpp"
#include "taxonomy/functions/taxonomy.hpp"
#include "taxonomy/functions/taxopath.hpp"
#include "taxonomy/taxon.hpp"
#include "taxonomy/taxonomy.hpp"
#include "taxonomy/taxopath.hpp"

#include "utils/core/fs.hpp"
#include "utils/core/std.hpp"
#include "utils/formats/csv/reader.hpp"
#include "utils/io/input_stream.hpp"
#include "utils/text/string.hpp"
#include "utils/tools/date_time.hpp"

#include "utils/math/histogram.hpp"
#include "utils/math/histogram/accumulator.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace genesis;
using namespace genesis::sequence;
using namespace genesis::taxonomy;

int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();

    // -------------------------------------------------------------------------
    //     Prepare.
    // -------------------------------------------------------------------------

    // Load the taxonomy.
    LOG_TIME << "Reading taxonomy...";
    std::string tax_file = "/home/lucas/Projects/data/silva/tax_slv_ssu_123.1.txt";
    auto tax = Taxonomy();
    TaxonomyReader().from_file( tax_file, tax );
    sort_by_name( tax );

    // tax.remove_child( "Bacteria" );
    // tax.remove_child( "Eukaryota" );

    LOG_TIME << "done";

    LOG_DBG << "total children:";
    LOG_DBG1 << "Archaea " << total_taxa_count( tax.get_child("Archaea") );
    LOG_DBG1 << "Bacteria " << total_taxa_count( tax.get_child("Bacteria") );
    LOG_DBG1 << "Eukaryota " << total_taxa_count( tax.get_child("Eukaryota") );
    // return 0;

    // Create a Sequence Count objeect for each taxon.
    LOG_TIME << "Preparing counts map...";
    auto counts_map = std::unordered_map< Taxon*, SequenceCounts >();
    auto add_counts_to_map = [&] ( Taxon& t ) {
        auto gen = TaxopathGenerator();
        auto name = gen(t);

        // if( ! utils::starts_with( name, "Archaea" )) {
        //     return;
        // }

        counts_map[ &t ] = SequenceCounts( "ACGU", 50000 );
    };
    preorder_for_each( tax, add_counts_to_map );
    LOG_TIME << "done";
    LOG_INFO << "counts map size " << counts_map.size();

    // -------------------------------------------------------------------------
    //     Read all sequences.
    // -------------------------------------------------------------------------

    // Prepare sequence input.
    // std::string aln_file = "/home/lucas/Projects/data/silva/archaea.fasta";
    std::string aln_file = "/home/lucas/Projects/data/silva/SILVA_123.1_SSURef_Nr99_tax_silva_full_align_trunc.fasta";
    utils::InputStream aln_is { utils::make_unique<utils::FileInputSource>( aln_file ) };
    auto reader = FastaReader();
    auto taxopath_parser = TaxopathParser();
    Sequence seq;

    // Add to count objects along their taxonomic path.
    LOG_TIME << "Start reading at " << utils::current_time();
    size_t seq_cnt = 0;
    while( reader.parse_sequence( aln_is, seq )) {
        if( seq_cnt % 100000 == 0 ) {
            LOG_TIME << "At sequence " << seq_cnt;
        }
        ++seq_cnt;

        if( seq.metadata() == "" ) {
            throw std::runtime_error( "Empty metadata." );
        }

        auto taxopath = taxopath_parser( seq.metadata() );
        taxopath.pop_back();
        auto taxp = find_taxon_by_taxopath( tax, taxopath );
        if( taxp == nullptr ) {
            throw std::runtime_error( "Sequence taxon not in taxonomy: " + seq.metadata() );
        }

        // if( taxopath[0] != "Archaea" ) {
        //     continue;
        // }

        auto cur_tax = taxp;
        do {
            if( counts_map.count( cur_tax ) == 0 ) {
                throw std::runtime_error( "Sequence taxon not in map: " + seq.metadata() );
            }

            counts_map[cur_tax].add_sequence( seq );
            cur_tax = cur_tax->parent();
        } while( cur_tax != nullptr );
    }
    LOG_TIME << "Done reading at " << utils::current_time();

    // -------------------------------------------------------------------------
    //     Entropy splitting Options.
    // -------------------------------------------------------------------------

    for( size_t i = 0; i < 4; ++i ) {
        auto opt = static_cast< SiteEntropyOptions >(i);

        LOG_DBG << "Option " << i;

        // Make a tmp copy so that we can remove stuff.
        auto& tax_cpy = tax;

        // -------------------------------------------------------------------------
        //     Splitting
        // -------------------------------------------------------------------------

        LOG_TIME << "entropy calculations";
        std::unordered_map< Taxon const*, double > entropies;
        auto calc_entropies = [&]( Taxon& t ) {
            auto gen = TaxopathGenerator();
            auto name = gen(t);

            // if( ! utils::starts_with( name, "Archaea" )) {
            //     return;
            // }

            if( counts_map.count( &t ) == 0 ) {
                throw std::runtime_error( "Taxon not in map: " + name );
            }
            if( entropies.count( &t ) != 0 ) {
                throw std::runtime_error( "Taxon already in entropy map: " + name );
            }

            auto& counts = counts_map[&t];
            entropies[ &t ] = averaged_entropy( counts, true, opt );
        };
        preorder_for_each( tax_cpy, calc_entropies );
        LOG_TIME << "finished entropy calculations";

        auto crop_list = split_taxonomy_by_entropy_with_target_size( tax_cpy, entropies, 1000 );
        LOG_DBG1 << "valid " << ( validated_splitted_taxonomy( tax_cpy, crop_list ) ? "1" : "0" );
        LOG_DBG1 << "total taxa leaf count: " << taxa_count_lowest_levels( tax_cpy );

        // -------------------------------------------------------------------------
        //     Output
        // -------------------------------------------------------------------------

        // Create a list of all (including parents) taxa that we want to keep.
        LOG_DBG1 << "split list size " << crop_list.size();
        LOG_DBG1 << "adding inner taxa";
        auto full_crop_list = fill_splitted_entropy_parents( crop_list );
        LOG_DBG1 << "full split list size " << full_crop_list.size();

        LOG_TIME << "write tax_entr_list_file";
        std::ofstream tax_entr_crop_file;
        std::ofstream tax_entr_all_file;
        tax_entr_crop_file.open( "/home/lucas/tax_entr_crop_det_" + std::to_string(i) );
        tax_entr_all_file.open( "/home/lucas/tax_entr_all_det_" + std::to_string(i) );

        auto print_crop_list = [&] ( Taxon const& t ) {
            tax_entr_all_file << std::string( taxon_level(t) * 4, ' ' );
            if( full_crop_list.count( &t ) > 0 ) {
                tax_entr_all_file << "X ";
            }
            tax_entr_all_file << t.name();
            if( entropies.count( &t ) > 0 ) {
                tax_entr_all_file << " (" + std::to_string( entropies.at( &t )) + ")";
            }
            tax_entr_all_file << "\n";

            if( full_crop_list.count( &t ) > 0 ) {
                tax_entr_crop_file << std::string( taxon_level(t) * 4, ' ' );
                tax_entr_crop_file << t.name();
                if( entropies.count( &t ) > 0 ) {
                    tax_entr_crop_file << " (" + std::to_string( entropies.at( &t )) + ")";
                }
                if( crop_list.count( &t ) > 0 ) {
                    tax_entr_crop_file << " subtree: " << total_taxa_count( t );
                    tax_entr_crop_file << " leaves: " << taxa_count_lowest_levels( t );
                }
                tax_entr_crop_file << "\n";
            }
        };
        preorder_for_each( tax_cpy, print_crop_list );

        tax_entr_crop_file.close();
        tax_entr_all_file.close();
        LOG_TIME << "finished tax_entr_list_file";

        auto total_taxa_cnt = count_splitted_taxonomy_total_size( tax_cpy, full_crop_list );
        LOG_DBG1 << "total count " << total_taxa_cnt;
        LOG_DBG1 << "total children:";
        for( auto const& base_taxon : tax_cpy ) {
            auto total_cnt = count_splitted_taxonomy_total_size( base_taxon, full_crop_list );
            LOG_DBG2 << base_taxon.name() << "\t" << total_cnt << "\t"
                     << static_cast<double>( total_cnt ) / static_cast<double>( total_taxa_cnt );
        }

        // remove_splitted_taxonomy_children( tax_cpy, crop_list );
        // auto total_taxa_cnt = total_taxa_count( tax_cpy );
        // LOG_DBG1 << "total count " << total_taxa_cnt;
        // LOG_DBG1 << "total children:";
        // for( auto const& base_taxon : tax_cpy ) {
        //     auto total_cnt = total_taxa_count( base_taxon );
        //     LOG_DBG2 << base_taxon.name() << "\t" << total_cnt << "\t"
        //              << static_cast<double>( total_cnt ) / static_cast<double>( total_taxa_cnt );
        // }
        //
        // LOG_DBG2 << "Archaea\t" << total_taxa_count( tax_cpy.get_child("Archaea") );
        // LOG_DBG2 << "Bacteria\t" << total_taxa_count( tax_cpy.get_child("Bacteria") );
        // LOG_DBG2 << "Eukaryota\t" << total_taxa_count( tax_cpy.get_child("Eukaryota") );

        // LOG_INFO << print_splitted_taxonomy( tax, crop_list, entropies );
    }

    return 0;

    // -------------------------------------------------------------------------
    //     Write entropy file.
    // -------------------------------------------------------------------------

    LOG_TIME << "write tax_entr_file";
    std::ofstream tax_entr_file;
    tax_entr_file.open( "/home/lucas/tax_entr" );
    auto print_entropy = [&]( Taxon& t ) {
        auto gen = TaxopathGenerator();
        auto name = gen(t);

        if( ! utils::starts_with( name, "Archaea" )) {
            return;
        }

        if( counts_map.count( &t ) == 0 ) {
            throw std::runtime_error( "Taxon not in map: " + name );
        }
        auto& counts = counts_map[&t];

        tax_entr_file << name;
        tax_entr_file << "\t" << taxon_level(t);
        tax_entr_file << "\t" << counts.added_sequences_count();

        tax_entr_file << "\t|";
        for( size_t i = 0; i < 4; ++i ) {
            auto opt = static_cast< SiteEntropyOptions >(i);
            tax_entr_file << "\t" << absolute_entropy( counts, opt );
        }
        tax_entr_file << "\t|";
        for( size_t i = 0; i < 4; ++i ) {
            auto opt = static_cast< SiteEntropyOptions >(i);
            tax_entr_file << "\t" << averaged_entropy( counts, false, opt );
        }
        tax_entr_file << "\t|";
        for( size_t i = 0; i < 4; ++i ) {
            auto opt = static_cast< SiteEntropyOptions >(i);
            tax_entr_file << "\t" << averaged_entropy( counts, true,  opt );
        }

        tax_entr_file << "\n";
    };
    preorder_for_each( tax, print_entropy );
    tax_entr_file.close();
    LOG_TIME << "finished tax_entr_file";

    return 0;

    // -------------------------------------------------------------------------
    //     Write consensus sequence file.
    // -------------------------------------------------------------------------

    LOG_TIME << "write tax_cons_file";
    std::ofstream tax_cons_file;
    tax_cons_file.open( "/home/lucas/tax_cons.fasta" );
    auto print_consensus = [&]( Taxon& t ) {
        auto gen = TaxopathGenerator();
        auto name = gen(t);

        if( counts_map.count( &t ) == 0 ) {
            throw std::runtime_error( "Taxon not in map: " + name );
        }
        auto const& counts = counts_map[&t];

        tax_cons_file << ">" << name << "\n";

        auto const sites = consensus_sequence_with_majorities( counts );
        for (size_t i = 0; i < sites.length(); i += 80) {
            tax_cons_file << sites.substr(i, 80) << "\n";
        }
        // tax_cons_file << "\n";
    };
    preorder_for_each( tax, print_consensus );
    tax_cons_file.close();
    LOG_TIME << "finished tax_cons_file";

    // -------------------------------------------------------------------------
    //     Write count matrices file.
    // -------------------------------------------------------------------------

    LOG_TIME << "write count_matrices_file";
    std::ofstream count_matrices_file;
    count_matrices_file.open( "/home/lucas/count_matrices" );
    auto print_count_matrices = [&]( Taxon& t ) {
        auto gen = TaxopathGenerator();
        auto name = gen(t);

        if( counts_map.count( &t ) == 0 ) {
            throw std::runtime_error( "Taxon not in map: " + name );
        }
        auto const& counts = counts_map[&t];

        count_matrices_file << name;
        count_matrices_file << "\n";
        for( size_t i = 0; i < counts.length(); ++i ) {
            for( size_t j = 0; j < counts.characters().size(); ++j ) {
                count_matrices_file << counts.count_at(i, j) << "\t";
            }
            count_matrices_file << "\n";
        }
        count_matrices_file << "\n";
    };
    preorder_for_each( tax, print_count_matrices );
    count_matrices_file.close();
    LOG_TIME << "finished count_matrices_file";

    LOG_INFO << "Finished.";
    return 0;
}
