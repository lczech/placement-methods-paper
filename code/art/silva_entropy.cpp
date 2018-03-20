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

#include "genesis/genesis.hpp"

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

// =================================================================================================
//     Write Entropy
// =================================================================================================

void write_entropy_files( Taxonomy const& tax, std::string outdir )
{
    LOG_TIME << "Writing entropy files.";

    std::ofstream tax_entr_prune_file;
    std::ofstream tax_entr_all_file;
    std::ofstream tax_entr_subtree_list_file;
    tax_entr_prune_file.open( outdir + "/tax_entr_prune" );
    tax_entr_all_file.open( outdir + "/tax_entr_all" );
    tax_entr_subtree_list_file.open( outdir + "/tax_entr_subtree_list_file" );

    auto print_prune_list = [&] ( Taxon const& t ) {
        tax_entr_all_file << std::string( taxon_level(t) * 4, ' ' );
        auto tot_tax_cnt = total_taxa_count(t);

        // if( t.data<EntropyTaxonData>().status != EntropyTaxonData::PruneStatus::kOutside ) {
        //     tax_entr_all_file << "X ";
        // }

        if( t.data<EntropyTaxonData>().status == EntropyTaxonData::PruneStatus::kOutside ) {
            tax_entr_all_file << "O ";
            if( tot_tax_cnt > 0 ) {
                tax_entr_subtree_list_file << "status O";
            }
        }
        if( t.data<EntropyTaxonData>().status == EntropyTaxonData::PruneStatus::kBorder ) {
            tax_entr_all_file << "B ";
            if( tot_tax_cnt > 0 ) {
                tax_entr_subtree_list_file << "status B";
            }
        }
        if( t.data<EntropyTaxonData>().status == EntropyTaxonData::PruneStatus::kInside ) {
            tax_entr_all_file << "I ";
            if( tot_tax_cnt > 0 ) {
                tax_entr_subtree_list_file << "status I";
            }
        }

        tax_entr_all_file << t.name();
        tax_entr_all_file << " (" + std::to_string( t.data<EntropyTaxonData>().entropy ) + ", "
                           + std::to_string(gap_site_count(t.data<EntropyTaxonData>().counts)) + ")";
        tax_entr_all_file << "\n";

        if( t.data<EntropyTaxonData>().status != EntropyTaxonData::PruneStatus::kOutside ) {
            tax_entr_prune_file << std::string( taxon_level(t) * 4, ' ' );
            tax_entr_prune_file << t.name();
            tax_entr_prune_file << " (" + std::to_string( t.data<EntropyTaxonData>().entropy) + ")";
            if( t.data<EntropyTaxonData>().status == EntropyTaxonData::PruneStatus::kBorder ) {
                tax_entr_prune_file << " subtree: " << tot_tax_cnt;
                tax_entr_prune_file << " leaves: " << taxa_count_lowest_levels( t );
                if( tot_tax_cnt > 0 ) {
                    tax_entr_prune_file << " <========";
                }
            }
            tax_entr_prune_file << "\n";
        }

        if( tot_tax_cnt > 0 ) {
            tax_entr_subtree_list_file
            << " subtree " << tot_tax_cnt
            << " leaves " << taxa_count_lowest_levels( t ) << "\n";
        }
    };
    preorder_for_each( tax, print_prune_list );

    tax_entr_prune_file.close();
    tax_entr_all_file.close();
    tax_entr_subtree_list_file.close();
    LOG_TIME << "finished tax_entr_list_file";

    auto total_taxa_cnt = count_taxa_with_prune_status( tax, EntropyTaxonData::PruneStatus::kBorder )
                        + count_taxa_with_prune_status( tax, EntropyTaxonData::PruneStatus::kInside );

    LOG_DBG1 << "total pruned inner/border count " << total_taxa_cnt;
    LOG_DBG1 << "of total taxa count " << total_taxa_count( tax );
    LOG_DBG1 << "total children:";
    for( auto const& base_taxon : tax ) {
        auto total_cnt = count_taxa_with_prune_status( base_taxon, EntropyTaxonData::PruneStatus::kBorder )
                       + count_taxa_with_prune_status( base_taxon, EntropyTaxonData::PruneStatus::kInside );

        LOG_DBG2 << base_taxon.name() << "\t" << total_cnt << "\t"
                 << static_cast<double>( total_cnt ) / static_cast<double>( total_taxa_cnt );
    }

    // */

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

    // return 0;

    // -------------------------------------------------------------------------
    //     Tabluated Entropy information file.
    // -------------------------------------------------------------------------

    LOG_TIME << "Writing tabulated tax entropy files.";
    std::ofstream tab_all_file;
    std::ofstream tab_pruned_file;
    tab_all_file.open( outdir + "/tax_all_entr.csv" );
    tab_pruned_file.open( outdir + "/tax_pruned_entr.csv" );

    tab_all_file    << "Name\tStatus\tLevel\tTotal Count\tChildren\tLeaves\tSequences\tEntropy\n";
    tab_pruned_file << "Name\tLevel\tTotal Count\tChildren\tLeaves\tSequences\tEntropy\n";

    auto print_entropy = [&]( Taxon const& t ) {
        auto gen = TaxopathGenerator();
        auto name = gen(t);

        auto const& counts = t.data<EntropyTaxonData>().counts;

        tab_all_file << name;
        tab_all_file << "\t" << EntropyTaxonData::status_abbreviation( t.data<EntropyTaxonData>().status );
        tab_all_file << "\t" << taxon_level(t);
        tab_all_file << "\t" << total_taxa_count(t);
        tab_all_file << "\t" << t.size();
        tab_all_file << "\t" << taxa_count_lowest_levels(t);
        tab_all_file << "\t" << counts.added_sequences_count();
        tab_all_file << "\t" << t.data<EntropyTaxonData>().entropy;
        tab_all_file << "\n";

        if( t.data<EntropyTaxonData>().status == EntropyTaxonData::PruneStatus::kBorder && t.size() > 0 ) {
            tab_pruned_file << name;
            tab_pruned_file << "\t" <<  EntropyTaxonData::status_abbreviation( t.data<EntropyTaxonData>().status );
            tab_pruned_file << "\t" << taxon_level(t);
            tab_pruned_file << "\t" << total_taxa_count(t);
            tab_pruned_file << "\t" << t.size();
            tab_pruned_file << "\t" << taxa_count_lowest_levels(t);
            tab_pruned_file << "\t" << counts.added_sequences_count();
            tab_pruned_file << "\t" << t.data<EntropyTaxonData>().entropy;
            tab_pruned_file << "\n";
        }
    };
    preorder_for_each( tax, print_entropy );
    tab_all_file.close();
    tab_pruned_file.close();
    LOG_TIME << "finished tax_all_entr";

}

// =================================================================================================
//     Write Taxonomy
// =================================================================================================

void write_taxonomy_files( Taxonomy const& tax, std::string outdir )
{
    LOG_DBG << "Writing pruned taxonomy and taxonomic assignment files.";

    std::ofstream tax_pruned_tax;
    std::ofstream tax_assign;
    tax_pruned_tax.open( outdir + "/tax_pruned.txt" );
    tax_assign.open( outdir + "/tax_assign.txt" );

    // std::ofstream tax_leaves_tax;
    std::ofstream tax_leaves_assign;
    // tax_leaves_tax.open( outdir + "/tax_tax_leaves.txt" );
    tax_leaves_assign.open( outdir + "/tax_assign_leaves.txt" );

    std::string subtax_dir = outdir + "/sub_taxonomies";
    utils::dir_create( subtax_dir );

    auto print_pruned_tax = [&]( Taxon const& t ) {
        auto gen = TaxopathGenerator();
        auto name = gen(t);
        auto san_name = sanitize_label( name );
        auto san_file = utils::sanitize_filname( san_name );

        if( t.data<EntropyTaxonData>().status != EntropyTaxonData::PruneStatus::kOutside ) {
            tax_pruned_tax << name << "\n";
        }

        if( t.data<EntropyTaxonData>().status == EntropyTaxonData::PruneStatus::kBorder ) {
            tax_assign << san_name << "\t" << name << "\n";
        }

        if( t.size() == 0 && t.data<EntropyTaxonData>().counts.added_sequences_count() > 0 ) {
            // tax_leaves_tax << name << "\n";
            tax_leaves_assign << san_name << "\t" << name << "\n";
        }

        // subtree taxonomies
        if( t.data<EntropyTaxonData>().status == EntropyTaxonData::PruneStatus::kBorder && total_taxa_count(t) > 0 ) {

            std::ofstream sub_pruned_tax;
            std::ofstream sub_assign;
            sub_pruned_tax.open( subtax_dir + "/" + san_file + "_tax.txt" );
            sub_assign.open( subtax_dir + "/" + san_file + "_assign.txt" );

            auto print_subtree_taxonomy = [&]( Taxon const& subt ) {

                // Only use those who actually have sequence data.
                if( subt.data<EntropyTaxonData>().counts.added_sequences_count() == 0 ) {
                    return;
                }

                // only use leaves.
                if( subt.size() > 0 ) {
                    return;
                }

                auto sub_name = gen(subt);
                auto sub_san_name = sanitize_label( sub_name );

                sub_pruned_tax << sub_name << "\n";
                sub_assign << sub_san_name << "\t" << sub_name << "\n";
            };
            preorder_for_each( t, print_subtree_taxonomy );

            sub_pruned_tax.close();
            sub_assign.close();
        }
    };
    preorder_for_each( tax, print_pruned_tax );
    tax_pruned_tax.close();
    tax_assign.close();
    // tax_leaves_tax.close();
    tax_leaves_assign.close();
    LOG_DBG << "done with writing pruned taxonomy.";
}

// =================================================================================================
//     Write Sequences
// =================================================================================================

void write_sequence_files( Taxonomy const& tax, std::string outdir )
{
    LOG_TIME << "Write Sequence files.";

    std::string subcons_dir = outdir + "/sub_alignments";
    utils::dir_create( subcons_dir );

    std::ofstream tax_cons_file_all;
    std::ofstream tax_cons_file_leaves;
    std::ofstream tax_cons_file_border;
    std::ofstream tax_cons_file_border_amb;
    std::ofstream tax_cons_file_selected;
    tax_cons_file_all.open( outdir + "/tax_cons_all.fasta" );
    tax_cons_file_leaves.open( outdir + "/tax_cons_leaves.fasta" );
    tax_cons_file_border.open( outdir + "/tax_cons_border.fasta" );
    tax_cons_file_border_amb.open( outdir + "/tax_cons_border_amb.fasta" );
    tax_cons_file_selected.open( outdir + "/tax_cons_selected.fasta" );

    auto write_fasta_sequence = [] ( std::ofstream& out, std::string name, std::string sites ) {
        out << ">" << name << "\n";
        for (size_t i = 0; i < sites.length(); i += 80) {
            out << sites.substr(i, 80) << "\n";
        }
    };

    auto print_consensus = [&]( Taxon const& t ) {
        auto gen = TaxopathGenerator();
        auto name = sanitize_label( gen(t) );
        auto san_file = utils::sanitize_filname( name );

        auto const& counts = t.data<EntropyTaxonData>().counts;
        auto const sites = consensus_sequence_with_threshold( counts, 0.90 );

        write_fasta_sequence( tax_cons_file_all, name, sites );

        if( t.data<EntropyTaxonData>().status == EntropyTaxonData::PruneStatus::kBorder ) {
            write_fasta_sequence( tax_cons_file_border, name, sites );
        }
        if( t.data<EntropyTaxonData>().status != EntropyTaxonData::PruneStatus::kOutside ) {
            write_fasta_sequence( tax_cons_file_selected, name, sites );
        }
        if( t.size() == 0 && counts.added_sequences_count() > 0 ) {
            write_fasta_sequence( tax_cons_file_leaves, name, sites );
        }

        if( t.data<EntropyTaxonData>().status == EntropyTaxonData::PruneStatus::kBorder ) {
            auto const cons_amb = consensus_sequence_with_ambiguities( counts, 0.75 );
            write_fasta_sequence( tax_cons_file_border_amb, name, cons_amb );
        }

        // subtrees
        if( t.data<EntropyTaxonData>().status == EntropyTaxonData::PruneStatus::kBorder && total_taxa_count(t) > 0 ) {

            std::ofstream sub_cons_file;
            sub_cons_file.open( subcons_dir + "/" + san_file + ".fasta" );

            auto print_subtree_cons = [&]( Taxon const& subt ) {
                auto const& sub_counts = subt.data<EntropyTaxonData>().counts;

                // only use those which have actually sequence data.
                if( sub_counts.added_sequences_count() == 0 ) {
                    return;
                }

                // only use leaves.
                if( subt.size() > 0 ) {
                    return;
                }

                auto sub_name = gen(subt);
                auto sub_san_name = sanitize_label( sub_name );
                auto const sub_sites = consensus_sequence_with_threshold( sub_counts, 0.90 );

                write_fasta_sequence( sub_cons_file, sub_san_name, sub_sites );
            };
            preorder_for_each( t, print_subtree_cons );

            sub_cons_file.close();
        }

        // tax_cons_file << "\n";
    };
    preorder_for_each( tax, print_consensus );
    tax_cons_file_all.close();
    tax_cons_file_leaves.close();
    tax_cons_file_border.close();
    tax_cons_file_selected.close();

    LOG_TIME << "finished sequence files";

    // */

    // -------------------------------------------------------------------------
    //     Write count matrices file.
    // -------------------------------------------------------------------------

    /*

    LOG_TIME << "write count_matrices_file";
    std::ofstream count_matrices_file;
    count_matrices_file.open( outdir + "/count_matrices" );
    auto print_count_matrices = [&]( Taxon& t ) {
        auto gen = TaxopathGenerator();
        auto name = gen(t);

        auto const& counts = t.data<EntropyTaxonData>().counts;

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
    */
}

// =================================================================================================
//     Main
// =================================================================================================

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    LOG_INFO << "Started " << utils::current_time();

    // -------------------------------------------------------------------------
    //     Input Output
    // -------------------------------------------------------------------------

    // Check if the command line contains the right number of arguments.
    if( argc != 4 ) {
        throw std::runtime_error(
            "Need to provide three arguments: taxonomy file, alignment file and output dir."
        );
    }

    // std::string tax_file = "data/silva/tax_slv_ssu_123.1.txt";
    // std::string aln_file = "data/silva/SILVA_123.1_SSURef_Nr99_tax_silva_full_align_trunc.fasta";
    // std::string outdir = "tax_entropy";

    auto const tax_file = std::string( argv[1] );
    auto const aln_file = std::string( argv[2] );
    auto outdir = utils::dir_normalize_path( std::string( argv[3] ));
    utils::dir_create(outdir);

    // -------------------------------------------------------------------------
    //     Read and Prepare Taxonomy.
    // -------------------------------------------------------------------------


    // Load the taxonomy.
    LOG_TIME << "Reading taxonomy...";
    auto tax = Taxonomy();
    TaxonomyReader().from_file( tax_file, tax );
    sort_by_name( tax );

    LOG_TIME << "done";

    LOG_DBG << "total children:";
    auto tot_tax_cnt = total_taxa_count( tax );
    for( auto const& base_taxon : tax ) {
        auto total_cnt = total_taxa_count( base_taxon );

        LOG_DBG1 << base_taxon.name() << "\t" << total_cnt << "\t"
                 << static_cast<double>( total_cnt ) / static_cast<double>( tot_tax_cnt );
    }
    LOG_DBG1 << "Archaea " << total_taxa_count( tax.get_child("Archaea") );
    LOG_DBG1 << "Bacteria " << total_taxa_count( tax.get_child("Bacteria") );
    LOG_DBG1 << "Eukaryota " << total_taxa_count( tax.get_child("Eukaryota") );

    // auto print_leaf_node_counts = [] ( Taxon const& t ) {
    //     if( taxon_level(t) == 1 ) {
    //         LOG_INFO << TaxopathGenerator()(t) << "\t" << taxa_count_lowest_levels(t);
    //     }
    // };
    // preorder_for_each( tax, print_leaf_node_counts );
    // return 0;

    // Create a Sequence Count objeect for each taxon.
    LOG_TIME << "Preparing counts map...";
    auto add_sequence_counts_to_taxonomy = [&] ( Taxon& taxon ) {
        auto gen = TaxopathGenerator();
        auto name = gen(taxon);

        // if( ! utils::starts_with( name, "Archaea" )) {
        //     return;
        // }

        taxon.reset_data( EntropyTaxonData::create() );
        taxon.data<EntropyTaxonData>().counts = SiteCounts( "ACGT", 50000 );
    };
    preorder_for_each( tax, add_sequence_counts_to_taxonomy );
    LOG_TIME << "done";

    // -------------------------------------------------------------------------
    //     Read all sequences.
    // -------------------------------------------------------------------------

    // Prepare sequence input.
    utils::InputStream aln_is { utils::make_unique<utils::FileInputSource>( aln_file ) };
    auto reader = FastaReader();
    auto taxopath_parser = TaxopathParser();
    Sequence seq;

    // Add to count objects along their taxonomic path.
    LOG_TIME << "Start reading at " << utils::current_time();
    size_t seq_cnt = 0;
    while( reader.parse_sequence( aln_is, seq )) {
        replace_u_with_t( seq );

        if( seq_cnt % 100000 == 0 ) {
            LOG_TIME << "At sequence " << seq_cnt;
        }
        ++seq_cnt;

        // if( seq.metadata() == "" ) {
        //     throw std::runtime_error( "Empty metadata." );
        // }

        std::string taxopath_str;
        auto const delim = seq.label().find_first_of( " \t" );
        if( delim == std::string::npos ) {
            taxopath_str = seq.label();
        } else {
            taxopath_str = seq.label().substr( delim + 1 );
        }

        auto taxopath = taxopath_parser( taxopath_str );
        taxopath.pop_back();
        auto taxp = find_taxon_by_taxopath( tax, taxopath );
        if( taxp == nullptr ) {
            throw std::runtime_error( "Sequence taxon not in taxonomy: " + taxopath_str );
        }

        // if( taxopath[0] != "Archaea" ) {
        //     continue;
        // }

        auto cur_tax = taxp;
        do {
            cur_tax->data<EntropyTaxonData>().counts.add_sequence( seq );
            cur_tax = cur_tax->parent();
        } while( cur_tax != nullptr );
    }
    LOG_TIME << "Done reading at " << utils::current_time();

    // -------------------------------------------------------------------------
    //     Entropy Calculations.
    // -------------------------------------------------------------------------

    // auto opt = SiteEntropyOptions::kWeighted;
    auto opt = SiteEntropyOptions::kIncludeGaps;

    LOG_TIME << "entropy calculations";
    auto calc_entropies = [&]( Taxon& t ) {
        // auto gen = TaxopathGenerator();
        // auto name = gen(t);

        // if( ! utils::starts_with( name, "Archaea" )) {
        //     return;
        // }

        auto const& counts = t.data<EntropyTaxonData>().counts;
        t.data<EntropyTaxonData>().entropy = averaged_entropy( counts, false, opt );
    };
    preorder_for_each( tax, calc_entropies );
    LOG_TIME << "finished entropy calculations";

    // -------------------------------------------------------------------------
    //     Pruning
    // -------------------------------------------------------------------------

    LOG_INFO << "=============================================================";
    LOG_INFO << "General";

    PruneByEntropySettings prune_settings;
    // prune_settings.min_subtaxonomy_size = 25;
    prune_settings.max_subtaxonomy_size = 2000;
    prune_settings.min_border_level = 2;

    prune_by_entropy( tax, 2000, prune_settings );
    auto valid_tax = validate_pruned_taxonomy( tax );
    LOG_DBG1 << "valid " << ( valid_tax ? "1" : "nooooooooooooooo" );
    if( ! valid_tax ) {
        return 0;
    }

    // User output.
    LOG_INFO << "taxonomy size: " << total_taxa_count( tax );
    LOG_INFO << "leaf count:    " << taxa_count_lowest_levels( tax );
    LOG_INFO << "inside count:  " << count_taxa_with_prune_status(
        tax, EntropyTaxonData::PruneStatus::kInside
    );
    LOG_INFO << "border count:  " << count_taxa_with_prune_status(
        tax, EntropyTaxonData::PruneStatus::kBorder
    );
    LOG_INFO << "outside count: " << count_taxa_with_prune_status(
        tax, EntropyTaxonData::PruneStatus::kOutside
    );

    // -------------------------------------------------------------------------
    //     Write Output Files
    // -------------------------------------------------------------------------

    utils::dir_create( outdir + "/General" );

    LOG_INFO << "Writing entropy output " << utils::current_time();
    write_entropy_files( tax, outdir + "/General"  );

    LOG_INFO << "Writing taxonomy output " << utils::current_time();
    write_taxonomy_files( tax, outdir + "/General"  );

    LOG_INFO << "Writing sequence output " << utils::current_time();
    write_sequence_files( tax, outdir + "/General"  );

    // -------------------------------------------------------------------------
    //     Domains
    // -------------------------------------------------------------------------

    for( auto& domain : tax ) {
        LOG_INFO << "=============================================================";
        LOG_INFO << "Domain " << domain.name();

        auto& tax_domain = tax;

        PruneByEntropySettings prune_settings;
        // prune_settings.min_subtaxonomy_size = 25;
        // prune_settings.max_subtaxonomy_size = 2000;
        prune_settings.min_border_level = 2;

        prune_by_entropy( tax_domain, 3, prune_settings );
        prune_by_entropy( tax_domain.get_child( domain.name() ), 1800, prune_settings );

        auto valid_tax = validate_pruned_taxonomy( tax_domain );
        LOG_DBG1 << "valid " << ( valid_tax ? "1" : "nooooooooooooooo" );
        if( ! valid_tax ) {
            return 0;
        }

        // User output.
        LOG_INFO << "taxonomy size: " << total_taxa_count( tax_domain );
        LOG_INFO << "leaf count:    " << taxa_count_lowest_levels( tax_domain );
        LOG_INFO << "inside count:  " << count_taxa_with_prune_status(
            tax_domain, EntropyTaxonData::PruneStatus::kInside
        );
        LOG_INFO << "border count:  " << count_taxa_with_prune_status(
            tax_domain, EntropyTaxonData::PruneStatus::kBorder
        );
        LOG_INFO << "outside count: " << count_taxa_with_prune_status(
            tax_domain, EntropyTaxonData::PruneStatus::kOutside
        );

        // -------------------------------------------------------------------------
        //     Write Output Files
        // -------------------------------------------------------------------------

        utils::dir_create( outdir + "/" + domain.name() );

        LOG_INFO << "Writing entropy output " << utils::current_time();
        write_entropy_files( tax_domain, outdir + "/" + domain.name()  );

        LOG_INFO << "Writing taxonomy output " << utils::current_time();
        write_taxonomy_files( tax_domain, outdir + "/" + domain.name()  );

        LOG_INFO << "Writing sequence output " << utils::current_time();
        write_sequence_files( tax_domain, outdir + "/" + domain.name()  );
    }

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
