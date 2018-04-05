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

/**
 * @brief Write the alignments of consensus sequences and the taxonomies for the five clades
 * that we want to evaluate for the Russian doll approach.
 */
int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;
    LOG_INFO << "Started " << utils::current_time();

    std::string basedir = "data/silva/";

    /*
        Cyanobacteria
        94 taxa
        11,379 sequences

        Proteobacteria
        1,569 taxa
        200,083 sequences

        Firmicutes
        652 taxa
        138,517 sequences

        Bacteroidetes
        440 taxa
        49,174 sequences

        Actinobacteria
        494 taxa
        51,160 sequence
    */

    // -------------------------------------------------------------------------
    //     Read Taxonomy.
    // -------------------------------------------------------------------------

    // Load the taxonomy.
    LOG_INFO << "Reading taxonomy...";
    std::string tax_file = basedir + "tax_slv_ssu_123.1.txt";
    auto tax = Taxonomy();
    TaxonomyReader().from_file( tax_file, tax );
    sort_by_name( tax );
    LOG_INFO << "done";

    // -------------------------------------------------------------------------
    //     Prepare Sequence Counts
    // -------------------------------------------------------------------------

    auto const phyla = std::unordered_set<std::string>({
        "Cyanobacteria",
        "Proteobacteria",
        "Firmicutes",
        "Bacteroidetes",
        "Actinobacteria"
    });

    // Create a Sequence Count objeect for each taxon.
    LOG_INFO << "Preparing counts map...";
    auto add_sequence_counts_to_taxonomy = [&] ( Taxon& taxon ) {
        taxon.reset_data( EntropyTaxonData::create() );
        taxon.data<EntropyTaxonData>().counts = SiteCounts( "ACGT", 50000 );
    };
    auto& tax_bact = tax.get_child( "Bacteria" );
    for( auto const& phylum : phyla ) {
        if( ! tax_bact.has_child( phylum ) ) {
            LOG_ERR << "no phylum " << phylum;
            return 1;
        }
        add_sequence_counts_to_taxonomy( tax_bact.get_child( phylum ) );
        preorder_for_each( tax_bact.get_child( phylum ), add_sequence_counts_to_taxonomy );
    }
    LOG_INFO << "done";

    // -------------------------------------------------------------------------
    //     Write Taxonomies.
    // -------------------------------------------------------------------------

    LOG_INFO << "Writing taxonomies";
    for( auto const& phylum : phyla ) {
        std::ofstream tax_tips;
        std::ofstream tax_all;
        tax_tips.open( basedir + "subtree_alignments/Bacteria_" + phylum + ".tips.tax" );
        tax_all.open(  basedir + "subtree_alignments/Bacteria_" + phylum + ".all.tax" );

        auto print_taxa = [&]( Taxon const& taxon ){
            auto const name     = TaxopathGenerator()(taxon);
            auto const san_name = sanitize_label( name );

            if( taxon.size() == 0 ) {
                tax_tips << san_name << "\t" << name << "\n";
            }
            tax_all << san_name << "\t" << name << "\n";
        };
        preorder_for_each( tax_bact.get_child( phylum ), print_taxa );
    }
    LOG_INFO << "Done writing taxonomies";

    // -------------------------------------------------------------------------
    //     Read all sequences.
    // -------------------------------------------------------------------------

    // Prepare sequence input.
    std::string aln_file = basedir + "SILVA_123.1_SSURef_Nr99_tax_silva_full_align_trunc.fasta";
    utils::InputStream aln_is { utils::make_unique<utils::FileInputSource>( aln_file ) };
    auto reader = FastaReader();
    auto taxopath_parser = TaxopathParser();
    Sequence seq;

    // auto write_fasta_sequence = [] ( std::ofstream& out, std::string name, std::string sites ) {
    //     out << ">" << name << "\n";
    //     for (size_t i = 0; i < sites.length(); i += 80) {
    //         out << sites.substr(i, 80) << "\n";
    //     }
    // };

    // Count how many sequences are used for each phylum
    std::unordered_map<std::string, size_t> added_seq_cnt;
    std::unordered_map<std::string, size_t> inner_seq_cnt;

    // Add to count objects along their taxonomic path.
    LOG_INFO << "Start reading";
    size_t seq_cnt = 0;
    while( reader.parse_sequence( aln_is, seq )) {
        replace_u_with_t( seq );

        if( seq_cnt % 100000 == 0 ) {
            LOG_INFO << "At sequence " << seq_cnt;
        }
        ++seq_cnt;

        if( seq.label() == "" ) {
            throw std::runtime_error( "Empty label." );
        }

        // Only process Bacteria, and only the wanted phyla
        auto taxopath = taxopath_parser( seq.label() );
        taxopath.pop_back();
        if( taxopath.size() < 2 || taxopath[0] != "Bacteria" || phyla.count( taxopath[1] ) == 0 ) {
            continue;
        }

        // Get taxon for this sequence.
        auto taxp = find_taxon_by_taxopath( tax, taxopath );
        if( taxp == nullptr ) {
            throw std::runtime_error( "Sequence taxon not in taxonomy: " + seq.label() );
        }

        // Sanity check
        if( taxp->data_cast<EntropyTaxonData>() == nullptr ) {
            LOG_ERR << "no entropy data at " << TaxopathGenerator()( *taxp );
            return 1;
        }

        // Add the sequence to the counts for its taxon.
        taxp->data<EntropyTaxonData>().counts.add_sequence( seq );
        ++added_seq_cnt[taxopath[1]];
        if( taxp->size() > 0 ) {
            ++inner_seq_cnt[taxopath[1]];
        }
    }
    LOG_INFO << "Done reading";
    // LOG_INFO << "Processed useful sequences " << add_cnt;

    for( auto const& phylum : phyla ) {
        LOG_INFO << "Phylum " << phylum;
        LOG_INFO << "Added sequences " << added_seq_cnt[phylum];
        LOG_INFO << "Inner sequences " << inner_seq_cnt[phylum];
    }

    // -------------------------------------------------------------------------
    //     Write Consensus Sequences
    // -------------------------------------------------------------------------

    LOG_INFO << "Making consensus sequences";

    std::vector<std::string> const methods = {
        "majorities", "cavener", "threshold_0.5", "threshold_0.6",
        "threshold_0.7", "threshold_0.8", "threshold_0.9", "threshold_0.95"
    };

    for( auto const& method : methods ) {
        LOG_INFO << "===========================================";
        LOG_INFO << "Method " << method;

        for( auto const& phylum : phyla ) {
            LOG_INFO << "----------------------------------";
            LOG_INFO << "At phylum " << phylum;

            if( ! tax_bact.has_child( phylum ) ) {
                throw std::runtime_error("! tax_bact.has_child( phylum )");
            }

            // Add seqs to set
            SequenceSet set_tips;
            SequenceSet set_all;
            auto gen = TaxopathGenerator();
            size_t inner_cnt = 0;
            auto make_consensus_sequences = [&] ( Taxon& taxon ) {
                auto const& counts = taxon.data<EntropyTaxonData>().counts;
                auto const name = sanitize_label( gen(taxon) );

                // Some sanity checks
                if( taxon.size() == 0 && counts.added_sequences_count() == 0 ) {
                    LOG_WARN << "no counts for terminal taxon " << name;
                }
                if( taxon.size() > 0 && counts.added_sequences_count() > 0 ) {
                    // LOG_WARN << "seq count " << counts.added_sequences_count() << " for inner taxon " << name;
                    inner_cnt += counts.added_sequences_count();
                }

                // auto const sites = consensus_sequence_with_threshold( counts, 0.90 );
                std::string sites;
                if( method == "majorities" ) {
                    sites = consensus_sequence_with_majorities( counts );
                } else if( method == "cavener" ) {
                    sites = consensus_sequence_cavener( counts );
                } else if( method == "threshold_0.5" ) {
                    sites = consensus_sequence_with_threshold( counts, 0.5 );
                } else if( method == "threshold_0.6" ) {
                    sites = consensus_sequence_with_threshold( counts, 0.6 );
                } else if( method == "threshold_0.7" ) {
                    sites = consensus_sequence_with_threshold( counts, 0.7 );
                } else if( method == "threshold_0.8" ) {
                    sites = consensus_sequence_with_threshold( counts, 0.8 );
                } else if( method == "threshold_0.9" ) {
                    sites = consensus_sequence_with_threshold( counts, 0.9 );
                } else if( method == "threshold_0.95" ) {
                    sites = consensus_sequence_with_threshold( counts, 0.95 );
                } else {
                    throw std::runtime_error( "unknown method " + method );
                }

                if( taxon.size() == 0 ) {
                    set_tips.add({ name, sites });
                }
                set_all.add({ name, sites });
            };
            preorder_for_each( tax_bact.get_child( phylum ), make_consensus_sequences );

            LOG_DBG << "inner_cnt " << inner_cnt;
            LOG_DBG << "taxon path " << sanitize_label( gen(tax_bact.get_child( phylum )) );

            std::string const outdir = basedir + method + "/";
            utils::dir_create( outdir );

            // Write tip seqs
            std::string const filename = outdir + "Bacteria_" + phylum + ".tips.fasta";
            LOG_INFO << "Writing " << set_tips.size() << " seqs to " << filename;
            FastaWriter().to_file( set_tips, filename );

            // Write reduced tip seqs
            remove_gap_sites( set_tips );
            std::string const filename2 = outdir + "Bacteria_" + phylum + ".tips.reduced.fasta";
            LOG_INFO << "Writing " << set_tips.size() << " seqs to " << filename2;
            FastaWriter().to_file( set_tips, filename2 );

            // Write all seqs
            std::string const filename3 = outdir + "Bacteria_" + phylum + ".all.fasta";
            LOG_INFO << "Writing " << set_all.size() << " seqs to " << filename3;
            FastaWriter().to_file( set_all, filename3 );

            // Write reduced all seqs
            remove_gap_sites( set_all );
            std::string const filename4 = outdir + "Bacteria_" + phylum + ".all.reduced.fasta";
            LOG_INFO << "Writing " << set_all.size() << " seqs to " << filename4;
            FastaWriter().to_file( set_all, filename4 );
        }
    }


    LOG_INFO << "Done with consensus sequences";
    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
