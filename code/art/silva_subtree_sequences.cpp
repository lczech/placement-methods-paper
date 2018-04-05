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
 * @brief Write fasta files that contain the aligned sequences of the five clades that
 * we want to evaluate for the Russian doll approach.
 */
int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;
    LOG_INFO << "Started " << utils::current_time();

    std::string basedir = "path/to/data/silva/";

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

    // -------------------------------------------------------------------------
    //     Process all sequences.
    // -------------------------------------------------------------------------

    // Prepare output streams for the phyla
    std::unordered_map<std::string, std::ofstream> outfiles;
    for( auto p : phyla ) {
        outfiles[ p ].open( basedir + "subtree_sequences/Bacteria_" + p + ".fasta" );
    }

    // Count how many seqs were written to each phylum
    std::unordered_map<std::string, size_t> seq_counts;

    // Prepare sequence input.
    std::string aln_file = basedir + "600k_taxa.fasta";
    utils::InputStream aln_is { utils::make_unique<utils::FileInputSource>( aln_file ) };
    auto reader = FastaReader();
    Sequence seq;

    auto write_fasta_sequence = [] ( std::ofstream& out, std::string name, std::string sites ) {
        out << ">" << name << "\n";
        for (size_t i = 0; i < sites.length(); i += 80) {
            out << sites.substr(i, 80) << "\n";
        }
    };

    // Add to count objects along their taxonomic path.
    LOG_INFO << "Start reading";
    size_t seq_cnt = 0;
    while( reader.parse_sequence( aln_is, seq )) {
        replace_u_with_t( seq );

        if( seq_cnt % 100000 == 0 ) {
            LOG_INFO << "At sequence " << seq_cnt;
        }
        ++seq_cnt;

        // Check if the sequence is one of our processed ones.
        if( ! utils::starts_with( seq.label(), "SEQ_" )) {
            LOG_ERR << "sequence name invalid: " << seq.label();
            return 1;
        }

        // Get the name, excluding the prefix `SEQ_000000_`.
        // SEQ_000000_Bacteria_xyz
        // ^         ^          ^
        auto seq_name = seq.label().substr( 11 );
        if( ! utils::starts_with( seq_name, "Bacteria_" ) ) {
            continue;
        }
        // LOG_DBG << seq_name;

        // Get the part after bacteria
        // Bacteria_xyz_abc
        // ^         ^
        size_t phylum_end_underscore = seq_name.find_first_of( "_", 9 );
        if( phylum_end_underscore == std::string::npos ) {
            // LOG_DBG << "no underscore";
            continue;
        }

        auto const phylum = seq_name.substr( 9, phylum_end_underscore - 9 );
        if( outfiles.count( phylum ) > 0 ) {
            write_fasta_sequence( outfiles[phylum], seq.label(), seq.sites() );
            ++seq_counts[phylum];
        }
    }
    LOG_INFO << "Done reading";

    for( auto const& e : seq_counts ) {
        LOG_INFO << e.first << ": " << e.second;
    }

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
