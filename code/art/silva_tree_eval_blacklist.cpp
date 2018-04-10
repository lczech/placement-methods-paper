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
//      Main
// =================================================================================================

/**
 * @brief Write a blacklist file that lists all sequence names (SEQ_xxxxxx_) that we want to exclude
 * from evaluation.
 *
 * Those are:
 *
 *   * all sequences form the Sativa mislabel list
 *   * all sequences that contain "incertae", "unclassified" or "unknown" in their taxopath
 *
 * This list has 25,910 entries out of 598,470 sequences, or 4.3% of the data.
 */
int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;
    LOG_INFO << "Started " << utils::current_time();

    // Prepare sequence input.
    std::string aln_file = "data/silva/SILVA_123.1_SSURef_Nr99_tax_silva_full_align_trunc.fasta";
    utils::InputStream aln_is { utils::make_unique<utils::FileInputSource>( aln_file ) };
    auto reader = FastaReader();
    Sequence seq;

    // Mislabel list
    auto mislabel_table = utils::CsvReader().from_file( "data/silva/mislabels.txt" );
    std::unordered_map<std::string, size_t> mislabel_map;
    for( auto const& line : mislabel_table ) {
        if( line.size() != 1 ) {
            LOG_ERR << "line size " << line.size();
            return 1;
        }
        if( mislabel_map.count( line[0] ) > 0 ) {
            LOG_WARN << "already in mislable map: " << line[0];
        }
        mislabel_map[ line[0] ] = 0;
    }

    // Output
    std::string outdir = "data/silva/";
    std::ofstream outfile;
    utils::file_output_stream( outdir + "blacklist.txt", outfile );

    size_t total_blacklist_count = 0;
    size_t cnt_mislabel = 0;
    size_t cnt_incertae = 0;
    size_t cnt_unclassified = 0;
    size_t cnt_unknown = 0;

    // Add to count objects along their taxonomic path.
    LOG_DBG << "Start reading at " << utils::current_time();
    size_t seq_cnt = 0;
    while( reader.parse_sequence( aln_is, seq )) {
        if( seq_cnt % 100000 == 0 ) {
            LOG_DBG << "At sequence " << seq_cnt;
        }
        ++seq_cnt;

        // if( seq.metadata() == "" ) {
        //     throw std::runtime_error( "Empty metadata." );
        // }

        bool bad = false;

        // Check if the sequence name is in the mislabel list.
        auto const ncbi = seq.label().substr(0, seq.label().find_first_of('.'));
        if( mislabel_map.count( ncbi ) > 0 ) {
            ++mislabel_map[ncbi];
            ++cnt_mislabel;
            ++total_blacklist_count;
            outfile << "SEQ_" << utils::to_string_leading_zeros( seq_cnt, 6 ) << "\t" << ncbi << "\n";
            bad = true;
        }

        // Check if the sequence metadata contains any of the bad words.
        auto const lower_label = utils::to_lower( seq.label() );
        if( lower_label.find( "incertae" ) != std::string::npos ) {
            ++cnt_incertae;
            if( !bad ) {
                bad = true;
                ++total_blacklist_count;
                outfile << "SEQ_" << utils::to_string_leading_zeros( seq_cnt, 6 ) << "\t" << "incertae" << "\n";
            }
        }
        if( lower_label.find( "unclassified" ) != std::string::npos ) {
            ++cnt_unclassified;
            if( !bad ) {
                bad = true;
                ++total_blacklist_count;
                outfile << "SEQ_" << utils::to_string_leading_zeros( seq_cnt, 6 ) << "\t" << "unclassified" << "\n";
            }
        }
        if( lower_label.find( "unknown" ) != std::string::npos ) {
            ++cnt_unknown;
            if( !bad ) {
                bad = true;
                ++total_blacklist_count;
                outfile << "SEQ_" << utils::to_string_leading_zeros( seq_cnt, 6 ) << "\t" << "unknown" << "\n";
            }
        }
    }
    LOG_DBG << "Done reading at " << utils::current_time();

    for( auto const& e : mislabel_map ) {
        if( e.second != 1 ) {
            LOG_DBG << "mislabel " << e.first << " " << e.second;
        }
    }

    LOG_INFO << "cnt_mislabel     " << cnt_mislabel;
    LOG_INFO << "cnt_incertae     " << cnt_incertae;
    LOG_INFO << "cnt_unclassified " << cnt_unclassified;
    LOG_INFO << "cnt_unknown      " << cnt_unknown;

    LOG_INFO << "total_blacklist_count " << total_blacklist_count;

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
