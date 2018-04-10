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
//     Main
// =================================================================================================

int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    LOG_INFO << "Started " << utils::current_time();

    // Check if the command line contains the right number of arguments.
    if (argc != 3) {
        throw std::runtime_error(
            "Need aln file."
        );
    }

    auto const aln_file      = std::string( argv[1] );
    auto const out_file      = std::string( argv[2] );

    // -------------------------------------------------------------------------
    //     Read all sequences.
    // -------------------------------------------------------------------------

    utils::InputStream aln_is { utils::make_unique<utils::FileInputSource>( aln_file ) };
    auto reader = FastaReader();
    Sequence seq;

    auto count_vec = std::vector<size_t>();

    size_t seq_cnt = 0;
    while( reader.parse_sequence( aln_is, seq )) {

        if( count_vec.empty() ) {
            count_vec = std::vector<size_t>( seq.size(), 0 );
            LOG_INFO << "alignment length " << seq.size();
        } else if( count_vec.size() != seq.size() ) {
            throw std::runtime_error( "Wrong seq size" );
        }

        // if( seq_cnt % 100000 == 0 ) {
        //     LOG_TIME << "At sequence " << seq_cnt;
        // }
        ++seq_cnt;

        for( size_t i = 0; i < seq.size(); ++i ) {
            if( seq[i] != '-' && seq[i] != '.' ) {
                ++count_vec[i];
            }
        }
    }
    LOG_TIME << "number of sequences " << seq_cnt;

    size_t min = *std::min_element( count_vec.begin(), count_vec.end() );
    size_t max = *std::max_element( count_vec.begin(), count_vec.end() );

    LOG_INFO << "number of non-gaps per site: min " << min << " max " << max;
    // LOG_INFO << "size " << count_vec.size();

    auto hist = std::vector<size_t>();

    size_t line_cnt = 0;
    std::ofstream ostr_all( out_file + "all.txt" );
    for( auto e : count_vec ) {
        ostr_all << e;

        ++line_cnt;
        if( line_cnt % 25 == 0 ) {
            ostr_all << "\n";
        } else {
            ostr_all << "\t";
        }

        if( hist.size() <= e ) {
            hist.resize( e+1, 0 );
        }
        ++hist[e];
    }

    // LOG_INFO << "hist size " << hist.size();

    line_cnt = 0;
    std::ofstream ostr_hist( out_file + "hist.txt" );
    for( auto e : hist ) {
        ostr_hist << e;

        ++line_cnt;
        if( line_cnt % 25 == 0 ) {
            ostr_hist << "\n";
        } else {
            ostr_hist << "\t";
        }
    }

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
