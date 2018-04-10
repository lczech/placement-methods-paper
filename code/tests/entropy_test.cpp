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

#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>

using namespace genesis;
using namespace genesis::sequence;
using namespace genesis::taxonomy;

// =================================================================================================
//     Main
// =================================================================================================

/**
 * Calculate the entropy values for the test case shown in the art paper.
 */
int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.date = true;
    utils::Logging::details.time = true;
    LOG_INFO << "Started " << utils::current_time();

    (void) argc;
    (void) argv;

    SequenceSet set;
    // set.add( Sequence( "", "ACGT----A" ));
    // set.add( Sequence( "", "AAAA--AAA" ));
    // set.add( Sequence( "", "AACC-AACC" ));
    // set.add( Sequence( "", "AACC-AAGG" ));
    // set.add( Sequence( "", "AACC-AATT" ));

    FastaReader().from_file( "manuscripts/art/svg/tax_entropy_selection.fasta", set );
    // FastaReader().from_file( "manuscripts/art/svg/tax_entropy_seqs.fasta", set );
    replace_u_with_t( set );

    LOG_INFO << "sites consensus_sequence_with_threshold " << consensus_sequence_with_threshold( set, 0.80, true );
    LOG_INFO << "sites consensus_sequence_with_majorities " << consensus_sequence_with_majorities( set );
    // LOG_INFO << "sites " << consensus_sequence_with_threshold( set, 0.90, false );;

    // auto counts = SequenceCounts( "ACGT", set[0].length() );
    // counts.add_sequences( set );
    //
    // // LOG_INFO << "sites " << consensus_sequence_with_threshold( counts, 0.90, true );;
    // // LOG_INFO << "sites " << consensus_sequence_with_threshold( counts, 0.90, false );;
    //
    // auto opt = SiteEntropyOptions::kIncludeGaps;
    // auto entr = averaged_entropy( counts, false, opt );
    // LOG_INFO << "entr " << entr;
    //
    // for( size_t i = 0; i < counts.length(); ++i ) {
    //     LOG_DBG << "site " << i << ": " << site_entropy( counts, i, opt );
    // }

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
