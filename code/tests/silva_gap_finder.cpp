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

#include "genesis/sequence.hpp"

#include "genesis/utils/core/logging.hpp"
#include "genesis/utils/core/std.hpp"
#include "genesis/utils/io/input_stream.hpp"
#include "genesis/utils/tools/date_time.hpp"
#include "genesis/utils/math/bitvector.hpp"

#include <fstream>
#include <string>

using namespace genesis;
using namespace genesis::sequence;

/**
 * @brief Take the 50k sites 600k sequences SILVA alignment, find the all-gap sites in it,
 * and write them out as a 50k bitvector file, `1` for gap sites, `0` for sites containing chars.
 */
int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();

    std::string const silva_path = "path/to/silva_files";
    std::string infile = silva_path + "/SILVA_123.1_SSURef_Nr99_tax_silva_full_align_trunc.fasta";

    // auto it = FastaInputIterator( ifs );
    auto reader = FastaReader();
    reader.site_casing( FastaReader::SiteCasing::kToUpper );

    utils::InputStream cit { utils::make_unique<utils::FileInputSource>(infile) };
    Sequence seq;

    auto bv = utils::Bitvector( 50000, true );

    auto lookup = utils::CharLookup<bool>( false );
    lookup.set_selection_upper_lower( nucleic_acid_codes_undetermined(), true );

    LOG_TIME << "Start reading at " << utils::current_time();
    size_t cnt = 0;
    while( reader.parse_sequence( cit, seq ) ) {
        if( cnt % 100000 == 0 ) {
            LOG_TIME << "At sequence " << cnt;
        }
        ++cnt;

        // replace_u_with_t( seq );
        auto const gaps = find_sites( seq, lookup );
        bv &= gaps;
    }

    LOG_TIME << "Done reading at " << utils::current_time();
    LOG_DBG << "seq cnt: " << cnt;
    LOG_DBG << "gap cnt: " << bv.count();

    std::ofstream out_stream;
    out_stream.open ( silva_path + "/Silva_full_align_reduced.gaps");
    for( size_t i = 0; i < bv.size(); ++i ) {
        if( bv[i] ) {
            out_stream.put( 1 );
        } else {
            out_stream.put( 0 );
        }
    }
    out_stream.close();

    LOG_INFO << "Finished.";
    return 0;
}
