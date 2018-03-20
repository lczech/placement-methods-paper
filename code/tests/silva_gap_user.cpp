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
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/core/std.hpp"
#include "genesis/utils/io/input_stream.hpp"
#include "genesis/utils/tools/date_time.hpp"
#include "genesis/utils/math/bitvector.hpp"

#include <fstream>
#include <string>
#include <stdexcept>

using namespace genesis;
using namespace genesis::sequence;

int main( int argc, char** argv )
{
    if( argc != 2 ) {
        throw std::runtime_error( "argc != 2" );
    }
    // Activate logging.
    utils::Logging::log_to_stdout();

    std::string const silva_path = "path/to/silva_files";
    std::string infile = std::string( argv[1] );

    auto const bv_str = utils::file_read( silva_path + "/Silva_full_align_reduced.gaps" );
    if( bv_str.size() != 50000 ) {
        throw std::runtime_error( "bv_str.size() != 50000" );
    }
    auto bv = utils::Bitvector( 50000 );
    for( size_t i = 0; i < bv_str.size(); ++i ) {
        bv.set( i, ( bv_str[i] == 1 ));
    }
    LOG_DBG << "bv cnt " << bv.count();

    // auto it = FastaInputIterator( ifs );
    auto reader = FastaReader();
    reader.site_casing( FastaReader::SiteCasing::kToUpper );
    auto writer = FastaWriter();
    // writer.enable_metadata(true);

    std::ofstream out_stream;
    out_stream.open ( utils::file_filename(infile) + "_no_silva_gaps.fasta" );

    utils::InputStream cit { utils::make_unique<utils::FileInputSource>(infile) };
    Sequence seq;

    LOG_TIME << "Start reading at " << utils::current_time();
    size_t cnt = 0;
    while( reader.parse_sequence( cit, seq ) ) {
        if( cnt % 100000 == 0 ) {
            LOG_TIME << "At sequence " << cnt;
        }
        ++cnt;

        replace_u_with_t( seq );
        remove_sites( seq, bv );
        for( auto& c : seq ) {
            if( c == '.' ) {
                c = '-';
            }
        }
        writer.write_sequence( seq, out_stream );
    }
    LOG_TIME << "Done reading at " << utils::current_time();
    LOG_DBG << "seq cnt: " << cnt;
    out_stream.close();

    LOG_INFO << "Finished.";
    return 0;
}
