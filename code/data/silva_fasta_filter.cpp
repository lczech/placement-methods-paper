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

#include "sequence/functions/codes.hpp"
#include "sequence/formats/fasta_input_iterator.hpp"
#include "sequence/formats/fasta_reader.hpp"
#include "sequence/formats/fasta_writer.hpp"
#include "sequence/functions/functions.hpp"
#include "sequence/sequence_set.hpp"

#include "taxonomy/formats/taxopath_generator.hpp"
#include "taxonomy/formats/taxopath_parser.hpp"
#include "taxonomy/functions/taxonomy.hpp"
#include "taxonomy/functions/taxopath.hpp"
#include "taxonomy/taxon.hpp"
#include "taxonomy/taxonomy.hpp"
#include "taxonomy/taxopath.hpp"

#include "utils/core/logging.hpp"
#include "utils/core/std.hpp"
#include "utils/io/input_stream.hpp"
#include "utils/tools/date_time.hpp"

#include <fstream>
#include <string>

using namespace genesis;
using namespace genesis::sequence;
using namespace genesis::taxonomy;

/**
 * Simple tool for extracting a part of the Silva sequences belonging to one taxonomic path.
 */
int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();

    std::string infile = "/path/to/silva/SILVA_123.1_SSURef_Nr99_tax_silva_full_align_trunc.fasta";
    // size_t ccnt = 0;

    std::ofstream out_stream;
    out_stream.open ("/path/to/output/archaea.fasta");
    auto writer = FastaWriter();

    // auto it = FastaInputIterator( ifs );
    auto reader = FastaReader();
    reader.to_upper( false );

    utils::InputStream cit { utils::make_unique<utils::FileInputSource>(infile) };

    Sequence seq;
    auto taxopath_parser = TaxopathParser();

    LOG_TIME << "Start reading at " << utils::current_time();
    size_t len = 0;
    size_t cnt = 0;
    // while( it ) {
    //     if( cnt % 10000 == 0 ) {
    //         LOG_TIME << "At sequence " << cnt;
    //     }
    //
    //     len = std::max( len, it->length() );
    //     ++cnt;
    //     ++it;
    // }

    // auto char_histogram = std::vector<size_t>( 128, 0 );

    while( reader.parse_sequence( cit, seq ) ) {
        if( cnt % 100000 == 0 ) {
        // if( cnt % 10000 == 0 ) {
            LOG_TIME << "At sequence " << cnt;
        }
        ++cnt;
        len = std::max( len, seq.length() );

        // for( auto const& site : seq ) {
        //     ++char_histogram[ site ];
        // }

        // if( seq.metadata().size() == 0 ) {
        //     continue;
        // }

        auto taxopath = taxopath_parser.from_string( seq.metadata() );
        if( taxopath.size() == 0 ) {
            LOG_DBG << "empty metadata at " << seq.label();
            continue;
        }

        if( taxopath[ 0 ] == "Archaea" ) {
            writer.write_sequence( seq, out_stream );
        }
    }

    // out_stream.close();
    LOG_TIME << "Done reading at " << utils::current_time();
    LOG_DBG << "count: " << cnt << ", longest: " << len;

    // for( size_t i = 0; i < char_histogram.size(); ++i ) {
    //     if( char_histogram[i] == 0 ) {
    //         continue;
    //     }
    //     LOG_INFO << static_cast<char>(i) << " " << char_histogram[i];
    // }

    LOG_INFO << "Finished.";
    return 0;
}
