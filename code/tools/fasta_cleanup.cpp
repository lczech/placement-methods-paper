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

using namespace genesis;
using namespace genesis::sequence;

int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    LOG_INFO << "Started " << utils::current_time();

    // Check if the command line contains the right number of arguments.
    if (argc != 2) {
        throw std::runtime_error(
            "Need to provide a fsa/fasta file.\n"
        );
    }

    auto infile  = std::string( argv[1] );
    auto outfile = utils::file_filename(infile) + ".upper";

    auto writer = FastaWriter();
    // writer.enable_metadata(false);

    std::ofstream out_stream;
    out_stream.open ( utils::file_filename(infile) + "_clean.fasta" );

    auto reader = FastaReader();
    reader.site_casing( FastaReader::SiteCasing::kToUpper );

    utils::InputStream cit { utils::make_unique<utils::FileInputSource>(infile) };
    Sequence seq;

    while( reader.parse_sequence( cit, seq ) ) {
        replace_u_with_t( seq );
        remove_characters( seq, "BDEFHIJKLMNOPQRSUVWXYZ" );
        writer.write_sequence( seq, out_stream );
    }
    out_stream.close();

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
