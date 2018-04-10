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

using namespace genesis;
using namespace genesis::sequence;

/**
 * @brief Remove all gaps from a given set of fasta files.
 */
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
            "Need to provide a path to a dir with fsa/fasta files.\n"
        );
    }

    // In out dirs.
    auto indir = utils::trim_right( std::string( argv[1] ), "/") + "/";

    // Reading and writing fasta.
    auto reader = FastaReader();
    // reader.to_upper( true );

    auto writer = FastaWriter();
    // writer.enable_metadata(false);

    // Process all fasta files, filter them, chunk them.
    auto fasta_filenames = utils::dir_list_files( indir, ".*\\.(fsa|fasta)" );
    for( auto const& fasta_filename : fasta_filenames ) {

        auto filename = utils::file_filename(fasta_filename);

        // Read.
        auto set = SequenceSet();
        reader.from_file( indir + fasta_filename, set );

        remove_all_gaps( set );

        // write.
        writer.to_file( set, indir + filename + "_unaligned.fasta" );
    }

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
