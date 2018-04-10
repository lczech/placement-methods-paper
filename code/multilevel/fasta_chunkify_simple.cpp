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

// =================================================================================================
//     Main
// =================================================================================================

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.date = true;
    utils::Logging::details.time = true;
    LOG_INFO << "Started " << utils::current_time();

    // Check if the command line contains the right number of arguments.
    if (argc != 3) {
        throw std::runtime_error(
            "Need to provide a path to a dir with fsa/fasta files and a path to an output dir.\n"
        );
    }

    // In out dirs.
    auto const infile = std::string( argv[1] );
    auto const outdir = utils::trim_right( std::string( argv[2] ), "/") + "/";

    // Infile
    utils::InputStream instr( utils::make_unique< utils::FileInputSource >( infile ));
    auto fasta_it = FastaInputIterator( instr );

    // Output
    auto writer = FastaWriter();
    writer.enable_metadata(false);

    // Collect sequences for a chunk here.
    SequenceSet chunk;
    size_t chunk_count = 0;

    while( fasta_it ) {

        // Relabel: remove count at the end.
        // auto parts = utils::split( fasta_it->label(), "_", true);
        // if( parts.size() != 2 ) {
        //     LOG_WARN << "bad sequence label. has size " << parts.size();
        //     return 0;
        // }

        // Make a sequence with the correct label, and store it in the chunk.
        // auto seq = *fasta_it;
        // seq.label( parts[0] );
        // chunk.add( seq );

        chunk.add( *fasta_it );

        // If a chunk is full, flush it.
        if( chunk.size() == 50000 ) {
            writer.to_file( chunk, outdir + "chunk_" + utils::to_string(chunk_count) + ".fasta" );
            ++chunk_count;
            chunk.clear();
        }

        ++fasta_it;
    }

    // Flush the remaining chunk.
    writer.to_file( chunk, outdir + "chunk_" + utils::to_string(chunk_count) + ".fasta" );

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
