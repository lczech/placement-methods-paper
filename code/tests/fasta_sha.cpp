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

#include "genesis.hpp"

#include <fstream>
#include <string>
#include <unordered_map>

using namespace genesis;
using namespace genesis::sequence;

/**
 * Test if there is a sha hash collision in a dataset.
 */
int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    LOG_INFO << "Started " << utils::current_time();
    size_t count = 0;

    // Check if the command line contains the right number of arguments.
    if (argc != 2) {
        throw std::runtime_error(
            "Need to provide a path to a dir with fsa/fasta files.\n"
        );
    }

    // In out dirs.
    // std::string indir = "data/hmp/links";
    auto indir = std::string( argv[1] );

    // Reading and writing fasta.
    auto reader = FastaReader();
    reader.to_upper( true );

    std::unordered_map< std::string, std::string > hashmap;

    size_t seqcnt = 0;

    // Process all fasta files, filter them, chunk them.
    auto fasta_filenames = utils::dir_list_files(indir);
    for( auto const& fasta_filename : fasta_filenames ) {
        if( count > 0 && count % 500 == 0 ) {
            LOG_DBG << "at " << count;
            LOG_INFO << "read " << seqcnt << " seqs, " << hashmap.size() << " uniq";
        }
        ++count;

        // Check file type.
        if( ! utils::ends_with( fasta_filename, ".fsa" ) && ! utils::ends_with( fasta_filename, ".fasta" ) ) {
            LOG_WARN << "skip invalid filename " << fasta_filename;
            continue;
        }

        auto filename = utils::file_filename(fasta_filename);

        // Read.
        auto set = SequenceSet();
        reader.from_file( indir + "/" + fasta_filename, set );

        // sha
        for( auto const& seq : set ) {
            auto hash = utils::SHA1::from_string_hex( seq.sites() );

            if( hashmap.count(hash) == 0 ) {
                hashmap[hash] = seq.sites();
            } else if( hashmap[hash] != seq.sites() ) {
                LOG_WARN << "collision at hash " << hash;
                LOG_DBG1 << hashmap[hash];
                LOG_DBG1 << seq.sites();
            }

            ++seqcnt;
        }
    }

    LOG_INFO << "at " << count;
    LOG_INFO << "read " << seqcnt << " seqs, " << hashmap.size() << " uniq";

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
