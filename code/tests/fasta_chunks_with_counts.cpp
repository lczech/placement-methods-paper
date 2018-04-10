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
 * @brief Take a chunk dir and a map dir (as produced for our data pipeline),
 * and prepend each sequence in the chunks with the total number of counts for that sequence
 * from all map files. Write the result into one huuuuge fasta file.
 */
int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;

    utils::Options::get().number_of_threads( 4 );
    // LOG_BOLD << utils::Options::get().info();
    // LOG_BOLD;

    LOG_INFO << "Started";

    // Check if the command line contains the right number of arguments.
    if (argc != 4) {
        throw std::runtime_error(
            "Need to provide three arguments.\n"
        );
    }

    // In out dirs.
    auto map_dir   = utils::dir_normalize_path( std::string( argv[1] ));
    auto chunk_dir = utils::dir_normalize_path( std::string( argv[2] ));
    auto out_file  = std::string( argv[3] );

    // Reader for map files.
    utils::CsvReader csv_reader;
    csv_reader.separator_chars("\t");

    // Find map files.
    auto map_filenames = utils::dir_list_files( map_dir );
    std::sort( map_filenames.begin(), map_filenames.end() );
    LOG_INFO << "processing " << map_filenames.size() << " maps";

    // Store map counts for all sequences.
    std::unordered_map<std::string, size_t> sequence_counts;

    // Process all map files.
    size_t cnt = 0;
    for( auto const& map_filename : map_filenames ) {

        if( cnt % 1000 == 0 ) {
            LOG_INFO << "at " << cnt;
        }
        ++cnt;

        // Read mapfile, add counts.
        auto const table = csv_reader.from_file( map_dir + map_filename );
        for( auto const& row : table ) {
            if( row.size() != 3 ) {
                LOG_WARN << "invalid row in map file " << map_filename;
                continue;
            }

            sequence_counts[ row[0] ] += stoi( row[1] );
        }
    }
    LOG_INFO << "at " << cnt;

    // Prepare fasta in and out
    auto fasta_reader = FastaReader();
    std::ofstream huge_fasta_of( out_file );
    auto fasta_out = FastaOutputIterator( huge_fasta_of );

    // Find seq files.
    auto seq_filenames = utils::dir_list_files( chunk_dir );
    std::sort( seq_filenames.begin(), seq_filenames.end() );
    LOG_INFO << "processing " << seq_filenames.size() << " chunks";

    // Process all chunk files.
    cnt = 0;
    for( auto const& seq_filename : seq_filenames ) {

        if( cnt % 200 == 0 ) {
            LOG_INFO << "at " << cnt;
        }
        ++cnt;

        // Read file
        auto seqs = fasta_reader.from_file( chunk_dir + seq_filename );
        remove_characters( seqs, "N" );

        for( auto& seq : seqs ) {
            if( sequence_counts.count( seq.label() ) == 0 ) {
                LOG_WARN << "cannot find sequence " << seq.label() << " in maps";
                continue;
            }
            auto const new_label = seq.label() + "_" + std::to_string( sequence_counts[ seq.label() ] );
            seq.label( new_label );
            fasta_out = seq;
        }
    }
    LOG_INFO << "at " << cnt;

    LOG_INFO << "Finished";
    return 0;
}
