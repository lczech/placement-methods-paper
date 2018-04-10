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
#include <unordered_set>

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::sequence;

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    LOG_INFO << "Started " << utils::current_time();

    // Check if the command line contains the right number of arguments.
    if (argc != 4) {
        throw std::runtime_error(
            "Need to provide three arguments.\n"
        );
    }

    // In out dirs.
    auto epadir = utils::dir_normalize_path( std::string( argv[1] ));
    // auto seqdir = utils::dir_normalize_path( std::string( argv[2] ));
    // auto outdir = utils::dir_normalize_path( std::string( argv[3] ));
    // utils::dir_create(outdir);

    // Prepare reader and writer
    JplaceReader jplace_reader;
    FastaWriter  fasta_writer;

    // Find fsa/fasta files.
    auto jplace_filenames = utils::dir_list_files( epadir, ".*/.*\\.jplace" );
    std::sort( jplace_filenames.begin(), jplace_filenames.end() );
    LOG_INFO << "processing " << jplace_filenames.size() << " files";

    for( auto const& jplace_filename : jplace_filenames ) {
        LOG_DBG << jplace_filename;
    }

    /*
    // Process all fasta files.
    size_t file_count = 0;
    for( auto const& jplace_filename : jplace_filenames ) {

        // Get the sample name from the map file name.
        auto map_filename_parts = utils::split( map_filename, "." );
        if( map_filename_parts.size() != 3 ) {
            LOG_WARN << "Weird map filename " << map_filename;
        }
        auto sample_name = map_filename_parts[0];

        // Progress output.
        if( file_count > 0 && file_count % 50 == 0 ) {
            LOG_DBG << "at map file " << file_count << " at " << utils::current_time();
        }
        ++file_count;

        // Read mapfile, get list of chunks that this sample needs, and a list of sequence names with multiplicities.
        utils::SortedVector<std::string> chunk_ids;
        std::unordered_map<std::string, size_t> sequence_names;
        auto table = csv_reader.from_file( mapdir + map_filename );
        for( auto const& row : table ) {
            if( row.size() != 3 ) {
                LOG_WARN << "invalid row in map file " << map_filename;
                continue;
            }
            chunk_ids.insert( row[2] );

            if( sequence_names.count( row[0] ) > 0 ) {
                LOG_WARN << "dup seq name in map file " << map_filename << ": " << row[0];
            }
            sequence_names[ row[0] ] = stoi( row[1] );
        }
        table.clear();

        // Prepare a sample.
        auto sample = Sample();

        // Read all chunks that have placements of the current sample.
        for( auto const& chunk_id : chunk_ids ) {
            auto chunk_sample = Sample();
            reader.from_file( epadir + "chunk_" + chunk_id + ".jplace" , chunk_sample );

            // If this is the first chunk for that sample, we use its tree. Otherwise, check
            // if the trees are compatible, just for safety.
            if( sample.empty() ) {
                sample = Sample( chunk_sample.tree() );
            } else if( ! compatible_trees( sample.tree(), chunk_sample.tree() )) {
                LOG_ERR << "incompatible trees in map " << map_filename << " and chunk " << chunk_id;
            }

            // Iterate all queries of the chunk and if needed, at them to the sample.
            for( auto const& pquery : chunk_sample ) {
                if( pquery.name_size() != 1 ) {
                    LOG_WARN << "name with size " << pquery.name_size() << " in map " << map_filename << " and chunk " << chunk_id;
                    continue;
                }

                // if the sample has a sequence of that name, add the pquery to it.
                auto entry = sequence_names.find( pquery.name_at(0).name );
                if( entry != sequence_names.end() ) {
                    auto& added_pquery = sample.add_pquery( pquery );
                    added_pquery.name_at(0).multiplicity = entry->second;
                }
            }
        }

        writer.to_file( sample, outdir + sample_name + ".jplace" );
    }
    */

    // Final output.
    // LOG_INFO << "at file " << file_count << " at " << utils::current_time();
    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
