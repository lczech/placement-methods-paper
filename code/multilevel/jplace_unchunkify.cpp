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
using namespace genesis::placement;

// =================================================================================================
//     unchunkify_with_named_chunks
// =================================================================================================

void unchunkify_with_named_chunks(
    std::string const& mapdir,
    std::string const& epadir,
    std::string const& outdir
) {
        // Prepare reader and writer
        JplaceReader reader;
        JplaceWriter writer;

        // Reader for map files.
        utils::CsvReader csv_reader;
        csv_reader.separator_chars("\t");

        // Find map files.
        auto map_filenames = utils::dir_list_files( mapdir );
        std::sort( map_filenames.begin(), map_filenames.end() );
        LOG_INFO << "processing " << map_filenames.size() << " files";

        // Process all map files.
        size_t file_count = 0;
        for( auto const& map_filename : map_filenames ) {

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
                auto chunk_sample = reader.from_file( epadir + "chunk_" + chunk_id + ".jplace" );

                // If this is the first chunk for that sample, we use its tree. Otherwise, check
                // if the trees are compatible, just for safety.
                if( sample.empty() ) {
                    sample = Sample( chunk_sample.tree() );
                } else if( ! compatible_trees( sample.tree(), chunk_sample.tree() )) {
                    LOG_ERR << "incompatible trees in map " << map_filename << " and chunk " << chunk_id;
                }

                // Iterate all queries of the chunk and if needed, add them to the sample.
                for( auto const& pquery : chunk_sample ) {
                    if( pquery.name_size() != 1 ) {
                        LOG_WARN << "name with size " << pquery.name_size() << " in map " << map_filename << " and chunk " << chunk_id;
                        continue;
                    }

                    // if the sample has a sequence of that name, add the pquery to it.
                    auto entry = sequence_names.find( pquery.name_at(0).name );
                    if( entry != sequence_names.end() ) {
                        auto& added_pquery = sample.add( pquery );
                        added_pquery.name_at(0).multiplicity = entry->second;
                    }
                }
            }

            writer.to_file( sample, outdir + sample_name + ".jplace" );
        }

        // Final output.
        LOG_INFO << "at file " << file_count << " at " << utils::current_time();
}

// =================================================================================================
//     unchunkify_with_wild_chunks
// =================================================================================================

struct MappedSample
{
    // Name of the sample, derived from map file name.
    std::string name;

    struct SeqInfo
    {
        size_t count;
        std::string chunk_name;
    };

    // Which sequence/placement has how many counts/multiplicity.
    std::unordered_map<std::string, SeqInfo> sequences;

    // Sample where placements are collected.
    Sample sample;
};

void unchunkify_with_wild_chunks(
    std::string const& mapdir,
    std::string const& epadir,
    std::string const& outdir
) {
    // Prepare reader and writer
    JplaceReader reader;
    JplaceWriter writer;

    // Reader for map files.
    utils::CsvReader csv_reader;
    csv_reader.separator_chars("\t");

    // Find map files.
    auto map_filenames = utils::dir_list_files( mapdir );
    std::sort( map_filenames.begin(), map_filenames.end() );
    LOG_INFO << "reading " << map_filenames.size() << " map files";

    // Prepare map and sample storage.
    std::vector<MappedSample> sample_maps;

    // Process all map files.
    size_t file_count = 0;
    for( auto const& map_filename : map_filenames ) {

        // Progress output.
        if( file_count > 0 && file_count % 50 == 0 ) {
            LOG_INFO << "at map file " << file_count << " (" << map_filename << ") at " << utils::current_time();
        }
        ++file_count;

        // Get the sample name from the map file name.
        auto map_filename_parts = utils::split( map_filename, "." );
        if( map_filename_parts.size() != 3 ) {
            LOG_WARN << "Weird map filename " << map_filename;
        }

        // Prepare sample map.
        MappedSample sample_map;
        sample_map.name = map_filename_parts[0];

        // Read mapfile, get list of sequence names with multiplicities.
        auto table = csv_reader.from_file( mapdir + map_filename );
        for( auto const& row : table ) {
            if( row.size() != 3 ) {
                LOG_WARN << "invalid row in map file " << map_filename;
                continue;
            }

            if( sample_map.sequences.count( row[0] ) > 0 ) {
                LOG_WARN << "dup seq name in map file " << map_filename << ": " << row[0];
            }
            auto freq = stoi( row[1] );
            sample_map.sequences[ row[0] ].count      = freq;
            sample_map.sequences[ row[0] ].chunk_name = row[2];

            if( freq <= 0 ) {
                LOG_WARN << "weird count " << freq << " in " << map_filename;
            }
        }

        sample_maps.push_back( std::move( sample_map ));
    }

    // Final output for map reading
    LOG_INFO << "at map file " << file_count << " at " << utils::current_time();
    LOG_INFO << "finished reading map files";
    LOG_INFO;

    // Get all jplace files, either in the epa dir, or in its subdirs
    std::vector<std::string> jplace_filenames;
    // jplace_filenames = utils::dir_list_files( epadir, ".*\\.jplace" );
    for( auto const& epa_subdir : utils::dir_list_directories( epadir ) ) {
        for( auto const& jplace_filename : utils::dir_list_files( epadir + epa_subdir, ".*\\.jplace" ) ) {
            jplace_filenames.push_back( utils::dir_normalize_path( epa_subdir ) + jplace_filename );
        }
    }
    std::sort( jplace_filenames.begin(), jplace_filenames.end() );
    LOG_INFO << "reading " << jplace_filenames.size() << " jplace chunk files";

    // Process all jplace files.
    file_count = 0;
    for( auto const& jplace_filename : jplace_filenames ) {

        // Progress output.
        if( file_count > 0 && file_count % 10 == 0 ) {
            LOG_INFO << "at jplace file " << file_count << " (" << jplace_filename << ") at " << utils::current_time();
        }
        ++file_count;

        // Read chunk.
        auto chunk_sample = reader.from_file( epadir + jplace_filename );
        std::string chunk_name_hypo;

        for( auto const& pquery : chunk_sample ) {
            // safety
            if( pquery.name_size() != 1 || pquery.name_at(0).multiplicity != 1.0 ) {
                LOG_WARN << "weird pquery in " << jplace_filename;
                continue;
            }

            auto const& name = pquery.name_at(0).name;

            // add the pquery to all samples where it occured.
            for( auto& sample_map : sample_maps ) {
                auto key = sample_map.sequences.find( name );
                if( key == sample_map.sequences.end() ) {
                    continue;
                }

                // Reverse-engineed the chunk name.
                if( chunk_name_hypo.empty() ) {
                    chunk_name_hypo = key->second.chunk_name;
                } else if( chunk_name_hypo != key->second.chunk_name ) {
                    LOG_WARN << "weird: " << jplace_filename << " has a pquery (" << name << ") "
                             << "with a chunk name conflict (" << chunk_name_hypo << " vs " << key->second.chunk_name << ")";
                }

                // If this is the first chunk for that sample, we use its tree. Otherwise, check
                // if the trees are compatible, just for safety.
                if( sample_map.sample.empty() ) {
                    sample_map.sample = Sample( chunk_sample.tree() );
                } else if( ! compatible_trees( sample_map.sample.tree(), chunk_sample.tree() )) {
                    LOG_ERR << "incompatible trees in sample " << sample_map.name << " and jplace chunk " << jplace_filename;
                    continue;
                }

                // Add the pquery to the sample, adjust multiplicity
                auto& added_pquery = sample_map.sample.add( pquery );
                added_pquery.name_at(0).multiplicity = key->second.count;

                // The sequence should only occur in one chunk jplace file. in order to verify this,
                // we set the count for the sequence to 0 as a marker.
                // If it then already is zero, the sequences occured before, which is wrong!
                if( key->second.count == 0 ) {
                    LOG_WARN << "in " << jplace_filename << " sequence " << name
                             << " occured, but was already zero in " << sample_map.name;
                }
                key->second.count = 0;
            }
        }

        utils::file_append(
            jplace_filename + "\t" + chunk_name_hypo + "\n",
            outdir + "chunknames.csv"
        );
    }

    // Final output for jplace reading
    LOG_INFO << "at jplace file " << file_count << " at " << utils::current_time();
    LOG_INFO << "finished reading chunk jplace files";
    LOG_INFO;

    // Process all samples.
    LOG_INFO << "writing sample jplace files";
    file_count = 0;
    std::unordered_map<std::string, size_t> missing_seqs;
    for( auto const& sample_map : sample_maps ) {

        // Progress output.
        if( file_count > 0 && file_count % 10 == 0 ) {
            LOG_INFO << "at sample file " << file_count << " (" << sample_map.name << ") at " << utils::current_time();
        }
        ++file_count;

        // write sample jplace file.
        writer.to_file( sample_map.sample, outdir + "samples/" + sample_map.name + ".jplace" );

        // write sample summary (how many pqueries it has)
        utils::file_append(
            sample_map.name + "\t" + std::to_string( sample_map.sample.size() ) + "\n",
            outdir + "overview.csv"
        );

        // process missing pquries for which there was nothing in the chunks.
        for( auto const& seqs : sample_map.sequences ) {
            if( seqs.second.count > 0 ) {
                continue;
            }
            utils::file_append(
                seqs.first + "\t" + std::to_string( seqs.second.count ) + "\t" + seqs.second.chunk_name + "\n",
                outdir + "missing/" + sample_map.name + ".csv"
            );
            missing_seqs[ seqs.first ] += seqs.second.count;
        }
    }

    // Final output for jplace reading
    LOG_INFO << "at sample file " << file_count << " at " << utils::current_time();
    LOG_INFO << "finished writing sample jplace files";
    LOG_INFO;

    // write missing seqs summary.
    for( auto const& mis : missing_seqs ) {
        utils::file_append(
            mis.first + "\t" + std::to_string( mis.second ) + "\n",
            outdir + "missing.csv"
        );
    }
}

// =================================================================================================
//     unchunkify_with_wild_chunks_fasterish
// =================================================================================================

void unchunkify_with_wild_chunks_fasterish(
    std::string const& mapdir,
    std::string const& epadir,
    std::string const& outdir
) {
    // Prepare reader and writer
    JplaceReader reader;
    JplaceWriter writer;
    SampleSerializer set_writer;

    // -------------------------------------------------------------------------
    //     Find and read chunk jplace files.
    // -------------------------------------------------------------------------

    // Get all jplace files, either in the epa dir, or in its subdirs
    std::vector<std::string> jplace_filenames;
    // jplace_filenames = utils::dir_list_files( epadir, ".*\\.jplace" );
    for( auto const& epa_subdir : utils::dir_list_directories( epadir ) ) {
        for( auto const& jplace_filename : utils::dir_list_files( epadir + epa_subdir, ".*\\.jplace" ) ) {
            jplace_filenames.push_back( utils::dir_normalize_path( epa_subdir ) + jplace_filename );
        }
    }
    std::sort( jplace_filenames.begin(), jplace_filenames.end() );
    LOG_INFO << "reading " << jplace_filenames.size() << " jplace chunk files";

    // Keep all chunk samples in memory. Hopefully you have enough ram!
    std::vector< Sample > chunk_samples;

    // Process all jplace files.
    size_t file_count = 0;
    for( auto const& jplace_filename : jplace_filenames ) {

        // Progress output.
        // if( file_count > 0 && file_count % 10 == 0 ) {
            LOG_INFO << "at jplace file " << file_count << " ("
                     << jplace_filename << ")";
        // }
        ++file_count;

        chunk_samples.push_back(
            reader.from_file( epadir + jplace_filename )
        );
    }

    // Final output for jplace reading
    assert( chunk_samples.size() == jplace_filenames.size() );
    LOG_INFO << "finished reading chunk jplace files";
    LOG_INFO;

    // Build THE tree used for all samples.
    LOG_INFO << "Building avg tree: 1";
    tree::TreeSet avg_tree_set;
    for( auto const& chunk_sample : chunk_samples ) {
        avg_tree_set.add( "", chunk_sample.tree() );
    }
    LOG_INFO << "Building avg tree: 2";
    auto avg_tree = average_branch_length_tree( avg_tree_set );
    avg_tree_set.clear();
    LOG_INFO << "Building avg tree: done";

    // -------------------------------------------------------------------------
    //     Construct lookup maps.
    // -------------------------------------------------------------------------

    // lookup to later find which chunk contains a certain sequence, plus a direct pointer to the pquery.
    // map from seq name to index in chunk_samples plus pointer to the pquery.
    std::unordered_map< std::string, std::pair< size_t, Pquery const* > > seqence_name_to_pquery;

    LOG_INFO << "Constructing sequence name to chunk lookup...";
    for( size_t i = 0; i < chunk_samples.size(); ++i ) {
        auto const& chunk_sample = chunk_samples[i];

        for( auto const& pquery : chunk_sample ) {

            // safety
            if( pquery.name_size() != 1 || pquery.name_at(0).multiplicity != 1.0 ) {
                LOG_WARN << "weird pquery in " << jplace_filenames[i];
                continue;
            }

            // Add pquery to lookup
            auto const& name = pquery.name_at(0).name;
            auto check = seqence_name_to_pquery.emplace( name, std::make_pair( i, &pquery ) );

            if( !check.second ) {
                LOG_WARN << "sequence " << name << " already in lookup at chunk"
                         << check.first->second.first;
            }
        }
    }
    LOG_INFO << "...done";
    LOG_INFO;

    // Lookup for chunk id's as produced by the "chunkify" app, using the same indices as
    // jplace_filenames.
    auto chunk_id_map = std::vector<std::string>( jplace_filenames.size(), "" );

    // -------------------------------------------------------------------------
    //     Process samples.
    // -------------------------------------------------------------------------

    // Reader for map files.
    utils::CsvReader csv_reader;
    csv_reader.separator_chars("\t");

    // Find map files.
    auto map_filenames = utils::dir_list_files( mapdir );
    std::sort( map_filenames.begin(), map_filenames.end() );
    LOG_INFO << "reading " << map_filenames.size() << " map files";

    std::unordered_map<std::string, size_t> missing_seqs;

    // Process all map files.
    file_count = 0;
    for( auto const& map_filename : map_filenames ) {

        // Progress output.
        // if( file_count > 0 && file_count % 50 == 0 ) {
            LOG_INFO << "at map file " << file_count << " (" << map_filename << ")";
        // }
        ++file_count;

        // Get the sample name from the map file name.
        auto map_filename_parts = utils::split( map_filename, "." );
        if( map_filename_parts.size() != 3 ) {
            LOG_WARN << "Weird map filename " << map_filename;
        }

        // Create a sample using the avg tree.
        auto map_sample = Sample( avg_tree );

        std::ofstream missing_seqs_of( outdir + "missing/" + map_filename_parts[0] + ".csv" );

        // Read mapfile, get list of sequence names with multiplicities.
        auto table = csv_reader.from_file( mapdir + map_filename );
        for( auto const& row : table ) {
            if( row.size() != 3 ) {
                LOG_WARN << "invalid row in map file " << map_filename;
                continue;
            }

            auto seq_name = row[0];
            auto freq = stoi( row[1] );
            auto chunk_id = row[2];
            if( freq <= 0 ) {
                LOG_WARN << "weird count " << freq << " in " << map_filename;
            }

            auto chunk_id_it = seqence_name_to_pquery.find( seq_name );
            if( chunk_id_it == seqence_name_to_pquery.end() ) {
                missing_seqs_of << seq_name << "\t"
                                << std::to_string( freq ) << "\t" << chunk_id << "\n";
                missing_seqs[ seq_name ] += freq;
                continue;
            }

            // Add the pquery to the sample
            auto& added_pquery = map_sample.add( *chunk_id_it->second.second );
            added_pquery.name_at(0).multiplicity = freq;

            // Back-map the chunk names.
            auto chunk_i = chunk_id_it->second.first;
            if( chunk_id_map[chunk_i] != "" && chunk_id_map[chunk_i] != chunk_id ) {
                LOG_WARN << "inconsistent chunk id. sample " << map_filename_parts[0]
                         << ", seq " << seq_name << ", chunk name " << jplace_filenames[chunk_i]
                         << ", chunk_id" << chunk_id;
            }
            chunk_id_map[chunk_i] = chunk_id;
        }

        writer.to_file(  map_sample, outdir + "samples/" + map_filename_parts[0] + ".jplace" );
        set_writer.save( map_sample, outdir + "samples/" + map_filename_parts[0] + ".bplace" );

        // write sample summary (how many pqueries it has)
        utils::file_append(
            map_filename_parts[0] + "\t" +
            std::to_string( map_sample.size() ) + "\t" +
            std::to_string( total_multiplicity( map_sample )) + "\n",
            outdir + "overview.csv"
        );
    }

    // Final output for map reading
    LOG_INFO << "finished reading map files";
    LOG_INFO;

    // -------------------------------------------------------------------------
    //     Final bookkeeping stuff.
    // -------------------------------------------------------------------------

    LOG_INFO << "Writing backmapping for chunk names";
    std::ofstream backmap_of( outdir + "chunknames.csv" );
    for( size_t i = 0; i < chunk_id_map.size(); ++i ) {
        backmap_of << jplace_filenames[i] << "\t" << chunk_id_map[i] << "\n";
    }
    backmap_of.close();
    LOG_INFO << "done";
    LOG_INFO;

    // write missing seqs summary.
    LOG_INFO << "Writing missing seqs overview";
    std::ofstream missing_of( outdir + "missing.csv" );
    for( auto const& mis : missing_seqs ) {
        missing_of << mis.first << "\t" << mis.second << "\n";
    }
    missing_of.close();
    LOG_INFO << "done";
    LOG_INFO;
}

// =================================================================================================
//     main
// =================================================================================================

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;
    LOG_INFO << "Started";

    // Check if the command line contains the right number of arguments.
    if (argc != 4) {
        throw std::runtime_error(
            "Need to provide three arguments.\n"
        );
    }

    // In out dirs.
    auto mapdir = utils::dir_normalize_path( std::string( argv[1] ));
    auto epadir = utils::dir_normalize_path( std::string( argv[2] ));
    auto outdir = utils::dir_normalize_path( std::string( argv[3] ));
    utils::dir_create(outdir);
    utils::dir_create(outdir + "samples");
    utils::dir_create(outdir + "missing");

    // unchunkify_with_named_chunks( mapdir, epadir, outdir );
    unchunkify_with_wild_chunks_fasterish( mapdir, epadir, outdir );

    LOG_INFO << "Finished";
    return 0;
}
