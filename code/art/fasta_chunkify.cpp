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

// namespace std
// {
//     /**
//      * @brief Hash function for SHA1 digestes.
//      *
//      * Basically, we re-hash from 160 bit to 64 bit. This is ugly, but currently faster to implement
//      * than a custom container that uses the full hash width. Might work on this in the future.
//      */
//     template<>
//     struct hash<genesis::utils::SHA1::DigestType>
//     {
//         using argument_type = genesis::utils::SHA1::DigestType;
//         using result_type   = std::size_t;
//
//         result_type operator()( argument_type const& s ) const {
//             result_type hash = 0;
//             hash ^= s[0] ^ ( static_cast<result_type>( s[1] ) << 32 );
//             hash ^= s[2] ^ ( static_cast<result_type>( s[3] ) << 32 );
//             hash ^= s[4];
//             return hash;
//         }
//     };
// }

// =================================================================================================
//     Dataset Specific Stuff
// =================================================================================================

/**
 * @brief Filter HMP files. Return whether the file should be kept.
 */
bool filter_hmp_sequences( SequenceSet& set, std::string const& fasta_filename )
{
    // // HMP: Remove unnecessary label parts, add filename to each label (need for chunking).
    // for( auto& seq : set ) {
    //     // s = Clear Span, nbp = No Barcodes or Primers, rc = Reverse Complemented
    //     auto label = seq.label();
    //     auto start_position_to_erase = label.find("_cs_nbp_rc");
    //     if( start_position_to_erase == std::string::npos ) {
    //         LOG_WARN << "unclean sequence: " << fasta_filename << " : " << label;
    //     }
    //     label.erase(start_position_to_erase);
    //     seq.label( utils::file_filename(fasta_filename) + "_" + label );
    // }

    // HMP: Remove too short and too long sequences.
    size_t index = 0;
    while( index < set.size() ) {
        auto const& seq = set[index];

        // Get primer.
        int primer = 0;
        if( seq.label().find( "primer=V1-V3" ) != std::string::npos ) {
            primer = 1;
        }
        if( seq.label().find( "primer=V3-V5" ) != std::string::npos ) {
            primer = 3;
        }
        if( seq.label().find( "primer=V6-V9" ) != std::string::npos ) {
            primer = 6;
        }
        if( primer == 0 ) {
            LOG_WARN << "invalid primer: " << fasta_filename << " : " << seq.label();
        }

        // Filter by length, depending on primer.
        // We set the max to the position where the last declining number appears.
        auto len = set.at(index).length();
        bool remove = false;
        if( len < 150 ) {
            remove = true;
        }
        if( primer == 1 && len > 540 ) {
            remove = true;
        }
        if( primer == 3 && len > 575 ) {
            remove = true;
        }
        if( primer == 6 && len > 540 ) {
            remove = true;
        }
        if( remove ) {
            set.remove(index);
        } else {
            ++index;
        }
    }

    // HMP: We don't want samples with too few sequences. Skip all output.
    if( set.size() < 1500 ) {
        return false;
    }

    // // Some last checks.
    // if( ! has_unique_labels(set, false) ) {
    //     LOG_WARN << "non unique labels in " << fasta_filename;
    //     return false;
    // }
    // if( ! has_valid_labels( set ) ) {
    //     LOG_WARN << "invalid labels in " << fasta_filename;
    //     sanitize_labels(set);
    // }

    return true;
}

bool filter_tara_sequences( SequenceSet& set, std::string const& fasta_filename )
{
    (void) fasta_filename;

    // Tara: Remove too short and too long sequences.
    size_t index = 0;
    while( index < set.size() ) {
        auto const& seq = set[index];

        if( seq.length() < 95 || seq.length() > 150 ) {
            set.remove(index);
        } else {
            ++index;
        }
    }

    // always keep the sample, never return false
    return true;
}

size_t process_tara_sequence( Sequence const& seq, std::string const& sha1_hex )
{
    auto parts = utils::split( seq.label(), ";", true);

    if( parts.size() != 2 ) {
        LOG_WARN << "bad sequence label. has size " << parts.size();
        return 0;
    }
    if( sha1_hex != parts[0] ) {
        LOG_WARN << "wrong sha 1 hex " << sha1_hex << " vs " << parts[0];
    }

    auto sizes = utils::split( parts[1], "=", false );
    if( sizes.size() != 2 || sizes[0] != "size" ) {
        LOG_WARN << "bad label size attrib";
        return 0;
    }

    return std::stoi( sizes[1] );
}

size_t process_nts_sequence( Sequence const& seq, std::string const& sha1_hex )
{
    auto parts = utils::split( seq.label(), "_" );

    if( parts.size() != 2 ) {
        LOG_WARN << "bad sequence label. has size " << parts.size();
        return 0;
    }
    if( sha1_hex != parts[0] ) {
        LOG_WARN << "wrong sha 1 hex " << sha1_hex << " vs " << parts[0];
    }

    return std::stoi( parts[1] );
}

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
    auto indir  = utils::trim_right( std::string( argv[1] ), "/") + "/";
    auto outdir = utils::trim_right( std::string( argv[2] ), "/") + "/";
    utils::dir_create(outdir);
    utils::dir_create(outdir + "maps");
    utils::dir_create(outdir + "chunks");
    utils::dir_create(outdir + "filtered");

    // Prepare reading and writing fasta.
    auto reader = FastaReader();
    // reader.to_upper( true );
    auto writer = FastaWriter();
    // writer.enable_metadata(false);

    // Collect sequences for a chunk here.
    SequenceSet chunk;
    size_t chunk_count = 0;

    // Sequences hashes, mapping to the chunk number where we store them.
    std::unordered_map< std::string, size_t > hashes;

    // Histogram
    std::vector<size_t> hist;
    hist.resize(5000);
    size_t hist_max = 0;

    // Count how many files and sequences were processed.
    size_t file_count = 0;
    size_t seqs_count = 0;

    // Find fsa/fasta files.
    auto fasta_filenames = utils::dir_list_files( indir, ".*\\.(fsa|fasta)" );
    std::sort( fasta_filenames.begin(), fasta_filenames.end() );
    LOG_INFO << "processing " << fasta_filenames.size() << " files";

    // Process all fasta files.
    for( auto const& fasta_filename : fasta_filenames ) {

        // Progress output.
        if( file_count > 0 && file_count % 20 == 0 ) {
            LOG_DBG << "at file " << file_count;
            LOG_INFO << "read " << seqs_count << " seqs, " << hashes.size() << " uniq";
        }
        ++file_count;

        // Check file type (again, safety first).
        if( ! utils::ends_with( fasta_filename, ".fsa" ) && ! utils::ends_with( fasta_filename, ".fasta" ) ) {
            LOG_WARN << "skip invalid filename " << fasta_filename;
            continue;
        }

        // Read.
        auto set = SequenceSet();
        reader.from_file( indir + fasta_filename, set );

        // -------------------------------------------------------------------------
        //     Dataset Specific Stuff >>>

        // // HMP
        // if( ! filter_hmp_sequences( set, fasta_filename )) {
        //     continue;
        // }

        // Tara
        // if( ! filter_tara_sequences( set, fasta_filename )) {
        //     continue;
        // }

        //     <<<
        // -------------------------------------------------------------------------

        // Count identical sequences, accesses via their hash.
        std::unordered_map< std::string, size_t > seq_freqs;

        // Precess all sequences: count how often they occur in the file, and if new, add to chunks.
        for( auto const& seq : set ) {
            auto hash_hex = utils::SHA1::from_string_hex( seq.sites() );

            // -------------------------------------------------------------------------
            //     Dataset Specific Stuff >>>

            // Increment seq counters.
            // ++seq_freqs[ hash_hex ];
            // seq_freqs[ hash_hex ] += process_tara_sequence( seq, hash_hex );
            seq_freqs[ hash_hex ] += process_nts_sequence( seq, hash_hex );

            //     <<<
            // -------------------------------------------------------------------------

            ++seqs_count;

            // Histogram.
            ++hist[ seq.length() ];
            if( seq.length() > hist_max ) {
                hist_max = seq.length();
            }

            // We saw that sequence before. Don't need to add it to the chunk.
            if( hashes.count(hash_hex) > 0 ) {
                continue;
            }

            // New sequence: never saw that hash before. Add it to the chunk, store chunk num.
            hashes[hash_hex] = chunk_count;
            chunk.add( Sequence( hash_hex, seq.sites() ) );

            // If a chunk is full, flush it.
            if( chunk.size() == 50000 ) {
                writer.to_file( chunk, outdir + "chunks/chunk_" + utils::to_string(chunk_count) + ".fasta" );
                ++chunk_count;
                chunk.clear();
            }
        }

        // Write filtered files with sha labels.
        relabel_sha1(set);
        writer.to_file( set, outdir + "filtered/" + fasta_filename );

        // Write a file that tells us which sequence (by its hash) occured how often and in which chunk.
        std::ofstream seq_freq_file( outdir + "maps/" + fasta_filename + ".map" );
        for( auto const& elem : seq_freqs ) {
            seq_freq_file << elem.first << "\t" << elem.second << "\t" << hashes.at(elem.first) << "\n";
        }
        seq_freq_file.close();
    }

    // Flush the remaining chunk.
    writer.to_file( chunk, outdir + "chunks/chunk_" + utils::to_string(chunk_count) + ".fasta" );

    // Final output.
    LOG_INFO << "at file " << file_count;
    LOG_INFO << "read " << seqs_count << " seqs, " << hashes.size() << " uniq";

    // Write histogram.
    std::ofstream hist_a_f;
    hist_a_f.open( outdir + "length_histogram" );
    for( size_t i = 0; i <= hist_max; ++i ) {
        hist_a_f << i << "\t" << hist[i] << "\n";
    }
    hist_a_f.close();

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
