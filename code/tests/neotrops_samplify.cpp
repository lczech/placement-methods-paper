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
#include <vector>

using namespace genesis;
using namespace genesis::sequence;

// =================================================================================================
//     Main
// =================================================================================================

/**
 * @brief Take the neotrop amplicons fasta file and abundance table, and make single fasta
 * sample files, each with abundance information in the sequence name.
 */
int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.date = true;
    utils::Logging::details.time = true;
    LOG_INFO << "Started " << utils::current_time();

    // In out dirs.
    std::string fasta_in   = "data/nts/sequences/neotropical_soil_175_samples_protist_amplicons.fas";
    std::string tabfile = "data/nts/sequences/neotropical_soil_175_samples_protist_amplicons.table";
    std::string outdir  = "data/nts/samples/";

    // Prepare reading and writing fasta.
    auto reader = FastaReader();
    reader.to_upper( true );
    auto writer = FastaWriter();
    writer.enable_metadata(false);
    auto csv_reader = utils::CsvReader();
    csv_reader.separator_chars("\t");

    LOG_INFO << "Reading fasta";
    SequenceSet all_seqs = reader.from_file( fasta_in );
    LOG_INFO << "done. found " << all_seqs.size() << " sequences.";

    LOG_INFO << "building lookup";
    std::unordered_map< std::string, size_t > lookup;
    for( size_t i = 0; i < all_seqs.size(); ++i ) {
        auto splitted = utils::split( all_seqs[i].label(), "_" );
        if( splitted.size() != 2 ) {
            LOG_ERR << "wrong split size";
            return 1;
        }
        auto seq_name = splitted[0];
        if( lookup.count( seq_name) > 0 ) {
            LOG_ERR << "duplicate sequence: " << seq_name;
            return 1;
        }
        lookup[ seq_name ] = i;
    }
    LOG_INFO << "done";

    auto sample_seqs = std::vector<SequenceSet>( 154 );

    LOG_INFO << "Reading table";
    size_t count = 0;
    utils::InputStream tab_in { utils::make_unique<utils::FileInputSource>(tabfile) };
    while( tab_in ) {
        if( count % 100000 == 0 ) {
            LOG_DBG1 << "reading at " << count;
        }

        // Get csv line.
        auto line = csv_reader.parse_line( tab_in );
        if( line.size() != 156 ) {
            LOG_ERR << "invalid line length " << line.size();
            return 1;
        }
        auto name = line[0];

        // get sequence from all seqs.
        if( lookup.count(name) == 0 ) {
            LOG_ERR << "unknown seq " << name;
            return 1;
        }
        Sequence cur_seq = all_seqs[ lookup[name] ];
        auto splitted = utils::split( cur_seq.label(), "_" );
        if( splitted.size() != 2 ) {
            LOG_ERR << "wrong split size";
            return 1;
        }
        auto seq_name = splitted[0];
        size_t seq_abun = stoi( splitted[1] );

        if( name != seq_name ) {
            LOG_ERR << "wrong name " << name << " at count " << count;
            return 1;
        }

        size_t sum = 0;
        for( size_t i = 0; i < 154; ++i ) {
            size_t cell_count = stoi( line[ i+1 ] );
            sum += cell_count;

            if( cell_count > 0 ) {
                cur_seq.label( seq_name + "_" + std::to_string( cell_count ) );
                sample_seqs[ i ].add( cur_seq );
            }
        }
        if( sum != static_cast<size_t>(stoi( line[155] )) ) {
            LOG_ERR << "wrong line sum at " << count << " with " << sum << " instead of " << line[155];
            return 1;
        }
        if( seq_abun != sum ) {
            LOG_ERR << "wrong abun at " << count << " with " << seq_abun << " instead of " << sum;
            return 1;
        }

        all_seqs[ lookup[name] ].clear();
        ++count;
    }
    LOG_INFO << "done. found " << count << " lines.";

    LOG_INFO << "writing sample fasta files";
    for( size_t i = 0; i < 154; ++i ) {
        writer.to_file( sample_seqs[i], outdir + "sample_" + std::to_string(i) + ".fasta" );
    }
    LOG_INFO << "done";

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
