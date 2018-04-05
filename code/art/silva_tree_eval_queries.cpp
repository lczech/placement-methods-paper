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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace genesis;
using namespace genesis::sequence;
using namespace genesis::taxonomy;

// =================================================================================================
//      Main
// =================================================================================================

/**
 * @brief Write out the 600k silva sequences in to chunks, and rename the sequences so that
 * their names reflect the taxo path they should be placed on.
 *
 * The renaming is done so that the evaulation is easy: simply look at the name to figure out where
 * it should go in the placement. We also prepend a unique counter prefix to each sequence, so that
 * epa and other tools don't complain about duplicate sequence names (e.g., from the same genus).
 * This fixed length id is easy to remove in the eval program.
 */
int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;
    LOG_INFO << "Started " << utils::current_time();

    // Load the taxonomy.
    LOG_DBG << "Reading taxonomy...";
    std::string tax_file = "/path/to/data/silva/tax_slv_ssu_123.1.txt";
    auto tax = Taxonomy();
    TaxonomyReader().from_file( tax_file, tax );
    sort_by_name( tax );
    LOG_DBG << "done";

    // Prepare sequence input.
    std::string aln_file = "/path/to/data/silva/SILVA_123.1_SSURef_Nr99_tax_silva_full_align_trunc.fasta";
    utils::InputStream aln_is { utils::make_unique<utils::FileInputSource>( aln_file ) };
    auto reader = FastaReader();
    auto taxopath_parser = TaxopathParser();
    Sequence seq;

    // Output
    auto writer = FastaWriter();
    // writer.enable_metadata(false);
    std::string outdir = "/path/to/data/silva/600k_taxa/";

    // Collect sequences for a chunk here.
    SequenceSet chunk;
    size_t chunk_count = 0;

    // Add to count objects along their taxonomic path.
    LOG_DBG << "Start reading at " << utils::current_time();
    size_t seq_cnt = 0;
    while( reader.parse_sequence( aln_is, seq )) {
        replace_u_with_t( seq );

        if( seq_cnt % 100000 == 0 ) {
            LOG_DBG << "At sequence " << seq_cnt;
        }
        ++seq_cnt;

        auto const meta = seq.label().substr( seq.label().find( " " ));
        if( meta == "" ) {
            throw std::runtime_error( "Empty metadata." );
        }

        auto taxopath = taxopath_parser( meta );
        // taxopath.pop_back();
        // if( find_taxon_by_taxopath( tax, taxopath ) == nullptr ) {
        //     LOG_DBG << "Sequence taxon not in taxonomy: " << seq.metadata();
        //     // throw std::runtime_error( "Sequence taxon not in taxonomy: " + seq.metadata() );
        // }

        // Set label to counter plus taxopath.
        auto gen = TaxopathGenerator();
        auto name = sanitize_label( gen(taxopath) );
        seq.label( "SEQ_" + utils::to_string_leading_zeros( seq_cnt, 6 ) + "_" + name );

        // Add to chunk
        chunk.add( seq );

        // If a chunk is full, flush it.
        if( chunk.size() == 100000 ) {
            writer.to_file( chunk, outdir + "chunk_" + utils::to_string(chunk_count) + ".fasta" );
            ++chunk_count;
            chunk.clear();
        }
    }
    LOG_DBG << "Done reading at " << utils::current_time();

    // Flush the remaining chunk.
    writer.to_file( chunk, outdir + "chunk_" + utils::to_string(chunk_count) + ".fasta" );

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
