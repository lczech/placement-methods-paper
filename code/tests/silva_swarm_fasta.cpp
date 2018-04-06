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

#include "genesis/sequence.hpp"
#include "genesis/taxonomy.hpp"
#include "genesis/utils.hpp"

#include <fstream>
#include <string>

using namespace genesis;
using namespace genesis::sequence;
using namespace genesis::taxonomy;

/**
 * @brief Prepare the Silva 600k sequences for running swarm on them.
 */
int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;

    std::string basedir = "data/silva/";

    // Read sequences
    LOG_INFO << "Start reading";
    std::string infile = basedir + "SILVA_123.1_SSURef_Nr99_tax_silva.fasta";
    auto reader = FastaReader();
    reader.to_upper( true );
    auto seqs = reader.from_file( infile );
    LOG_INFO << "Read " << seqs.size() << " sequences";

    // Prepare stuff
    std::string const acgt = "ACGT";
    auto writer = FastaWriter();
    writer.enable_metadata(true);

    // Process
    LOG_INFO << "Start processing";
    auto taxopath_parser = TaxopathParser();
    for( size_t i = 0; i < seqs.size(); ++i ) {
        auto& seq = seqs[i];
        replace_u_with_t(seq);

        // Check taxonomy
        auto taxopath = taxopath_parser.from_string( seq.metadata() );
        if( taxopath.size() == 0 ) {
            LOG_DBG << "empty metadata at " << seq.label();
            continue;
        }

        // Remove everything that is not ACGT
        auto is_search_char = [&] ( char c ) {
            return acgt.find( c ) == std::string::npos;
        };
        auto& str = seq.sites();
        str.erase( std::remove_if( str.begin(), str.end(), is_search_char ), str.end() );
    }
    LOG_INFO << "Done processing";

    // Write raw seqs
    LOG_INFO << "Writing raw seqs";
    writer.to_file( seqs, basedir + "swarm/raw.fasta" );
    LOG_INFO << "Done writing";

    // Relabel seqs with their sha1, and make a mapping file
    LOG_INFO << "SHA";
    std::ofstream sha_map_file;
    utils::file_output_stream( basedir + "swarm/shamap.txt",  sha_map_file );
    for( auto& seq : seqs ) {
        auto const digest = utils::SHA1::from_string_hex( seq.sites() );
        sha_map_file << digest << "\t" << seq.label() << "\t" << seq.metadata() << "\n";
        seq.label( digest );
    }
    sha_map_file.close();

    // Check if there are duplicate sequence labels.
    LOG_INFO << "Dups?";
    std::unordered_set<std::string> hashes;
    bool has_dups = false;
    for( auto const& seq : seqs ) {
        if( hashes.count(seq.label()) > 0 ) {
            has_dups = true;
            break;
        } else {
            hashes.insert(seq.label());
        }
    }
    hashes.clear();

    // If so, remove and append their count.
    // If not, simply add a count of 1 to each sequence.
    if( has_dups ) {
        LOG_INFO << "HAS DUPS!";
        merge_duplicate_sequences( seqs, MergeDuplicateSequencesCountPolicy::kAppendToLabel );
    } else {
        LOG_INFO << "No dups! Phew.";
        for( auto& seq : seqs ) {
            seq.label( seq.label() + "_1" );
        }
    }

    // Write final fasta file
    LOG_INFO << "Writing final seqs";
    writer.enable_metadata(false);
    writer.to_file( seqs, basedir + "swarm/sequences.fasta" );
    LOG_INFO << "Wrote " << seqs.size() << " sequences";

    LOG_INFO << "Finished.";
    return 0;
}
