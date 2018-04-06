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

size_t match_score( std::string const& consensus, std::string const& candidate )
{
    if( consensus.size() != candidate.size() ) {
        throw std::runtime_error( "Comparing apples and oranges." );
    }

    size_t matches = 0;
    for( size_t i = 0; i < consensus.size(); ++i ) {
        auto const amb = nucleic_acid_ambiguities( consensus[i] );
        if( amb.find( candidate[i] ) != std::string::npos ) {
            ++matches;
        }
    }
    return matches;
}

/**
 * @brief Find the Silva 600k sequences that are closest to the enotrpy-chosen consensus seqs.
 * This was a test: instead of using consensus sequences, use actual species sequences
 * as representatives of a clde. in order to have them still represent the calde well,
 * select the ones that are closest to a consenes seq of the clade...
 * kind of a weird approchac, so we did not contiue with it.
 */
int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;

    if (argc != 2) {
        throw std::runtime_error(
            "Need to provide a reference name."
        );
    }
    auto const reference = std::string( argv[1] );
    LOG_INFO << "=================================================================================";
    LOG_INFO << "Processing " << reference;

    std::string silva_dir = "data/silva/";
    std::string ref_dir = "00_reference/";
    std::string out_dir = "11_single_seqs/00_reference/";

    // Read ref seq labels
    LOG_INFO << "Reading ref seqs";
    std::string ref_file = ref_dir + reference + "/tax_cons_border.fasta";
    auto ref_seqs = FastaReader().from_file( ref_file );

    // get map of seq name to seq for fast lookup
    LOG_INFO << "Processing ref";
    std::unordered_map<std::string, std::string> ref_map;
    for( auto const& seq : ref_seqs ) {
        if( ref_map.count(seq.label()) > 0 ) {
            throw std::runtime_error("ref seq duplicate: " + seq.label());
        }
        ref_map[ seq.label() ] = seq.sites();
    }
    LOG_INFO << "Found " << ref_map.size() << " ref seqs";

    // Prepare result collection.
    struct BestSequence
    {
        std::string sequence;
        size_t      score;
    };
    std::unordered_map< std::string, BestSequence> matches;

    // Prepare sequence input.
    std::string aln_file = silva_dir + "SILVA_123.1_SSURef_Nr99_tax_silva_full_align_trunc.fasta";
    utils::InputStream aln_is { utils::make_unique<utils::FileInputSource>( aln_file ) };
    auto reader = FastaReader();
    auto taxopath_parser = TaxopathParser();
    auto taxopath_gen = TaxopathGenerator();
    Sequence seq;

    LOG_INFO << "Start reading 600k sequences";
    size_t seq_cnt = 0;
    size_t no_taxon_count = 0;
    while( reader.parse_sequence( aln_is, seq )) {
        replace_u_with_t( seq );
        for( auto& c : seq ) {
            if( c == '.' ) {
                c = '-';
            }
        }

        if( seq_cnt % 50000 == 0 ) {
            LOG_INFO << "At sequence " << seq_cnt;
        }
        ++seq_cnt;

        if( seq.metadata() == "" ) {
            throw std::runtime_error( "Empty metadata." );
        }

        // auto taxopath = taxopath_parser( seq.metadata() );
        // auto name = sanitize_label( taxopath_gen(taxopath) );
        // if( taxopath_gen(taxopath) != seq.metadata() ) {
        //     LOG_DBG << taxopath_gen(taxopath) << " != " << seq.metadata();
        // }
        auto name = sanitize_label( seq.metadata() );

        while( name.size() > 0 && ref_map.count( name ) == 0 ) {
            size_t lastindex = name.find_last_of("_");
            name = name.substr(0, lastindex);

            if( lastindex == std::string::npos ) {
                break;
            }
        }
        if( name.size() == 0 || ref_map.count( name ) == 0 ) {
            ++no_taxon_count;
            continue;
        }

        auto const score = match_score( ref_map[ name ], seq.sites() );
        if( score > matches[name].score ) {
            matches[name].score = score;
            matches[name].sequence = seq.sites();
        }
    }
    LOG_INFO << "Done reading";
    LOG_INFO << "no_taxon_count " << no_taxon_count;

    LOG_INFO << "writing full file";
    std::ofstream full_file;
    SequenceSet all_set;
    utils::file_output_stream( out_dir + reference + "_detail.txt",  full_file );
    for( auto const& match : matches ) {
        full_file << match.first << ": " << match.second.score << "\n";
        full_file << match.second.sequence << "\n";
        full_file << ref_map[ match.first ] << "\n\n";

        all_set.add({ match.first, match.second.sequence });
    }
    full_file.close();

    // Write final fasta file
    LOG_INFO << "Writing final seqs";
    auto writer = FastaWriter();
    writer.enable_metadata(false);
    writer.to_file( all_set, out_dir + reference + "_sequences.fasta" );
    LOG_INFO << "Wrote " << all_set.size() << " sequences";

    LOG_INFO << "Finished.";
    return 0;
}
