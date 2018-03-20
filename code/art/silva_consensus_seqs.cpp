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
 * @brief Calculate consensus sequences with different methods,
 * using the output data from the `silva_entropy.cpp` program as input.
 *
 * This is mainly a cheap trick in order to not have to re-run the entropy calculations,
 * but instead use the sequence names as identifiers for which sequences to make new consensus sequences
 * with different other consensus methods.
 */
int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;

    if( argc != 5 ) {
        throw std::runtime_error(
            "Need to provide a taxon name, the silve data dir, the taxon_name dir, and an out dir."
        );
    }

    // auto const silva_dir = "data/silva/";
    // auto const ref_dir = "01_backbone/00_reference/";
    // auto const out_dir = "11_consensus_seqs/00_reference/";

    auto const taxon_name = std::string( argv[1] );
    auto const silva_dir = utils::dir_normalize_path( std::string( argv[2] ));
    auto const ref_dir   = utils::dir_normalize_path( std::string( argv[3] ));
    auto const out_dir   = utils::dir_normalize_path( std::string( argv[4] ));


    LOG_INFO << "=================================================================================";
    LOG_INFO << "Processing " << taxon_name;


    // Read ref seq labels
    LOG_INFO << "Reading ref seqs";
    std::string ref_file = ref_dir + taxon_name + "/tax_cons_border.fasta";
    auto ref_seqs = FastaReader().from_file( ref_file );

    // get map of seq name to counts
    LOG_INFO << "Processing ref";
    std::unordered_map<std::string, SiteCounts> ref_map;
    for( auto const& seq : ref_seqs ) {
        if( ref_map.count(seq.label()) > 0 ) {
            throw std::runtime_error("ref seq duplicate: " + seq.label());
        }
        ref_map[ seq.label() ] = SiteCounts( "ACGT", 50000 );
    }
    LOG_INFO << "Found " << ref_map.size() << " ref seqs";

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

        std::string taxopath_str;
        auto const delim = seq.label().find_first_of( " \t" );
        if( delim == std::string::npos ) {
            taxopath_str = seq.label();
        } else {
            taxopath_str = seq.label().substr( delim + 1 );
        }

        if( taxopath_str == "" ) {
            throw std::runtime_error( "Empty metadata." );
        }

        auto name = sanitize_label( taxopath_str );

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

        ref_map[ name ].add_sequence( seq );
    }
    LOG_INFO << "Done reading";
    LOG_INFO << "no_taxon_count " << no_taxon_count;

    auto writer = FastaWriter();

    // Majorities
    utils::dir_create(out_dir + "majorities");
    {
        LOG_INFO << "Making consensus_sequence_with_majorities";
        SequenceSet cons_seqs;
        for( auto const& ref : ref_map ) {
            cons_seqs.add({ ref.first, consensus_sequence_with_majorities( ref.second ) });
        }

        LOG_INFO << "Writing consensus_sequence_with_majorities";
        writer.to_file( cons_seqs, out_dir + "majorities/" + taxon_name + "_sequences.fasta" );
    }

    // Cavener
    utils::dir_create(out_dir + "cavener");
    {
        LOG_INFO << "Making consensus_sequence_cavener";
        SequenceSet cons_seqs;
        for( auto const& ref : ref_map ) {
            cons_seqs.add({ ref.first, consensus_sequence_cavener( ref.second ) });
        }

        LOG_INFO << "Writing consensus_sequence_cavener";
        writer.to_file( cons_seqs, out_dir + "cavener/" + taxon_name + "_sequences.fasta" );
    }

    // Threshold
    std::vector<double> thresholds = { 0.5, 0.6, 0.7, 0.8, 0.9, 0.95 };
    for( auto threshold : thresholds ) {
        auto const od = out_dir + "threshold_" + utils::to_string(threshold);
        utils::dir_create(od);

        LOG_INFO << "Making consensus_sequence_with_threshold " << threshold;
        SequenceSet cons_seqs;
        for( auto const& ref : ref_map ) {
            cons_seqs.add({ ref.first, consensus_sequence_with_threshold( ref.second, threshold ) });
        }

        LOG_INFO << "Writing consensus_sequence_with_threshold " << threshold;
        writer.to_file( cons_seqs, od + "/" + taxon_name + "_sequences.fasta" );
    }

    LOG_INFO << "Finished.";
    return 0;
}
