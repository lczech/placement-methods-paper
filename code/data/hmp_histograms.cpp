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

#include "genesis/sequence/functions/codes.hpp"
#include "genesis/sequence/formats/fasta_input_iterator.hpp"
#include "genesis/sequence/formats/fasta_reader.hpp"
#include "genesis/sequence/formats/fasta_writer.hpp"
#include "genesis/sequence/functions/functions.hpp"
#include "genesis/sequence/sequence_set.hpp"
#include "genesis/sequence/sequence.hpp"

#include "genesis/utils/core/logging.hpp"
#include "genesis/utils/core/std.hpp"
#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/io/input_stream.hpp"
#include "genesis/utils/tools/date_time.hpp"
#include "genesis/utils/text/string.hpp"

#include <fstream>
#include <string>

using namespace genesis;
using namespace genesis::sequence;

int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    LOG_INFO << "Started " << utils::current_time();

    // Modify this to get the code to work:
    std::string basedir = "path/to/hmp/data/";

    std::string indir1 = basedir + "by_sample";
    std::string indir2 = basedir + "SRP002860";
    std::string outdir = basedir;

    auto reader = FastaReader();
    // reader.to_upper( false );

    auto ff1  = utils::dir_list_files(indir1);
    auto ff2 = utils::dir_list_files(indir2);
    std::vector<std::string> fasta_filenames;
    fasta_filenames.reserve( ff1.size() + ff2.size() );
    for( auto const& f : ff1 ) {
        fasta_filenames.push_back( indir1 + "/" + f );
    }
    for( auto const& f : ff2 ) {
        fasta_filenames.push_back( indir2 + "/" + f );
    }

    std::ofstream seqnames;
    seqnames.open( outdir + "/seqnames" );

    std::vector<size_t> hist_1;
    hist_1.resize(5000);
    auto hist_3 = hist_1;
    auto hist_6 = hist_1;
    auto hist_a = hist_1;
    size_t max = 0;
    size_t count = 0;

    size_t cnt_1 = 0;
    size_t cnt_3 = 0;
    size_t cnt_6 = 0;
    size_t cnt_a = 0;

    for( auto const& fasta_filename : fasta_filenames ) {
        if( ! utils::ends_with( fasta_filename, ".fsa" ) ) {
            continue;
        }

        if( count > 0 && count % 500 == 0 ) {
            LOG_DBG << "at " << count;
        }
        ++count;

        auto set = SequenceSet();
        reader.from_file( fasta_filename, set );

        for( auto const& seq : set ) {
            seqnames << seq.label() << "\n";

            bool found = false;
            if( seq.label().find( "V1-V3" ) != std::string::npos ) {
                ++hist_1[ seq.length() ];
                found = true;
                ++cnt_1;
            }
            if( seq.label().find( "V3-V5" ) != std::string::npos ) {
                ++hist_3[ seq.length() ];
                found = true;
                ++cnt_3;
            }
            if( seq.label().find( "V6-V9" ) != std::string::npos ) {
                ++hist_6[ seq.length() ];
                found = true;
                ++cnt_6;
            }

            if( ! found ) {
                LOG_DBG << "no region found for " << fasta_filename << " : " << seq.label();
            }

            if( seq.length() > max ) {
                max = seq.length();
            }

            ++hist_a[ seq.length() ];
            ++cnt_a;
        }
    }

    seqnames.close();

    LOG_INFO << "count V1-V3 " << cnt_1;
    LOG_INFO << "count V3-V5 " << cnt_3;
    LOG_INFO << "count V6-V9 " << cnt_6;
    LOG_INFO << "count all   " << cnt_a;

    std::ofstream hist_1_f;
    std::ofstream hist_3_f;
    std::ofstream hist_6_f;
    std::ofstream hist_a_f;
    hist_1_f.open( outdir + "/hist_1" );
    hist_3_f.open( outdir + "/hist_3" );
    hist_6_f.open( outdir + "/hist_6" );
    hist_a_f.open( outdir + "/hist_a" );

    for( size_t i = 0; i <= max; ++i ) {
        hist_1_f << i << "\t" << hist_1[i] << "\n";
        hist_3_f << i << "\t" << hist_3[i] << "\n";
        hist_6_f << i << "\t" << hist_6[i] << "\n";
        hist_a_f << i << "\t" << hist_a[i] << "\n";
    }

    hist_1_f.close();
    hist_3_f.close();
    hist_6_f.close();
    hist_a_f.close();

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
