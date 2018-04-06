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

#include "sequence/counts.hpp"
#include "sequence/formats/fasta_input_iterator.hpp"
#include "sequence/formats/fasta_reader.hpp"
#include "sequence/formats/fasta_writer.hpp"
#include "sequence/functions/consensus.hpp"
#include "sequence/functions/entropy.hpp"
#include "sequence/functions/functions.hpp"
#include "sequence/sequence_set.hpp"
#include "sequence/sequence.hpp"

#include "taxonomy/formats/taxonomy_reader.hpp"
#include "taxonomy/formats/taxopath_generator.hpp"
#include "taxonomy/formats/taxopath_parser.hpp"
#include "taxonomy/functions/entropy.hpp"
#include "taxonomy/functions/operators.hpp"
#include "taxonomy/functions/taxonomy.hpp"
#include "taxonomy/functions/taxopath.hpp"
#include "taxonomy/iterator/postorder.hpp"
#include "taxonomy/iterator/preorder.hpp"
#include "taxonomy/taxon.hpp"
#include "taxonomy/taxonomy.hpp"
#include "taxonomy/taxopath.hpp"

#include "utils/core/fs.hpp"
#include "utils/core/logging.hpp"
#include "utils/core/std.hpp"
#include "utils/formats/csv/reader.hpp"
#include "utils/io/input_stream.hpp"
#include "utils/math/histogram.hpp"
#include "utils/math/histogram/accumulator.hpp"
#include "utils/text/string.hpp"
#include "utils/tools/date_time.hpp"

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

int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    LOG_INFO << "Started " << utils::current_time();

    // -------------------------------------------------------------------------
    //     Read all sequences.
    // -------------------------------------------------------------------------

    std::string indir = "tax/";

    std::vector<std::string> files = {
        "sub_alignments/Archaea.fasta",
        "sub_alignments/Bacteria_Actinobacteria.fasta",
        "sub_alignments/Bacteria_Bacteroidetes.fasta",
        "sub_alignments/Bacteria_Firmicutes.fasta",
        "sub_alignments/Bacteria_Proteobacteria.fasta",
        "sub_alignments/Eukaryota_Archaeplastida_Chloroplastida_Charophyta.fasta",
        "sub_alignments/Eukaryota_Archaeplastida_Chloroplastida_Chlorophyta.fasta",
        "sub_alignments/Eukaryota_Archaeplastida_Rhodophyceae.fasta",
        "sub_alignments/Eukaryota_Opisthokonta_Nucletmycea.fasta",
        "sub_alignments/Eukaryota_SAR_Alveolata.fasta",
        "sub_alignments/Eukaryota_SAR_Rhizaria_Cercozoa.fasta",
        "sub_alignments/Eukaryota_SAR_Stramenopiles.fasta",
        "tax_cons_border.fasta"
    };

    // sequence input.
    auto reader = FastaReader();
    LOG_DBG << "reading seqs";
    std::vector<SequenceSet> setset;
    for( auto const& filename : files ) {
        LOG_DBG1 << filename;
        setset.push_back({});
        reader.from_file( filename, setset.back() );
    }

    auto gapsites = utils::Bitvector( setset[0][0].length() );
    for( auto const& set : setset ) {
        auto setgapsites = gap_sites( set );

        if( setgapsites.size() != gapsites.size() ) {
            LOG_WARN << "diff size: " << setgapsites.size() << " and " << gapsites.size();
            return 0;
        }
    }

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
