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

using namespace genesis;

/**
 * @brief Take a set of jplace files and merge them into one.
 *
 * The command line arguments are:
 *
 * ./merge_jplace_files [<tree_file>] <input_jplace_file>... <output_jplace_file>
 *
 * That is, the program takes an optional reference tree file, then a list of input jplace files,
 * and finally an output jplace file.
 *
 * If the first file does not end in "jplace", it is considered to be a tree in Newick format.
 * It is used as reference tree for merging the other jplace files, which is particularly relevant
 * for the branch lengths.
 *
 * If no tree file is given as first arugment (that is, all file names end in "jplace"),
 * instead the tree with average branch lengths of all input jplace files is used as the
 * reference tree for merging the samples.
 *
 * Then, all input jplace files are merged onto the reference tree. It is important to note that
 * in this process, the position of the placements on the branches are adjusted so that they
 * are relatively the same as before. So, a placement that is located at a third of its branch
 * in the input file is still located at a third in the merged output, even if the branch
 * has a different length. For most downstream analysis, this slight shift does not matter too much.
 *
 * The merged output is then written to the output jplace file, which must not exist before.
 */
int main( int argc, const char* argv[] )
{
    using namespace ::genesis::placement;

    // Activate logging, print genesis header.
    utils::Logging::log_to_stdout();
    LOG_BOLD << genesis_header();

    LOG_INFO << "Started 'merge jplace files' program.";

    // Check if the command line contains the right number of arguments.
    if (argc < 4) {
        LOG_ERR << "Invalid number of command line arguments.";
        return 1;
    }

    // Prepare a Jplace reader that reports wrong values.
    auto jplace_reader = JplaceReader();
    jplace_reader.invalid_number_behaviour(
        genesis::placement::JplaceReader::InvalidNumberBehaviour::kLogAndCorrect
    );

    // Prepare an empty sample to collect all pqueries.
    Sample out_sample;

    // If the first arg is a jplace file, we read all of them into a set,
    // and then merge them onto the average branch length tree.
    if( utils::ends_with( utils::to_lower( argv[1] ), "jplace" )) {

        // Create list of file names.
        std::vector<std::string> files;
        std::string file_list;
        for( size_t i = 1; i < static_cast<size_t>( argc - 1 ); ++i ) {
            files.push_back( argv[i] );
            file_list += "\n    " + std::string( argv[i] );
        }
        LOG_INFO << "Reading " << (argc - 2) << " input Jplace files:" << file_list;

        // Read the files into a set, and merge it to get the final sample.
        auto sample_set = jplace_reader.from_files( files );
        out_sample = merge_all( sample_set );
        LOG_INFO << "Found " << out_sample.size() << " Pqueries in total.";

    // If the first arg is not a jplace file, we take it as a Newick tree,
    // which is used are reference tree.
    // All other files are then merged onto that tree.
    } else {
        LOG_INFO << "First argument is not a Jplace file. "
                 << "Assuming that it is a reference tree in Newick format.";
        LOG_INFO << "Reading input Newick file:\n    " << argv[1];

        // Read the tree and create an empty sample out of it.
        auto newick_reader = tree::DefaultTreeNewickReader();
        newick_reader.enable_tags( true );
        out_sample = Sample(
            convert_default_tree_to_placement_tree( newick_reader.from_file( argv[1] ))
        );

        std::string file_list;
        for( size_t i = 2; i < static_cast<size_t>( argc - 1 ); ++i ) {
            file_list += "\n    " + std::string( argv[i] );
        }
        LOG_INFO << "Reading " << (argc - 3) << " input Jplace files:" << file_list;

        // Read all jplace files and merge them into the sample.
        for( size_t i = 2; i < static_cast<size_t>( argc - 1 ); ++i ) {
            copy_pqueries( jplace_reader.from_file( argv[i] ), out_sample );
        }
        LOG_INFO << "Found " << out_sample.size() << " Pqueries in total.";
    }

    // Finally, write out the merged sample.
    LOG_INFO << "Writing output jplace file:\n    " << argv[argc - 1];
    JplaceWriter().to_file( out_sample,  argv[argc - 1] );
    LOG_INFO << "Finished.";
    return 0;
}
