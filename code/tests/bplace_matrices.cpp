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

using namespace genesis;
using namespace genesis::placement;

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;

    utils::Options::get().number_of_threads( 12 );
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    LOG_INFO << "Started";

    // Check if the command line contains the right number of arguments.
    if (argc != 3) {
        throw std::runtime_error(
            "Need to provide three arguments.\n"
        );
    }

    // In out dirs.
    auto epadir = utils::dir_normalize_path( std::string( argv[1] ));
    auto outdir = utils::dir_normalize_path( std::string( argv[2] ));
    utils::dir_create(outdir);

    // -------------------------------------------------------------------------
    //     Find and read sample bplace files.
    // -------------------------------------------------------------------------

    // Get all bplace files, either in the epa dir, or in its subdirs
    auto bplace_filenames = utils::dir_list_files( epadir, ".*\\.bplace" );
    std::sort( bplace_filenames.begin(), bplace_filenames.end() );

    // Output list of samples in the order that we use them for the matrix
    LOG_INFO << "Writing bplace sample order list";
    std::ofstream bplace_files_of( outdir + "bplace_order_mats.list" );
    for( auto const& elem : bplace_filenames ) {
        bplace_files_of << elem << "\n";
    }
    bplace_files_of.close();
    LOG_INFO << "done";
    LOG_INFO;

    SampleSet sample_set;
    SampleSerializer bplace_loader;

    // Process all jplace files.
    LOG_INFO << "reading " << bplace_filenames.size() << " bplace sample files";
    size_t file_count = 0;
    for( auto const& bplace_filename : bplace_filenames ) {

        // Progress output.
        // LOG_INFO << "at jplace file " << file_count << " (" << bplace_filename << ")";
        ++file_count;

        sample_set.add( bplace_loader.load( epadir + bplace_filename ), bplace_filename );
    }

    // Final output for jplace reading
    assert( sample_set.size() == bplace_filenames.size() );
    LOG_INFO << "finished reading sample bplace files";
    LOG_INFO;

    // Write all kinds of matrices

    // Average tree.
    tree::TreeSet avg_tree_set;
    for( auto const& smp : sample_set ) {
        avg_tree_set.add( "", smp.sample.tree() );
    }
    auto const avg_tree = tree::average_branch_length_tree( avg_tree_set );
    tree::DefaultTreeNewickWriter().to_file( avg_tree, outdir + "avg_tree.newick" );
    avg_tree_set.clear();

    // Calculate the pairwise distance between all pairs of nodes.
    {
        auto const node_dists = node_branch_length_distance_matrix( avg_tree );
        utils::file_write( utils::to_string( node_dists ), outdir + "node_dists.mat" );
    }

    // Imbalance matrix.
    {
        auto const imbalance_matrix = epca_imbalance_matrix( sample_set );
        utils::file_write( utils::to_string( imbalance_matrix ), outdir + "imbalance.mat" );
    }

    // Edge weight matrix.
    {
        auto const ew_matrix = placement_weight_per_edge( sample_set );
        utils::file_write( utils::to_string( ew_matrix ), outdir + "edge_weights.mat" );
    }

    LOG_INFO << "Finished";
    return 0;
}
