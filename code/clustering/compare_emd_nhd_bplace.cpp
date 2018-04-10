/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2018 Lucas Czech and HITS gGmbH

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
using namespace genesis::tree;

void merge_order_rec( SquashClustering const& sc, size_t index, std::vector<size_t>& result )
{
    auto const& merger = sc.mergers()[ index ];

    // Get how many trees were clustered with squash clustering.
    auto const num_elems = sc.clusters().size() - sc.mergers().size();

    // Now, if the merge indices are lower, these are tips of the cluster tree,
    // i.e., they represent actual trees. If not, they are mergers() and we need to recurse.
    if( merger.index_a < num_elems ) {
        result.push_back( merger.index_a );
    } else {
        merge_order_rec( sc, merger.index_a - num_elems, result );
    }
    if( merger.index_b < num_elems ) {
        result.push_back( merger.index_b );
    } else {
        merge_order_rec( sc, merger.index_b - num_elems, result );
    }
}

std::vector<size_t> merge_order( SquashClustering const& sc )
{
    auto result = std::vector<size_t>();
    merge_order_rec( sc, sc.mergers().size() - 1, result );
    return result;
}

utils::Matrix<double> reorder_matrix( utils::Matrix<double> const& mat, std::vector<size_t> const& order )
{
    if( mat.rows() != mat.cols() ) {
        throw std::runtime_error( "mat not symmetrical" );
    }

    auto res = mat;
    for( size_t i = 0; i < mat.rows(); ++i ) {
        for( size_t j = 0; j < mat.cols(); ++j ) {
            res( i, j ) = mat( order[i], order[j] );
        }
    }
    return res;
}

// =================================================================================================
//      Main
// =================================================================================================

/**
 * Compare NHD and EMD to each other and to Squash Clustering,
 * by sorting the distance matrices using the squash cluster merge order.
 * This yields distance matrices sorted in a way that brings similar rows/colors close to each other,
 * resulting in a nice visualization of the distnace matrix a as a heat map.
 */
int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.date = true;
    utils::Logging::details.time = true;

    // utils::Options::get().number_of_threads( 4 );
    // LOG_BOLD << utils::Options::get().info();
    // LOG_BOLD;

    LOG_INFO << "Started";

    // Check if the command line contains the right number of arguments.
    if (argc != 4) {
        throw std::runtime_error(
            "Need to provide three arguments.\n"
        );
    }

    // In out dirs.
    auto const threads = std::stoi( argv[1] );
    auto epadir = utils::dir_normalize_path( std::string( argv[2] ));
    auto outdir = utils::dir_normalize_path( std::string( argv[3] ));
    utils::dir_create(outdir);

    utils::Options::get().command_line( argc, argv );
    utils::Options::get().number_of_threads( threads );
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    // -------------------------------------------------------------------------
    //     Find and read sample bplace files.
    // -------------------------------------------------------------------------

    // Get all bplace files, either in the epa dir, or in its subdirs
    auto bplace_filenames = utils::dir_list_files( epadir, ".*\\.bplace" );
    // std::sort( bplace_filenames.begin(), bplace_filenames.end() );

    // Output list of samples in the order that we use them for the matrix
    // LOG_INFO << "Writing bplace sample order list";
    // std::ofstream bplace_files_of( outdir + "bplace_order.list" );
    // for( auto const& elem : bplace_filenames ) {
    //     bplace_files_of << elem << "\n";
    // }
    // bplace_files_of.close();
    // LOG_INFO << "done";
    // LOG_INFO;

    SampleSet sset;
    SampleSerializer bplace_loader;

    // Process all bplace files.
    LOG_INFO << "reading " << bplace_filenames.size() << " bplace sample files";
    for( auto const& bplace_filename : bplace_filenames ) {
        sset.add( bplace_loader.load( epadir + bplace_filename ), bplace_filename );
    }

    // Final output for bplace reading
    assert( sset.size() == bplace_filenames.size() );
    LOG_INFO << "finished reading sample bplace files";
    LOG_INFO;

    // Output list of samples in the order that we use them for the matrix
    LOG_INFO << "Writing bplace sample order list";
    std::ofstream bplace_files_of( outdir + "bplace_order.list" );
    for( auto const& elem : bplace_filenames ) {
        bplace_files_of << elem << "\n";
    }
    bplace_files_of.close();
    LOG_INFO << "done";
    LOG_INFO;

    LOG_INFO << "Using " << utils::Options::get().number_of_threads() << " threads.";

    LOG_INFO << "adjust_to_average_branch_lengths";
    adjust_to_average_branch_lengths( sset );

    LOG_INFO << "Converting Trees";
    auto mass_trees = convert_sample_set_to_mass_trees( sset );
    // sset.clear();

    // -------------------------------------------------------------------------
    //     Do the calculations.
    // -------------------------------------------------------------------------

    LOG_INFO << "NHD Matrix calculation started";
    auto const nhd_matrix = node_histogram_distance( sset );
    utils::file_write( utils::to_string( nhd_matrix ), outdir + "nhd_unordered.mat" );
    LOG_INFO << "Done";


    LOG_INFO << "EMD Matrix calculation started";
    auto const emd_matrix = earth_movers_distance( sset );
    utils::file_write( utils::to_string( emd_matrix ), outdir + "emd_unordered.mat" );
    LOG_INFO << "Done";


    LOG_INFO << "Starting squash clustering";
    auto sc = tree::SquashClustering();
    sc.run(  std::move( mass_trees.first ) );
    LOG_INFO << "Finished squash clustering";

    LOG_INFO << "Writing cluster results";
    std::ofstream file_clust_results;
    utils::file_output_stream( outdir + "/cluster.info",  file_clust_results );

    file_clust_results << "Clusters\n";
    for( size_t i = 0; i < sc.clusters().size(); ++i ) {
        file_clust_results << i << ": " << ( sc.clusters()[i].active ? "act " : "dea " )
                 << sc.clusters()[i].count << " "
                 << sc.clusters()[i].distances.size() << " "
                 << sc.clusters()[i].tree.node_count()
                 << "\n";
    }
    file_clust_results << "\nMergers\n";
    for( size_t i = 0; i < sc.mergers().size(); ++i ) {
        file_clust_results << i << ": "
                 << sc.mergers()[i].index_a << " " << sc.mergers()[i].distance_a << " | "
                 << sc.mergers()[i].index_b << " " << sc.mergers()[i].distance_b
                 << "\n";
    }

    // Remove bplace part from file names
    for( auto& fn : bplace_filenames ) {
        fn = utils::replace_all( fn, ".bplace", "" );
    }

    LOG_INFO << "Writing cluster tree";
    std::ofstream file_clust_tree;
    utils::file_output_stream( outdir + "/cluster.newick",  file_clust_tree );
    file_clust_tree << sc.tree_string( bplace_filenames );

    LOG_INFO << "Reordering matrices";
    auto const order = merge_order( sc );
    utils::file_write( utils::to_string( reorder_matrix( nhd_matrix, order ) ), outdir + "nhd_sc-order.mat" );
    utils::file_write( utils::to_string( reorder_matrix( emd_matrix, order ) ), outdir + "emd_sc-order.mat" );

    file_clust_results << "\nReorder\n";
    for( size_t i = 0; i < order.size(); ++i ) {
        file_clust_results << order[i] << " ";
    }
    file_clust_results << "\n";

    LOG_INFO << "Finished";
    return 0;
}
