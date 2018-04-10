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

#include <chrono>
#include <fstream>
#include <string>
#include <unordered_map>

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::utils;

/**
 * Testing how binning affects the accuracy and speed of EMD calculations.
 */
int main( int argc, char** argv )
{

    // -------------------------------------------------------------------------
    //     Input shit.
    // -------------------------------------------------------------------------

    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;

    utils::Options::get().number_of_threads( 16 );
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
    //     Read shit.
    // -------------------------------------------------------------------------

    // Get all bplace files, either in the epa dir, or in its subdirs
    auto bplace_filenames = utils::dir_list_files( epadir, ".*\\.bplace" );
    std::sort( bplace_filenames.begin(), bplace_filenames.end() );

    // Output list of samples in the order that we use them for the matrix
    std::ofstream bplace_files_of( outdir + "bplace_order.list" );
    for( auto const& elem : bplace_filenames ) {
        bplace_files_of << elem << "\n";
    }
    bplace_files_of.close();

    SampleSet sample_set;
    SampleSerializer bplace_loader;

    // Process all jplace files.
    LOG_INFO << "reading " << bplace_filenames.size() << " bplace sample files";
    for( auto const& bplace_filename : bplace_filenames ) {
        sample_set.add( bplace_loader.load( epadir + bplace_filename ), bplace_filename );
    }

    // Final output for jplace reading
    assert( sample_set.size() == bplace_filenames.size() );
    LOG_INFO << "reading done";
    LOG_INFO;

    // -------------------------------------------------------------------------
    //     Helper shit.
    // -------------------------------------------------------------------------

    auto matrix_stats = []( utils::Matrix<double> const& mat ){
        auto const minmax = matrix_minmax(mat);
        LOG_DBG1 << "min    " << minmax.min;
        LOG_DBG1 << "max    " << minmax.max;
        LOG_DBG1 << "avg    " << ( matrix_sum( mat ) / static_cast<double>( mat.size() ) );

        auto const meanstddev = matrix_mean_stddev( mat );
        LOG_DBG1 << "mean   " << meanstddev.mean;
        LOG_DBG1 << "stddev " << meanstddev.stddev;

        auto const quarts = matrix_quartiles(mat);
        LOG_DBG1 << "q0     " << quarts.q0;
        LOG_DBG1 << "q1     " << quarts.q1;
        LOG_DBG1 << "q2     " << quarts.q2;
        LOG_DBG1 << "q3     " << quarts.q3;
        LOG_DBG1 << "q4     " << quarts.q4;
    };

    // -------------------------------------------------------------------------
    //     Calcualte shit.
    // -------------------------------------------------------------------------

    LOG_INFO << "Using " << utils::Options::get().number_of_threads() << " threads.";

    LOG_INFO << "Convert to mass trees";
    auto const mass_trees = convert_sample_set_to_mass_trees( sample_set );
    LOG_INFO << "Convert done.";

    LOG_INFO << "Full Matrix calculation";
    auto const c_start = std::chrono::steady_clock::now();
    auto const emd_mat = tree::earth_movers_distance( mass_trees.first );
    auto const c_duration = std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::steady_clock::now() - c_start
    );
    LOG_INFO << "Full Matrix calculation done";
    LOG_INFO << "Time: " << c_duration.count() << "s";
    utils::file_write( utils::to_string( emd_mat ), outdir + "emd_full.mat" );
    LOG_INFO;

    matrix_stats( emd_mat );
    LOG_INFO;

    std::vector<size_t> bins { 1, 2, 4, 8, 16, 32, 64, 128, 256 };
    for( auto bin : bins ) {
        LOG_INFO << "===============================================================";
        LOG_INFO << "Binning to " << bin;

        // Make copy and binify it.
        auto binned_trees = mass_trees;
        for( auto& tree : binned_trees.first ) {
            mass_tree_binify_masses( tree, bin );
        }

        LOG_INFO << "Matrix calculation";
        auto const c_start = std::chrono::steady_clock::now();
        auto const emd_mat_binned = tree::earth_movers_distance( binned_trees.first );
        auto const c_duration = std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::steady_clock::now() - c_start
        );
        LOG_INFO << "Matrix calculation done";
        LOG_INFO << "Time: " << c_duration.count() << "s";
        utils::file_write( utils::to_string( emd_mat_binned ), outdir + "emd_" + std::to_string(bin) + ".mat" );
        LOG_INFO;

        matrix_stats(emd_mat_binned);
        LOG_INFO;

        LOG_INFO << "Difference matrix:";
        matrix_stats( matrix_subtraction( emd_mat, emd_mat_binned ));
        LOG_INFO;
    }

    LOG_INFO << "Finished";
    return 0;
}
