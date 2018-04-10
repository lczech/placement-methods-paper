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

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

void run_nhd( std::vector<std::string> const& bplace_filenames, std::string const& outdir )
{
    size_t const bins = 50;
    auto const set_size = bplace_filenames.size();

    auto hist_vecs = std::vector<NodeDistanceHistogramSet>( bplace_filenames.size() );
    SampleSerializer bplace_loader;

    Matrix<double> node_distances;
    Matrix<signed char> node_sides;
    // size_t node_count;

    LOG_INFO << "init";
    {
        // Load the first one twice... can live with that for now
        auto smp = bplace_loader.load( bplace_filenames[0] );
        node_distances = node_branch_length_distance_matrix( smp.tree() );
        node_sides = node_root_direction_matrix( smp.tree() );
        // node_count = smp.tree().node_count();
    }

    auto nhd_matrix = utils::Matrix<double>( set_size, set_size, 0.0 );

    LOG_INFO << "fill histogams";
    #pragma omp parallel for
    for( size_t fi = 0; fi < bplace_filenames.size(); ++fi ) {
        // if( fi > 0 ) {
        //     if( ! compatible_trees( sample_set[ fi - 1 ].sample, sample_set[ fi ].sample )) {
        //         throw std::invalid_argument(
        //             "Trees in SampleSet not compatible for calculating Node Histogram Distance."
        //         );
        //     }
        // }

        // TODO check compatibility!!!

        #pragma omp critical(PRINT)
        {
            LOG_DBG1 << "load " << fi;
        }
        auto smp = bplace_loader.load( bplace_filenames[fi] );

        #pragma omp critical(PRINT)
        {
            LOG_DBG1 << "make " << fi;
        }
        hist_vecs[fi] = node_distance_histogram_set( smp, node_distances, node_sides, bins );

        // #pragma omp critical(PRINT)
        // {
        //     LOG_DBG1 << "fill " << fi;
        // }
        // fill_node_distance_histograms( smp, node_distances, node_sides, hist_vecs[fi] );

        #pragma omp critical(PRINT)
        {
            LOG_DBG1 << "done " << fi;
        }
    }

    LOG_INFO << "calc dists";

    // We only need to calculate the upper triangle. Get the number of indices needed
    // to describe this triangle.
    size_t const max_k = utils::triangular_size( set_size );

    // Calculate distance matrix for every pair of samples.
    #pragma omp parallel for
    for( size_t k = 0; k < max_k; ++k ) {

        // For the given linear index, get the actual position in the Matrix.
        auto const ij = utils::triangular_indices( k, set_size );
        auto const i = ij.first;
        auto const j = ij.second;

        // Calculate and store distance.
        auto const dist = node_histogram_distance( hist_vecs[ i ], hist_vecs[ j ] );
        nhd_matrix(i, j) = dist;
        nhd_matrix(j, i) = dist;
    }
    LOG_INFO << "finished";

    utils::file_write( utils::to_string( nhd_matrix ), outdir + "nhd_unordered.mat" );
    LOG_INFO << "written";

}

void run_emd( std::vector<std::string> const& bplace_filenames, std::string const& outdir )
{
    SampleSerializer bplace_loader;

    // Process all bplace files.
    LOG_INFO << "reading " << bplace_filenames.size() << " bplace sample files";
    SampleSet sset = bplace_loader.load( bplace_filenames );

    // Final output for bplace reading
    assert( sset.size() == bplace_filenames.size() );
    LOG_INFO << "finished reading sample bplace files";
    LOG_INFO;

    LOG_INFO << "adjust_to_average_branch_lengths";
    adjust_to_average_branch_lengths( sset );

    LOG_INFO << "EMD Matrix calculation started";
    auto const emd_matrix = earth_movers_distance( sset );
    utils::file_write( utils::to_string( emd_matrix ), outdir + "emd_unordered.mat" );
    LOG_INFO << "finished";

    utils::file_write( utils::to_string( emd_matrix ), outdir + "emd_sc-order.mat" );
    LOG_INFO << "written";
}

// =================================================================================================
//      Main
// =================================================================================================

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
    if (argc != 5) {
        throw std::runtime_error(
            "Need to provide four arguments.\n"
        );
    }

    // In out dirs.
    auto const mode = std::string( argv[1] );
    auto const threads = std::stoi( argv[2] );
    auto epadir = utils::dir_normalize_path( std::string( argv[3] ));
    auto outdir = utils::dir_normalize_path( std::string( argv[4] ));
    utils::dir_create(outdir);

    utils::Options::get().command_line( argc, argv );
    utils::Options::get().number_of_threads( threads );
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    LOG_INFO << "Using " << utils::Options::get().number_of_threads() << " threads.";

    // -------------------------------------------------------------------------
    //     Find and read sample bplace files.
    // -------------------------------------------------------------------------

    // Get all bplace files, either in the epa dir, or in its subdirs
    auto bplace_filenames = utils::dir_list_files( epadir, true, ".*\\.bplace" );
    std::sort( bplace_filenames.begin(), bplace_filenames.end() );

    // -------------------------------------------------------------------------
    //     Do the calculations.
    // -------------------------------------------------------------------------

    if( mode == "nhd" ) {
        run_nhd( bplace_filenames, outdir );
    } else if( mode == "emd" ) {
        run_emd( bplace_filenames, outdir );
    } else {
        LOG_ERR << "nope";
    }

    LOG_INFO << "Finished";
    return 0;
}
