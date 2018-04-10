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

/**
 * Some tests while developing the NH distance.
 */
int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.date = true;
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
    std::ofstream bplace_files_of( outdir + "bplace_order_nhd.list" );
    for( auto const& elem : bplace_filenames ) {
        bplace_files_of << elem << "\n";
    }
    bplace_files_of.close();
    LOG_INFO << "done";
    LOG_INFO;

    SampleSet sset;
    SampleSerializer bplace_loader;

    // Process all jplace files.
    LOG_INFO << "reading " << bplace_filenames.size() << " bplace sample files";
    // size_t file_count = 0;
    for( auto const& bplace_filename : bplace_filenames ) {

        // Progress output.
        // LOG_INFO << "at jplace file " << file_count << " (" << bplace_filename << ")";
        // ++file_count;

        sset.add( bplace_loader.load( epadir + bplace_filename ), bplace_filename );
    }

    // Final output for jplace reading
    assert( sset.size() == bplace_filenames.size() );
    LOG_INFO << "finished reading sample bplace files";
    LOG_INFO;

    LOG_INFO << "Using " << utils::Options::get().number_of_threads() << " threads.";

    LOG_INFO << "Old NHD Matrix calculation started";
    auto nhd_matrix = node_histogram_distance( sset );
    LOG_INFO << "Done";
    LOG_INFO << "Writing...";
    utils::file_write( utils::to_string( nhd_matrix ), outdir + "nhd_old.mat" );
    LOG_INFO << "Done";

    LOG_INFO << "New NHD Matrix calculation started, positive only";
    nhd_matrix = node_histogram_distance_2( sset, 25, false );
    LOG_INFO << "Done";
    LOG_INFO << "Writing...";
    utils::file_write( utils::to_string( nhd_matrix ), outdir + "nhd_pos.mat" );
    LOG_INFO << "Done";

    LOG_INFO << "New NHD Matrix calculation started, negative axis";
    nhd_matrix = node_histogram_distance_2( sset, 25, true );
    LOG_INFO << "Done";
    LOG_INFO << "Writing...";
    utils::file_write( utils::to_string( nhd_matrix ), outdir + "nhd_neg.mat" );
    LOG_INFO << "Done";

    LOG_INFO << "New NHD single calculation started, positive only";
    nhd_matrix = utils::Matrix<double>( 10, 10 );
    for( size_t i = 0; i < 10; ++i ) {
        for( size_t j = 0; j < 10; ++j ) {
            if( i < sset.size() && j < sset.size() ) {
                nhd_matrix( i, j ) = node_histogram_distance_2( sset[i].sample, sset[j].sample, 25, false );
            }
        }
    }
    LOG_INFO << "Done";
    LOG_INFO << "Writing...";
    utils::file_write( utils::to_string( nhd_matrix ), outdir + "nhd_single_pos.mat" );
    LOG_INFO << "Done";

    LOG_INFO << "New NHD single calculation started, negative axis";
    nhd_matrix = utils::Matrix<double>( 10, 10 );
    for( size_t i = 0; i < 10; ++i ) {
        for( size_t j = 0; j < 10; ++j ) {
            if( i < sset.size() && j < sset.size() ) {
                nhd_matrix( i, j ) = node_histogram_distance_2( sset[i].sample, sset[j].sample, 25, true );
            }
        }
    }
    LOG_INFO << "Done";
    LOG_INFO << "Writing...";
    utils::file_write( utils::to_string( nhd_matrix ), outdir + "nhd_single_neg.mat" );
    LOG_INFO << "Done";

    LOG_INFO << "Finished";
    return 0;
}
