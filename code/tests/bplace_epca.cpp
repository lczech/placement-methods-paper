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
    //     Find and read sample bplace files.
    // -------------------------------------------------------------------------

    // Get all bplace files, either in the epa dir, or in its subdirs
    auto bplace_filenames = utils::dir_list_files( epadir, ".*\\.bplace" );
    std::sort( bplace_filenames.begin(), bplace_filenames.end() );

    // Output list of samples in the order that we use them for the matrix
    LOG_INFO << "Writing bplace sample order list";
    std::ofstream bplace_files_of( outdir + "bplace_order.list" );
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
    size_t file_count = 0;
    for( auto const& bplace_filename : bplace_filenames ) {

        // Progress output.
        // LOG_INFO << "at jplace file " << file_count << " (" << bplace_filename << ")";
        ++file_count;

        sset.add( bplace_loader.load( epadir + bplace_filename ), bplace_filename );
    }

    // Final output for jplace reading
    assert( sset.size() == bplace_filenames.size() );
    LOG_INFO << "finished reading sample bplace files";
    LOG_INFO;

    LOG_INFO << "Using " << utils::Options::get().number_of_threads() << " threads.";

    LOG_INFO << "Edpge PCA calculation started";
    auto epca_mat = epca( sset );
    LOG_INFO << "Edpge PCA calculation finished";

    utils::file_write( utils::to_string( epca_mat.projection ), outdir + "epca_projection.mat" );
    utils::file_write( utils::to_string( epca_mat.eigenvectors ), outdir + "epca_eigenvectors.mat" );

    std::ostringstream evout;
    for( auto const& v : epca_mat.eigenvalues ) {
        evout << v << "\n";
    }
    utils::file_write( evout.str(), outdir + "epca_eigenvalues.mat" );

    LOG_INFO << "Finished";
    return 0;
}
