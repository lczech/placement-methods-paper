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
 * @brief Test the speed of our EMD implementation, using different numbers of threads.
 */
int main( int argc, char** argv )
{

    // -------------------------------------------------------------------------
    //     Input shit.
    // -------------------------------------------------------------------------

    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;

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
    //     Read shit.
    // -------------------------------------------------------------------------

    // Get all jplace files, either in the epa dir, or in its subdirs
    auto jplace_filenames = utils::dir_list_files( epadir, true, ".*\\.jplace" );

    // SampleSet sset;
    JplaceReader jplace_reader;

    // Process all jplace files.
    LOG_INFO << "reading " << jplace_filenames.size() << " jplace sample files";
    // for( auto const& jplace_filename : jplace_filenames ) {
    //     sset.add( jplace_reader.from_file( epadir + jplace_filename ), jplace_filename );
    // }
    auto const r_start = std::chrono::steady_clock::now();
    auto const sset = jplace_reader.from_files( jplace_filenames );
    auto const r_duration = std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::steady_clock::now() - r_start
    );

    // Final output for jplace reading
    assert( sset.size() == jplace_filenames.size() );
    LOG_INFO << "finished reading sample jplace files";
    LOG_INFO << "Time: " << r_duration.count() << " s";
    LOG_INFO;

    // -------------------------------------------------------------------------
    //     Calcualte shit.
    // -------------------------------------------------------------------------

    LOG_INFO << "Using " << utils::Options::get().number_of_threads() << " threads.";

    LOG_INFO << "Matrix calculation started";
    auto const c_start = std::chrono::steady_clock::now();
    auto const emd_matrix = earth_movers_distance( sset );
    auto const c_duration = std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::steady_clock::now() - c_start
    );
    LOG_INFO << "Matrix calculation finished";
    LOG_INFO << "Time: " << c_duration.count() << " s";

    utils::file_write( utils::to_string( emd_matrix ), outdir + "emd_" + std::to_string(threads) + ".mat" );

    LOG_INFO << "Finished";
    return 0;
}
