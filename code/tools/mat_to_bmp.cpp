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

#include <algorithm>

using namespace genesis;
using namespace genesis::utils;

int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    // Check if the command line contains the right number of arguments.
    if (argc != 2) {
        throw std::runtime_error(
            "Need to provide an EMD matrix file.\n"
        );
    }

    // In out dirs.
    auto emd_mat_file = std::string( argv[1] );
    auto emd_bmp_file = utils::file_filename( emd_mat_file ) + ".bmp";

    // Read table.
    LOG_INFO << "Reading file " << emd_mat_file;
    auto reader = CsvReader();
    reader.separator_chars( " " );
    auto table  = reader.from_file( emd_mat_file );
    if( table.size() == 0 ) {
        LOG_INFO << "Empty";
        return 0;
    }

    LOG_INFO << "Processing Matrix";

    // Get and check dimensions.
    auto const row_num = table.size();
    auto const col_num = table[0].size();
    if( row_num != col_num ) {
        throw std::invalid_argument( "Input is not a symmetric matrix." );
    }
    for( auto const& row : table ) {
        if( row.size() != col_num ) {
            throw std::invalid_argument( "Input is not a symmetric matrix." );
        }
    }

    // Convert to double.
    auto dmat = Matrix<double>( row_num, col_num );
    for( size_t i = 0; i < row_num; ++i ) {
        for( size_t j = 0; j < col_num; ++j ) {
            dmat( i, j ) = stod( table[i][j] );
        }
    }
    table.clear();

    // Find max value.
    auto maxy = *std::max_element( dmat.begin(), dmat.end() );

    // // Get row sums.
    // auto row_sums = std::vector<double>( dmat.rows(), 0 );
    // for( size_t i = 0; i < row_num; ++i ) {
    //     for( size_t j = 0; j < col_num; ++j ) {
    //         row_sums[i] += dmat( i, j );
    //     }
    // }
    //
    // // Sort by row sum.
    // auto si = sort_indices( row_sums.begin(), row_sums.end() );
    // auto dmat2 = Matrix<double>( row_num, col_num );
    // for( size_t i = 0; i < row_num; ++i ) {
    //     for( size_t j = 0; j < col_num; ++j ) {
    //         dmat2( i, j ) = dmat( si[i], si[j] );
    //     }
    // }

    LOG_INFO << "Creating Bitmap";

    // Convert to grayscale.
    auto bmat = Matrix<unsigned char>( row_num, col_num );
    for( size_t i = 0; i < row_num; ++i ) {
        for( size_t j = 0; j < col_num; ++j ) {
            bmat( i, j ) = 255.0 * dmat( i, j ) / maxy;
        }
    }

    // // Nice color palette.
    // auto gradient   = std::map<double, utils::Color>();
    // gradient[ 0.0 ] = utils::color_from_hex("#81bfff");
    // gradient[ 0.5 ] = utils::color_from_hex("#c040be");
    // gradient[ 1.0 ] = utils::color_from_hex("#000000");
    // auto palette = std::vector<Color>( 256 );
    // for( size_t i = 0; i < 256; ++i ) {
    //     auto c = utils::gradient( gradient, static_cast<double>(i) / 255.0 );
    //     palette[i].r( c.r() );
    //     palette[i].g( c.g() );
    //     palette[i].b( c.b() );
    // }

    LOG_INFO << "Writing Bitmap";

    // Write to bitmap file.
    // BmpWriter().to_file( bmat, palette, emd_bmp_file );
    BmpWriter().to_file( bmat, emd_bmp_file );

    LOG_INFO << "Finished";
    return 0;
}
