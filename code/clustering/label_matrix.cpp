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
#include <string>
#include <unordered_map>
#include <unordered_set>

using namespace genesis;
using namespace genesis::utils;

// =================================================================================================
//     Main
// =================================================================================================

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.date = true;
    utils::Logging::details.time = true;
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    // Check if the command line contains the right number of arguments.
    if (argc != 5 && argc != 6) {
        throw std::runtime_error(
            "wrong num args"
        );
    }

    // In out dirs.

    // Actual labels that we want for each sample.
    auto const label_file      = std::string( argv[1] );

    // List of the possible labels, in the order that we want for the matrix rows.
    auto const label_list_file = std::string( argv[2] );

    // Cluster assignment for each sample.
    auto const clust_file      = std::string( argv[3] );

    // Output bmp file.
    auto const bmp_file        = std::string( argv[4] );

    bool const normalize = ( argc == 6 );
    if( argc == 6 && std::string(argv[5]) != "normalize" ) {
        throw std::runtime_error(
            "Last arg needs to be 'normalize'"
        );
    }

    auto reader = CsvReader();
    reader.separator_chars( "\t" );

    // Read files.
    LOG_INFO << "Reading file " << label_file;
    auto label_table  = reader.from_file( label_file );
    LOG_INFO << "Lines: " << label_table.size();

    LOG_INFO << "Reading file " << label_list_file;
    auto label_list_table  = reader.from_file( label_list_file );
    LOG_INFO << "Lines: " << label_list_table.size();

    LOG_INFO << "Reading file " << clust_file;
    auto clusts_table  = reader.from_file( clust_file );
    LOG_INFO << "Lines: " << clusts_table.size();

    // Map from sample name to metadata name.
    std::unordered_map<std::string, std::string> label_lookup;

    // Sort labels into lookup table and label set.
    for( auto const& line : label_table ) {
        if( line.size() != 2 ) {
            LOG_WARN << "Label line of size " << line.size();
            continue;
        }

        label_lookup[ line[0] ] = line[ 1 ];

        // if( label_list.count( line [1] ) == 0 ) {
        //     label_list[ line[1] ] = label_list.size();
        // }
    }
    // LOG_INFO << "Distinct labels: " << label_list.size();

    // Map from metadata name to a unique index.
    std::unordered_map<std::string, size_t> label_list;
    size_t max_label_list_index = 0;

    for( auto const& line : label_list_table ) {
        if( line.size() != 1 && line.size() != 2 ) {
            LOG_WARN << "Label list line of size " << line.size();
            continue;
        }

        if( label_list.count( line [0] ) == 0 ) {
            if( line.size() == 1 ) {
                label_list[ line[0] ] = label_list.size();
            } else if( line.size() == 2 ) {
                label_list[ line[0] ] = stoi( line[1] );
            }
        }
        max_label_list_index = std::max( max_label_list_index, label_list[ line[0] ] );
    }

    // Debug output
    LOG_DBG << "max_label_list_index " << max_label_list_index;
    auto label_list_print = std::vector<std::string>( max_label_list_index + 1 );
    LOG_DBG << "metadata name to row number";
    for( auto const& elem : label_list ) {
        label_list_print[ elem.second ] += elem.first + " ";
        // LOG_DBG1 << "label " << elem.first << " -> row " << elem.second;
    }
    for( size_t i = 0; i < label_list_print.size(); ++i ) {
        LOG_DBG1 << "row " << i <<  " label " << label_list_print[i];
    }

    // Map from cluster index to unique index.
    std::unordered_map<std::string, size_t> clust_list;

    // Get cluster names
    for( auto const& line : clusts_table ) {
        if( line.size() != 2 ) {
            LOG_WARN << "Clust line of size " << line.size();
            continue;
        }

        if( clust_list.count( line [1] ) == 0 ) {
            clust_list[ line[1] ] = clust_list.size();
        }
    }
    LOG_INFO << "Distinct clusts: " << clust_list.size();

    // Debug output
    // LOG_DBG << "cluster index to unique id";
    // for( auto const& elem : clust_list ) {
    //     LOG_DBG1 << elem.first << " -> " << elem.second;
    // }

    // Matrix: rows are metadata indices, colums are cluster indices.
    auto count_matrix = Matrix<size_t>( max_label_list_index + 1, clust_list.size(), 0 );
    LOG_INFO << "count_matrix: " << count_matrix.rows() << " x " << count_matrix.cols();

    // Fill matrix
    for( auto const& line : clusts_table ) {
        if( line.size() != 2 ) {
            LOG_WARN << "Clust line of size " << line.size();
            continue;
        }

        auto const sample_name = file_basename( line[0] );
        auto const meta_name = label_lookup[ sample_name ];

        if( label_list.count( meta_name ) == 0 ) {
            LOG_DBG << "skipping " << sample_name << " with " << meta_name;
            continue;
        }
        if( clust_list.count( line[1] ) == 0 ) {
            LOG_DBG << "clust_list misses " << line[1];
        }

        auto const clust_idx = clust_list[ line[1] ];
        auto const label_idx = label_list[ meta_name ];

        ++count_matrix.at( label_idx, clust_idx );
    }


    // Sort by row sum.
    // auto dmat2 = matrix_sort_by_row_sum( count_matrix );

    // // Get row sums.
    // auto row_sums = matrix_row_sums( count_matrix );
    // // Sort by row sum.
    // auto si = sort_indices( row_sums.begin(), row_sums.end() );
    // auto dmat2 = Matrix<double>( label_list.size(), clust_list.size() );
    // for( size_t i = 0; i < label_list.size(); ++i ) {
    //     for( size_t j = 0; j < clust_list.size(); ++j ) {
    //         dmat2( i, j ) = count_matrix( si[i], si[j] );
    //     }
    // }

    // weighted column averages.
    auto col_avgs = std::vector<double>( count_matrix.cols(), 0 );
    for( size_t j = 0; j < count_matrix.cols(); ++j ) {
        double sum = 0;
        for( size_t i = 0; i < count_matrix.rows(); ++i ) {
            col_avgs[ j ] += i * count_matrix( i, j );
            sum += count_matrix( i, j );
        }
        col_avgs[ j ] /= sum;
    }

    // Sort by col avg.
    auto si = sort_indices( col_avgs.begin(), col_avgs.end() );
    auto dmat2 = Matrix<double>( count_matrix.rows(), count_matrix.cols() );
    for( size_t i = 0; i < count_matrix.rows(); ++i ) {
        for( size_t j = 0; j < count_matrix.cols(); ++j ) {
            dmat2( i, j ) = static_cast<double>( count_matrix( i, si[j] ));
        }
    }

    // Normalize columns.
    if( normalize ) {
        auto col_minmax = matrix_col_minmax( dmat2 );
        for( size_t r = 0; r < dmat2.rows(); ++r ) {
            for( size_t c = 0; c < dmat2.cols(); ++c ) {
                dmat2( r, c ) /= col_minmax[c].max;
            }
        }
    }

    // dmat2 = matrix_sort_diagonal_symmetric(dmat2);
    // auto dmat2 = count_matrix;

    // Find max value.
    auto maxy = *std::max_element( dmat2.begin(), dmat2.end() );

    // Convert to grayscale.
    auto bmat = Matrix<unsigned char>( dmat2.rows(), dmat2.cols() );
    for( size_t i = 0; i < dmat2.rows(); ++i ) {
        for( size_t j = 0; j < dmat2.cols(); ++j ) {
            bmat( i, j ) = 255.0 * static_cast<double>( dmat2( i, j )) / static_cast<double>( maxy );
        }
    }

    // auto colors = color_list_viridis();
    // auto cm = ColorMap( color_list_orrd() );
    auto cm = ColorMap( color_list_orrd() );
    auto colors = cm.color_list( 256 );
    // std::reverse( colors.begin(), colors.end() );

    // Write to bitmap file.
    LOG_INFO << "Writing Bitmap";
    BmpWriter().to_file( bmat, colors, bmp_file );

    LOG_INFO << "Finished";

    return 0;
}
