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
#include <fstream>
#include <limits>
#include <string>
#include <unordered_map>

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

std::vector<Color> vec_to_col( std::vector<double> const& values )
{
    auto res = std::vector<Color>( values.size() );

    auto const invalid_value = 99999;

    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::lowest();

    for( size_t i = 0; i < values.size(); ++i ) {
        if( values[i] == invalid_value ) {
            continue;
        }
        if( values[i] < min ) {
            min = values[i];
        }
        if( values[i] > max ) {
            max = values[i];
        }
    }

    // auto const minmax = std::minmax_element( values.begin(), values.end() );
    // auto const min = *minmax.first;
    // auto const max = *minmax.second;

    for( size_t i = 0; i < values.size(); ++i ) {
        if( values[i] == invalid_value ) {
            res[i] = Color( 255, 0, 0 );
            continue;
        }

        auto const rel = ( values[i] - min ) / ( max - min );
        res[i] = color_list_viridis()[ rel * 255 ];
    }

    return res;
}

// =================================================================================================
//      Main
// =================================================================================================

/**
 * Take a cluster tree (from Squash clustering), as well as a metadata file, and annotate the
 * leaves of the tree with colors according to metadta feature values.
 * Because each leaf of the cluster tree corresponds to a samples,
 * there is one leaf per metadata table row. The result then visualizes the metadata distribution
 * of samples on the clsuter tree, as for example shown in the ART evaluation of the BV data,
 * using their Nugent score metadata.
 */
int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
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

    // Newick tree.
    auto const tree_file   = std::string( argv[1] );

    // Meta data table file, tab separted.
    // First row is the header, which contains the names of the meta data columns.
    // First column is the sample name. THose do not need to be in the same order as the
    // matrix files above. In order to sort them so that the lines are fitting, we use the order
    // file.
    auto const meta_file   = std::string( argv[2] );

    // Output prefix, that is, path plus file name prefix.
    auto const out_pref    = utils::dir_normalize_path( std::string( argv[3] ));

    // utils::Options::get().command_line( argc, argv );
    // utils::Options::get().number_of_threads( threads );
    // LOG_BOLD << utils::Options::get().info();
    // LOG_BOLD;

    // -------------------------------------------------------------------------
    //     Read stuff
    // -------------------------------------------------------------------------

    // Cluster tree
    LOG_INFO << "Reading tree";
    auto clust_tree = DefaultTreeNewickReader().from_file( tree_file );
    ladderize(clust_tree);

    // Read metadata table
    LOG_INFO << "Reading meta data";
    auto reader = utils::CsvReader();
    reader.separator_chars( "\t" );
    auto meta_table = reader.from_file( meta_file );

    // Get column names from first row of meta table
    auto meta_names = std::vector<std::string>();
    for( size_t i = 1; i < meta_table[0].size(); ++i ) {
        meta_names.push_back( meta_table[0][i] );
    }

    // Get sample names from first col of meta table
    auto sample_names = std::vector<std::string>();
    for( size_t i = 1; i < meta_table.size(); ++i ) {
        sample_names.push_back( meta_table[i][0] );
    }

    // Go through all meta data.
    for( size_t mi = 0; mi < meta_names.size(); ++mi ) {

        LOG_INFO << "Meta: " << meta_names[mi];

        try {
            stod( meta_table.at(1).at(mi+1) );
        } catch(...) {
            LOG_WARN << "Skipping";
            continue;
        }

        // Fill vec with medatadata for this col
        std::vector<double> data;
        for( size_t r = 0; r < sample_names.size(); ++r ) {
            data.push_back( stod( meta_table.at(r+1).at(mi+1) ));
        }

        // Make nice colors for the metadata
        auto const colors = vec_to_col( data );
        auto layout = RectangularLayout( clust_tree );

        // Set colourful node shapes.
        std::vector<utils::SvgGroup> node_shapes;
        node_shapes.resize( clust_tree.node_count() );
        for( size_t i = 0; i < clust_tree.node_count(); ++i ) {
            if( ! clust_tree.node_at(i).is_leaf() ) {
                continue;
            }

            // Get the sample position in the data matrix.
            auto const& node_name = clust_tree.node_at( i ).data<DefaultNodeData>().name;
            auto const pos_it = std::find( sample_names.begin(), sample_names.end(), node_name );
            if( pos_it == sample_names.end() ) {
                LOG_WARN << "No sample name " << node_name;
                continue;
            }
            auto const pos = pos_it - sample_names.begin();

            // Add a shape according to the sample color.
            node_shapes[i].add( utils::SvgCircle(
                utils::SvgPoint( 0, 0 ),
                5,
                utils::SvgStroke(),
                utils::SvgFill( colors[pos] )
            ));
        }
        layout.set_node_shapes( node_shapes );

        std::ostringstream out;
        layout.to_svg_document().write( out );
        utils::file_write( out.str(), out_pref + "clust_tree_" + meta_names[mi] + ".svg" );
    }

    // // Make a table of meta data for all columns, and fill it with the data
    // auto meta_data = utils::Matrix<double>( sample_names.size(), meta_names.size() );
    // for( size_t r = 0; r < meta_data.rows(); ++r ) {
    //
    //     // Now, fill the data columns.
    //     for( size_t c = 0; c < meta_data.cols(); ++c ) {
    //         // LOG_DBG1 << meta_table[meta_index][c+1];
    //         meta_data( r, c ) = stod( meta_table.at(r+1).at(c+1) );
    //     }
    // }


    LOG_INFO << "Finished";
    return 0;
}
