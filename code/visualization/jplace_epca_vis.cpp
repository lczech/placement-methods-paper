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
using namespace genesis::utils;

// =================================================================================================
//      EdgeOCA To Colors
// =================================================================================================

struct ColPalPair
{
    std::vector<utils::Color> colors;
    std::pair<SvgGradientLinear, SvgGroup> palette;
};

ColPalPair epca_to_colors(
    tree::Tree const& tree,
    EpcaData const& data,
    size_t eigenvector_index
) {
    // Create the result vector. Init to grey version of color at 0.0 of spectral palette.
    std::vector<utils::Color> color_vector( tree.edge_count(), utils::Color( 0.874509804, 0.874509804, 0.874509804 ) );

    auto const eigen_vec = data.eigenvectors.col( eigenvector_index );

    // Use "spectral" color palette.
    auto norm = utils::ColorNormalizationDiverging( eigen_vec );
    norm.make_centric();
    // auto pal = utils::ColorPalette( utils::color_list_spectral(), norm );
    auto map = ColorMap( utils::color_list_spectral() );

    // Get the column we are interested in.
    auto const eigen_color_vector = map( norm, eigen_vec );
    assert( data.edge_indices.size() == eigen_color_vector.size() );

    // For each edge that has an eigenvector, get its color and store it.
    for( size_t i = 0; i < data.edge_indices.size(); ++i ) {
        auto const edge_index = data.edge_indices[i];

        color_vector[ edge_index ] = eigen_color_vector[i];
    }

    auto sp = SvgColorBarSettings();
    // sp.height = svg_doc.bounding_box().height() / 2.0;
    // sp.width = sp.height / 10.0;
    auto spg = make_svg_color_bar( sp, map, norm );

    return { color_vector, spg };
}

// =================================================================================================
//     Write Color Tree To SVG
// =================================================================================================

/**
 * @brief Write an SVG file containing a tree with colored branches.
 */
void write_color_tree_to_svg(
    placement::PlacementTree const&  tree,
    ColPalPair const& colors_per_branch,
    std::string const&               svg_filename
) {
    if( colors_per_branch.colors.size() != tree.edge_count() ) {
        throw std::runtime_error( "Branch colors vec of wrong size" );
    }

    auto copy = tree;
    ladderize(copy);

    // Make a layout tree.
    // auto layout = tree::RectangularLayout( tree );
    auto layout = tree::CircularLayout( copy, tree::LayoutType::kCladogram );
    // layout.scaler_r( 500.0 );
    // layout.text_template().font.size /= 3.5;

    // Set edge colors.
    std::vector<utils::SvgStroke> strokes;
    for( auto color : colors_per_branch.colors ) {
        strokes.push_back( utils::SvgStroke() );
        strokes.back().color = color;
        strokes.back().line_cap = utils::SvgStroke::LineCap::kRound;
        // strokes.back().width = 3;
        strokes.back().width = 8;
    }
    layout.set_edge_strokes( strokes );

    // Prepare svg doc.
    std::ostringstream out;
    auto svg_doc = layout.to_svg_document();

    // Add palette
    auto svgpal = colors_per_branch.palette;
    svgpal.second.transform.append( SvgTransform::Translate(
        svg_doc.bounding_box().width() / 2.0, svg_doc.bounding_box().height() / 2.0
    ));
    svg_doc.defs.push_back( svgpal.first );
    svg_doc.add( svgpal.second );

    // Write to file
    svg_doc.margin.left = svg_doc.margin.top = svg_doc.margin.bottom = svg_doc.margin.right = 200;
    svg_doc.write( out );
    utils::file_write( out.str(), svg_filename );
}

// =================================================================================================
//     Write Svg Plot of PCA Matrix
// =================================================================================================

// void write_pca_mat_to_svg(
//     utils::Matrix<double> const&     mat,
//     // std::vector<utils::Color> const& colors_per_row,
//     std::string const&               svg_filename
// )  {
//     using namespace genesis::utils;
//     if( mat.cols() < 2 ) {
//         throw std::runtime_error( "mat too small" );
//     }
//
//     double min_x = std::numeric_limits<double>::max();
//     double max_x = std::numeric_limits<double>::lowest();
//     double min_y = std::numeric_limits<double>::max();
//     double max_y = std::numeric_limits<double>::lowest();
//     for( size_t r = 0; r < mat.rows(); ++r ) {
//         min_x = std::min( min_x, mat(r,0) );
//         max_x = std::max( max_x, mat(r,0) );
//         min_y = std::min( min_y, mat(r,1) );
//         max_y = std::max( max_y, mat(r,1) );
//     }
//     double rx = 20.0 / ( max_x - min_x );
//     double ry = 20.0 / ( max_y - min_y );
//
//     SvgDocument doc;
//     for( size_t r = 0; r < mat.rows(); ++r ) {
//         doc << SvgCircle(
//            utils::SvgPoint( rx * mat(r,0), ry * mat(r,1) ),
//            1.0,
//            utils::SvgStroke( SvgStroke::Type::kNone ),
//            utils::SvgFill( Color() )
//            // utils::SvgFill( colors_per_row[r] )
//        );
//     }
//
//     std::ostringstream out;
//     doc.write( out );
//     utils::file_write( out.str(), svg_filename );
// }

// =================================================================================================
//     Meta Data Prep
// =================================================================================================

// struct MetaData
// {
//     std::vector<std::string> names;
//     utils::Matrix<double>    data;
// };
//
// MetaData meta_data_prep(
//     std::string const& meta_file,
//     size_t edge_w_rows,
//     std::vector<std::string> const& order,
//     // std::string const& out_pref
// ) {
//     // Read metadata table
//     LOG_INFO << "Reading meta data";
//     auto reader = CsvReader();
//     reader.separator_chars( "\t" );
//     auto meta_table = reader.from_file( meta_file );
//
//     // Meta data needs to have one more row (header line) than the matrices.
//     if( meta_table.size() != edge_w_rows + 1 ) {
//         throw std::runtime_error( "Meta table has not length of matrix rows." );
//     }
//
//     // Get column names of the meta data
//     auto meta_names = std::vector<std::string>();
//     for( size_t i = 1; i < meta_table[0].size(); ++i ) {
//         meta_names.push_back( meta_table[0][i] );
//     }
//
//     // Make a table of meta data for all columns, and fill it with the data, sorted
//     // by the sample names, so that each row corresponds to the other matrix rows.
//     auto meta_data = utils::Matrix<double>( edge_w_rows, meta_names.size() );
//     if( meta_data.rows() != order.size() ) {
//         throw std::runtime_error( "Impossibru." );
//     }
//     for( size_t r = 0; r < meta_data.rows(); ++r ) {
//
//         // For the given matrix row, find the line in the meta table that corresponds to this sample.
//         auto const& sample_name = order[r];
//         size_t meta_index = meta_table.size();
//         for( size_t j = 1; j < meta_table.size(); ++j ) {
//             if( meta_table[j][0] == sample_name ) {
//                 meta_index = j;
//                 break;
//             }
//         }
//         if( meta_index == meta_table.size() ) {
//             throw std::runtime_error( "Cannot find meta data table index for sample " + sample_name );
//         }
//         if(  meta_table[meta_index].size() != meta_data.cols() + 1 ) {
//             throw std::runtime_error( "Meta data table has wrong line length at " + std::to_string(r) );
//         }
//
//         // Now, fill the data columns.
//         for( size_t c = 0; c < meta_data.cols(); ++c ) {
//             // LOG_DBG1 << meta_table[meta_index][c+1];
//             meta_data( r, c ) = stod( meta_table[meta_index][c+1] );
//         }
//     }
//
//     // Write meta matrix.
//     // LOG_INFO << "Writing meta data matrix";
//     // utils::file_write( utils::to_string( meta_data ), out_pref + "metadata.mat" );
//
//     return { meta_names, meta_data };
// }

// =================================================================================================
//     Get Meta
// =================================================================================================

// std::vector<utils::Color> meta_colors( std::string const& meta_file, std::vector<std::string> const& order )
// {
//     auto const meta_data = meta_data_prep(  meta_file, edge_w.rows(), order, out_pref );
//
// }

// =================================================================================================
//     Main
// =================================================================================================

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;

    utils::Options::get().number_of_threads( 4 );
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    LOG_INFO << "Started";

    // Check if the command line contains the right number of arguments.
    if (argc != 3 && argc != 4) {
        throw std::runtime_error(
            "Need to provide two or three arguments.\n"
        );
    }

    // In out dirs.
    auto epadir = utils::dir_normalize_path( std::string( argv[1] ));
    auto outdir = utils::dir_normalize_path( std::string( argv[2] ));
    utils::dir_create(outdir);

    std::string metafile;
    if( argc == 4 ) {
        metafile = std::string( argv[2] );
    }

    // -------------------------------------------------------------------------
    //     Find and read sample jplace files.
    // -------------------------------------------------------------------------

    // Get all jplace files, either in the epa dir, or in its subdirs
    auto jplace_filenames = utils::dir_list_files( epadir, ".*\\.jplace" );
    std::sort( jplace_filenames.begin(), jplace_filenames.end() );

    LOG_INFO << "reading " << jplace_filenames.size() << " jplace sample files";
    auto const sample_set = JplaceReader().from_files( jplace_filenames );
    assert( sample_set.size() == jplace_filenames.size() );
    LOG_INFO << "finished reading sample jplace files";
    LOG_INFO;

    // -------------------------------------------------------------------------
    //     Calc and write edge pca
    // -------------------------------------------------------------------------

    LOG_INFO << "Using " << utils::Options::get().number_of_threads() << " threads.";

    size_t const components = 5;

    LOG_INFO << "Edge PCA calculation started";
    auto const epca_data = epca( sample_set, 1.0, -1.0, components );
    LOG_INFO << "Edge PCA calculation finished";

    utils::file_write( utils::to_string( epca_data.projection ), outdir + "epca_projection.mat" );
    utils::file_write( utils::to_string( epca_data.eigenvectors ), outdir + "epca_eigenvectors.mat" );

    LOG_DBG << "projection   " << epca_data.projection.rows() << " x " << epca_data.projection.cols();

    std::ostringstream proj_out;
    for( size_t r = 0; r < epca_data.projection.rows(); ++r ) {
        proj_out << sample_set[r].name << "," << epca_data.projection( r, 0 ) << "," << epca_data.projection( r, 1 ) << "\n";
    }
    utils::file_write( proj_out.str(), outdir + "proj.csv" );

    LOG_DBG << "eigenvectors " << epca_data.eigenvectors.rows() << " x " << epca_data.eigenvectors.cols();
    auto const ev_mm = utils::matrix_col_minmax( epca_data.eigenvectors );
    for( size_t c = 0; c < ev_mm.size(); ++c ) {
        LOG_DBG1 << "col " << c << " min " << ev_mm[c].min << " max " << ev_mm[c].max;
    }

    LOG_DBG << "eigenvalues  " << epca_data.eigenvalues.size();

    size_t inner_branch_count = 0;
    for( auto const& e : sample_set.at(0).sample.tree().edges() ) {
        if( e->secondary_node().is_inner() ) {
            ++inner_branch_count;
        }
    }
    LOG_DBG << "inner_branch_count " << inner_branch_count;

    std::ostringstream evout;
    for( auto const& v : epca_data.eigenvalues ) {
        evout << v << "\n";
    }
    utils::file_write( evout.str(), outdir + "epca_eigenvalues.mat" );

    // -------------------------------------------------------------------------
    //     write vis files
    // -------------------------------------------------------------------------

    auto const tree = average_branch_length_tree( sample_set );

    for( size_t i = 0; i < components; ++i ) {
        auto const color_vec = epca_to_colors( tree, epca_data, i );

        auto const filepref = outdir + "epca_tree_" + std::to_string( i+1 );
        tree::write_color_tree_to_nexus_file( tree, color_vec.colors, filepref + ".nexus" );
        write_color_tree_to_svg( tree, color_vec, filepref + ".svg" );
    }

    // write_pca_mat_to_svg( epca_data.projection, outdir + "epca_vis.svg" );

    // -------------------------------------------------------------------------
    //     PCA on edge weights
    // -------------------------------------------------------------------------

    // Bullshit. pca on weights is not useful.

    // LOG_INFO << "PCA calculation started";
    // auto const ew_matrix = placement_weight_per_edge( sample_set );
    // auto const pca = utils::principal_component_analysis( ew_matrix, components );
    // LOG_INFO << "PCA calculation finished";
    //
    // utils::file_write( utils::to_string( pca.projection ), outdir + "pca_projection.mat" );
    // utils::file_write( utils::to_string( pca.eigenvectors ), outdir + "pca_eigenvectors.mat" );
    //
    // std::ostringstream pca_proj_out;
    // for( size_t r = 0; r < pca.projection.rows(); ++r ) {
    //     pca_proj_out << sample_set[r].name << "," << pca.projection( r, 0 ) << "," << pca.projection( r, 1 ) << "\n";
    // }
    // utils::file_write( pca_proj_out.str(), outdir + "pca_proj.csv" );

    LOG_INFO << "Finished";
    return 0;
}
