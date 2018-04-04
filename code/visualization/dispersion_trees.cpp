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
#include <cmath>
#include <string>
#include <unordered_map>
#include <unordered_set>

using namespace genesis;
using namespace genesis::utils;

// =================================================================================================
//      Counts To Colorss
// =================================================================================================

struct ColPalPair
{
    std::vector<utils::Color> colors;
    // ColorPalette palette;
    std::pair<SvgGradientLinear, SvgGroup> palette;
};

/**
 * @brief Given a vector of doubles, return a vector of Colors representing the distribution
 * of the double values
 */
ColPalPair cov_to_colors(
    std::vector<double> const& disp_vector
) {
    // Create the result vector.
    std::vector<utils::Color> color_vector( disp_vector.size(), Color( 0, 1, 0 ) );

    // Find min max
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::lowest();
    size_t has_nonfinite = 0;
    for( auto e : disp_vector ) {
        if( ! std::isfinite( e ) ) {
            ++has_nonfinite;
            continue;
        }
        if( e < min ) {
            min = e;
        }
        if( e > max ) {
            max = e;
        }
    }

    LOG_DBG1 << "size    " << disp_vector.size();
    LOG_DBG1 << "non fin " << has_nonfinite;
    LOG_DBG1 << "min     " << min;
    LOG_DBG1 << "max     " << max;

    // Create a color gradient in "blue pink black".
    // auto color_range   = std::map<double, utils::Color>();
    // auto pal = ColorPalette( color_list_spectral() );
    // auto pal = ColorPalette( color_list_bupubk() );
    auto map = ColorMap( color_list_viridis() );
    // auto pal = ColorPalette( color_list_spectral() );
    // auto pal = ColorPalette( color_list_bupubk() );
    map.reverse(true);
    auto norm = utils::ColorNormalizationLinear();

    LOG_DBG1 << "gradient: 0.0, " << ", " << max;
    norm.scale( 0.0, max );

    // Calculate the resulting colors.
    for( size_t i = 0; i < disp_vector.size(); ++i ) {
        if( ! std::isfinite( disp_vector[i] ) ) {
            color_vector[i] = Color( 0.6, 0.6, 0 );
            continue;
        }

        // color_vector[i] = gradient( color_range, disp_vector[i] );
        color_vector[i] = map( norm, disp_vector[i] );
    }

    auto sp = SvgColorBarSettings();
    // sp.height = svg_doc.bounding_box().height() / 2.0;
    // sp.width = sp.height / 10.0;
    auto spg = make_svg_color_bar( sp, map, norm );

    return { color_vector, spg };
}

ColPalPair vmr_to_colors(
    std::vector<double> const& vmrs
) {
    // Create the result vector.
    std::vector<utils::Color> color_vector( vmrs.size(), Color( 0, 1, 0 ) );

    // Find min max
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::lowest();
    size_t has_nonfinite = 0;
    for( auto e : vmrs ) {
        if( ! std::isfinite( e ) ) {
            ++has_nonfinite;
            continue;
        }
        if( e < min ) {
            min = e;
        }
        if( e > max ) {
            max = e;
        }
    }

    LOG_DBG1 << "size    " << vmrs.size();
    LOG_DBG1 << "non fin " << has_nonfinite;
    LOG_DBG1 << "min     " << min;
    LOG_DBG1 << "max     " << max;

    // double center = 1.0;
    // LOG_DBG1 << "gradient: 0.0, " << center << ", " << max;
    // if( max <= center ) {
    //     center = max / 2.0;
    //     LOG_DBG1 << "now: gradient: 0.0, " << center << ", " << max;
    // }

    // Create a color gradient in "blue pink black".
    // auto color_range   = std::map<double, utils::Color>();
    // auto norm = utils::ColorNormalizationDiverging( 0.0, center, max );
    auto norm = utils::ColorNormalizationLogarithmic( 1.0, max );
    auto map = ColorMap( color_list_viridis() );
    // auto pal = ColorPalette( color_list_spectral() );
    // auto pal = ColorPalette( color_list_bupubk() );
    map.reverse(true);

    map.clip_under(true);
    // map.clip_over(true);
    // if( max > 10 ) {
    //     max = 10;
    // }

    // Calculate the resulting colors.
    for( size_t i = 0; i < vmrs.size(); ++i ) {
        if( ! std::isfinite( vmrs[i] ) ) {
            color_vector[i] = Color( 0.6, 0.6, 0 );
            continue;
        }

        // color_vector[i] = gradient( color_range, vmrs[i] );
        color_vector[i] = map( norm, vmrs[i] );
    }

    auto sp = SvgColorBarSettings();
    // sp.height = svg_doc.bounding_box().height() / 2.0;
    // sp.width = sp.height / 10.0;
    auto spg = make_svg_color_bar( sp, map, norm );

    return { color_vector, spg };
}

ColPalPair qcod_to_colors(
    std::vector<double> const& disp_vector
) {
    // Create the result vector.
    std::vector<utils::Color> color_vector( disp_vector.size(), Color( 0, 1, 0 ) );

    // Find min max
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::lowest();
    size_t has_nonfinite = 0;
    for( auto e : disp_vector ) {
        if( ! std::isfinite( e ) ) {
            ++has_nonfinite;
            continue;
        }
        if( e < min ) {
            min = e;
        }
        if( e > max ) {
            max = e;
        }
    }

    LOG_DBG1 << "size    " << disp_vector.size();
    LOG_DBG1 << "non fin " << has_nonfinite;
    LOG_DBG1 << "min     " << min;
    LOG_DBG1 << "max     " << max;

    // Create a color gradient in "blue pink black".
    // auto color_range   = std::map<double, utils::Color>();
    // auto pal = ColorPalette( color_list_spectral() );
    // auto pal = ColorPalette( color_list_bupubk() );
    auto map = ColorMap( color_list_viridis() );
    // auto pal = ColorPalette( color_list_spectral() );
    // auto pal = ColorPalette( color_list_bupubk() );
    map.reverse(true);


    LOG_DBG1 << "gradient: 0.0, " << ", " << max;
    // pal.normalization().scale( min, max );
    auto norm = utils::ColorNormalizationLinear( min, max );

    // Calculate the resulting colors.
    for( size_t i = 0; i < disp_vector.size(); ++i ) {
        if( ! std::isfinite( disp_vector[i] ) ) {
            color_vector[i] = Color( 0.6, 0.6, 0 );
            continue;
        }

        // color_vector[i] = gradient( color_range, disp_vector[i] );
        color_vector[i] = map( norm, disp_vector[i] );
    }

    auto sp = SvgColorBarSettings();
    // sp.height = svg_doc.bounding_box().height() / 2.0;
    // sp.width = sp.height / 10.0;
    auto spg = make_svg_color_bar( sp, map, norm );

    return { color_vector, spg };
}

ColPalPair variances_to_colors(
    std::vector<double> const& variances
) {
    // auto pal = ColorPalette( color_list_viridis() );
    // // auto pal = ColorPalette( color_list_spectral() );
    // // auto pal = ColorPalette( color_list_bupubk() );
    // pal.reverse(true);
    //
    // pal.range( variances );
    // pal.min( 0.0 );
    // if( pal.max() > 1000 ) {
    //     pal.max( 1000.0 );
    //     pal.mid( 1.0 );
    // }
    // pal.clip(true);
    //
    // LOG_DBG1 << "size    " << variances.size();
    // LOG_DBG1 << "min     " << pal.min();
    // LOG_DBG1 << "max     " << pal.max();
    // pal.min( 0.0 );
    //
    // return { pal.sequential_colors( variances ), pal };

    auto map = ColorMap( color_list_viridis() );
    // auto pal = ColorPalette( color_list_spectral() );
    // auto pal = ColorPalette( color_list_bupubk() );
    map.reverse(true);

    // Find min max
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::lowest();
    size_t has_nonfinite = 0;
    for( auto e : variances ) {
        if( ! std::isfinite( e ) ) {
            ++has_nonfinite;
            continue;
        }
        if( e < min ) {
            min = e;
        }
        if( e > max ) {
            max = e;
        }
    }

    LOG_DBG1 << "size    " << variances.size();
    LOG_DBG1 << "non fin " << has_nonfinite;
    LOG_DBG1 << "min     " << min;
    LOG_DBG1 << "max     " << max;

    // Only use log if we have a large scale
    if( max < 1.0 ) {
        auto norm = utils::ColorNormalizationLinear( variances );
        // pal.normalization().autoscale( variances );
        norm.min_value( 0.0 );
        // pal.clip(true);

        auto sp = SvgColorBarSettings();
        // sp.height = svg_doc.bounding_box().height() / 2.0;
        // sp.width = sp.height / 10.0;
        auto spg = make_svg_color_bar( sp, map, norm );

        return { map( norm, variances ), spg };
    }

    // pal.range( variances );
    // pal.min( 0.0 );
    // if( pal.max() > 10 ) {
    //     pal.max( 10.0 );
    //     pal.mid( 1.0 );
    // }
    map.clip_under(true);
    // auto norm = utils::ColorNormalization( variances );
    auto norm = utils::ColorNormalizationLogarithmic( variances );
    norm.min_value( 1.0 );
    // pal.mid( 1.0001 );
    // norm.max( log(max) );
    norm.max_value( max );

    LOG_DBG1 << "using log:";
    LOG_DBG1 << "min     " << norm.min_value();
    LOG_DBG1 << "max     " << norm.max_value();

    // Calculate the resulting colors.
    std::vector<utils::Color> color_vector( variances.size(), Color( 0, 1, 0 ) );
    for( size_t i = 0; i < variances.size(); ++i ) {
        if( ! std::isfinite( variances[i] ) ) {
            color_vector[i] = Color( 0.6, 0.6, 0 );
            continue;
        }

        // auto pos = log(variances[i]); // / log(max);
        // color_vector[i] = map( norm, pos );
        color_vector[i] = map( norm, variances[i] );
    }

    auto sp = SvgColorBarSettings();
    // sp.height = svg_doc.bounding_box().height() / 2.0;
    // sp.width = sp.height / 10.0;
    auto spg = make_svg_color_bar( sp, map, norm );
    return { color_vector, spg };
}

// =================================================================================================
//     Write Color Tree To Nexus
// =================================================================================================

/**
 * @brief Write a nexus file containing a tree with colored branches.
 *
 * The file format can be read and visualized by, e.g., FigTree.
 *
 * The nexus classes of genesis are currently only rudimentary. They do their job, but are not
 * particularly nice to use. This might change in the future.
 */
void write_color_tree_to_nexus(
    placement::PlacementTree const&  tree,
    ColPalPair const& colors_per_branch,
    std::string const&               nexus_filename
) {
    // We use a normal Newick writer for PlacementTrees, but also wrap it in a Color Mixin
    // in order to allow for color annotated branches.
    auto writer = placement::PlacementTreeNewickWriter();
    auto color_plugin = tree::NewickColorWriterPlugin();
    color_plugin.register_with( writer );

    // Get the Newick representation of the tree, with color annotated branches.
    writer.enable_edge_nums(false);
    color_plugin.edge_colors(colors_per_branch.colors);
    std::string newick_tree = writer.to_string(tree);

    // Create an (empty) Nexus document.
    auto nexus_doc = utils::NexusDocument();

    // Add the taxa of the tree to the document.
    auto taxa = utils::make_unique<utils::NexusTaxa>();
    taxa->add_taxa(node_names(tree));
    nexus_doc.set_block( std::move(taxa) );

    // Add the tree itself to the document.
    auto trees = utils::make_unique<utils::NexusTrees>();
    trees->add_tree( "tree1", newick_tree );
    nexus_doc.set_block( std::move(trees) );

    // Write the document to a Nexus file.
    auto nexus_writer = utils::NexusWriter();
    nexus_writer.to_file( nexus_doc, nexus_filename );
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
    auto copy = tree;
    ladderize(copy);

    // Make a layout tree.
    // auto layout = tree::RectangularLayout( tree );
    auto layout = tree::CircularLayout( copy );
    // layout.scaler_r( 500.0 );
    // layout.text_template().font.size /= 3.0;

    // Set edge colors.
    std::vector<utils::SvgStroke> strokes;
    for( auto color : colors_per_branch.colors ) {
        strokes.push_back( utils::SvgStroke() );
        strokes.back().color = color;
        strokes.back().line_cap = utils::SvgStroke::LineCap::kRound;
        // strokes.back().width = 3;
        strokes.back().width = 5;
    }
    layout.set_edge_strokes( strokes );

    // Prepare svg doc.
    std::ostringstream out;
    auto svg_doc = layout.to_svg_document();

    // Add palette
    // auto sp = SvgPalette( colors_per_branch.palette );
    // sp.height = svg_doc.bounding_box().height() / 2.0;
    // sp.width = sp.height / 10.0;
    // auto spg = sp.make();
    auto spg = colors_per_branch.palette;
    spg.second.transform.append( SvgTransform::Translate(
        svg_doc.bounding_box().width() / 2.0, svg_doc.bounding_box().height() / 2.0
    ));
    svg_doc.defs.push_back( spg.first );
    svg_doc.add( spg.second );

    // Write to file
    svg_doc.margin.left = svg_doc.margin.top = svg_doc.margin.bottom = svg_doc.margin.right = 200;
    svg_doc.write( out );
    utils::file_write( out.str(), svg_filename );
}

// =================================================================================================
//     Read Matrix
// =================================================================================================

Matrix<double> read_double_matrix( std::string const& filename )
{
    auto reader = CsvReader();
    reader.separator_chars( " " );
    auto table  = reader.from_file( filename );

    if( table.size() == 0 ) {
        return {};
    }

    // Get and check dimensions.
    auto const row_num = table.size();
    auto const col_num = table[0].size();
    for( auto const& row : table ) {
        if( row.size() != col_num ) {
            throw std::invalid_argument( "Input is not a matrix." );
        }
    }

    // Convert to double.
    auto dmat = Matrix<double>( row_num, col_num );
    for( size_t i = 0; i < row_num; ++i ) {
        for( size_t j = 0; j < col_num; ++j ) {
            dmat( i, j ) = stod( table[i][j] );
        }
    }

    return dmat;
}

// =================================================================================================
//     Process
// =================================================================================================

void process(
    tree::DefaultTree const& tree,
    Matrix<double> const& edge_w,
    Matrix<double> edge_i,
    std::string const& out_pref
) {
    LOG_INFO << "tree: " << tree.edge_count() << " edges";
    LOG_INFO << "edge_w: " << edge_w.cols() << " cols and " << edge_w.rows() << " rows";
    LOG_INFO << "edge_i: " << edge_i.cols() << " cols and " << edge_i.rows() << " rows";

    if( edge_w.cols() != tree.edge_count() ) {
        throw std::runtime_error( "Edge weights does not have corrent length" );
    }
    if( edge_i.cols() != tree.edge_count() ) {
        // throw std::runtime_error( "Edge imbalance does not have corrent length" );

        // Get the indices of all edges that do not lead to a tip.
        std::vector<size_t> inner_edge_indices;
        for( auto const& edge_it : tree.edges() ) {
            if( edge_it->secondary_node().is_inner() ) {
                inner_edge_indices.push_back( edge_it->index() );
            }
        }

        if( inner_edge_indices.size() != edge_i.cols() ) {
            throw std::runtime_error( "Edge imbalance does not have corrent length" );
        }

        // Reverse copy the values to a bigger matrix.
        auto imbalance_matrix = utils::Matrix<double>( edge_i.rows(), tree.edge_count(), 0.0 );
        for( size_t s = 0; s < edge_i.rows(); ++s ) {
            for( size_t i = 0; i < inner_edge_indices.size(); ++i ) {
                auto idx = inner_edge_indices[i];
                imbalance_matrix( s, idx ) = edge_i( s, i );
            }
        }
        edge_i = imbalance_matrix;
    }

    // Last checks
    if( edge_i.cols() != tree.edge_count() ) {
        throw std::runtime_error( "Edge imbalance does not have corrent length" );
    }
    if( edge_w.rows() != edge_i.rows() ) {
        throw std::runtime_error( "Edge matrices do not have correct size" );
    }

    LOG_INFO << "Calculating stuff";

    // LOG_DBG1 << "edge_w\n" << edge_w;
    // LOG_DBG1 << "edge_i\n" << edge_i;

    // Make imbalances absolute.
    // for( auto& e : edge_i ) {
    //     if( e < 0.0 ) {
    //         e = -e;
    //     }
    // }

    auto w_mean_stddev = matrix_col_mean_stddev( edge_w );
    auto i_mean_stddev = matrix_col_mean_stddev( edge_i );
    auto w_quarties = matrix_col_quartiles( edge_w );
    auto i_quarties = matrix_col_quartiles( edge_i );

    auto w_cv   = std::vector<double>( w_mean_stddev.size(), 0.0 );
    auto i_cv   = std::vector<double>( i_mean_stddev.size(), 0.0 );
    auto w_vmr  = std::vector<double>( w_mean_stddev.size(), 0.0 );
    auto i_vmr  = std::vector<double>( i_mean_stddev.size(), 0.0 );
    auto w_qcod = std::vector<double>( w_mean_stddev.size(), 0.0 );
    auto i_qcod = std::vector<double>( i_mean_stddev.size(), 0.0 );
    auto w_var  = std::vector<double>( w_mean_stddev.size(), 0.0 );
    auto i_var  = std::vector<double>( i_mean_stddev.size(), 0.0 );

    for( size_t i = 0; i < w_mean_stddev.size(); ++i ) {
        w_cv[ i ]  = w_mean_stddev[ i ].stddev / w_mean_stddev[ i ].mean;
        i_cv[ i ]  = std::abs( i_mean_stddev[ i ].stddev / i_mean_stddev[ i ].mean );
        w_vmr[ i ] = w_mean_stddev[ i ].stddev * w_mean_stddev[ i ].stddev / w_mean_stddev[ i ].mean;
        i_vmr[ i ] = std::abs(i_mean_stddev[ i ].stddev * i_mean_stddev[ i ].stddev / i_mean_stddev[ i ].mean );
        w_var[ i ] = w_mean_stddev[ i ].stddev * w_mean_stddev[ i ].stddev;
        i_var[ i ] = i_mean_stddev[ i ].stddev * i_mean_stddev[ i ].stddev;

        // if( w_vmr[ i ] > 1000 ) {
        //     LOG_DBG << "w_vmr[ i ] = " << w_vmr[ i ];
        //     for( auto e : edge_w.col(i) ) {
        //         LOG_DBG1 << e;
        //     }
        // }

        // if( w_var[ i ] > 29490 ) {
        //     LOG_DBG << "w_var[ i ] = " << w_var[ i ];
        //     for( auto e : edge_w.col(i) ) {
        //         LOG_DBG1 << e;
        //     }
        // }

        // if( i_vmr[ i ] < 0.0 ) {
            // LOG_DBG << "i_vmr[ " << i << " ] " << i_vmr[ i ] << " " << i_mean_stddev[ i ].stddev << " " << i_mean_stddev[ i ].mean;
        // }

        // auto qw = matrix_col_quartiles( edge_w, i );
        // auto qi = matrix_col_quartiles( edge_i, i );
        // LOG_DBG1 << "qw i " << i << " q0 " << qw.q0 << " q1 " << qw.q1 << " q2 " << qw.q2 << " q3 " << qw.q3 << " q4 " << qw.q4;
        // LOG_DBG1 << "qi i " << i << " q0 " << qi.q0 << " q1 " << qi.q1 << " q2 " << qi.q2 << " q3 " << qi.q3 << " q4 " << qi.q4;

        w_qcod[i] = quartile_coefficient_of_dispersion( w_quarties[i] );
        i_qcod[i] = std::abs( quartile_coefficient_of_dispersion( i_quarties[i] ));
    }

    LOG_DBG << "w cv";
    auto color_w_cv = cov_to_colors( w_cv );

    LOG_DBG << "i cv";
    auto color_i_cv = cov_to_colors( i_cv );

    LOG_DBG << "w vmr";
    auto color_w_vmr = vmr_to_colors( w_vmr );

    LOG_DBG << "i vmr";
    auto color_i_vmr = vmr_to_colors( i_vmr );

    // LOG_DBG << "w qcod";
    // auto color_w_qcod = qcod_to_colors( w_qcod );
    //
    // LOG_DBG << "i qcod";
    // auto color_i_qcod = qcod_to_colors( i_qcod );

    LOG_DBG << "w var";
    auto color_w_var = variances_to_colors( w_var );

    LOG_DBG << "i var";
    auto color_i_var = variances_to_colors( i_var );

    LOG_INFO << "Writing stuff";

    write_color_tree_to_nexus( tree, color_w_cv,  out_pref + "disp_edge_weights_cv.nexus" );
    write_color_tree_to_svg(   tree, color_w_cv,  out_pref + "disp_edge_weights_cv.svg" );
    write_color_tree_to_nexus( tree, color_i_cv,  out_pref + "disp_edge_imbalance_cv.nexus" );
    write_color_tree_to_svg(   tree, color_i_cv,  out_pref + "disp_edge_imbalance_cv.svg" );

    write_color_tree_to_nexus( tree, color_w_vmr, out_pref + "disp_edge_weights_vmr.nexus" );
    write_color_tree_to_svg(   tree, color_w_vmr, out_pref + "disp_edge_weights_vmr.svg" );
    write_color_tree_to_nexus( tree, color_i_vmr, out_pref + "disp_edge_imbalance_vmr.nexus" );
    write_color_tree_to_svg(   tree, color_i_vmr, out_pref + "disp_edge_imbalance_vmr.svg" );

    // write_color_tree_to_nexus( tree, color_w_qcod, out_pref + "disp_edge_weights_qcod.nexus" );
    // write_color_tree_to_svg(   tree, color_w_qcod, out_pref + "disp_edge_weights_qcod.svg" );
    // write_color_tree_to_nexus( tree, color_i_qcod, out_pref + "disp_edge_imbalance_qcod.nexus" );
    // write_color_tree_to_svg(   tree, color_i_qcod, out_pref + "disp_edge_imbalance_qcod.svg" );

    write_color_tree_to_svg(   tree, color_w_var, out_pref + "disp_edge_weights_variance.svg" );
    write_color_tree_to_svg(   tree, color_i_var, out_pref + "disp_edge_imbalance_variance.svg" );
}

// =================================================================================================
//     Process 4
// =================================================================================================

void process_4( int argc, char** argv ) {
    // Check if the command line contains the right number of arguments.
    if (argc != 5) {
        throw std::runtime_error(
            "Need to provide 4 arguments.\n"
        );
    }

    // In out dirs.
    auto const tree_file   = std::string( argv[1] );
    auto const edge_w_file = std::string( argv[2] );
    auto const edge_i_file = std::string( argv[3] );
    auto const out_pref    = std::string( argv[4] );

    LOG_INFO << "Reading stuff";

    auto tree = tree::DefaultTreeNewickReader().from_file( tree_file );

    auto edge_w = read_double_matrix( edge_w_file );
    auto edge_i = read_double_matrix( edge_i_file );

    process( tree, edge_w, edge_i, out_pref );
}

// =================================================================================================
//     Process 2
// =================================================================================================

void process_2( int argc, char** argv ) {
    using namespace genesis::placement;

    // Check if the command line contains the right number of arguments.
    if (argc != 3) {
        throw std::runtime_error(
            "Need to provide 2 arguments.\n"
        );
    }

    auto const epadir   = std::string( argv[1] );
    auto const out_pref = std::string( argv[2] );

    SampleSet sample_set;

    // Get all bplace files, either in the epa dir, or in its subdirs
    auto jplace_filenames = utils::dir_list_files( epadir, false, ".*\\.jplace" );
    std::sort( jplace_filenames.begin(), jplace_filenames.end() );

    if( jplace_filenames.size() > 0 ) {
        // order.clear();
        for( auto& e : jplace_filenames ) {
            // order.push_back( utils::replace_all( e, ".jplace", "" ) );
            e = epadir + e;
        }

        LOG_INFO << "reading " << jplace_filenames.size() << " jplace sample files";
        auto jplace_reader = JplaceReader();
        sample_set = jplace_reader.from_files( jplace_filenames );
        assert( sample_set.size() == jplace_filenames.size() );
        LOG_INFO << "finished reading sample jplace files";
        LOG_INFO;
    }

    auto const tree = average_branch_length_tree( sample_set );
    auto const edge_i = epca_imbalance_matrix( sample_set, true );
    auto const edge_w = placement_weight_per_edge( sample_set );

    process( tree, edge_w, edge_i, out_pref );
}

// =================================================================================================
//     Main
// =================================================================================================

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.date = true;
    utils::Logging::details.time = true;
    utils::Options::get().number_of_threads( 0 );
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    if( argc == 3 ) {
        process_2( argc, argv );
    } else if( argc == 5 ) {
        process_4( argc, argv );
    } else {
        throw std::runtime_error(
            "Need to provide 2 or 4 arguments: (tree edge_weights edge_imbalance out_prefix) "
            "or (jplace out_prefix)\n"
        );
    }

    LOG_INFO << "Finished";

    return 0;
}
