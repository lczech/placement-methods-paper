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

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

#ifdef GENESIS_PTHREADS
#    include <thread>
#endif

using namespace genesis;
using namespace genesis::utils;

// =================================================================================================
//      Counts To Colors
// =================================================================================================

/**
 * @brief Given a vector of doubles, return a vector of Colors representing the PCC values as colors
 */
std::vector<Color> cc_to_colors(
    std::vector<double> const& srcc_vector
) {
    // Create the result vector.
    std::vector<Color> color_vector( srcc_vector.size(), Color( 1, 0, 1 ) );

    // // Create a color gradient in "red grey blue".
    // auto color_range    = std::map<double, Color>();
    // color_range[ -1.0 ] = color_from_hex("#ff0000");
    // color_range[  0.0 ] = color_from_hex("#f0f0f0");
    // color_range[  1.0 ] = color_from_hex("#0000ff");

    // Use "spectral" color palette.
    auto norm = utils::ColorNormalizationDiverging( -1.0, 0.0, 1.0 );
    auto map = ColorMap( color_list_spectral() );

    // Use the grey equivalent of the central (0.0) color of the palette.
    map.mask_color( Color( 0.874509804, 0.874509804, 0.874509804 ) );

    color_vector = map( norm, srcc_vector );

    // // Calculate the resulting colors.
    // for( size_t i = 0; i < srcc_vector.size(); ++i ) {
    //     // if( ! std::isfinite( srcc_vector[i] ) ) {
    //     //     color_vector[i] = Color( 255, 255, 0 );
    //     //     continue;
    //     // }
    //
    //     // color_vector[i] = gradient( color_range, srcc_vector[i] );
    //
    //     color_vector[i] = pal.diverging_color( srcc_vector[i] );
    // }

    return color_vector;
}

std::vector<Color> vec_to_col( std::vector<double> const& values )
{
    // auto pal = ColorPalette( color_list_viridis() );

    std::vector<Color> list = { Color(0,0,0,0), Color(0,0,0,1) };
    auto map = ColorMap( list );
    auto norm = utils::ColorNormalizationLinear();

    norm.mask_value( 99999 );
    norm.autoscale( values.begin(), values.end() );

    return map( norm, values );
    // return pal( values );

    // auto res = std::vector<Color>( values.size() );
    // auto const invalid_value = 99999;
    //
    // double min = std::numeric_limits<double>::max();
    // double max = std::numeric_limits<double>::lowest();
    //
    // for( size_t i = 0; i < values.size(); ++i ) {
    //     if( values[i] == invalid_value ) {
    //         continue;
    //     }
    //     if( values[i] < min ) {
    //         min = values[i];
    //     }
    //     if( values[i] > max ) {
    //         max = values[i];
    //     }
    // }
    //
    // // auto const minmax = std::minmax_element( values.begin(), values.end() );
    // // auto const min = *minmax.first;
    // // auto const max = *minmax.second;
    //
    // for( size_t i = 0; i < values.size(); ++i ) {
    //     if( values[i] == invalid_value ) {
    //         res[i] = Color( 255, 0, 0 );
    //         continue;
    //     }
    //
    //     auto const rel = ( values[i] - min ) / ( max - min );
    //     assert( 0.0 <= rel && rel <= 1.0 );
    //     res[i] = color_list_viridis()[ 255 - rel * 255 ];
    // }
    //
    // return res;
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
    std::vector<utils::Color> const& colors_per_branch,
    std::string const&               nexus_filename
) {
    // We use a normal Newick writer for PlacementTrees, but also wrap it in a Color Mixin
    // in order to allow for color annotated branches.
    auto writer = placement::PlacementTreeNewickWriter();
    auto color_plugin = tree::NewickColorWriterPlugin();
    color_plugin.register_with( writer );

    // Get the Newick representation of the tree, with color annotated branches.
    writer.enable_edge_nums(false);
    color_plugin.edge_colors(colors_per_branch);
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
    std::vector<utils::Color> const& colors_per_branch,
    std::vector<utils::Color> const& colors_per_node,
    std::string const&               svg_filename
) {
    if( colors_per_branch.size() != tree.edge_count() ) {
        throw std::runtime_error( "Branch colors vec of wrong size" );
    }
    if( !colors_per_node.empty() && colors_per_node.size() != tree.edge_count() ) {
        throw std::runtime_error( "Node colors vec of wrong size" );
    }

    auto copy = tree;
    ladderize(copy);

    // Make a layout tree.
    // auto layout = tree::RectangularLayout( tree );
    auto layout = tree::CircularLayout( copy, tree::LayoutType::kCladogram );
    // layout.scaler_r( 500.0 );
    layout.text_template().font.size /= 3.5;

    // Set edge colors.
    std::vector<utils::SvgStroke> strokes;
    for( auto color : colors_per_branch ) {
        strokes.push_back( utils::SvgStroke() );
        strokes.back().color = color;
        strokes.back().line_cap = utils::SvgStroke::LineCap::kRound;
        // strokes.back().width = 3;
        strokes.back().width = 2.3;
    }
    layout.set_edge_strokes( strokes );

    // Set colourful node shapes.
    // std::vector<utils::SvgGroup> node_shapes;
    // node_shapes.resize( tree.node_count() );
    // for( size_t i = 0; i < tree.edge_count(); ++i ) {
    //     // Add a shape according to the sample color.
    //     auto index = tree.edge_at(i).secondary_node().index();
    //     node_shapes[index].add( utils::SvgCircle(
    //         utils::SvgPoint( 0, 0 ),
    //         3.0,
    //         utils::SvgStroke( SvgStroke::Type::kNone ),
    //         utils::SvgFill( colors_per_node[i] )
    //     ));
    // }
    // layout.set_node_shapes( node_shapes );

    // Write to file
    std::ostringstream out;
    auto svg_doc = layout.to_svg_document();
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
//     Meta Data Prep
// =================================================================================================

struct MetaData
{
    std::vector<std::string> names;
    utils::Matrix<double>    data;
};

MetaData meta_data_prep(
    std::string const& meta_file,
    size_t edge_w_rows,
    std::vector<std::string> const& order,
    std::string const& out_pref
) {
    // Read metadata table
    LOG_INFO << "Reading meta data";
    auto reader = CsvReader();
    reader.separator_chars( "\t" );
    auto meta_table = reader.from_file( meta_file );

    // Meta data needs to have one more row (header line) than the matrices.
    if( meta_table.size() != edge_w_rows + 1 ) {
        throw std::runtime_error( "Meta table has not length of matrix rows." );
    }

    // Get column names of the meta data
    auto meta_names = std::vector<std::string>();
    for( size_t i = 1; i < meta_table[0].size(); ++i ) {
        meta_names.push_back( meta_table[0][i] );
    }

    // Make a table of meta data for all columns, and fill it with the data, sorted
    // by the sample names, so that each row corresponds to the other matrix rows.
    auto meta_data = utils::Matrix<double>( edge_w_rows, meta_names.size() );
    if( meta_data.rows() != order.size() ) {
        throw std::runtime_error( "Impossibru." );
    }
    for( size_t r = 0; r < meta_data.rows(); ++r ) {

        // For the given matrix row, find the line in the meta table that corresponds to this sample.
        auto const& sample_name = order[r];
        size_t meta_index = meta_table.size();
        for( size_t j = 1; j < meta_table.size(); ++j ) {
            if( meta_table[j][0] == sample_name ) {
                meta_index = j;
                break;
            }
        }
        if( meta_index == meta_table.size() ) {
            throw std::runtime_error( "Cannot find meta data table index for sample " + sample_name );
        }
        if(  meta_table[meta_index].size() != meta_data.cols() + 1 ) {
            throw std::runtime_error( "Meta data table has wrong line length at " + std::to_string(r) );
        }

        // Now, fill the data columns.
        for( size_t c = 0; c < meta_data.cols(); ++c ) {
            // LOG_DBG1 << meta_table[meta_index][c+1];
            meta_data( r, c ) = stod( meta_table[meta_index][c+1] );
        }
    }

    // Write meta matrix.
    LOG_INFO << "Writing meta data matrix";
    utils::file_write( utils::to_string( meta_data ), out_pref + "metadata.mat" );

    return { meta_names, meta_data };
}

// =================================================================================================
//     Process
// =================================================================================================

void process(
    tree::DefaultTree const& tree,
    Matrix<double> const& edge_w,
    Matrix<double> const& edge_i,
    std::vector<std::string> const& order,
    MetaData const& meta,
    std::string const& out_pref
) {


    // User info
    LOG_INFO << "tree: " << tree.edge_count() << " edges";
    LOG_INFO << "edge_w: " << edge_w.cols() << " cols and " << edge_w.rows() << " rows";
    LOG_INFO << "edge_i: " << edge_i.cols() << " cols and " << edge_i.rows() << " rows";
    LOG_INFO << "order: " << order.size() << " rows";
    LOG_INFO << "meta_data: " << meta.data.cols() << " cols and " << meta.data.rows() << " rows";

    // currently ignored: test to visualize total sum of weights on each edge
    // LOG_INFO << "Calculating tree edge weights";
    // auto const edge_w_vec = matrix_col_sums( edge_w );
    // auto const edge_side_mat = tree::edge_sides( tree );
    // auto const edge_i_vec = matrix_multiplication( edge_side_mat, edge_w_vec );
    //
    // auto const edge_w_colors = vec_to_col( edge_w_vec );
    // auto const edge_i_colors = vec_to_col( edge_i_vec );
    auto const edge_w_colors = std::vector<utils::Color>();
    auto const edge_i_colors = std::vector<utils::Color>();

    LOG_INFO << "Calculating meta";

    // For each meta datum, make a pcc and srcc coloured tree with both edge weight and imbalance.
    for( size_t c = 0; c < meta.data.cols(); ++c ) {

        auto name = utils::sanitize_filname( meta.names[c] );
        LOG_INFO << "Meta data: " << name;

        // How many row were filtered due to bad metadata
        size_t filtered = 0;

        auto pcc_w = std::vector<double>( edge_w.cols(), 0.0 );
        auto pcc_i = std::vector<double>( edge_w.cols(), 0.0 );
        auto srcc_w = std::vector<double>( edge_w.cols(), 0.0 );
        auto srcc_i = std::vector<double>( edge_w.cols(), 0.0 );
        for( size_t j = 0; j < edge_w.cols(); ++j ) {

            // Get only rows that have good metadata!
            std::vector<double> edge_w_col;
            std::vector<double> edge_i_col;
            std::vector<double> meta_col;
            for( size_t i = 0; i < meta.data.rows(); ++i ) {
                if( meta.data( i, c ) != 99999 ) {
                    edge_w_col.push_back( edge_w( i,j ) );
                    edge_i_col.push_back( edge_i( i,j ) );
                    meta_col.push_back( meta.data( i,c ) );
                }
            }

            if( filtered == 0 ) {
                filtered = meta.data.rows() - meta_col.size();
            } else if( filtered != meta.data.rows() - meta_col.size() ) {
                LOG_ERR << filtered << " = filtered != meta.data.rows() - meta_col.size() = " << (meta.data.rows() - meta_col.size());
            }

            pcc_w[j] = pearson_correlation_coefficient( edge_w_col, meta_col );
            pcc_i[j] = pearson_correlation_coefficient( edge_i_col, meta_col );
            srcc_w[j] = spearmans_rank_correlation_coefficient( edge_w_col, meta_col );
            srcc_i[j] = spearmans_rank_correlation_coefficient( edge_i_col, meta_col );

            // Simple version that uses the whole matrix.
            // srcc_w[j] = matrix_col_spearmans_rank_correlation_coefficient( edge_w, j, meta.data, c );
            // srcc_i[j] = matrix_col_spearmans_rank_correlation_coefficient( edge_i, j, meta.data, c );
        }

        LOG_INFO << "Filtered out " << filtered << " metadata rows";

        auto pcc_col_w = cc_to_colors( pcc_w );
        auto pcc_col_i = cc_to_colors( pcc_i );
        auto srcc_col_w = cc_to_colors( srcc_w );
        auto srcc_col_i = cc_to_colors( srcc_i );

        write_color_tree_to_nexus( tree, pcc_col_w, out_pref + name + "_pcc_weights.nexus" );
        write_color_tree_to_svg(   tree, pcc_col_w, edge_w_colors, out_pref + name + "_pcc_weights.svg" );
        write_color_tree_to_nexus( tree, pcc_col_i, out_pref + name + "_pcc_imbalance.nexus" );
        write_color_tree_to_svg(   tree, pcc_col_i, edge_i_colors,  out_pref + name + "_pcc_imbalance.svg" );

        write_color_tree_to_nexus( tree, srcc_col_w, out_pref + name + "_srcc_weights.nexus" );
        write_color_tree_to_svg(   tree, srcc_col_w, edge_w_colors, out_pref + name + "_srcc_weights.svg" );
        write_color_tree_to_nexus( tree, srcc_col_i, out_pref + name + "_srcc_imbalance.nexus" );
        write_color_tree_to_svg(   tree, srcc_col_i, edge_i_colors,  out_pref + name + "_srcc_imbalance.svg" );
    }

    LOG_INFO << "Finished";
}

// =================================================================================================
//     Process 6
// =================================================================================================

/**
 * @brief Version of the programm that takes 6 command line args.
 */
void process_6( int argc, char** argv )
{
    // Check if the command line contains the right number of arguments.
    if (argc != 7) {
        throw std::runtime_error(
            "Need to provide 6 arguments.\n"
        );
    }

    // Newick tree.
    auto const tree_file   = std::string( argv[1] );

    // Weight of placement mass per sample per edge.
    // Rows correspond to samples, cols to edges of the tree.
    auto const edge_w_file = std::string( argv[2] );

    // Imbalance of placement mass per sample per edge.
    // Rows correspond to samples, cols to edges of the tree.
    // Edges can either be all, or just the non-tip ones.
    auto const edge_i_file = std::string( argv[3] );

    // Simple list file, where each row corresponds to the sample name that was used for the
    // previous to matrix files.
    auto const order_file  = std::string( argv[4] );

    // Meta data table file, tab separted.
    // First row is the header, which contains the names of the meta data columns.
    // First column is the sample name. THose do not need to be in the same order as the
    // matrix files above. In order to sort them so that the lines are fitting, we use the order
    // file.
    auto const meta_file   = std::string( argv[5] );

    // Output prefix, that is, path plus file name prefix.
    auto const out_pref    = std::string( argv[6] );


    LOG_INFO << "Reading stuff";

    // Read tree
    LOG_INFO << "Reading tree";
    auto tree = tree::DefaultTreeNewickReader().from_file( tree_file );

    // Read edge matrices
    LOG_INFO << "Reading matrices";
    auto edge_w = read_double_matrix( edge_w_file );
    auto edge_i = read_double_matrix( edge_i_file );

    // Check and correct matrix sizes
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

    // Last checks of matrix sizes
    if( edge_i.cols() != tree.edge_count() ) {
        throw std::runtime_error( "Edge imbalance does not have corrent length" );
    }
    if( edge_w.rows() != edge_i.rows() ) {
        throw std::runtime_error( "Edge matrices do not have correct size" );
    }

    // Read order list. This tells us which matrix rows are coming from which sample.
    LOG_INFO << "Reading order";
    auto reader = CsvReader();
    auto order_table  = reader.from_file( order_file );
    std::vector<std::string> order;
    for( auto const& line : order_table ) {
        if( line.size() != 1 ) {
            throw std::runtime_error( "Order list entry with more than 1 column." );
        }
        if( ! utils::ends_with( line[0], ".bplace" ) ) {
            throw std::runtime_error( "Order list entry without bplace" );
        }
        order.push_back( utils::replace_all( line[0], ".bplace", "" ) );
    }
    if( order.size() != edge_w.rows() ) {
        throw std::runtime_error( "Order has not length of matrix rows." );
    }

    auto const meta_data = meta_data_prep(  meta_file, edge_w.rows(), order, out_pref );

    process( tree, edge_w, edge_i, order, meta_data, out_pref );
}

// =================================================================================================
//     Process 3
// =================================================================================================

/**
 * @brief Versoin that takes two args
 */
void process_3( int argc, char** argv )
{
    // Check if the command line contains the right number of arguments.
    if (argc != 4 && argc != 5) {
        throw std::runtime_error(
            "Need to provide three or four arguments.\n"
        );
    }

    using namespace placement;

    // In out dirs.
    auto epadir   = utils::dir_normalize_path( std::string( argv[1] ));
    auto const meta_file   = std::string( argv[2] );
    auto out_pref = utils::dir_normalize_path( std::string( argv[3] ));

    std::string sequence_name_file;
    if( argc == 5 ) {
        sequence_name_file = std::string( argv[4] );
    }

    // -------------------------------------------------------------------------
    //     Find and read sample [jb]place files.
    // -------------------------------------------------------------------------

    SampleSet sample_set;
    std::vector<std::string> order;

    // Get all bplace files, either in the epa dir, or in its subdirs
    auto bplace_filenames = utils::dir_list_files( epadir, false, ".*\\.bplace" );
    std::sort( bplace_filenames.begin(), bplace_filenames.end() );

    if( bplace_filenames.size() > 0 ) {
        order.clear();
        for( auto const& e : bplace_filenames ) {
            order.push_back( utils::replace_all( e, ".bplace", "" ) );
        }

        // Output list of samples in the order that we use them for the matrix
        // LOG_INFO << "Writing bplace sample order list";
        // std::ofstream bplace_files_of( outdir + "bplace_order_mats.list" );
        // for( auto const& elem : bplace_filenames ) {
        //     bplace_files_of << elem << "\n";
        // }
        // bplace_files_of.close();
        // LOG_INFO << "done";
        // LOG_INFO;

        SampleSerializer bplace_loader;

        // Process all jplace files.
        LOG_INFO << "reading " << bplace_filenames.size() << " bplace sample files";
        // for( auto const& bplace_filename : bplace_filenames ) {
        //     sample_set.add( bplace_loader.load( epadir + bplace_filename ), bplace_filename );
        // }

        auto tmp = std::vector<Sample>( bplace_filenames.size() );

        // Parallel parsing.
        #pragma omp parallel for
        for( size_t i = 0; i < bplace_filenames.size(); ++i ) {
            tmp[ i ] = bplace_loader.load( epadir + bplace_filenames[i] );
        }

        // Move to target SampleSet.
        for( size_t i = 0; i < bplace_filenames.size(); ++i ) {
            sample_set.add( std::move( tmp[i] ), bplace_filenames[i] );
        }


        // Final output for jplace reading
        assert( sample_set.size() == bplace_filenames.size() );
        LOG_INFO << "finished reading sample bplace files";
        LOG_INFO;
    }

    // Get all bplace files, either in the epa dir, or in its subdirs
    auto jplace_filenames = utils::dir_list_files( epadir, false, ".*\\.jplace" );
    std::sort( jplace_filenames.begin(), jplace_filenames.end() );

    if( jplace_filenames.size() > 0 ) {
        order.clear();
        for( auto& e : jplace_filenames ) {
            order.push_back( utils::replace_all( e, ".jplace", "" ) );
            e = epadir + e;
        }

        LOG_INFO << "reading " << jplace_filenames.size() << " jplace sample files";
        auto jplace_reader = JplaceReader();
        sample_set = jplace_reader.from_files( jplace_filenames );
        assert( sample_set.size() == jplace_filenames.size() );
        LOG_INFO << "finished reading sample jplace files";
        LOG_INFO;
    }

    // -------------------------------------------------------------------------
    //     Get matrices etc
    // -------------------------------------------------------------------------

    // Average tree.
    tree::TreeSet avg_tree_set;
    for( auto const& smp : sample_set ) {
        avg_tree_set.add( "", smp.sample.tree() );
    }
    auto const tree = tree::average_branch_length_tree( avg_tree_set );
    // tree::DefaultTreeNewickWriter().to_file( avg_tree, outdir + "avg_tree.newick" );
    avg_tree_set.clear();

    // Rename tips of the tree to more readable proper names
    if( ! sequence_name_file.empty() ) {
        auto reader = CsvReader();
        reader.separator_chars( "\t" );
        auto seq_table = reader.from_file( sequence_name_file );
        LOG_INFO << "found seq table with " << seq_table.size() << " entries, tree has " << tree::leaf_node_count(tree) << " nodes";

        std::unordered_map<std::string, std::string> seq_map;
        for( auto const& line : seq_table ) {
            if( line.size() != 2  ) {
                throw std::runtime_error( "seq table line length error" );
            }
            if( seq_map.count(line[0]) > 0 ) {
                throw std::runtime_error( "seq table insertion error" );
            }

            seq_map[line[0]] = line[1];
        }

        for( auto& node : tree.nodes() ) {
            auto& nn = node->data<tree::DefaultNodeData>().name;
            if( nn.empty() ) {
                continue;
            }
            if( seq_map.count(nn) == 0 ) {
                LOG_WARN << "no seq name for " << nn;
            } else {
                nn = seq_map[nn];
            }
        }
    }

    auto const edge_i = epca_imbalance_matrix( sample_set, true );
    auto const edge_w = placement_weight_per_edge( sample_set );

    // utils::file_write( utils::to_string( edge_i ), outdir + "imbalance.mat" );
    // utils::file_write( utils::to_string( edge_w ), outdir + "edge_weights.mat" );

    auto const meta_data = meta_data_prep(  meta_file, edge_w.rows(), order, out_pref );

    process( tree, edge_w, edge_i, order, meta_data, out_pref );
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

    #ifdef GENESIS_PTHREADS
    utils::Options::get().number_of_threads( std::thread::hardware_concurrency() );
    #endif

    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    if (argc == 7) {
        process_6( argc, argv );
    } else if( argc == 4 || argc == 5 ) {
        process_3( argc, argv );
    } else {
        LOG_ERR << "wrong num args";
    }

    return 0;
}
