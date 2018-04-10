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
using namespace genesis::tree;
using namespace genesis::utils;

// =================================================================================================
//      Counts To Colors
// =================================================================================================

/**
 * @brief Given a vector of doubles, return a vector of Colors representing the distribution
 * of the double values.
 *
 * The resulting vector contains a color indicating how high the value of each input double is,
 * compared to the other values. This means: First, we find the highest value in the vector. Then,
 * for all values in the vector, we calculate the log-scaled relative value compared to the highest
 * value. This relative value is then turned into a color gradient.
 *
 * This way, the resulting vector has ligh blue colors for lower numbers, purple for medium numbers,
 * and darker colors up to black for higher numbers (when using the given color gradient scheme).
 */
// std::vector<utils::Color> counts_to_colors(
//     std::vector<double> const& count_vector
// ) {
//     // Create a color gradient in "blue pink black".
//     auto gradient   = std::map<double, utils::Color>();
//     gradient[ 0.0 ] = utils::color_from_hex("#81bfff");
//     gradient[ 0.5 ] = utils::color_from_hex("#c040be");
//     gradient[ 1.0 ] = utils::color_from_hex("#000000");
//     auto base_color = utils::color_from_hex("#81bfff");
//
//     // Create the result vector.
//     std::vector<utils::Color> color_vector( count_vector.size(), base_color );
//
//     // Find the highest value in the input vector.
//     auto mass_max = *std::max_element (count_vector.begin(), count_vector.end());
//
//     // Output for building the color gradient: show the color values.
//     LOG_INFO << "In order to build a color gradient, use these color stops:\n"
//              << "    0.0 -> #81bfff (light blue, no placement mass)\n"
//              << "    0.5 -> #c040be (purple, medium placement mass)\n"
//              << "    1.0 -> #000000 (black, maximum placement mass)";
//     LOG_INFO << "The edge with the maximum placement mass has a mass of " << mass_max;
//
//     // More output for the color gradient: show the label positions.
//     std::string labels = "The following list shows position on the color gradient ";
//     labels += "with their correspinding placement mass. Use this to label your color gradient.\n";
//     size_t pow_i = 0;
//     double pow_v = 0.0;
//     do {
//         pow_v = log( pow( 10, pow_i )) / log( mass_max );
//         if( pow_v <= 1.0 ) {
//             labels += "    " + utils::to_string_precise( pow_v, 6 ) + " -> ";
//             labels += utils::to_string( pow( 10, pow_i )) + "\n";
//         }
//         ++pow_i;
//     } while( pow_v <= 1.0 );
//     labels += "    1.000000 -> " + utils::to_string( mass_max ) + "\n";
//     LOG_INFO << labels;
//
//     // Calculate the resulting colors.
//     for( size_t i = 0; i < count_vector.size(); ++i ) {
//         if( count_vector[i] > 0.0 ) {
//             double log_val = log( static_cast<double>(count_vector[i]) ) / log( mass_max );
//             color_vector[i] = utils::gradient( gradient, log_val );
//         }
//     }
//
//     return color_vector;
// }

// std::vector<utils::Color> counts_to_colors(
//     std::vector<double> const& count_vector,
//     bool use_log = true
// ) {
//     auto map = ColorMap( color_list_bupubk() );
//     if( use_log ) {
//
//         auto norm = ColorNormalizationLogarithmic( count_vector );
//         return map( norm, count_vector );
//
//     } else {
//
//         auto norm = ColorNormalization( count_vector );
//         return map( norm, count_vector );
//     }


    // // pal.reverse( true );
    // // auto base_color = utils::color_from_hex("#FF00FF");
    // auto base_color = map( 0.0 );
    //
    // // Create the result vector.
    // std::vector<utils::Color> color_vector( count_vector.size(), base_color );
    //
    // // Find the highest value in the input vector.
    // auto mass_max = *std::max_element (count_vector.begin(), count_vector.end());
    //
    // if( use_log ) {
    //     // More output for the color gradient: show the label positions.
    //     std::string labels = "The following list shows the log positions on the color gradient ";
    //     labels += "with their correspinding placement mass. Use this to label your color gradient.\n";
    //     double scaler = 1000.0;
    //     double pow_i  = -1.0;
    //     double pow_v  = 0.0;
    //     do {
    //         pow_v = log( scaler * pow( 10, pow_i )) / log( scaler * mass_max );
    //         if( pow_v <= 1.0 ) {
    //             labels += "    " + utils::to_string_precise( pow_v, 6 ) + " -> ";
    //             labels += utils::to_string( scaler * pow( 10, pow_i ) ) + " with ";
    //             labels += utils::to_string(pow_i) + "\n";
    //         }
    //         pow_i += 0.25;
    //     } while( pow_v <= 1.0 );
    //     labels += "    1.000000 -> " + utils::to_string( scaler * mass_max ) + "\n";
    //     // LOG_INFO << labels;
    //
    //     // Calculate the resulting colors.
    //     for( size_t i = 0; i < count_vector.size(); ++i ) {
    //         if( count_vector[i] > 0.0 ) {
    //             double log_val = log( scaler * static_cast<double>(count_vector[i]) ) / log( scaler * mass_max );
    //             if( log_val < 0.0 ) {
    //                 log_val = 0.0;
    //             }
    //             color_vector[i] = map( log_val );
    //         }
    //     }
    //
    // } else {
    //     // More output for the color gradient: show the label positions.
    //     std::string labels = "The following list shows the linear positions on the color gradient ";
    //     labels += "with their correspinding placement mass. Use this to label your color gradient.\n";
    //
    //     double oom = floor( log10( mass_max ));
    //     double lin_inc = pow( 10, oom );
    //     double lin_cur = 0.0;
    //
    //     do {
    //         double lin_pos = lin_cur / mass_max;
    //         labels += "    " + utils::to_string_precise( lin_pos, 6 ) + " -> ";
    //         labels += utils::to_string( lin_cur ) + "\n";
    //         lin_cur += lin_inc;
    //     } while( lin_cur <= mass_max );
    //     labels += "    1.000000 -> " + utils::to_string( mass_max ) + "\n";
    //     // LOG_INFO << labels;
    //
    //     // Calculate the resulting colors.
    //     for( size_t i = 0; i < count_vector.size(); ++i ) {
    //         if( count_vector[i] > 0.0 ) {
    //             // double log_val = log( static_cast<double>(count_vector[i]) ) / log( mass_max );
    //             double log_val = static_cast<double>(count_vector[i]) / mass_max;
    //             color_vector[i] = map( log_val );
    //         }
    //     }
    // }
    //
    // return color_vector;
// }

// std::vector<utils::Color> counts_to_colors_ranked(
//     std::vector<double> const& count_vector
// ) {
//     size_t const bins = 5;
//
//     auto const ranks = ranking_ordinal( count_vector );
//     auto const order = sort_indices( ranks.begin(), ranks.end() );
//
//     auto pal = ColorPalette( color_list_bupubk() );
//     auto base_color = pal.sequential_color( 0.0 );
//     pal.max( bins - 1 );
//     std::vector<utils::Color> color_vector( count_vector.size(), base_color );
//
//     double const frac = static_cast<double>( bins ) / static_cast<double>( order.size());
//     for( size_t i = 0; i < order.size(); ++i ) {
//         auto const bin = static_cast<size_t>( std::floor( static_cast<double>(i) * frac ));
//         assert( bin < bins );
//
//         color_vector[ order[i] ] = pal.sequential_color( bin );
//     }
//
//     return color_vector;
// }

// =================================================================================================
//     Write Squash Tree To SVG
// =================================================================================================

/**
 * @brief Write an SVG file containing a tree with colored branches.
 */
void write_squash_tree_to_svg(
    placement::PlacementTree const&  tree,
    std::vector<utils::Color> const& colors_per_node,
    std::string const&               svg_filename
) {
    auto copy = tree;
    ladderize(copy);
    // LOG_DBG << "Ladderize tree: " << ( validate_topology(copy) ? "valid!" : "INVALID!" );

    // Make a layout tree.
    // auto layout = tree::RectangularLayout( tree );
    auto layout = tree::RectangularLayout( copy, tree::LayoutType::kPhylogram );
    // layout.scaler_x( 300 );

    // Set colourful node shapes.
    if( colors_per_node.size() != tree.node_count() ) {
        throw std::runtime_error( "Node colors vec of wrong size" );
    }

    std::vector<utils::SvgGroup> node_shapes;
    node_shapes.resize( tree.node_count() );
    for( size_t i = 0; i < tree.node_count(); ++i ) {
        if( colors_per_node[i] == Color(1,0,1) ) {
            continue;
        }
        // Add a shape according to the sample color.
        node_shapes[i].add( utils::SvgCircle(
            utils::SvgPoint( 0, 0 ),
            6.0,
            utils::SvgStroke( SvgStroke::Type::kNone ),
            utils::SvgFill( colors_per_node[i] )
        ));
    }
    layout.set_node_shapes( node_shapes );

    // Write to file
    std::ostringstream out;
    layout.to_svg_document().write( out );
    utils::file_write( out.str(), svg_filename + ".svg" );
}

void make_squash_tree(
    std::vector<size_t> const& ref_assignments,
    DefaultTree const& sc_tree,
    std::vector<std::string> const& sample_names,
    std::string const& outfile
)
{
    // auto const& ref_assignments = massignments;
    // auto const& ref_assignments = iassignments;

    if( ref_assignments.size() != leaf_node_count(sc_tree) ) {
        LOG_WARN << ref_assignments.size() << " == ref_assignments.size() != leaf_node_count(sc_tree) == " << leaf_node_count(sc_tree);
        return;
    }

    if( ref_assignments.size() != sample_names.size() ) {
        LOG_WARN << ref_assignments.size() << " == ref_assignments.size() != sample_names.size() == " << sample_names.size();
        return;
    }

    auto colors_per_node = std::vector<utils::Color>( sc_tree.node_count(), Color(1,0,1) );
    auto const& node_pal = color_list_set1();
    size_t tmp_sum = 0;
    for( size_t i = 0; i < sc_tree.node_count(); ++i ) {
        auto const& name = sc_tree.node_at(i).data<DefaultNodeData>().name;
        if( name.empty() ) {
            continue;
        }

        // Fake test to see if it is a num
        bool is_num = false;
        try{
            // LOG_DBG << stoi(name);
            tmp_sum += stoi(name);
            is_num = true;
        } catch(...){}
        if( is_num ) {
            continue;
        }

        // find index in sample name list
        auto match = std::find(sample_names.begin(), sample_names.end(), name);
        size_t index;
        if(match != sample_names.end()) {
            index = match - sample_names.begin();
        } else {
            LOG_WARN << "cannot find sample name " << name;
            continue;
        }

        colors_per_node[ i ] = node_pal[ ref_assignments[index] ];
    }

    //  do not write node names
    auto cpy = sc_tree;
    for( auto& node : cpy.nodes() ) {
        node->data<DefaultNodeData>().name = "";
    }

    // write_color_tree_to_nexus( cpy, {}, outdir + "/cluster_kmeans.nexus" );
    write_squash_tree_to_svg( cpy, colors_per_node, outfile );
}

// =================================================================================================
//      Main
// =================================================================================================

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
    if (argc != 5) {
        throw std::runtime_error(
            "Need to provide four arguments.\n"
        );
    }

    // In out dirs.
    auto const threads = std::stoi( argv[1] );
    auto const k       = std::stoi( argv[2] );
    auto epadir = utils::dir_normalize_path( std::string( argv[3] ));
    auto outdir = utils::dir_normalize_path( std::string( argv[4] ));
    utils::dir_create(outdir);

    utils::Options::get().command_line( argc, argv );
    utils::Options::get().number_of_threads( threads );
    utils::Options::get().random_seed( 1455697401 );
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    // -------------------------------------------------------------------------
    //     Find and read sample jplace files.
    // -------------------------------------------------------------------------

    // Get all jplace files, either in the epa dir, or in its subdirs
    auto jplace_filenames = utils::dir_list_files( epadir, ".*\\.jplace" );
    std::sort( jplace_filenames.begin(), jplace_filenames.end() );

    // Output list of samples in the order that we use them for the matrix
    // LOG_INFO << "Writing jplace sample order list";
    // std::ofstream jplace_files_of( outdir + "jplace_order.list" );
    // for( auto const& elem : jplace_filenames ) {
    //     jplace_files_of << elem << "\n";
    // }
    // jplace_files_of.close();
    // LOG_INFO << "done";
    // LOG_INFO;

    SampleSet sset;
    JplaceReader jplace_reader;

    // Process all jplace files.
    LOG_INFO << "reading " << jplace_filenames.size() << " jplace sample files";
    sset = jplace_reader.from_files( jplace_filenames );

    // size_t file_count = 0;
    // for( auto const& jplace_filename : jplace_filenames ) {
    //
    //     // Progress output.
    //     // LOG_INFO << "at jplace file " << file_count << " (" << jplace_filename << ")";
    //     ++file_count;
    //
    //     sset.add( jplace_reader.from_file( epadir + jplace_filename ), jplace_filename );
    // }

    // Final output for jplace reading
    assert( sset.size() == jplace_filenames.size() );
    LOG_INFO << "finished reading sample jplace files";
    LOG_INFO;

    // -------------------------------------------------------------------------
    //     Tree Kmeans
    // -------------------------------------------------------------------------

    // Options::get().random_seed( 1623390399 );

    std::vector<std::string> sample_names;
    for( auto const& smp : sset ) {
        sample_names.push_back( smp.name );
        // LOG_DBG << smp.name;
    }

    LOG_INFO << "Converting Trees";
    auto mass_trees = convert_sample_set_to_mass_trees( sset );
    auto const avg_tree =average_branch_length_tree( sset );
    // sset.clear();

    auto mt2 = mass_trees;

    // size_t bins = 1;
    // LOG_INFO << "Binning";
    // for( auto& tree : mt2.first ) {
    //     mass_tree_binify_masses( tree, bins );
    // }

    LOG_INFO << "Kmeans started";
    auto mkmeans = tree::MassTreeKmeans();
    // mkmeans.accumulate_centroid_masses(bins);
    mkmeans.run( mt2.first, k );
    LOG_INFO << "Kmeans finished";
    LOG_INFO << "assignments().size() " << mkmeans.assignments().size();
    LOG_INFO << "centroids().size()   " << mkmeans.centroids().size();

    // Write assignments
    LOG_INFO << "Write assignments";
    std::ofstream file_mkmeans_ass;
    file_output_stream( outdir + "emd_assignments.csv",  file_mkmeans_ass );
    auto const& massignments = mkmeans.assignments();
    auto pqry_cnts = std::vector<size_t>( mkmeans.centroids().size(), 0 );
    for( size_t i = 0; i < massignments.size(); ++i ) {
        file_mkmeans_ass << file_filename( sample_names[i] );
        file_mkmeans_ass << "\t" << massignments[i];
        file_mkmeans_ass << "\n";

        pqry_cnts[ massignments[i] ] += sset[i].sample.size();
    }
    file_mkmeans_ass.close();

    // Write centroids
    auto param = tree::LayoutParameters();
    param.stroke.width = 6.0;

    LOG_INFO << "Write centroids";
    std::ofstream file_mkmeans_cent;
    file_output_stream( outdir + "emd_centroids.csv",  file_mkmeans_cent );
    auto const& mcentroids = mkmeans.centroids();
    for( size_t i = 0; i < mcentroids.size(); ++i ) {
        auto masses = mass_tree_mass_per_edge( mcentroids[i] );

        // Rescale to what the original masses per sample were.
        for( auto& mass : masses ) {
            mass *= pqry_cnts[i];
        }

        // auto cent = mcentroids[i];
        // mass_tree_scale_masses( cent, list[i].size() );
        // auto colors_per_branch = counts_to_colors( mass_tree_mass_per_edge( cent ));

        // auto colors_per_branch = counts_to_colors( masses, false );
        // auto colors_per_branch = counts_to_colors_ranked( masses );
        // write_color_tree_to_nexus( avg_tree, colors_per_branch, outdir + "tree_emd_" + std::to_string(i) + ".nexus" );
        // write_color_tree_to_svg( avg_tree, colors_per_branch, outdir + "tree_emd_" + std::to_string(i) );

        auto map = ColorMap( color_list_bupubk() );
        auto norm_seq = ColorNormalizationLinear( masses );
        auto col_seq = map( norm_seq, masses );

        tree::write_color_tree_to_nexus_file( avg_tree, col_seq, outdir + "tree_emd_" + std::to_string(i) + "_lin.nexus" );
        tree::write_color_tree_to_svg_file( avg_tree, param, col_seq, map, norm_seq, outdir + "tree_emd_" + std::to_string(i) + "_lin.svg" );

        auto norm_log = ColorNormalizationLogarithmic( masses );
        norm_log.min_value( 1.0 );
        map.clip_under( true );

        auto col_log = map( norm_log, masses );

        tree::write_color_tree_to_nexus_file( avg_tree, col_log, outdir + "tree_emd_" + std::to_string(i) + "_log.nexus" );
        tree::write_color_tree_to_svg_file( avg_tree, param, col_log, map, norm_log, outdir + "tree_emd_" + std::to_string(i) + "_log.svg" );

        // file_mkmeans_cent << i;
        for( auto const& mass : masses ) {
            if( &mass != &masses[0] ) {
                file_mkmeans_cent << "\t";
            }
            file_mkmeans_cent << mass;
        }
        file_mkmeans_cent << "\n";
    }
    file_mkmeans_cent.close();

    // -------------------------------------------------------------------------
    //     Edge Imbalance Kmeans
    // -------------------------------------------------------------------------

    // Prepare
    LOG_INFO << "Calculating Matrix";
    auto edge_imb_mat = epca_imbalance_matrix( sset, true );
    auto const columns = epca_filter_constant_columns( edge_imb_mat, 0.001 );

    auto edge_imb_vec = std::vector<std::vector<double>>();
    edge_imb_vec.resize( sset.size() );

    for( size_t i = 0; i < edge_imb_mat.rows(); ++i ) {
        edge_imb_vec[i].resize( edge_imb_mat.cols() );
        for( size_t j = 0; j < edge_imb_mat.cols(); ++j ) {
            edge_imb_vec[i][j] = edge_imb_mat(i, j);
        }
    }

    // Write imbalance matrix
    std::ofstream file_imbmat;
    file_output_stream( outdir + "/imbalance_matrix.csv",  file_imbmat );
    file_imbmat << edge_imb_mat;
    file_imbmat.close();

    // Kmeans
    LOG_INFO << "Kmeans started";
    auto ikmeans = EuclideanKmeans( edge_imb_mat.cols() );
    ikmeans.run( edge_imb_vec, k );
    LOG_INFO << "Kmeans finished";

    // Write assignments
    LOG_INFO << "Write assignments";
    std::ofstream file_ikmeans_ass;
    file_output_stream( outdir + "/imbalance_assignments.csv",  file_ikmeans_ass );
    auto const& iassignments = ikmeans.assignments();
    for( size_t i = 0; i < iassignments.size(); ++i ) {
        file_ikmeans_ass << file_filename( sset[i].name );
        file_ikmeans_ass << "\t" << iassignments[i];
        file_ikmeans_ass << "\n";
    }
    file_ikmeans_ass.close();

    // Colors: positive and negative gradient, and base color
    auto p_gradient   = std::map<double, utils::Color>();
    p_gradient[ 0.0 ] = utils::color_from_hex("#000000");
    p_gradient[ 1.0 ] = utils::color_from_hex("#0000ff");
    auto n_gradient   = std::map<double, utils::Color>();
    n_gradient[ 0.0 ] = utils::color_from_hex("#000000");
    n_gradient[ 1.0 ] = utils::color_from_hex("#ff0000");
    auto base_color = utils::color_from_hex("#008800");

    // Write centroids
    LOG_INFO << "Write centroids";
    std::ofstream file_ikmeans_cent;
    file_output_stream( outdir + "/imbalance_centroids.csv",  file_ikmeans_cent );
    auto const& icentroids = ikmeans.centroids();
    for( size_t i = 0; i < icentroids.size(); ++i ) {
        auto const& cent = icentroids[i];

        auto col_vec = std::vector<utils::Color>( sset[0].sample.tree().edge_count(), base_color );
        assert( columns.size() == cent.size() );
        for( size_t j = 0; j < columns.size(); ++j ) {
            assert( -1.0 <= cent[ j ] && cent[ j ] <= 1.0 );
            if( cent[ j ] >= 0.0 ) {
                col_vec[ columns[j] ] = utils::gradient( p_gradient, cent[ j ] );
            } else {
                col_vec[ columns[j] ] = utils::gradient( n_gradient, - cent[ j ] );
            }
        }
        tree::write_color_tree_to_nexus_file( sset[0].sample.tree(), col_vec, outdir + "tree_imb_" + std::to_string(i) + ".nexus" );
        tree::write_color_tree_to_svg_file( sset[0].sample.tree(), param, col_vec, outdir + "tree_imb_" + std::to_string(i) + "_aweful.svg" );

        for( auto const& v : cent ) {
            if( &v != &cent[0] ) {
                file_ikmeans_cent << "\t";
            }
            file_ikmeans_cent << v;
        }
        file_ikmeans_cent << "\n";

        // Get the mass of all trees in this cluster.
        tree::MassTree comb_mass;
        for( size_t a = 0; a < ikmeans.assignments().size(); ++a ) {
            if( ikmeans.assignments()[a] != i ) {
                continue;
            }
            if( comb_mass.empty() ) {
                comb_mass = mass_trees.first[a];
            } else {
                mass_tree_merge_trees_inplace( comb_mass, mass_trees.first[a] );
            }
        }
        auto masses = mass_tree_mass_per_edge( comb_mass );

        // Rescale to what the original masses per sample were.
        for( auto& mass : masses ) {
            mass *= pqry_cnts[i];
        }

        auto map = ColorMap( color_list_bupubk() );
        auto norm_seq = ColorNormalizationLinear( masses );
        auto col_seq = map( norm_seq, masses );

        tree::write_color_tree_to_nexus_file( avg_tree, col_seq, outdir + "tree_imb_" + std::to_string(i) + "_lin.nexus" );
        tree::write_color_tree_to_svg_file( avg_tree, param, col_seq, map, norm_seq, outdir + "tree_imb_" + std::to_string(i) + "_lin.svg" );

        auto norm_log = ColorNormalizationLogarithmic( masses );
        norm_log.min_value( 1.0 );
        map.clip_under( true );

        auto col_log = map( norm_log, masses );

        tree::write_color_tree_to_nexus_file( avg_tree, col_log, outdir + "tree_imb_" + std::to_string(i) + "_log.nexus" );
        tree::write_color_tree_to_svg_file( avg_tree, param, col_log, map, norm_log, outdir + "tree_imb_" + std::to_string(i) + "_log.svg" );
    }
    file_ikmeans_cent.close();

    // -------------------------------------------------------------------------
    //     Squash
    // -------------------------------------------------------------------------

    LOG_INFO << "Starting squash clustering";
    auto sc = tree::SquashClustering();
    sc.run( std::move( mass_trees.first ) );
    LOG_INFO << "Finished squash clustering";

    // LOG_DBG << "Clusters";
    // for( size_t i = 0; i < sc.clusters.size(); ++i ) {
    //     LOG_DBG1 << i << ": " << ( sc.clusters[i].active ? "act " : "dea " )
    //              << sc.clusters[i].count << " "
    //              << sc.clusters[i].distances.size() << " "
    //              << sc.clusters[i].tree.node_count();
    // }
    //
    // LOG_DBG << "Mergers";
    // for( size_t i = 0; i < sc.mergers.size(); ++i ) {
    //     LOG_DBG1 << i << ": "
    //              << sc.mergers[i].index_a << " " << sc.mergers[i].distance_a << " | "
    //              << sc.mergers[i].index_b << " " << sc.mergers[i].distance_b;
    // }

    // Remove jplace part from file names
    // for( auto& fn : jplace_filenames ) {
    //     fn = utils::replace_all( fn, ".jplace", "" );
    // }

    LOG_INFO << "Writing cluster tree";
    std::ofstream file_clust_tree;
    utils::file_output_stream( outdir + "/cluster.newick",  file_clust_tree );
    auto const sct_str = sc.tree_string( sample_names );
    file_clust_tree << sct_str;

    // LOG_INFO << "Writing fat trees";
    // for( size_t i = 0; i < sc.clusters.size(); ++i ) {
    //     auto const& cc = sc.clusters[i];
    //
    //     auto const cv = mass_tree_mass_per_edge( cc.tree );
    //     auto const colors = counts_to_colors(cv, false);
    //
    //     write_color_tree_to_nexus( avg_tree, colors, outdir + "/tree_" + std::to_string(i) + ".nexus" );
    //     write_color_tree_to_svg( avg_tree, colors, outdir + "/tree_" + std::to_string(i) );
    // }

    // -------------------------------------------------------------------------
    //     Squash Tree with kmeans node colors
    // -------------------------------------------------------------------------

    auto sc_tree = DefaultTreeNewickReader().from_string( sct_str );

    // auto const& ref_assignments = massignments;
    // auto const& ref_assignments = iassignments;

    make_squash_tree( massignments, sc_tree, sample_names, outdir + "/cluster_kmeans_mass" );
    make_squash_tree( iassignments, sc_tree, sample_names, outdir + "/cluster_kmeans_imb" );

    LOG_INFO << "Finished ";
    return 0;
}
