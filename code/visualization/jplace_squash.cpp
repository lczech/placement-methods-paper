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

std::vector<utils::Color> counts_to_colors(
    std::vector<double> const& count_vector,
    bool use_log = true
) {
    // Create a color gradient in "blue pink black".
    auto gradient   = std::map<double, utils::Color>();
    gradient[ 0.0 ] = utils::color_from_hex("#81bfff");
    gradient[ 0.5 ] = utils::color_from_hex("#c040be");
    gradient[ 1.0 ] = utils::color_from_hex("#000000");
    auto base_color = utils::color_from_hex("#81bfff");

    // Create the result vector.
    std::vector<utils::Color> color_vector( count_vector.size(), base_color );

    // Find the highest value in the input vector.
    auto mass_max = *std::max_element (count_vector.begin(), count_vector.end());

    // Output for building the color gradient: show the color values.
    // LOG_INFO << "In order to build a color gradient, use these color stops:\n"
    //          << "    0.0 -> #81bfff (light blue, no placement mass)\n"
    //          << "    0.5 -> #c040be (purple, medium placement mass)\n"
    //          << "    1.0 -> #000000 (black, maximum placement mass)";
    // LOG_INFO << "The edge with the maximum placement mass has a mass of " << mass_max;

    if( use_log ) {
        // More output for the color gradient: show the label positions.
        std::string labels = "The following list shows the log positions on the color gradient ";
        labels += "with their correspinding placement mass. Use this to label your color gradient.\n";
        double scaler = 1000.0;
        double pow_i  = -1.0;
        double pow_v  = 0.0;
        do {
            pow_v = log( scaler * pow( 10, pow_i )) / log( scaler * mass_max );
            if( pow_v <= 1.0 ) {
                labels += "    " + utils::to_string_precise( pow_v, 6 ) + " -> ";
                labels += utils::to_string( scaler * pow( 10, pow_i ) ) + " with ";
                labels += utils::to_string(pow_i) + "\n";
            }
            pow_i += 0.25;
        } while( pow_v <= 1.0 );
        labels += "    1.000000 -> " + utils::to_string( scaler * mass_max ) + "\n";
        // LOG_INFO << labels;

        // Calculate the resulting colors.
        for( size_t i = 0; i < count_vector.size(); ++i ) {
            if( count_vector[i] > 0.0 ) {
                double log_val = log( scaler * static_cast<double>(count_vector[i]) ) / log( scaler * mass_max );
                color_vector[i] = utils::gradient( gradient, log_val );
            }
        }

    } else {
        // More output for the color gradient: show the label positions.
        std::string labels = "The following list shows the linear positions on the color gradient ";
        labels += "with their correspinding placement mass. Use this to label your color gradient.\n";

        double oom = floor( log10( mass_max ));
        double lin_inc = pow( 10, oom );
        double lin_cur = 0.0;

        do {
            double lin_pos = lin_cur / mass_max;
            labels += "    " + utils::to_string_precise( lin_pos, 6 ) + " -> ";
            labels += utils::to_string( lin_cur ) + "\n";
            lin_cur += lin_inc;
        } while( lin_cur <= mass_max );
        labels += "    1.000000 -> " + utils::to_string( mass_max ) + "\n";
        // LOG_INFO << labels;

        // Calculate the resulting colors.
        for( size_t i = 0; i < count_vector.size(); ++i ) {
            if( count_vector[i] > 0.0 ) {
                // double log_val = log( static_cast<double>(count_vector[i]) ) / log( mass_max );
                double log_val = static_cast<double>(count_vector[i]) / mass_max;
                color_vector[i] = utils::gradient( gradient, log_val );
            }
        }
    }

    return color_vector;
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

// /**
//  * @brief Write an SVG file containing a tree with colored branches.
//  */
// void write_color_tree_to_svg(
//     placement::PlacementTree const&  tree,
//     std::vector<utils::Color> const& colors_per_branch,
//     std::string const&               svg_filename
// ) {
//     auto copy = tree;
//     ladderize(copy);
//     // LOG_DBG << "Ladderize tree: " << ( validate_topology(copy) ? "valid!" : "INVALID!" );
//
//     {
//         // Make a layout tree.
//         // auto layout = tree::RectangularLayout( tree );
//         auto layout = tree::CircularLayout( copy, tree::CircularLayout::Type::kCladogram );
//
//         // Set edge colors.
//         std::vector<utils::SvgStroke> strokes;
//         for( auto color : colors_per_branch ) {
//             strokes.push_back( utils::SvgStroke() );
//             strokes.back().color = color;
//             strokes.back().line_cap = utils::SvgStroke::LineCap::kRound;
//             // strokes.back().width = 3;
//             strokes.back().width = 1;
//         }
//         layout.set_edge_strokes( strokes );
//
//         // Write to file
//         std::ostringstream out;
//         layout.to_svg_document().write( out );
//         utils::file_write( out.str(), svg_filename + "_c.svg" );
//     }
//     {
//         // Make a layout tree.
//         auto layout = tree::RectangularLayout( tree, tree::RectangularLayout::Type::kCladogram );
//
//         // Set edge colors.
//         std::vector<utils::SvgStroke> strokes;
//         for( auto color : colors_per_branch ) {
//             strokes.push_back( utils::SvgStroke() );
//             strokes.back().color = color;
//             strokes.back().line_cap = utils::SvgStroke::LineCap::kRound;
//             // strokes.back().width = 3;
//             strokes.back().width = 1;
//         }
//         layout.set_edge_strokes( strokes );
//
//         // Write to file
//         std::ostringstream out;
//         layout.to_svg_document().write( out );
//         utils::file_write( out.str(), svg_filename + "_r.svg" );
//     }
// }

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
    //     Find and read sample jplace files.
    // -------------------------------------------------------------------------

    // Get all jplace files, either in the epa dir, or in its subdirs
    auto jplace_filenames = utils::dir_list_files( epadir, ".*\\.jplace" );
    // std::sort( jplace_filenames.begin(), jplace_filenames.end() );

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
    size_t file_count = 0;
    for( auto const& jplace_filename : jplace_filenames ) {

        // Progress output.
        // LOG_INFO << "at jplace file " << file_count << " (" << jplace_filename << ")";
        ++file_count;

        sset.add( jplace_reader.from_file( jplace_filename ), jplace_filename );
    }

    // Final output for jplace reading
    assert( sset.size() == jplace_filenames.size() );
    LOG_INFO << "finished reading sample jplace files";
    LOG_INFO;

    LOG_INFO << "Using " << utils::Options::get().number_of_threads() << " threads.";

    // -------------------------------------------------------------------------
    //     Calculations
    // -------------------------------------------------------------------------

    tree::TreeSet avg_tree_set;
    for( auto const& smp : sset ) {
        avg_tree_set.add( "", smp.sample.tree() );
    }
    auto const avg_tree = tree::average_branch_length_tree( avg_tree_set );
    avg_tree_set.clear();

    LOG_INFO << "Converting Trees";
    auto mass_trees = convert_sample_set_to_mass_trees( sset );
    sset.clear();

    LOG_INFO << "Starting squash clustering";
    // auto sc = tree::squash_clustering( std::move( mass_trees.first ));
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
    for( auto& fn : jplace_filenames ) {
        fn = utils::replace_all( fn, ".jplace", "" );
    }

    LOG_INFO << "Writing cluster tree";
    std::ofstream file_clust_tree;
    utils::file_output_stream( outdir + "/cluster.newick",  file_clust_tree );
    file_clust_tree << sc.tree_string( jplace_filenames );

    LOG_INFO << "Writing fat trees";
    for( size_t i = 0; i < sc.clusters().size(); ++i ) {
        auto const& cc = sc.clusters()[i];

        auto const cv = mass_tree_mass_per_edge( cc.tree );
        auto const colors = counts_to_colors(cv, false);

        write_color_tree_to_nexus( avg_tree, colors, outdir + "/tree_" + std::to_string(i) + ".nexus" );
        // write_color_tree_to_svg( avg_tree, colors, outdir + "/tree_" + std::to_string(i) );
        write_color_tree_to_svg_file( avg_tree, {}, colors, outdir + "/tree_" + std::to_string(i) + ".svg" );
    }


    // utils::file_write( utils::to_string( emd_matrix ), outdir + "emd.mat" );
    // utils::file_write( utils::to_string( nhd_matrix ), outdir + "nhd.mat" );

    LOG_INFO << "Finished";
    return 0;
}
