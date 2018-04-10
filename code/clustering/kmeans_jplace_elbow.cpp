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
std::vector<utils::Color> counts_to_colors(
    std::vector<double> const& count_vector
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
    LOG_INFO << "In order to build a color gradient, use these color stops:\n"
             << "    0.0 -> #81bfff (light blue, no placement mass)\n"
             << "    0.5 -> #c040be (purple, medium placement mass)\n"
             << "    1.0 -> #000000 (black, maximum placement mass)";
    LOG_INFO << "The edge with the maximum placement mass has a mass of " << mass_max;

    // More output for the color gradient: show the label positions.
    std::string labels = "The following list shows position on the color gradient ";
    labels += "with their correspinding placement mass. Use this to label your color gradient.\n";
    int    pow_i = 0;
    double pow_v = 0.0;
    do {
        pow_v = std::abs( log( pow( 10, pow_i )) / log( mass_max ));
        if( pow_v <= 1.0 ) {
            labels += "    " + utils::to_string_precise( pow_v, 6 ) + " -> ";
            labels += utils::to_string( pow( 10, pow_i )) + "\n";
        }
        if( mass_max > 1.0 ) {
            ++pow_i;
        } else {
            --pow_i;
        }
    } while( pow_v <= 1.0 && mass_max > 1.0 );
    labels += "    1.000000 -> " + utils::to_string( mass_max ) + "\n";
    LOG_INFO << labels;

    // Calculate the resulting colors.
    for( size_t i = 0; i < count_vector.size(); ++i ) {
        if( count_vector[i] > 0.0 ) {
            double log_val = log( static_cast<double>(count_vector[i]) ) / log( mass_max );
            color_vector[i] = utils::gradient( gradient, log_val );
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
    auto writer = tree::DefaultTreeNewickWriter();
    auto color_plugin = tree::NewickColorWriterPlugin();
    color_plugin.register_with( writer );

    // Get the Newick representation of the tree, with color annotated branches.
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
    tree::DefaultTree const&  tree,
    std::vector<utils::Color> const& colors_per_branch,
    std::string const&               svg_filename
) {
    auto copy = tree;
    ladderize(copy);
    LOG_DBG << "Ladderize tree: " << ( validate_topology(copy) ? "valid!" : "INVALID!" );

    {
        // Make a layout tree.
        // auto layout = tree::RectangularLayout( tree );
        auto layout = tree::CircularLayout( copy );

        // Set edge colors.
        std::vector<utils::SvgStroke> strokes;
        for( auto color : colors_per_branch ) {
            strokes.push_back( utils::SvgStroke() );
            strokes.back().color = color;
            strokes.back().line_cap = utils::SvgStroke::LineCap::kRound;
            // strokes.back().width = 3;
            strokes.back().width = 1;
        }
        layout.set_edge_strokes( strokes );

        // Write to file
        std::ostringstream out;
        layout.to_svg_document().write( out );
        utils::file_write( out.str(), svg_filename + "_c.svg" );
    }
    {
        // Make a layout tree.
        auto layout = tree::RectangularLayout( tree );

        // Set edge colors.
        std::vector<utils::SvgStroke> strokes;
        for( auto color : colors_per_branch ) {
            strokes.push_back( utils::SvgStroke() );
            strokes.back().color = color;
            strokes.back().line_cap = utils::SvgStroke::LineCap::kRound;
            // strokes.back().width = 3;
            strokes.back().width = 1;
        }
        layout.set_edge_strokes( strokes );

        // Write to file
        std::ostringstream out;
        layout.to_svg_document().write( out );
        utils::file_write( out.str(), svg_filename + "_r.svg" );
    }
}

// =================================================================================================
//     Phylo kmeans
// =================================================================================================

std::pair<double, double> run_pkmeans(
    SampleSet const& sset,
    std::vector<tree::MassTree> const& trees,
    size_t const k,
    std::string outdir
) {

    outdir = outdir + "k_" + std::to_string(k) + "/";
    utils::dir_create( outdir );

    // -------------------------------------------------------------------------
    //     Tree Kmeans
    // -------------------------------------------------------------------------

    LOG_INFO << "Kmeans started with k=" << k;
    auto mkmeans = tree::MassTreeKmeans();
    mkmeans.run( trees, k );
    LOG_INFO << "Kmeans finished";

    // utils::file_write( utils::to_string( emd_matrix ), outdir + "emd.mat" );

    // -------------------------------------------------------------------------
    //     Info
    // -------------------------------------------------------------------------

    LOG_INFO << "cluster info";
    auto const clustinfo = mkmeans.cluster_info( trees );
    auto const& massignments = mkmeans.assignments();
    if( massignments.size() != trees.size() ){
        LOG_ERR << "massignments.size() != trees.size()";
    }

    double avg_dist = 0.0;
    double avg_var = 0.0;
    for( size_t ik = 0; ik < k; ++ik ) {

        double cavg_dist = 0.0;
        size_t cavg_cnt = 0;
        for( size_t ai = 0; ai < massignments.size(); ++ai ) {
            if( massignments[ai] != ik ) {
                continue;
            }
            avg_dist  += clustinfo.distances[ai];
            avg_var   += clustinfo.distances[ai] * clustinfo.distances[ai];
            cavg_dist += clustinfo.distances[ai];
            ++cavg_cnt;
        }
        cavg_dist /= static_cast<double>( cavg_cnt );

        LOG_DBG1 << "cluster " << ik << ": " << clustinfo.counts[ik] << " samples, with a variance of " << clustinfo.variances[ik] << " and average dist " << cavg_dist;
    }
    avg_dist /= static_cast<double>( massignments.size() );
    avg_var  /= static_cast<double>( massignments.size() );

    double avg_dist2 = 0.0;
    double avg_var2 = 0.0;
    for( size_t i = 0; i < clustinfo.distances.size(); ++i ) {
        avg_dist2 += clustinfo.distances[i];
        avg_var2  += clustinfo.distances[i] * clustinfo.distances[i];
    }
    avg_dist2 /= static_cast<double>( clustinfo.distances.size() );
    avg_var2  /= static_cast<double>( clustinfo.distances.size() );
    LOG_DBG << "avg dist " << avg_dist << " and " << avg_dist2;
    LOG_DBG << "avg var " << avg_var << " and " << avg_var2;

    // -------------------------------------------------------------------------
    //     Output
    // -------------------------------------------------------------------------

    // Write assignments
    LOG_INFO << "Write assignments";
    std::ofstream file_mkmeans_ass;
    file_output_stream( outdir + "emd_assignments.csv",  file_mkmeans_ass );

    auto pqry_cnts = std::vector<size_t>( mkmeans.centroids().size(), 0 );
    for( size_t i = 0; i < massignments.size(); ++i ) {
        file_mkmeans_ass << file_filename( sset[i].name );
        file_mkmeans_ass << "\t" << massignments[i];
        file_mkmeans_ass << "\n";

        pqry_cnts[ massignments[i] ] += sset[i].sample.size();
    }
    file_mkmeans_ass.close();

    // Write centroids
    LOG_INFO << "Write centroids";
    std::ofstream file_mkmeans_cent;
    file_output_stream( outdir + "emd_centroids.csv",  file_mkmeans_cent );
    auto const& mcentroids = mkmeans.centroids();
    for( size_t i = 0; i < mcentroids.size(); ++i ) {
        auto masses = mass_tree_mass_per_edge( mcentroids[i] );
        for( auto& mass : masses ) {
            mass *= pqry_cnts[i];
        }

        // auto cent = mcentroids[i];
        // mass_tree_scale_masses( cent, list[i].size() );
        // auto colors_per_branch = counts_to_colors( mass_tree_mass_per_edge( cent ));

        auto colors_per_branch = counts_to_colors( masses );
        write_color_tree_to_nexus( sset[0].sample.tree(), colors_per_branch, outdir + "tree_emd_" + std::to_string(i) + ".nexus" );
        write_color_tree_to_svg( sset[0].sample.tree(), colors_per_branch, outdir + "tree_emd_" + std::to_string(i) );

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

    return { avg_dist2, avg_var2 };
}

std::pair<double, double> run_ikmeans(
    SampleSet const& sset,
    utils::Matrix<double> const& edge_imb_mat,
    std::vector<std::vector<double>> const& edge_imb_vec,
    std::vector<size_t> const& columns,
    size_t const k,
    std::string outdir
) {
    outdir = outdir + "k_" + std::to_string(k) + "/";
    utils::dir_create( outdir );

    // Kmeans
    LOG_INFO << "Kmeans started";
    auto ikmeans = EuclideanKmeans( edge_imb_mat.cols() );
    ikmeans.run( edge_imb_vec, k );
    LOG_INFO << "Kmeans finished";

    // -------------------------------------------------------------------------
    //     Info
    // -------------------------------------------------------------------------

    LOG_INFO << "cluster info";
    auto const clustinfo = ikmeans.cluster_info( edge_imb_vec );
    auto const& massignments = ikmeans.assignments();
    // assert( massignments.size() == trees.size() );

    double avg_dist = 0.0;
    double avg_var = 0.0;
    for( size_t ik = 0; ik < k; ++ik ) {

        double cavg_dist = 0.0;
        size_t cavg_cnt = 0;
        for( size_t ai = 0; ai < massignments.size(); ++ai ) {
            if( massignments[ai] != ik ) {
                continue;
            }
            avg_dist  += clustinfo.distances[ai];
            avg_var   += clustinfo.distances[ai] * clustinfo.distances[ai];
            cavg_dist += clustinfo.distances[ai];
            ++cavg_cnt;
        }
        cavg_dist /= static_cast<double>( cavg_cnt );

        LOG_DBG1 << "cluster " << ik << ": " << clustinfo.counts[ik] << " samples, with a variance of " << clustinfo.variances[ik] << " and average dist " << cavg_dist;
    }
    avg_dist /= static_cast<double>( massignments.size() );
    avg_var  /= static_cast<double>( massignments.size() );

    double avg_dist2 = 0.0;
    double avg_var2 = 0.0;
    for( size_t i = 0; i < clustinfo.distances.size(); ++i ) {
        avg_dist2 += clustinfo.distances[i];
        avg_var2  += clustinfo.distances[i] * clustinfo.distances[i];
    }
    avg_dist2 /= static_cast<double>( clustinfo.distances.size() );
    avg_var2  /= static_cast<double>( clustinfo.distances.size() );
    LOG_DBG << "avg dist " << avg_dist << " and " << avg_dist2;
    LOG_DBG << "avg var " << avg_var << " and " << avg_var2;

    // Write assignments
    LOG_INFO << "Write assignments";
    std::ofstream file_ikmeans_ass;
    file_output_stream( outdir + "imbalance_assignments.csv",  file_ikmeans_ass );
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
    file_output_stream( outdir + "imbalance_centroids.csv",  file_ikmeans_cent );
    auto const& icentroids = ikmeans.centroids();
    for( size_t i = 0; i < icentroids.size(); ++i ) {
        auto const& cent = icentroids[i];

        auto col_vec = std::vector<utils::Color>( sset[0].sample.tree().edge_count(), base_color );
        if( columns.size() != cent.size() ){
            LOG_ERR << "columns.size() != cent.size() ";
        }
        for( size_t j = 0; j < columns.size(); ++j ) {
            if(!( -1.0 <= cent[ j ] && cent[ j ] <= 1.0 )) {
                LOG_ERR << "!( -1.0 <= cent[ j ] && cent[ j ] <= 1.0 )";
            }
            if( cent[ j ] >= 0.0 ) {
                col_vec[ columns[j] ] = utils::gradient( p_gradient, cent[ j ] );
            } else {
                col_vec[ columns[j] ] = utils::gradient( n_gradient, - cent[ j ] );
            }
        }
        write_color_tree_to_nexus( sset[0].sample.tree(), col_vec, outdir + "tree_imb_" + std::to_string(i) + ".nexus" );
        write_color_tree_to_svg( sset[0].sample.tree(), col_vec, outdir + "tree_imb_" + std::to_string(i) );

        for( auto const& v : cent ) {
            if( &v != &cent[0] ) {
                file_ikmeans_cent << "\t";
            }
            file_ikmeans_cent << v;
        }
        file_ikmeans_cent << "\n";
    }
    file_ikmeans_cent.close();

    return { avg_dist2, avg_var2 };
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

    // Check if the command line contains the right number of arguments.
    if (argc != 5) {
        throw std::runtime_error(
            "Need to provide 4 arguments: threads max-k epadir outdir\n"
        );
    }

    // In out dirs.
    auto const threads = std::stoi( argv[1] );
    auto const maxk    = static_cast<size_t>( std::stoi( argv[2] ));
    auto const epadir  = utils::dir_normalize_path( std::string( argv[3] ));
    auto const outdir  = utils::dir_normalize_path( std::string( argv[4] ));
    utils::dir_create(outdir);

    utils::Options::get().command_line( argc, argv );
    utils::Options::get().number_of_threads( threads );
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    LOG_INFO << "Started";

    // -------------------------------------------------------------------------
    //     Find and read sample jplace files.
    // -------------------------------------------------------------------------

    // Get all jplace files, either in the epa dir, or in its subdirs
    auto jplace_filenames = utils::dir_list_files( epadir, ".*\\.jplace" );
    std::sort( jplace_filenames.begin(), jplace_filenames.end() );

    LOG_INFO << "reading " << jplace_filenames.size() << " jplace sample files";
    auto const sset = JplaceReader().from_files( jplace_filenames );
    assert( sset.size() == jplace_filenames.size() );
    LOG_INFO << "finished reading sample jplace files";
    LOG_INFO;

    std::ofstream file_avg_dists;
    file_output_stream( outdir + "avg_dists.txt",  file_avg_dists );
    std::ofstream file_avg_vars;
    file_output_stream( outdir + "avg_vars.txt",  file_avg_vars );

    // -------------------------------------------------------------------------
    //     Tree Kmeans
    // -------------------------------------------------------------------------

    // Kmeans
    LOG_INFO << "Using " << utils::Options::get().number_of_threads() << " threads.";

    // Options::get().random_seed( 1623390399 );

    LOG_INFO << "Converting Trees";
    auto mass_trees = convert_sample_set_to_mass_trees( sset );

    file_avg_dists << "phylo kmeans\n";
    file_avg_vars << "phylo kmeans\n";
    for( size_t ik = 1; ik < maxk; ++ik ) {
        auto const dv = run_pkmeans(
            sset,
            mass_trees.first,
            ik,
            outdir
        );
        file_avg_dists << "k " << ik << "\td " << dv.first << "\n";
        file_avg_vars << "k " << ik << "\td " << dv.second << "\n";
    }

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
    file_output_stream( outdir + "imbalance_matrix.csv",  file_imbmat );
    file_imbmat << edge_imb_mat;
    file_imbmat.close();

    file_avg_dists << "imbal kmeans\n";
    file_avg_vars << "imbal kmeans\n";
    for( size_t ik = 1; ik < maxk; ++ik ) {
        auto const dv = run_ikmeans(
            sset,
            edge_imb_mat,
            edge_imb_vec,
            columns,
            ik,
            outdir
        );
        file_avg_dists << "k " << ik << "\td " << dv.first << "\n";
        file_avg_vars << "k " << ik << "\td " << dv.second << "\n";
    }

    LOG_INFO << "Finished";
    return 0;
}
