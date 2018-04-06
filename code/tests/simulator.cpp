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

/**
 * This is the demo "Visualize Placements". See the Manual for more information.
 */

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "placement/formats/jplace_reader.hpp"
#include "placement/formats/jplace_writer.hpp"
#include "placement/formats/newick_writer.hpp"
#include "placement/function/functions.hpp"
#include "placement/function/helper.hpp"
#include "placement/function/measures.hpp"
#include "placement/function/operators.hpp"
#include "placement/placement_tree.hpp"
#include "placement/sample.hpp"
#include "placement/sample_set.hpp"

#include "placement/simulator/functions.hpp"
#include "placement/simulator/simulator.hpp"

#include "tree/default/distances.hpp"
#include "tree/default/functions.hpp"
#include "tree/default/newick_reader.hpp"
#include "tree/formats/newick/color_writer_mixin.hpp"
#include "tree/formats/newick/reader.hpp"
#include "tree/function/functions.hpp"

#include "utils/core/fs.hpp"
#include "utils/core/logging.hpp"

#include "utils/formats/nexus/document.hpp"
#include "utils/formats/nexus/taxa.hpp"
#include "utils/formats/nexus/trees.hpp"
#include "utils/formats/nexus/writer.hpp"

#include "utils/text/string.hpp"
#include "utils/tools/color.hpp"
#include "utils/tools/color/gradient.hpp"
#include "utils/tools/color/operators.hpp"
#include "utils/tools/date_time.hpp"

#include "utils/math/matrix.hpp"
#include "utils/math/matrix/operators.hpp"

using namespace genesis;

// =================================================================================================
//      Count Placement Mass Per Edge
// =================================================================================================

/**
 * @brief Examine all Placements in a Sample and add their like_weight_ratio to the branch where
 * the placement is located.
 *
 * The function loops over all Placements of all Pqueries of the given Sample. For each Placement,
 * it adds the like_weight_ratio to the given vector at the index position of the branch where
 * the Placement is located.
 */
void count_placement_mass_per_edge(
    placement::Sample const& sample,
    std::vector<double>&     placement_mass
) {
    using namespace ::genesis::placement;

    // Check whether the provided vector has the same number of elements as the tree has edges.
    if( placement_mass.size() != sample.tree().edge_count() ) {
        throw std::runtime_error( "Placement mass vector has wrong size." );
    }

    // Loop over all placements of all pqueries of the sample and accumulate the mass.
    for( auto const& pquery : sample ) {
        for( auto const& placement : pquery.placements() ) {
            auto index = placement.edge().index();
            placement_mass[index] += placement.like_weight_ratio;
        }
    }
}

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
 */
void write_color_tree_to_nexus(
    placement::PlacementTree const&  tree,
    std::vector<utils::Color> const& colors_per_branch,
    std::string const&               nexus_filename
) {
    // We use a normal Newick writer for PlacementTrees, but also wrap it in a Color Mixin
    // in order to allow for color annotated branches.
    using ColorWriter = tree::NewickColorWriterMixin<placement::PlacementTreeNewickWriter>;

    // Get the Newick representation of the tree, with color annotated branches.
    auto tree_writer = ColorWriter();
    tree_writer.enable_edge_nums(false);
    tree_writer.edge_colors(colors_per_branch);
    std::string newick_tree = tree_writer.to_string(tree);

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
//     Simulate Stuff
// =================================================================================================

void simulate_stuff(
    placement::PlacementTree const& ptree,
    size_t random_placements,
    size_t subtree_placements,
    std::string outdir
) {
    using namespace ::genesis::placement;
    using namespace ::genesis::tree;

    utils::dir_create( outdir );
    utils::dir_create( outdir + "trees" );
    utils::dir_create( outdir + "samples" );

    SampleSet smps;
    for( size_t i = 0; i < ptree.link_count(); ++i ) {
        if( ptree.link_at(i).is_leaf() ) {
            continue;
        }

        auto smp = Sample( ptree );

        Simulator sim;
        set_uniform_weights( smp, sim.edge_distribution() );
        sim.generate( smp, random_placements );

        set_subtree_weights( smp, i, sim.edge_distribution() );
        sim.generate( smp, subtree_placements );

        smps.add( smp );
    }

    size_t cnt = 0;
    for( auto const& nsmp : smps ) {
        validate( nsmp.sample );

        auto placement_mass = std::vector<double>( nsmp.sample.tree().edge_count(), 0.0 );
        count_placement_mass_per_edge( nsmp.sample, placement_mass);
        auto col_vec = counts_to_colors( placement_mass );

        write_color_tree_to_nexus(
            nsmp.sample.tree(),
            col_vec,
            outdir + "trees/tree_" + std::to_string( cnt ) + ".nexus"
        );

        JplaceWriter().to_file(
            nsmp.sample,
            outdir + "samples/sample_" + std::to_string( cnt ) + ".jplace"
        );

        ++cnt;
    }

    LOG_DBG2 << "starting emd " << utils::current_time();
    auto emd_matrix = earth_movers_distance( smps );
    LOG_DBG2 << "finished emd " << utils::current_time();
    // LOG_DBG << emd_matrix;

    utils::file_write( utils::to_string( emd_matrix ), outdir + "emd.mat" );

    LOG_DBG2 << "starting nhd " << utils::current_time();
    auto nhd_matrix = node_histogram_distance( smps );
    LOG_DBG2 << "finished nhd " << utils::current_time();
    // LOG_DBG << nhd_matrix;

    utils::file_write( utils::to_string( nhd_matrix ), outdir + "nhd.mat" );
}

// =================================================================================================
//     Main Function
// =================================================================================================

int main( int argc, char** argv )
{
    using namespace ::genesis::placement;
    using namespace ::genesis::tree;

    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();

    // std::string tree_file = "tax_slv_ssu_123.1.tre";
    // std::string tree_file = "genesis/comb2.newick";
    // std::string tree_file = "balanced3.newick";

    std::string base_dir = "placement_simulation/";
    std::vector<std::string> tree_file_names = { "balanced2", "balanced3", "comb" };

    for( auto const& tree_file_name : tree_file_names ) {
        LOG_DBG << "=====================================";
        LOG_DBG << "Tree " << tree_file_name;

        // LOG_DBG << "Loading tree...";
        DefaultTree intree;
        DefaultTreeNewickReader().from_file( base_dir + tree_file_name + ".newick", intree );
        set_all_branch_lengths( intree, 1.0 );
        // LOG_DBG << "done";

        // LOG_DBG << "deepest dist " << deepest_distance( intree );
        // LOG_DBG << "dia          " << diameter( intree );
        // LOG_DBG << "length       " << length( intree );
        // LOG_DBG << "nodes        " << intree.node_count();
        // LOG_DBG << "leaves       " << leaf_node_count( intree );
        // LOG_DBG << "inner        " << inner_node_count( intree );

        PlacementTree ptree = convert_to_placement_tree( intree );

        std::vector<size_t> randomness = { 0, 5, 10, 25, 50, 75, 100, 150, 250, 500 };
        for( auto const& rand_frac : randomness ) {
            LOG_DBG1 << "Randomness " << rand_frac;

            size_t sn = 1000;
            size_t rn = sn * rand_frac / 100;

            simulate_stuff( ptree, rn, sn, base_dir + tree_file_name + "_" + std::to_string(rand_frac) +  "/" );
        }
    }

    LOG_INFO << "Finished.";
    return 0;
}
