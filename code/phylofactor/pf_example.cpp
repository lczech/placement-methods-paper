/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2019 Lucas Czech and HITS gGmbH

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

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

// =================================================================================================
//     write_heat_tree
// =================================================================================================

void write_heat_tree(
    Matrix<double> const& raw_edge_masses,
    Tree const& tree,
    std::string const& out_dir
) {
    LOG_DBG << "write_heat_tree";

    assert( raw_edge_masses.cols() == tree.edge_count() );

    // Transpose and sort the matrix by edges/nodes.
    std::vector<size_t> tmp;
    for( auto const& e : tree.edges() ) {
        tmp.push_back( e.secondary_link().node().index() );
    }
    auto const sorting = utils::sort_indices( tmp.begin(), tmp.end() );
    assert( sorting.size() == tree.node_count() - 1 );
    auto matrix = Matrix<double>( raw_edge_masses.cols(), raw_edge_masses.rows() );
    for( size_t i = 0; i < sorting.size(); ++i ) {
        assert( sorting[i] < matrix.rows() );
        matrix.row( sorting[i] ) = raw_edge_masses.col(i).to_vector();
    }

    // center log ratio transform per sample.
    for( size_t c = 0; c < matrix.cols(); ++c ) {
        auto log_vec = matrix.col(c).to_vector();
        for( auto& e : log_vec ) {
            e = std::log(e);
        }
        auto const avg = arithmetic_mean(log_vec);
        for( size_t r = 0; r < log_vec.size(); ++r ) {
            auto& e = log_vec[r];
            e -= avg;
        }
        matrix.col(c) = log_vec;
        // auto const mm = minimum_maximum( log_vec );
        // LOG_DBG1 << "heat mat col " << c << " min " << mm.min << " max " << mm.max << " avg " << avg;
    }

    // Calc tree masses
    auto const edge_masses = matrix_col_sums( raw_edge_masses );
    assert( edge_masses.size() == tree.edge_count() );

    // Make color maps and norms.
    auto const mm    = matrix_minmax(matrix);
    auto matrix_norm = ColorNormalizationLinear( mm.min, mm.max );
    auto tree_norm   = ColorNormalizationLogarithmic( edge_masses.begin(), edge_masses.end() );

    auto write_heattree_ = [&](
        ColorMap const& map_mat, ColorMap const& map_tree, std::string const& mname
    ){

        // Make color matrix
        auto mat_col = Matrix<Color>( matrix.rows(), matrix.cols() );
        for( size_t r = 0; r < matrix.rows(); ++r ) {
            for( size_t c = 0; c < matrix.cols(); ++c ) {
                mat_col( r, c ) = map_mat( matrix_norm, matrix(r,c) );
            }
        }

        // Prepare heat tree
        HeatTreeParameters params;
        params.tree = tree;
        params.ladderize = false;
        params.type = LayoutType::kPhylogram;
        params.color_per_branch = map_tree( tree_norm, edge_masses );
        params.matrix = mat_col;

        // params.column_labels =
        auto const doc = heat_tree( params, map_mat , matrix_norm, map_tree, tree_norm );
        std::ostringstream out;
        doc.write( out );
        utils::file_write( out.str(), out_dir + "heat_tree_" + mname + ".svg" );
    };

    // write_heattree_( ColorMap( color_list_viridis() ), ColorMap( color_list_viridis() ), "viridis" );
    // write_heattree_( ColorMap( color_list_heat() ), ColorMap( color_list_heat() ), "heat" );
    // write_heattree_( ColorMap( color_list_heat() ), ColorMap( color_list_heat() ), "heat" );
    // write_heattree_( ColorMap( color_list_heat() ), ColorMap( color_list_viridis() ), "viridis_heat" );
    write_heattree_( ColorMap( color_list_heat() ), ColorMap( color_list_bupubk() ), "bupubk_heat" );
    // write_heattree_( ColorMap( color_list_orrd() ), "orrd" );
    // write_heattree_( ColorMap( color_list_ylorrd() ), "ylorrd" );
}

// =================================================================================================
//     main
// =================================================================================================

int main()
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.date = true;
    utils::Logging::details.time = true;
    utils::Options::get().number_of_threads( 4 );
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    std::string base_dir = "path/to/svg/phylofactor/";
    std::string in_tree = base_dir + "phylofactor.newick";

    auto tree = convert_common_tree_to_placement_tree(
        CommonTreeNewickReader().read( from_file( in_tree ))
    );

    auto add_edge_weights = [&]( std::vector<std::string> const& names, double extra_weight, std::vector<double>& edge_weights ){
        std::vector<TreeNode const*> nodes;
        for( auto const& taxon : names ) {
            auto const np = find_node( tree, taxon );
            if( !np ) {
                LOG_ERR << "no node " << taxon;
            }
            nodes.push_back(np);
        }

        auto const bips = bipartition_set( tree );

        auto with_split_edge = find_monophyletic_subtree_edges( tree, bips, nodes, true );
        for( auto f : with_split_edge ) {
            edge_weights[f] += extra_weight;
        }

        auto without_split_edge = find_monophyletic_subtree_edges( tree, bips, nodes, false );
        for( auto f : without_split_edge ) {
            edge_weights[f] += extra_weight;
        }

        return edge_weights;
    };

    Simulator sim;
    SampleSet samples;
    auto edge_weights = std::vector<double>( tree.edge_count(), 2 );

    add_edge_weights({ "A", "B", "C", "D" }, 5, edge_weights );
    sim.edge_distribution().edge_weights = edge_weights;
    for( size_t i = 0; i < 5; ++i ) {
        auto smpl = Sample( tree );
        sim.generate( smpl, 1000 );
        samples.add( smpl );
    }

    add_edge_weights({ "A", "B", "C", "D" }, -3, edge_weights );
    add_edge_weights({ "F", "G" }, 3, edge_weights);
    sim.edge_distribution().edge_weights = edge_weights;
    for( size_t i = 0; i < 5; ++i ) {
        auto smpl = Sample( tree );
        sim.generate( smpl, 1000 );
        samples.add( smpl );
    }

    auto const masses = placement_mass_per_edges_with_multiplicities(samples);
    write_heat_tree( masses, tree, base_dir );

    // Matrix<double> masses

    // auto edge_vals = std::vector<size_t>( tree.edge_count(), 0 );
    // auto colors = std::vector<Color>( tree.edge_count() );
    // LayoutParameters lm;
    // lm.ladderize = true;
    // lm.stroke.width = 10;
    //
    // write_color_tree_to_svg_file(
    //     tree, lm, colors,
    //     base_dir + "multi_factors_tree.svg"
    // );

    return 0;
}
