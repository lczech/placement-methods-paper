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
//     main
// =================================================================================================

int main()
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.date = true;
    utils::Logging::details.time = true;
    utils::Options::get().number_of_threads( 40 );
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    std::string base_dir = "path/to/multi_tree/";
    std::string in_tree = base_dir + "best_tree.newick";
    std::string in_tw   = base_dir + "tw_taxa.txt";
    std::string in_no   = base_dir + "no_tw_taxa.txt";
    std::string in_pf   = base_dir + "pf_taxa.txt";

    auto tree = CommonTreeNewickReader().read( from_file( in_tree ));
    auto edge_vals = std::vector<size_t>( tree.edge_count(), 0 );

    auto set_vals = [&]( std::string const& file, size_t or_val ){
        auto const fc = file_read_lines( file );

        std::vector<TreeNode const*> nodes;
        for( auto const& taxon : fc ) {
            if( taxon.empty() ) {
                continue;
            }

            auto const np = find_node( tree, taxon );
            if( !np ) {
                LOG_ERR << "no node " << taxon;
            }
            nodes.push_back(np);
        }

        auto const bips = bipartition_set( tree );
        auto const edges = find_monophyletic_subtree_edges( tree, bips, nodes );

        for( auto e : edges ) {
            edge_vals[e] |= or_val;
        }
    };

    set_vals( in_tw, 1 );
    set_vals( in_no, 2 );
    set_vals( in_pf, 4 );

    auto color_set = std::vector<Color>({
        Color( 0.800000, 0.800000, 0.800000 ),
        Color( 0.301961, 0.686275, 0.290196 ), // green
        Color( 0.215686, 0.494118, 0.721569 ), // blue
        Color( 0.090196, 0.733333, 0.701960 ), // my teal
        Color( 0.894118, 0.101961, 0.109804 ), // red
        Color( 0.901961, 0.670588, 0.007843 ), // banana
        Color( 0.596078, 0.305882, 0.639216 ), // purple
        Color( 0.400000, 0.400000, 0.400000 ),
        // Color( 0, 0, 0 ),
        // Color( 0, 0, 1 ),
        // Color( 0, 1, 0 ),
        // Color( 0, 1, 1 ),
        // Color( 1, 0, 0 ),
        // Color( 1, 0, 1 ),
        // Color( 1, 1, 0 ),
        // Color( 1, 1, 1 ),
    });

    auto colors = std::vector<Color>( tree.edge_count() );
    for( size_t i = 0; i < tree.edge_count(); ++i ) {
        colors[i] = color_set[ edge_vals[i] ];
    }

    LayoutParameters lm;
    lm.ladderize = true;
    lm.stroke.width = 10;

    write_color_tree_to_svg_file(
        tree, lm, colors,
        base_dir + "multi_factors_tree.svg"
    );

    return 0;
}
