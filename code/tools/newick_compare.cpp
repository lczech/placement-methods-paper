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
#include <vector>

using namespace genesis;
using namespace genesis::tree;
using namespace genesis::utils;

/**
 * @brief Compare two newick files and visualize their confilcting branches.
 * This is baiscally a viz of the RF distances.
 */
int main( int argc, char** argv )
{
    // Check if the command line contains the right number of arguments.
    if (argc != 3) {
        throw std::runtime_error(
            "Need to provide two newick tree files.\n"
        );
    }

    // Activate logging.
    utils::Logging::log_to_stdout();
    // LOG_INFO << "Started " << utils::current_time();

    auto const newick_path_a = std::string( argv[1] );
    auto const newick_path_b = std::string( argv[2] );
    auto const tree_a = DefaultTreeNewickReader().from_file( newick_path_a );
    auto const tree_b = DefaultTreeNewickReader().from_file( newick_path_b );

    auto const names_a = node_names( tree_a, true );
    auto const names_b = node_names( tree_b, true );
    if( names_a != names_b ) {
        throw std::runtime_error( "Not same taxa set." );
    }

    auto const bips_a = bipartition_set( tree_a );
    auto const bips_b = bipartition_set( tree_b );

    auto make_edge_colors = [](
        Tree const& tree_lhs, std::vector<Bipartition> const& bips_lhs,
        Tree const& tree_rhs, std::vector<Bipartition> const& bips_rhs
    ){
        (void) tree_rhs;
        auto const bad = Color( 0.8, 0.117, 0.286 );
        auto colors = std::vector<Color>( tree_lhs.edge_count(), bad );

        for( size_t i = 0; i < bips_lhs.size(); ++i ) {
            bool found = false;
            for( size_t j = 0; j < bips_rhs.size(); ++j ) {
                if( bips_lhs[i].leaf_nodes() == bips_rhs[j].leaf_nodes() ) {
                    found = true;
                    break;
                }
            }

            // if( bips_lhs[i].link().node().is_leaf()  ) {
            //     LOG_DBG1 << i << " " << bips_lhs[i].link().node().data<DefaultNodeData>().name << " " << bips_lhs[i].leaf_nodes();
            // } else {
            //     LOG_DBG1 << i << "   " << bips_lhs[i].leaf_nodes();
            // }

            assert( colors[ bips_lhs[i].link().edge().index() ] == bad );
            if( found ) {
                colors[ bips_lhs[i].link().edge().index() ] = Color( 0.117, 0.8, 0.286 );
            }
            if( bips_lhs[i].link().node().is_leaf()  ) {
                colors[ bips_lhs[i].link().edge().index() ] = Color( 0.7, 0.7, 0.7 );
            }
        }
        return colors;
    };

    // LOG_DBG << PrinterDetailed().print( tree_a );
    // LOG_DBG << PrinterDetailed().print( tree_b );

    LOG_DBG << "A";
    auto const colors_a = make_edge_colors( tree_a, bips_a, tree_b, bips_b );

    LOG_DBG << "B";
    auto const colors_b = make_edge_colors( tree_b, bips_b, tree_a, bips_a );

    LayoutParameters params;
    params.stroke.width = 15;

    write_color_tree_to_svg_file( tree_a, params, colors_a, newick_path_a + ".svg" );
    write_color_tree_to_svg_file( tree_b, params, colors_b, newick_path_b + ".svg" );

    // LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
