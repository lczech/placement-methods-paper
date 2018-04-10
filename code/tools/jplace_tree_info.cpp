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

#include "genesis/placement.hpp"
#include "genesis/tree.hpp"
#include "genesis/utils.hpp"

#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>

using namespace genesis;
using namespace genesis::placement;

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    LOG_INFO << "Started " << utils::current_time();

    // Check if the command line contains the right number of arguments.
    if (argc != 2) {
        throw std::runtime_error(
            "Need to provide a jplace file path.\n"
        );
    }

    auto infile = std::string( argv[1] );
    JplaceReader reader;
    Sample smp = reader.from_file( infile );
    auto const& tree = smp.tree();

    LOG_INFO << "Bifurcating:  " << ( is_bifurcating( tree ) ? "yes" : "no");
    LOG_INFO << "Edges:        " << tree.edge_count();
    LOG_INFO << "Nodes:        " << tree.node_count();
    LOG_INFO << "Leaves:       " << leaf_node_count( tree );
    LOG_INFO << "Length:       " << length( tree ) << " (sum of all branch lengths)";
    LOG_INFO << "Height:       " << height( tree ) << " (longest distance from root to a leaf)";
    LOG_INFO << "Diameter:     " << diameter( tree ) << " (longest distance between any two nodes)";

    std::vector<double> all_bl;
    std::vector<double> tip_bl;
    std::vector<double> inner_bl;

    for( auto const& edge : tree.edges() ) {
        auto const bl = edge->data<tree::DefaultEdgeData>().branch_length;
        all_bl.push_back( bl );

        if( edge->secondary_node().is_leaf() ) {
            tip_bl.push_back(bl);
        } else {
            inner_bl.push_back(bl);
        }
    }

    LOG_INFO << "Avg BL:       " << utils::mean_stddev( all_bl ).mean << " (average branch length)";
    LOG_INFO << "Avg tip BL:   " << utils::mean_stddev( tip_bl ).mean << " (average branch length of just leaf branches)";
    LOG_INFO << "Avg inner BL: " << utils::mean_stddev( inner_bl ).mean << " (average branch length of just inner branches)";

    tree::DefaultTreeNewickWriter().to_file( tree, infile + ".newick" );

    // Final output.
    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
