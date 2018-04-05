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
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

// =================================================================================================
//     Typedefs
// =================================================================================================

/**
 * @brief Contains a list of clades, each itself containing a list of taxa belonging to that clade.
 */
using CladeTaxaList = std::unordered_map<std::string, std::vector<std::string>>;

/**
 * @brief A jplace reference tree.
 */
using TreeType      = placement::PlacementTree;

/**
 * @brief Contains a list of clades, each itself containing a list of edge indices belonging to
 * that clade.
 */
using CladeEdgeList = std::unordered_map<std::string, std::unordered_set<size_t>>;

// =================================================================================================
//     Get Clade Taxa Lists
// =================================================================================================

CladeTaxaList get_clade_taxa_lists( std::vector<std::string> const& clade_prefixes, TreeType const& tree )
{
    LOG_INFO << "tree has " << leaf_node_count( tree ) << " leaves";

    // Create a list of all clades and fill each clade with its taxa.
    CladeTaxaList clades;
    for( auto const& prefix : clade_prefixes ) {
        for( auto const& ni : tree.nodes() ) {
            auto const& name = ni->data<DefaultNodeData>().name;
            if( starts_with( name, prefix ) ) {
                clades[ prefix ].push_back( name );
            }
        }

        LOG_INFO << prefix << " has " << clades[ prefix ].size() << " taxa";
    }

    return clades;
}

// =================================================================================================
//     Get Clade Edges
// =================================================================================================

CladeEdgeList get_clade_edges( CladeTaxaList const& clades, TreeType const& tree )
{
    // Prepare the result map.
    CladeEdgeList clade_edges;

    // Process all clades.
    for( auto const& clade : clades ) {

        // Find the nodes that belong to the taxa of this clade.
        std::vector< tree::TreeNode const* > node_list;
        for( auto const& taxon : clade.second ) {
            auto node = find_node( tree, taxon );
            if( node == nullptr ) {
                LOG_WARN << "Cannot find taxon " << taxon;
                continue;
            }
            node_list.push_back( node );
        }

        // Find the edges that are part of the monophyletic subtrees of this clade.
        // This part is a bit messy and might be cleaned up in the future.
        auto bipartitions = bipartition_set( tree );
        auto subedges     = find_monophyletic_subtree_edges( tree, bipartitions, node_list );
        auto const subedge_map = std::unordered_set<size_t>(
            subedges.begin(), subedges.end()
        );

        // Add them to the clade edges list.
        clade_edges[ clade.first ] = subedge_map;
    }

    return clade_edges;
}

// =================================================================================================
//      Visualize Clades
// =================================================================================================

void visualize_clades(
    std::vector<std::string> const& clade_prefixes,
    CladeEdgeList const& clade_edges,
    TreeType const& tree,
    std::string const& out_file
) {
    auto const base_color = Color( 0.6, 0.6, 0.6 );
    std::vector<utils::Color> color_vector( tree.edge_count(), base_color );
    auto cm = ColorMap( color_list_set1() );

    // Count if there were any branches that are in more than once clade.
    // should not happen!
    size_t double_bookings = 0;

    size_t pos = 0;
    for( auto prefix : clade_prefixes ) {
        LOG_INFO << prefix << " gets color " << pos << " = " << cm.color( pos );

        auto const& cle = clade_edges.at( prefix );
        for( auto ei : cle ) {
            if( color_vector[ ei ] != base_color ) {
                ++double_bookings;
            }

            color_vector[ ei ] = cm.color( pos );
        }
        ++pos;
    }

    LOG_INFO << "had " << double_bookings << " double bookings";

    LayoutParameters params;
    params.stroke.width = 8;

    write_color_tree_to_svg_file( tree, params, color_vector, out_file );
}

// =================================================================================================
//      Eval
// =================================================================================================

void eval( CladeEdgeList const& clade_edges, Sample const& sample )
{
    size_t total_count = 0;

    using ListType = std::unordered_map< std::string, size_t >;

    // count how many seqs per clade we found and how many were correctly placed in the clade.
    ListType total;
    ListType best;
    ListType mass_50;
    ListType mass_80;
    ListType mass_90;
    ListType mass_95;

    for( auto const& pquery : sample ) {
        ++total_count;
        if( pquery.name_size() != 1 ) {
            LOG_ERR << "pquery.name_size() != 1";
            continue;
        }
        if( pquery.placement_size() < 1 ) {
            LOG_ERR << "pquery.placement_size() < 1";
            continue;
        }

        auto const name = pquery.name_at(0).name.substr( 11 );
        auto const ppos = name.find_first_of( "_", 10 );
        auto const prefix = name.substr( 0, ppos + 1 );

        // Skip sequences that are not in one of our clades
        if( clade_edges.count( prefix ) == 0 ) {
            continue;
        }

        ++total[ prefix ];
        if( clade_edges.at( prefix ).count( pquery.placement_at(0).edge().index() ) > 0 ) {
            ++best[ prefix ];
        }

        // count how much mass is correctly placed
        double mass = 0.0;
        for( auto const& pl : pquery.placements() ) {
            if( clade_edges.at( prefix ).count( pl.edge().index() ) > 0 ) {
                mass += pl.like_weight_ratio;
            }
        }
        if( mass > 0.5 ) {
            ++mass_50[ prefix ];
        }
        if( mass > 0.8 ) {
            ++mass_80[ prefix ];
        }
        if( mass > 0.9 ) {
            ++mass_90[ prefix ];
        }
        if( mass > 0.95 ) {
            ++mass_95[ prefix ];
        }
    }

    auto print = [&]( ListType const& amount, ListType const& total ){
        for( auto const& e : clade_edges ) {
            auto const prefix = e.first;
            auto ratio = static_cast<double>( amount.at( prefix ) ) / static_cast<double>( total.at( prefix ) );
            LOG_DBG1 << prefix << ": " << amount.at( prefix ) << " / " << total.at( prefix ) << " = " << ratio;
        }
    };

    LOG_INFO << "best placement in correct clade:";
    print( best, total );

    LOG_INFO << "mass_50 placement in correct clade:";
    print( mass_50, total );

    LOG_INFO << "mass_80 placement in correct clade:";
    print( mass_80, total );

    LOG_INFO << "mass_90 placement in correct clade:";
    print( mass_90, total );

    LOG_INFO << "mass_95 placement in correct clade:";
    print( mass_95, total );
}

// =================================================================================================
//      Main
// =================================================================================================

int main( int argc, char** argv )
{
    // Activate logging.
    Logging::log_to_stdout();
    Logging::details.time = true;
    LOG_INFO << "Started";

    // Load new base dir if given.
    if( argc != 3 ) {
        throw std::runtime_error( "expects in and out dir as params" );
    }

    std::string jplace_file = argv[1];
    std::string outdir  = utils::dir_normalize_path( argv[2] );

    auto const prefixes = std::vector<std::string>{
        "Bacteria_Actinobacteria_",
        "Bacteria_Firmicutes_",
        "Bacteria_Proteobacteria_",
        "Bacteria_Bacteroidetes_",
        "Bacteria_Cyanobacteria_"
    };

    LOG_INFO << "reading sample";
    auto sample = JplaceReader().from_file( jplace_file );
    sort_placements_by_weight( sample );

    LOG_INFO << "making clades";
    auto const taxa_list = get_clade_taxa_lists( prefixes, sample.tree() );
    auto const clade_edges = get_clade_edges( taxa_list, sample.tree() );

    LOG_INFO << "viz";
    visualize_clades( prefixes, clade_edges, sample.tree(), outdir + "clades.svg" );

    LOG_INFO << "eval";
    eval( clade_edges, sample );

    LOG_INFO << "Finished";
    return 0;
}
