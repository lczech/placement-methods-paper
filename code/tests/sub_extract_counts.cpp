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

#include "genesis.hpp"

#include <fstream>
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>

/**
	last used in 2017-01
	used for the old first subtree graphs for evaluating the accuarcy of the art trees.
    see 01_backbone/05_epa_sub_accuracy for the results.
*/

using namespace genesis;
using namespace genesis::placement;

struct SubtreeStats
{
    size_t total_placements   = 0;
    size_t correct_placements = 0;
    size_t edge_dist_sum      = 0;
    double bl_dist_sum        = 0.0;
    double avg_bl             = 0.0;
};

std::string process_jplace( std::string const& jplace_file, std::string const& target_taxon, SubtreeStats& stats )
{
    // Prepare reader and writer
    JplaceReader reader;
    Sample smp;

    // count how many placements there are on each branch, using the best placement only.
    std::unordered_map< size_t, size_t > counts;
    reader.from_file( jplace_file, smp );
    for( auto const& pquery : smp.pqueries() ) {
        ++counts[  pquery.placement_at(0).edge().index() ];

        // if we wanted, we could also use the weights, to be more granular.
        // but for our use case, this is not needed.
        // for( auto const& place : pquery.placements() ) {
        //     ++counts[ place.edge().index() ];
        // }
    }

    // LOG_INFO << target_taxon;
    // LOG_INFO << smp.pquery_size() << " pqueries";

    // count sorts, most populated branch first
    std::vector< std::pair<size_t, size_t> > counts_sort;
    for( auto const& entry : counts ) {
        counts_sort.push_back(entry);
    }
    std::sort(
        counts_sort.begin(),
        counts_sort.end(),
        [] ( std::pair<size_t, size_t> lhs, std::pair<size_t, size_t> rhs ) {
            return lhs.second > rhs.second;
        }
    );

    // prepare output for this sample
    std::string lines;

    // prepare lookups for how far the placements are from the target taxon:
    // how far in terms of edges between placement and target, and in terms of distance along branches
    auto target_node = tree::find_node( smp.tree(), target_taxon );
    if( target_node == nullptr ) {
        LOG_WARN << "No target node for " << target_taxon;
    }
    auto edge_dists = tree::edge_path_length_vector( smp.tree(), target_node->link().edge() );
    auto edge_bldis = tree::edge_branch_length_distance_vector( smp.tree(), &target_node->link().edge() );
    size_t dist_sum = 0;
    double bldi_sum = 0;
    size_t target_taxon_count = 0;

    // go through the sorted list of counts per branch, and accumulate distances
    for( auto const& entry : counts_sort ) {
        auto const& edge = smp.tree().edge_at( entry.first );

        // check whether we are at the target taxon
        auto tipname = edge.secondary_node().data<tree::DefaultNodeData>().name;
        if( tipname == "" ) {
            tipname = "inner";
        }
        if( tipname == target_taxon ) {
            target_taxon_count = entry.second;
        }

        // accumulate distances
        auto dist = edge_dists[ edge.index() ];
        auto bldi = edge_bldis[ edge.index() ];
        dist_sum += dist;
        bldi_sum += bldi;

        // write a line with detail information to the output for this subgroup
        std::ostringstream out;
        out << tipname << "\t"                                   // name of the taxon of the branch
            << edge.index() << "\t"                              // index of the branch
            << edge.data<PlacementEdgeData>().edge_num() << "\t" // edge num
            << entry.second << "\t"                              // number of placements on the branch
            << dist << "\t"                                      // distance in number of edges
            << bldi                                              // distance in branch length units
        ;
        lines += out.str() + "\n";
    }

    // write subgroup info file
    utils::file_write( lines, jplace_file + ".csv" );

    // prepare some stats
    auto correct_frac = static_cast<double>( target_taxon_count ) / static_cast<double>( smp.pquery_size() );
    auto avg_dist     = static_cast<double>(dist_sum) / static_cast<double>(smp.pquery_size());
    auto avg_bldi     = static_cast<double>(bldi_sum) / static_cast<double>(smp.pquery_size());
    auto avg_bl       = tree::length( smp.tree() ) / static_cast<double>( smp.tree().edge_count() );

    // accumulate stats
    stats.total_placements   += smp.pquery_size();
    stats.correct_placements += target_taxon_count;
    stats.edge_dist_sum      += dist_sum;
    stats.bl_dist_sum        += bldi_sum;
    stats.avg_bl              = avg_bl;

    // return a line for the overview info file
    std::ostringstream out;
    out << target_taxon << "\t"       // name of the target taxon
        << smp.pquery_size() << "\t"  // total number of placments
        << target_taxon_count << "\t" // number of correct placements
        << correct_frac << "\t"       // fraction of exaclty correct placements
        << dist_sum << "\t"           // total dist in edges
        << avg_dist << "\t"           // avg dist in edges
        << bldi_sum << "\t"           // total dist in bl units
        << avg_bldi << "\t"           // avg dist in bl units
        << avg_bl                     // avg branch length of the tree
    ;
    return out.str();
}

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    // LOG_INFO << "Started " << utils::current_time();

    // Check if the command line contains the right number of arguments.
    if (argc != 2) {
        throw std::runtime_error(
            "Need to provide dir path.\n"
        );
    }

    // Get list of jplace directories in input dir
    auto input_dir = std::string( argv[1] );
    LOG_INFO << "Input dir: " << input_dir;
    auto input_list = utils::dir_list_directories( input_dir, "job_.*" );
    LOG_INFO << "Found " << input_list.size() << " dirs.";

    // prepare output
    std::string lines;
    lines += "Target Name\tTotal Placements\tCorrect Placements\tCorrect Fraction\t";
    lines += "Edge Dist Total\tEdge Dist Avg\tBL Dist Total\tBL Dist Avg\tAvg BL\n";
    auto total_stats = SubtreeStats();

    // process all jplace files
    for( auto const& jplace_dir : input_list ) {
        LOG_DBG1 << "Processing " << jplace_dir.substr(0, 80) << ( jplace_dir.length() > 80 ? "..." : "" );

        // we assume there is only one jplace file per directory.
        auto epa_files = utils::dir_list_files( input_dir + "/" + jplace_dir, "epa_result\\..*\\.jplace" );
        if( epa_files.size() != 1 ) {
            LOG_WARN << "Dir has " << epa_files.size() << " jplace files.";
        }

        // get the name of the taxon of this sub-group
        auto target_taxon = utils::file_basename( jplace_dir ).substr( 4 );
        // LOG_DBG2 << "Target taxon: " << target_taxon;

        auto line = process_jplace( input_dir + "/" + jplace_dir + "/" + "epa_result.0.jplace", target_taxon, total_stats );
        lines += line + "\n";
    }

    // prepare stats.
    auto correct_frac  = static_cast<double>( total_stats.correct_placements ) / static_cast<double>( total_stats.total_placements );
    auto avg_edge_dist = static_cast<double>( total_stats.edge_dist_sum ) / static_cast<double>( total_stats.total_placements );
    auto avg_bl_dist   = static_cast<double>( total_stats.bl_dist_sum ) / static_cast<double>( total_stats.total_placements );

    // write summary stats to file.
    std::ostringstream out;
    out << "\t"
        << total_stats.total_placements << "\t"
        << total_stats.correct_placements << "\t"
        << correct_frac << "\t"
        << total_stats.edge_dist_sum << "\t"
        << avg_edge_dist << "\t"
        << total_stats.bl_dist_sum << "\t"
        << avg_bl_dist << "\t"
        << avg_bl_dist / total_stats.avg_bl
    ;
    lines += out.str();

    utils::file_write( lines, input_dir + "/overview.csv" );
    // LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
