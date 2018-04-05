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

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

// =================================================================================================
//      Preprocess
// =================================================================================================

struct Lookups
{
    // Blacklist of sequences that we do not want to use, because either they are considered as
    // mislabels by sativa, or because their name contains "incertae" "unclassified" or "unknown",
    // indicating that this sequence is not trustworthy.
    std::unordered_set<std::string> blacklist;

    // Map from tip names to edge indices in the tree.
    std::unordered_map<std::string, size_t> tip_name_to_edge_index;

    // Matrix of edge distances in the tree, that is, number of nodes between two edges.
    utils::Matrix<size_t> edge_dists;

    // Matrix of branch distances, that is, distance using the branch lengths.
    utils::Matrix<double> branch_dists;

    // Matrix of pairwise dists between all nodes of the tree. Needed for edpl.
    utils::Matrix<double> node_dists;
};

struct Results
{
    std::string name;
    std::vector<double> weighted_placements_edge_distances;
    std::vector<double> weighted_placements_branch_distances;
};

Lookups preprocess( Sample const& sample )
{
    // -------------------------------------------------------------------------
    //     Prepare Lookups
    // -------------------------------------------------------------------------

    Lookups lookup;
    LOG_INFO << "Prepare lookups.";

    // Blacklist of sequences that we do not want to use, because either they are considered as
    // mislabels by sativa, or because their name contains "incertae" "unclassified" or "unknown",
    // indicating that this sequence is not trustworthy.
    auto blacklist_table = utils::CsvReader().separator_chars("\t").from_file( "/path/to/data/silva/blacklist.txt" );
    for( auto const& line : blacklist_table ) {
        if( line.size() != 2 ) {
            LOG_ERR << "line size " << line.size();
            continue;
        }
        if( lookup.blacklist.count( line[0] ) > 0 ) {
            LOG_WARN << "already in blacklist: " << line[0];
        }
        lookup.blacklist.insert( line[0] );
    }

    // Map from tip names to edge indices in the tree.
    for( auto const& node : sample.tree().nodes() ) {
        if( ! node->is_leaf() ) {
            continue;
        }
        if( node->rank() > 0 ) {
            LOG_ERR << "Tip node with more than one child";
            continue;
        }

        auto const name = node->data<DefaultNodeData>().name;
        if( lookup.tip_name_to_edge_index.count( name ) > 0 ) {
            LOG_ERR << "Name already in tip node list: " << name;
        }
        lookup.tip_name_to_edge_index[ name ] = node->link().edge().index();
    }

    // Matrix of edge distances in the tree, that is, number of nodes between two edges.
    lookup.edge_dists = edge_path_length_matrix( sample.tree() );

    // Matrix of branch distances, that is, distance using the branch lengths.
    lookup.branch_dists = edge_branch_length_distance_matrix( sample.tree() );

    // Matrix of pairwise dists between all nodes of the tree. Needed for edpl.
    lookup.node_dists = node_branch_length_distance_matrix( sample.tree() );

    LOG_INFO << "Finished lookups.";
    return lookup;
}

// =================================================================================================
//      Eval
// =================================================================================================

Results evil_eval( std::string const& outfile, Sample const& sample, Lookups const& lookup, std::string const& taxon_prefix, bool use_blacklist )
{
    std::ofstream out;
    out.open( outfile );

    // -------------------------------------------------------------------------
    //     Get Input
    // -------------------------------------------------------------------------

    out << "\n";
    out << "===========================================================================\n";
    out << "\n";
    out << "    Taxon prefix: " << taxon_prefix << "\n";
    out << "    Blacklist: " << ( use_blacklist ? "yes" : "no" ) << "\n";
    out << "\n";
    out << "===========================================================================\n";
    out << "\n";

    // -------------------------------------------------------------------------
    //     Result Storage
    // -------------------------------------------------------------------------

    // Helper function for actual distance in branch length units that consideres
    // the distance from the placement to the beginning of the branch
    // (or, if it is on the same branch, just gives 0).
    auto actual_branch_dist = [&](
        size_t const correct_edge_index, PqueryPlacement const& placement
    ){
        auto const placement_edge_index = placement.edge().index();
        if( correct_edge_index == placement_edge_index ) {
            return 0.0;
        }

        // Distance between the mid points of the edges.
        auto dist = lookup.branch_dists( correct_edge_index, placement_edge_index );

        // Subtract half the correct edge branch length to get to its beginning.
        dist -= sample.tree().edge_at( correct_edge_index ).data<DefaultEdgeData>().branch_length / 2.0;

        // Get length of branch where the placement is.
        auto const placement_bl = sample.tree().edge_at( placement_edge_index ).data<DefaultEdgeData>().branch_length;

        // Now add the half of the placement branch length, to get its full length,
        // and then subtract the proximal length to finally get the position of the placement itself.
        dist += ( placement_bl / 2.0 ) - placement.proximal_length;

        if( dist < 0.0 ) {
            LOG_WARN << "actual_branch_dist == " << dist;
        }
        return dist;
    };

    out << "Prepare result storage.\n";

    // How many of the best placements to take into account.
    size_t const n_best = 3;

    // Collection: how many edges is the best placement away from its correct edge?
    auto best_placement_edge_distances_int = std::vector<std::vector<size_t>>( n_best );

    // Collection: how far is the best placement away from its correct edge, in branch length units?
    auto best_placement_branch_distances = std::vector<std::vector<double>>( n_best );

    // Collection: how many edges are the pqueries away from the correct one?
    // This uses LWR to caluclate a weighted distance for each pquery.
    auto weighted_placements_edge_distances = std::vector<double>();

    // Collection: how far are the pqueries away from the correct edge, in branch length units?
    // This uses LWR to caluclate a weighted distance for each pquery.
    auto weighted_placements_branch_distances = std::vector<double>();

    // Collection: how many edges are the placements away from the correct one?
    auto all_placements_edge_distances = std::vector<std::pair<double,double>>();

    // Collection: how far are the placements away from the correct edge, in branch length units?
    auto all_placements_branch_distances = std::vector<std::pair<double,double>>();

    // Integer histogram: how many placements does each pquery have?
    auto num_placements_hist = std::vector<size_t>( 10, 0 );

    // Count how many of the best placements are exaclty on the right edge.
    auto correct_placements = std::vector<size_t>( n_best, 0 );

    // Collection: How far is the best placement away from the correct branch,
    // normalized by the local branch length of the correct branch.
    auto best_placement_local_branch_distances = std::vector<double>();

    // Collection: How far is the pquery away from the correct branch, using all placements, and
    // normalized with the branch length of the correct branch.
    auto all_placements_local_correct_branch_distances = std::vector<double>();

    // Collection: How far is the pquery away from the correct branch, using weighted normalization
    // of the local branches on which it was placed.
    auto all_placements_weighted_local_branch_distances = std::vector<double>();

    // Collection: How far is the pquery away from the correct branch, using the average
    // of the local branches on which it was placed.
    auto all_placements_avg_local_branch_distances = std::vector<double>();

    // Sum up distances.
    double sum_weighted_edge_distances = 0.0;
    double sum_weighted_branch_distances = 0.0;

    // LWRs
    auto all_lwrs_accu  = HistogramAccumulator();
    auto best_lwrs_accu = std::vector<HistogramAccumulator>( n_best );

    // EDPL
    HistogramAccumulator edpl_accu;
    std::vector<double> edpls;

    out << "Finished result storage.\n";

    // -------------------------------------------------------------------------
    //     Process Sample
    // -------------------------------------------------------------------------

    out << "Processing sample.\n";

    // Count how many pqueries are already processed.
    size_t total_count = 0;
    size_t blacklist_count = 0;
    size_t taxon_skip_count = 0;
    size_t no_taxon_count = 0;
    size_t processed_count = 0;

    for( auto const& pquery : sample ) {
        // if( total_count % 50000 == 0 ) {
        //     out << "Processed " << total_count << " pqueries\n";
        // }
        ++total_count;
        if( pquery.name_size() != 1 ) {
            LOG_ERR << "pquery.name_size() != 1";
            continue;
        }
        if( pquery.placement_size() < 1 ) {
            LOG_ERR << "pquery.placement_size() < 1";
            continue;
        }

        // Check if the sequence is on the blacklist.
        if( use_blacklist ) {
            auto seq_id = pquery.name_at(0).name.substr( 0, 10 );
            if( lookup.blacklist.count( seq_id ) > 0 ) {
                ++blacklist_count;
                continue;
            }
        }

        // Get the name, excluding the prefix `SEQ_000000_`.
        auto seq_name = pquery.name_at(0).name.substr( 11 );

        // Skip if not the wanted taxon.
        if( taxon_prefix != "" && ! starts_with( seq_name, taxon_prefix ) ) {
            ++taxon_skip_count;
            continue;
        }

        // Then, remove suffixes `_xyz` until we find a name that is a taxon in the tree.
        while( seq_name.size() > 0 && lookup.tip_name_to_edge_index.count( seq_name ) == 0 ) {
            size_t lastindex = seq_name.find_last_of("_");
            seq_name = seq_name.substr(0, lastindex);

            if( lastindex == std::string::npos ) {
                break;
            }
        }
        if( seq_name.size() == 0 || lookup.tip_name_to_edge_index.count( seq_name ) == 0 ) {
            // LOG_ERR << "Cannot find any taxon for sequence " << pquery.name_at(0).name.substr( 11 );
            ++no_taxon_count;
            continue;
        }

        // Now we have everything and can savely do all the counting.
        ++processed_count;
        auto const correct_edge_index = lookup.tip_name_to_edge_index.at( seq_name );

        // Collect best placement values.
        auto const min_n = std::min( n_best, pquery.placement_size() );
        for( size_t i = 0; i < min_n; ++i ) {
            best_placement_edge_distances_int[i].push_back(
                lookup.edge_dists( correct_edge_index, pquery.placement_at(i).edge().index() )
            );
            best_placement_branch_distances[i].push_back(
                actual_branch_dist( correct_edge_index, pquery.placement_at(i) )
            );
            best_lwrs_accu[i].increment( pquery.placement_at(i).like_weight_ratio );
        }

        // Local dist of best placement.
        auto const best_bd = actual_branch_dist( correct_edge_index, pquery.placement_at(0) );
        auto const corr_bl = sample.tree().edge_at( correct_edge_index ).data<DefaultEdgeData>().branch_length;
        best_placement_local_branch_distances.push_back( best_bd / corr_bl );

        // Collect sums of distances for the pquery to calculate their average.
        double weighted_sum_ed = 0.0;
        double weighted_sum_bd = 0.0;

        // Collect local weighted and avg branch length.
        double local_weighted_bl = 0.0;
        double local_avg_bl = 0.0;

        // Collect all placement values.
        for( auto const& placement : pquery.placements() ) {
            auto const ed = lookup.edge_dists( correct_edge_index, placement.edge().index() );
            auto const bd = actual_branch_dist( correct_edge_index, placement );

            all_placements_edge_distances.push_back({   ed, placement.like_weight_ratio });
            all_placements_branch_distances.push_back({ bd, placement.like_weight_ratio });

            sum_weighted_edge_distances   += static_cast<double>( ed ) * placement.like_weight_ratio;
            sum_weighted_branch_distances += bd * placement.like_weight_ratio;

            weighted_sum_ed += static_cast<double>( ed ) * placement.like_weight_ratio;
            weighted_sum_bd += bd * placement.like_weight_ratio;

            auto const bl = sample.tree().edge_at( placement.edge().index() ).data<DefaultEdgeData>().branch_length;
            local_weighted_bl += placement.like_weight_ratio * bl;
            local_avg_bl      += bl;

            all_lwrs_accu.increment( placement.like_weight_ratio );
        }
        local_avg_bl /= static_cast<double>( pquery.placement_size() );

        // Increment histograms and counters.
        weighted_placements_edge_distances.push_back( weighted_sum_ed );
        weighted_placements_branch_distances.push_back( weighted_sum_bd );
        all_placements_local_correct_branch_distances.push_back( weighted_sum_bd / corr_bl );
        all_placements_weighted_local_branch_distances.push_back( weighted_sum_bd / local_weighted_bl );
        all_placements_avg_local_branch_distances.push_back( weighted_sum_bd / local_avg_bl );

        ++num_placements_hist[ pquery.placement_size() < 10 ? pquery.placement_size() : 9 ];
        for( size_t i = 0; i < min_n; ++i ) {
            if( pquery.placement_at(i).edge().index() == correct_edge_index ) {
                ++correct_placements[i];
            }
        }

        auto const edplv = edpl( pquery, lookup.node_dists );
        edpl_accu.increment( edplv );
        edpls.push_back( edplv );
    }

    out << "Finished processing " << total_count << " pqueries.\n";
    out << "Processed " << processed_count << " or "
             << ( 100.0 * static_cast<double>( processed_count ) / static_cast<double>( total_count )) << "% pqueries." << "\n";
    out << "Skipped " << blacklist_count << " or "
             << ( 100.0 * static_cast<double>( blacklist_count ) / static_cast<double>( total_count )) << "% blacklisted pqueries" << "\n";
    if( taxon_prefix != "" ) {
        out << "Skipped " << taxon_skip_count << " or "
                 << ( 100.0 * static_cast<double>( taxon_skip_count ) / static_cast<double>( total_count ))
                 << "% pqueries that are not " << taxon_prefix << "\n";
    }
    out << "Found " << no_taxon_count << " or "
             << ( 100.0 * static_cast<double>( no_taxon_count ) / static_cast<double>( total_count )) << "% no taxon pqueries." << "\n";
    out << "processed_count + blacklist_count + taxon_skip_count + no_taxon_count = "
             << ( processed_count + blacklist_count + taxon_skip_count + no_taxon_count )
             << " should be equal to total_count = " << total_count << "\n";
    if( processed_count + blacklist_count + taxon_skip_count + no_taxon_count != total_count ) {
        LOG_ERR << "counts went wrong!";
        return {};
    } else {
        out << "it is! phew\n";
    }
    out << "\n";

    // out << "Calculating EDPLs\n";
    // auto const exdpl_vec = edpl( sample );
    // for( auto v : exdpl_vec ) {
    //     edpl_accu.increment(v);
    // }
    // out << "Finished EDPLs\n";
    // out << "\n";

    // -------------------------------------------------------------------------
    //     Statistics
    // -------------------------------------------------------------------------

    size_t const hist_resolution = 20;

    // Sanity.
    if( total_count != sample.size() ) {
        LOG_ERR << "Count " << total_count << " != sample.size() " << sample.size();
    }

    out << "Prepare statistics.\n";

    // Make a double vector to be able to use statistics functions. Ugly, but okay for now.
    auto best_placement_edge_distances = std::vector<std::vector<double>>( n_best );
    auto best_placement_edge_distance_hist = std::vector<std::vector<size_t>>( n_best );
    for( size_t i = 0; i < n_best; ++i ) {
        best_placement_edge_distances[i] = std::vector<double>(
            best_placement_edge_distances_int[i].begin(), best_placement_edge_distances_int[i].end()
        );

        // Integer histogram: how many edges is the best placement away from its correct edge?
        best_placement_edge_distance_hist[i] = std::vector<size_t>( sample.tree().edge_count(), 0 );
        for( auto e : best_placement_edge_distances_int[i] ) {
            ++( best_placement_edge_distance_hist[i][e] );
        }
        while( best_placement_edge_distance_hist[i].size() > 0 && best_placement_edge_distance_hist[i].back() == 0 ) {
            best_placement_edge_distance_hist[i].pop_back();
        }
    }

    auto const avg_branch_length = tree::length( sample.tree() ) / static_cast<double>( sample.tree().edge_count() );

    auto const hist_lwr_all = all_lwrs_accu.build_uniform_ranges_histogram( hist_resolution, 0.0, 1.0 );
    auto const hist_edpl    = edpl_accu.build_uniform_ranges_histogram( hist_resolution, 0.0, hist_resolution * avg_branch_length );

    auto hist_lwrs    = std::vector<Histogram>();
    for( size_t i = 0; i < best_lwrs_accu.size(); ++i ) {
        hist_lwrs.push_back( best_lwrs_accu[i].build_uniform_ranges_histogram( hist_resolution, 0.0, 1.0 ) );
    }

    // count the number of tree taxa with taxon_prefix
    size_t tree_taxon_cnt = 0;
    for( auto const& node : sample.tree().nodes() ) {
        auto const name = node->data<tree::DefaultNodeData>().name;
        if( taxon_prefix != "" && ! starts_with( name, taxon_prefix ) ) {
            ++tree_taxon_cnt;
        }
    }

    out << "Finished preparing statistics.\n";
    out << "\n";
    out << "How often was the nth placement on the correct edge?\n";
    for( size_t i = 0; i < n_best; ++i ) {
        out << "Correct pqueries " << i << ": " << correct_placements[i] << " / " << processed_count << " = "
        << (static_cast<double>(correct_placements[i]) / static_cast<double>(processed_count)) << "\n";
    }
    out << "\n";

    out << "tree node_count " << sample.tree().node_count() << "\n";
    out << "tree edge_count " << sample.tree().edge_count() << "\n";
    out << "tree leaf count " << leaf_node_count( sample.tree() ) << "\n";
    out << "tree taxa with taxon_prefix names " << tree_taxon_cnt << "\n";
    out << "\n";

    out << "avg_branch_length of ref tree " << avg_branch_length << "\n";
    out << "avg_weighted_edge_distances   " << ( sum_weighted_edge_distances   / static_cast<double>(processed_count) ) << "\n";
    out << "avg_weighted_branch_distances " << ( sum_weighted_branch_distances / static_cast<double>(processed_count) ) << "\n";
    out << "These two numbers are the averages of the weighted distances from a ppquery to the correct edge.\n"
             << "That is, for all its placements, distance * lwr." << "\n";
    out << "\n";

    // -------------------------------------------------------------------------
    //     Generate Output
    // -------------------------------------------------------------------------

    auto print_stats = [&]( std::vector<double>& data )
    {
        std::sort( data.begin(), data.end() );
        auto const ms  = mean_stddev( data );
        auto const q   = quartiles(   data );
        auto const sum = std::accumulate(data.begin(), data.end(), 0.0);

        out << "entries\t" << data.size() << "\n";
        out << "sum\t" << sum << "\n";
        out << "mean\t" << ms.mean << "\n";
        out << "stddev\t" << ms.stddev << "\n";
        out << "q0\t" << q.q0 << "\n";
        out << "q1\t" << q.q1 << "\n";
        out << "q2\t" << q.q2 << "\n";
        out << "q3\t" << q.q3 << "\n";
        out << "q4\t" << q.q4 << "\n";
        out << "\n";
    };

    auto print_histogram = [&]( Histogram const& hist, double const total )
    {
        auto const hist_sum = sum(hist);
        double accumulated = 0.0;

        out << "total\t" << total << "\n";
        out << "sum\t" << hist_sum << "\n";
        out << "mean\t" << mean(hist) << "\n";
        out << "stddev\t" << sigma(hist) << "\n";

        out << "Bin" << "\t" << "Range" << "\t" << "Range Start" << "\t" << "Range End"
                 << "\t" << "Bin Name" << "\t" << "Value" << "\t" << "Percentage" << "\t" << "Acc Val" << "\t" << "Acc Perc" << "\n";
        for( size_t i = 0; i < hist.bins(); ++i ) {
            auto range = hist.bin_range( i );
            accumulated += hist[i];
            out << i
                     << "\t" << "\"[" << range.first << ", " << range.second << ")\""
                     << "\t" << range.first << "\t" << range.second
                     << "\t" << ">= " << range.first
                     << "\t" << hist[ i ]
                     << "\t" << ( hist[i] / hist_sum )
                     << "\t" << accumulated
                     << "\t" << ( accumulated/ hist_sum )
                     << "\n";
        }

        out << "\n";
    };

    /**
     * @brief Print Stats and Histogram for a vec of doubles, with a manual hist max value.
     */
    auto print_stats_best = [&]( std::vector<double>& data, double const hist_max )
    {
        auto const ha = HistogramAccumulator( data );
        auto const hi = ha.build_uniform_ranges_histogram( hist_resolution, 0.0, hist_max );

        print_stats( data );
        print_histogram( hi, data.size() );
    };

    /**
     * @brief Print Stats and Histogram for pairs of distance and LWR (in that order).
     */
    auto print_stats_all = [&]( std::vector<std::pair<double,double>> const& data, double const hist_max )
    {
        std::vector<double> values;
        double total = 0.0;
        for( auto d : data ) {
            values.push_back( d.first * d.second );
            total += d.second;
        }

        // Histogram of data
        auto const had = HistogramAccumulator( data );
        auto const hid = had.build_uniform_ranges_histogram( hist_resolution, 0.0, hist_max );

        // Histogram of values
        // auto const hav = HistogramAccumulator( values );
        // auto const hiv = had.build_uniform_ranges_histogram( hist_resolution, 0.0, hist_max );

        print_stats( values );
        print_histogram( hid, total );
        // print_histogram( hiv, total );
    };

    out << "--------------------------------------------------\n";
    out << "\n";

    out << "weighted_placements_edge_distances\n";
    out << "how many edges are the pqueries away from the correct one?\n";
    out << "This uses LWR to caluclate a weighted distance for each pquery.\n";
    out << "\n";
    print_stats_best( weighted_placements_edge_distances, hist_resolution );

    out << "weighted_placements_branch_distances\n";
    out << "how far are the pqueries away from the correct edge, in branch length units?\n";
    out << "This uses LWR to caluclate a weighted distance for each pquery.\n";
    out << "\n";
    print_stats_best( weighted_placements_branch_distances, hist_resolution * avg_branch_length );

    out << "--------------------------------------------------\n";
    out << "\n";

    out << "best_placement_edge_distances\n";
    out << "how many edges is the best placement away from its correct edge?\n";
    out << "\n";
    for( size_t i = 0; i < n_best; ++i ) {
        out << "best_placement_edge_distances " << i << "\n";
        print_stats_best( best_placement_edge_distances[i], hist_resolution );
    }

    out << "--------------------------------------------------\n";
    out << "\n";

    out << "best_placement_branch_distances\n";
    out << "how far is the best placement away from its correct edge, in branch length units?\n";
    out << "\n";
    for( size_t i = 0; i < n_best; ++i ) {
        out << "best_placement_branch_distances " << i << "\n";
        print_stats_best( best_placement_branch_distances[i], hist_resolution * avg_branch_length );
    }

    out << "--------------------------------------------------\n";
    out << "\n";

    out << "all_placements_edge_distances\n";
    out << "how many edges are the placements away from the correct one?\n";
    out << "\n";
    print_stats_all( all_placements_edge_distances, hist_resolution );

    out << "--------------------------------------------------\n";
    out << "\n";

    out << "all_placements_branch_distances\n";
    out << "how far are the placements away from the correct edge, in branch length units?\n";
    out << "\n";
    print_stats_all( all_placements_branch_distances, hist_resolution * avg_branch_length );

    // out << "--------------------------------------------------\n";
    // out << "\n";
    //
    // // Already covered in best_placement_edge_distances
    //
    // out << "best_placement_edge_distance_hist\n";
    // out << "how many edges is the best placement away from its correct edge?\n";
    // out << "\n";
    // for( size_t i = 0; i < n_best; ++i ) {
    //     out << "best_placement_edge_distance_hist " << i;
    //     for( size_t j = 0; j < best_placement_edge_distance_hist.size(); ++j ) {
    //         out << j << "\t" << best_placement_edge_distance_hist[i][j];
    //     }
    //     out << "\n";
    // }

    out << "--------------------------------------------------\n";
    out << "\n";

    out << "best_placement_local_branch_distances\n";
    out << "How far is the best placement away from the correct branch, normalized by the local branch length of the correct branch.?\n";
    out << "\n";
    print_stats_best( best_placement_local_branch_distances, hist_resolution );

    out << "--------------------------------------------------\n";
    out << "\n";

    out << "all_placements_local_correct_branch_distances\n";
    out << "How far is the pquery away from the correct branch, using all placements, and normalized with the branch length of the correct branch?\n";
    out << "\n";
    print_stats_best( all_placements_local_correct_branch_distances, hist_resolution );

    out << "--------------------------------------------------\n";
    out << "\n";

    out << "all_placements_weighted_local_branch_distances\n";
    out << "How far is the pquery away from the correct branch, using weighted normalization of the local branches on which it was placed?\n";
    out << "\n";
    print_stats_best( all_placements_weighted_local_branch_distances, hist_resolution / 2 );

    out << "--------------------------------------------------\n";
    out << "\n";

    out << "all_placements_avg_local_branch_distances\n";
    out << "How far is the pquery away from the correct branch, using the average of the local branches on which it was placed?\n";
    out << "\n";
    print_stats_best( all_placements_avg_local_branch_distances, hist_resolution / 2 );

    out << "--------------------------------------------------\n";
    out << "\n";

    out << "num_placements_hist\n";
    out << "how many placements does each pquery have?\n";
    out << "\n";
    for( size_t i = 0; i < num_placements_hist.size(); ++i ) {
        out << i << "\t" << num_placements_hist[i] << "\n";
    }
    out << "\n";

    out << "--------------------------------------------------\n";
    out << "\n";

    auto print_hist_line = [&]( std::ostringstream& out, Histogram const& hist, size_t i ) {
        auto range = hist.bin_range( i );
        out << i;
        out << "\t" << "\"[" << range.first << ", " << range.second << ")\"";
        out << "\t" << range.first << "\t" << range.second;
        out << "\t" << ">= " << range.first;
        out << "\t" << hist[ i ];
    };

    out << "hist_lwr_all and hist_lwrs\n";
    out << "what is the distrubtion of lwrs in total and per nth best placement?\n";
    out << "\n";
    auto print_lwrs_hist = [&](){
        std::ostringstream s;
        s << "Bin" << "\t" << "Range" << "\t" << "Range Start" << "\t" << "Range End"
          << "\t" << "Bin Name" << "\t" << "LWR Total";
        for( size_t i = 0; i < hist_lwrs.size(); ++i ) {
            s << "\t" << "LWR " << ( i + 1 );
        }
        s << "\n";

        for( size_t i = 0; i < hist_resolution; ++i ) {
            print_hist_line( s, hist_lwr_all, i );
            for( size_t h = 0; h < hist_lwrs.size(); ++h ) {
                s << "\t" << hist_lwrs[ h ][ i ];
            }
            s << "\n";
        }

        out << s.str();
    };
    print_lwrs_hist();

    out << "\n";
    out << "--------------------------------------------------\n";
    out << "\n";

    auto const avg_edpl = std::accumulate( edpls.begin(), edpls.end(), 0.0) / static_cast<double>( edpls.size() );
    out << "average edpl: " << avg_edpl << "\n\n";

    out << "hist_edpl\n";
    out << "what is the distribution of edpls for the whol sample?\n";
    out << "\n";

    out << "Bin" << "\t" << "Range" << "\t" << "Range Start" << "\t" << "Range End" << "\t"
             << "Bin Name" << "\t" << "EDPL\n";
    for( size_t i = 0; i < hist_edpl.bins(); ++i ) {
        std::ostringstream s;
        print_hist_line( s, hist_edpl, i );
        out << s.str() << "\n";
    }
    out << "\n";

    // Return the most important results.
    return {
        "",
        weighted_placements_edge_distances,
        weighted_placements_branch_distances
    };
}

// =================================================================================================
//      Process Sample
// =================================================================================================

std::pair<Results, Results> process_sample( std::string const& jplace_file, std::string const& taxon_prefix, std::string const& name )
{
    // Read input file
    LOG_INFO << "Reading " << jplace_file;
    auto jplace_reader = JplaceReader();
    jplace_reader.invalid_number_behaviour(
        genesis::placement::JplaceReader::InvalidNumberBehaviour::kCorrect
    );
    auto sample = jplace_reader.from_file( jplace_file );
    sort_placements_by_weight( sample );
    LOG_INFO << "Done reading.";

    // Prepare all importatn lookup tables and lists.
    auto const lookup = preprocess( sample );
    auto const path = utils::file_path( jplace_file );

    LOG_INFO << "Start eval.";
    auto res_no_bl = evil_eval( path + "/silva_tree_eval_no-tax_no-blacklist.log", sample, lookup, "", false );
    auto res_bl    = evil_eval( path + "/silva_tree_eval_no-tax_blacklist.log", sample, lookup, "", true );
    if( taxon_prefix != "" ) {
        res_no_bl = evil_eval( path + "/silva_tree_eval_tax_no-blacklist.log", sample, lookup, taxon_prefix, false );
        res_bl    = evil_eval( path + "/silva_tree_eval_tax_blacklist.log", sample, lookup, taxon_prefix, true );
    }
    LOG_INFO << "Finished eval.";

    res_no_bl.name = name;
    res_bl.name    = name;
    return { res_no_bl, res_bl };
}

// =================================================================================================
//      Collector Function
// =================================================================================================

/**
 * @brief Takes the results from some eval runs and makes combined csv tables out of them,
 * which are then used for visualization.
 */
void collector( std::vector<Results> const& results, std::string const& path_prefix, std::string const& path_suffix )
{
    if( results.size() == 0 ) {
        return;
    }

    // Resolution of the histogram for visualization.
    size_t const hist_resolution = 200;

    // How many edges and branch length units do we want to visualize?
    // That is, max x of the graph.
    double const edge_dist_max = 10.0;
    double const branch_dist_max = 1.0;

    // Histograms of all results.
    auto hist_edge   = std::vector<Histogram>();
    auto hist_branch = std::vector<Histogram>();

    // For each result, sum of all entries in the data.
    // As we truncate the histograms at the max values, the sum() function for histograms
    // would yield a too small value, so we have to store this.
    std::vector<double> sum_edge;
    std::vector<double> sum_branch;

    // Fill the histograms.
    for( size_t i = 0; i < results.size(); ++i ) {
        auto const& res = results[i];

        sum_edge.push_back( res.weighted_placements_edge_distances.size() );
        hist_edge.push_back( Histogram( hist_resolution, 0.0, edge_dist_max ));
        hist_edge[i].out_of_range_behaviour( Histogram::OutOfRangeBehaviour::kIgnore );
        for( auto const& elem : res.weighted_placements_edge_distances ) {
            hist_edge[i].increment( elem );
        }

        sum_branch.push_back( res.weighted_placements_branch_distances.size() );
        hist_branch.push_back( Histogram( hist_resolution, 0.0, branch_dist_max ));
        hist_branch[i].out_of_range_behaviour( Histogram::OutOfRangeBehaviour::kIgnore );
        for( auto const& elem : res.weighted_placements_branch_distances ) {
            hist_branch[i].increment( elem );
        }
    }

    // Write out results to tab separated files
    auto write_res_hists = [&]( std::string const& outfile, std::vector<Histogram> const& hists, std::vector<double> sums ){
        std::ofstream out;
        out.open( outfile );

        // Safety
        if( hists.size() != results.size() ) {
            throw std::runtime_error( "hists.size() != results.size()" );
        }
        if( hists.size() != sums.size() ) {
            throw std::runtime_error( "hists.size() != sums.size()" );
        }

        // Header
        out << "Bin Start\tBin End";
        for( auto const& res : results ) {
            out << "\t" << res.name;
        }
        out << "\n";

        // Lines
        auto accumulations = std::vector<double>( hists.size(), 0.0 );
        for( size_t b = 0; b < hist_resolution; ++b ) {
            auto range = hists[0].bin_range( b );
            out << utils::to_string_precise( range.first,  3 ) << "\t";
            out << utils::to_string_precise( range.second, 3 );

            for( size_t j = 0; j < hists.size(); ++j ) {
                accumulations[j] += hists[j][b];
                out << "\t" << ( accumulations[j] / sums[j] );
            }
            out << "\n";
        }
    };

    write_res_hists( path_prefix + "tables/weighted_placements_edge_distances"   + path_suffix + ".csv" ,   hist_edge,   sum_edge );
    write_res_hists( path_prefix + "tables/weighted_placements_branch_distances" + path_suffix + ".csv", hist_branch, sum_branch );

    auto write_res_lists = []( std::string const& outfile, std::vector<double> const& list ){
        std::ofstream out;
        out.open( outfile );

        for( auto const& elem : list ) {
            out << elem << "\n";
        }
    };

    for( size_t i = 0; i < results.size(); ++i ) {
        write_res_lists(
            path_prefix + "lists/weighted_placements_edge_distances" + path_suffix + "_" + results[i].name + ".csv",
            results[i].weighted_placements_edge_distances
        );
        write_res_lists(
            path_prefix + "lists/weighted_placements_branch_distances" + path_suffix + "_" + results[i].name + ".csv",
            results[i].weighted_placements_branch_distances
        );
    }
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

    // Options::get().allow_file_overwriting(true);

    // // Check if the command line contains the right number of arguments.
    // if (argc != 2 && argc != 3) {
    //     throw std::runtime_error(
    //         "Need to provide jplace file and optional taxonomic prefix for filtering.\n"
    //     );
    // }

    // auto const jplace_file  = std::string( argv[1] );
    // auto const taxon_prefix = ( argc == 3 ? std::string( argv[2] ) : std::string() );
    // process_sample( jplace_file, taxon_prefix, "" );

    // Load new base dir if given.
    if( argc != 3 ) {
        throw std::runtime_error( "expects in and out dir as params" );
    }

    // std::string basedir; // = "/path/to/results/01_backbone/11_tree_eval/";
    std::string basedir = utils::dir_normalize_path( argv[1] );
    std::string outdir  = utils::dir_normalize_path( argv[2] );


    auto res_unconstr_no_bl = std::vector<Results>();
    auto res_unconstr_bl    = std::vector<Results>();
    auto res_constr_no_bl   = std::vector<Results>();
    auto res_constr_bl      = std::vector<Results>();

    struct Job
    {
        std::string jplace_file;
        std::string taxon_prefix;
        std::string name;
    };

    // Old dir names
    // std::vector<Job> jobs_unconstr {
    //     { basedir + "Gen_Unconstr_backbone/epa_result.jplace", "", "General" },
    //     { basedir + "Arch_Unconstr_backbone/epa_result.jplace", "Arch", "Archaea" },
    //     { basedir + "Bact_Unconstr_backbone/epa_result.jplace", "Bact", "Bacteria" },
    //     { basedir + "Euks_Unconstr_backbone/epa_result.jplace", "Euk", "Eukaryota" }
    // };
    // std::vector<Job> jobs_constr {
    //     { basedir + "Gen_Constr_backbone/epa_result.jplace", "", "General" },
    //     { basedir + "Arch_Constr_backbone/epa_result.jplace", "Arch", "Archaea" },
    //     { basedir + "Bact_Constr_backbone/epa_result.jplace", "Bact", "Bacteria" },
    //     { basedir + "Euks_Constr_backbone/epa_result.jplace", "Euk", "Eukaryota" }
    // };

    // Names for /path/to/results/11_consensus_seqs/05_accuracy_viz
    std::vector<Job> jobs_unconstr {
        { basedir + "General_Unconstr/epa_result.jplace", "", "General" },
        { basedir + "Archaea_Unconstr/epa_result.jplace", "Arch", "Archaea" },
        { basedir + "Bacteria_Unconstr/epa_result.jplace", "Bact", "Bacteria" },
        { basedir + "Eukaryota_Unconstr/epa_result.jplace", "Euk", "Eukaryota" }
    };
    std::vector<Job> jobs_constr {
        { basedir + "General_Constr/epa_result.jplace", "", "General" },
        { basedir + "Archaea_Constr/epa_result.jplace", "Arch", "Archaea" },
        { basedir + "Bacteria_Constr/epa_result.jplace", "Bact", "Bacteria" },
        { basedir + "Eukaryota_Constr/epa_result.jplace", "Euk", "Eukaryota" }
    };

    // Names for just the opt bact tree, in order to test whether opt gives better results.
    // std::vector<Job> jobs_unconstr {
    //     { basedir + "epa_result.jplace", "Bact", "Bacteria" }
    // };
    // std::vector<Job> jobs_constr {
    // };

    // Names for /path/to/results/02_multilevel/04_epa
    // std::vector<Job> jobs_unconstr {
    //     { basedir + "Actinobacteria_Unconstr/epa_result.jplace", "Bacteria_Actinobacteria_", "Actinobacteria" },
    //     { basedir + "Cyanobacteria_Unconstr/epa_result.jplace", "Bacteria_Cyanobacteria_", "Cyanobacteria" },
    //     { basedir + "Proteobacteria_Unconstr/epa_result.jplace", "Bacteria_Proteobacteria_", "Proteobacteria" },
    //     { basedir + "Firmicutes_Unconstr/epa_result.jplace", "Bacteria_Firmicutes_", "Firmicutes" },
    //     { basedir + "Bacteroidetes_Unconstr/epa_result.jplace", "Bacteria_Bacteroidetes_", "Bacteroidetes" }
    // };
    // std::vector<Job> jobs_constr {
    //     { basedir + "Actinobacteria_Constr/epa_result.jplace", "Bacteria_Actinobacteria_", "Actinobacteria" },
    //     { basedir + "Cyanobacteria_Constr/epa_result.jplace", "Bacteria_Cyanobacteria_", "Cyanobacteria" },
    //     { basedir + "Proteobacteria_Constr/epa_result.jplace", "Bacteria_Proteobacteria_", "Proteobacteria" },
    //     { basedir + "Firmicutes_Constr/epa_result.jplace", "Bacteria_Firmicutes_", "Firmicutes" },
    //     { basedir + "Bacteroidetes_Constr/epa_result.jplace", "Bacteria_Bacteroidetes_", "Bacteroidetes" }
    // };

    LOG_INFO << "==================================================";
    LOG_INFO << "Collecting Unconstr Data";
    LOG_INFO << "==================================================";

    for( size_t i = 0; i < jobs_unconstr.size(); ++i ) {
        auto const& job = jobs_unconstr[i];

        if( ! utils::file_exists( job.jplace_file ) ) {
            LOG_WARN << "Skipping eval dir " << job.jplace_file;
            continue;
        }

        auto const res = process_sample( job.jplace_file, job.taxon_prefix, job.name );
        res_unconstr_no_bl.push_back( res.first );
        res_unconstr_bl.push_back( res.second );
    }

    // res_unconstr.push_back( process_sample( basedir + "Gen_Unconstr_backbone/epa_result.jplace", "", "General" ));
    // res_unconstr.push_back( process_sample( basedir + "Arch_Unconstr_backbone/epa_result.jplace", "Arch", "Archaea" ));
    // res_unconstr.push_back( process_sample( basedir + "Bact_Unconstr_backbone/epa_result.jplace", "Bact", "Bacteria" ));
    // res_unconstr.push_back( process_sample( basedir + "Euks_Unconstr_backbone/epa_result.jplace", "Euk", "Eukaryota" ));

    LOG_INFO << "==================================================";
    LOG_INFO << "Collecting Constr Data";
    LOG_INFO << "==================================================";

    for( size_t i = 0; i < jobs_constr.size(); ++i ) {
        auto const& job = jobs_constr[i];

        if( ! utils::file_exists( job.jplace_file ) ) {
            LOG_WARN << "Skipping eval dir " << job.jplace_file;
            continue;
        }

        auto const res = process_sample( job.jplace_file, job.taxon_prefix, job.name );
        res_constr_no_bl.push_back( res.first );
        res_constr_bl.push_back( res.second );
    }

    // res_constr.push_back( process_sample( basedir + "Gen_Constr_backbone/epa_result.jplace", "", "General" ));
    // res_constr.push_back( process_sample( basedir + "Arch_Constr_backbone/epa_result.jplace", "Arch", "Archaea" ));
    // res_constr.push_back( process_sample( basedir + "Bact_Constr_backbone/epa_result.jplace", "Bact", "Bacteria" ));
    // res_constr.push_back( process_sample( basedir + "Euks_Constr_backbone/epa_result.jplace", "Euk", "Eukaryota" ));

    LOG_INFO << "==================================================";
    LOG_INFO << "Writing Results";
    LOG_INFO << "==================================================";

    utils::dir_create( outdir + "viz" );
    utils::dir_create( outdir + "viz/lists" );
    collector( res_unconstr_no_bl, outdir + "viz/", "_unconstr_no-blacklist" );
    collector( res_unconstr_bl,    outdir + "viz/", "_unconstr_blacklist" );
    collector( res_constr_no_bl,   outdir + "viz/", "_constr_no-blacklist" );
    collector( res_constr_bl,      outdir + "viz/", "_constr_blacklist" );

    LOG_INFO << "Finished";
    return 0;
}
