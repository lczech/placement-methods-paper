/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2018 Lucas Czech and HITS gGmbH

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
#include <cmath>
#include <string>
#include <unordered_map>
#include <unordered_set>

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

// =================================================================================================
//     Main
// =================================================================================================

int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.date = true;
    utils::Logging::details.time = true;
    utils::Options::get().number_of_threads( 4 );
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    // File paths.
    // std::string in_jplace = "03_bv/03_epa/orig_queries_jplace_clean_15";
    std::string in_jplace = "03_bv/03_epa/orig_queries_jplace_clean";
    std::string in_meta = "data/bv/meta/meta_simple.csv";
    std::string out_dir = "philr_bv/";

    LOG_INFO << "Reading data";
    auto dfr = DataframeReader<double>();
    dfr.csv_reader().separator_chars( "\t" );
    auto meta = dfr.read( from_file( in_meta ));
    LOG_INFO << "meta: " << meta.rows() << " x " << meta.cols();
    auto samples = JplaceReader().read( from_files( utils::dir_list_files( in_jplace, ".*\\.jplace" )));
    LOG_INFO << "samples: " << samples.size();

    LOG_INFO << "Converting to Mass Trees";
    auto trees = convert_sample_set_to_mass_trees( samples, false ).first;

    for( size_t i = 0; i < trees.size(); ++i ) {
        // LOG_DBG << "tree " << i << " " << trees.name_at(i);
        make_rooted( trees[i] );
        ladderize( trees[i] );
    }

    LOG_INFO << "Computing Balance data";
    BalanceSettings settings;
    settings.pseudo_count_summand_all = 0.65;
    auto const bal_data = mass_balance_data( trees, settings );

    // {`
    //     // LOG_DBG << "taxon weights " << utils::join( bal_data.taxon_weights );
    //
    //     auto lm = LayoutParameters();
    //     lm.ladderize = false;
    //     lm.stroke.width = 10;
    //     auto cm = ColorMap( color_list_viridis() );
    //     cm.mask_color( utils::color_from_hex( "#dfdfdf" ));
    //     write_color_tree_to_svg_file(
    //         bal_data.tree, lm, bal_data.taxon_weights, cm, ColorNormalizationLinear(),
    //         out_dir + "taxon_weights.svg"
    //     );
    // }`

    LOG_INFO << "Computing Phylo ILR";
    auto const ilr = phylogenetic_ilr_transform( bal_data, true );
    LOG_INFO << "ilr: " << ilr.rows() << " x " << ilr.cols();

    LOG_INFO << "Computing Edge Balances";
    auto const bals = edge_balances( bal_data, true );
    LOG_INFO << "bals: " << bals.rows() << " x " << bals.cols();
    assert( bals.rows() == trees.size() );
    assert( bals.cols() == trees[0].edge_count() );

    for( auto const& col : meta ) {
        LOG_INFO << "meta " << col.name();

        // Copy the meta data in the correct sample order.
        std::vector<double> meta_vals( trees.size(), 0.0 );
        for( size_t i = 0; i < trees.size(); ++i ) {
            meta_vals[i] = col.as<double>()[ trees.name_at(i) ];
        }

        // Need to have one meta value per irl and bal row (ie per sample),
        // and a col for every node or edge.
        assert( meta_vals.size() == ilr.rows() );
        assert( meta_vals.size() == bals.rows() );
        assert( bal_data.tree.node_count() == ilr.cols() );
        assert( bal_data.tree.edge_count() == bals.cols() );

        std::vector<double> edge_vals_ilr( bal_data.tree.edge_count(), 0.0 );
        for( size_t i = 0; i < bal_data.tree.node_count(); ++i ) {
            auto corr = spearmans_rank_correlation_coefficient( meta_vals, ilr.col(i) );
            edge_vals_ilr[ bal_data.tree.node_at(i).primary_edge().index() ] = corr;
        }

        double max_corr = -1.0;
        size_t max_edge = 0;
        double min_corr = 1.0;
        size_t min_edge = 0;

        std::vector<double> edge_vals_bal( bal_data.tree.edge_count(), 0.0 );
        for( size_t i = 0; i < bal_data.tree.edge_count(); ++i ) {
            // auto corr = spearmans_rank_correlation_coefficient( meta_vals, bals.col(i) );
            auto corr = pearson_correlation_coefficient( meta_vals, bals.col(i) );
            edge_vals_bal[ i ] = corr;

            if( corr > max_corr ) {
                max_corr = corr;
                max_edge = i;
            }
            if( corr < min_corr ) {
                min_corr = corr;
                min_edge = i;
            }
        }

        auto lm = LayoutParameters();
        lm.ladderize = false;
        lm.stroke.width = 10;
        auto cm = ColorMap( color_list_spectral() );
        cm.mask_color( utils::color_from_hex( "#dfdfdf" ));

        auto colors = cm( ColorNormalizationDiverging(), edge_vals_bal );
        // colors[max_edge] = utils::Color();

        write_color_tree_to_svg_file(
            bal_data.tree, lm, edge_vals_ilr, cm, ColorNormalizationDiverging(),
            out_dir + "philr_corr_" + col.name() + ".svg"
        );
        write_color_tree_to_svg_file(
            bal_data.tree, lm, colors, cm, ColorNormalizationDiverging(),
            out_dir + "edgebals_corr_" + col.name() + ".svg"
        );

        // Min max correlation tests.
        auto max_colors = colors;
        max_colors[max_edge] = utils::Color();
        write_color_tree_to_svg_file(
            bal_data.tree, lm, max_colors, cm, ColorNormalizationDiverging(),
            out_dir + "edgebals_corr_" + col.name() + "_max.svg"
        );
        auto min_colors = colors;
        min_colors[min_edge] = utils::Color();
        write_color_tree_to_svg_file(
            bal_data.tree, lm, min_colors, cm, ColorNormalizationDiverging(),
            out_dir + "edgebals_corr_" + col.name() + "_min.svg"
        );

        auto minmax_cols = utils::Matrix<double>( meta_vals.size(), 3 );
        minmax_cols.col(0) = meta_vals;
        minmax_cols.col(1) = bals.col(min_edge);
        minmax_cols.col(2) = bals.col(max_edge);
        utils::MatrixWriter<double>().to_file(
            minmax_cols, out_dir + "minmax_corr_" + col.name() + ".csv",
            trees.names(), { "Nugent", "Min", "Max" }, "Sample"
        );
    }

    LOG_INFO << "balance pca";
    auto const pca = utils::principal_component_analysis(
        bals, 5, utils::PcaStandardization::kCovariance
    );
    assert( pca.eigenvalues.size()  == 5 );
    assert( pca.eigenvectors.rows() == trees[0].edge_count() );
    assert( pca.eigenvectors.cols() == 5 );
    assert( pca.projection.rows()   == trees.size() );
    assert( pca.projection.cols()   == 5 );

    utils::MatrixWriter<double>().to_file(
        pca.eigenvectors, out_dir + "eigenvectors.csv"
    );
    utils::MatrixWriter<double>().to_file(
        pca.projection, out_dir + "projection.csv", trees.names()
    );

    for( size_t i = 0; i < 5; ++i ) {
        auto lm = LayoutParameters();
        lm.ladderize = false;
        lm.stroke.width = 10;
        auto cm = ColorMap( color_list_spectral() );
        cm.mask_color( utils::color_from_hex( "#dfdfdf" ));
        auto cn = ColorNormalizationDiverging( pca.eigenvectors.col(i) );
        cn.make_centric();

        write_color_tree_to_svg_file(
            bal_data.tree, lm, pca.eigenvectors.col(i), cm, cn,
            out_dir + "eigenvector_" + std::to_string(i) + ".svg"
        );
    }

    LOG_INFO << "Finished";

    return 0;
}
