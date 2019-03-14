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
    utils::Options::get().number_of_threads( 40 );
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    // File paths.
    std::string in_jplace = "path/to/06_samples/hmp_Bact_Unconstr_backbone/samples_jplace_clean/";
    std::string in_meta = "path/to/17_balances/hmp_pca/meta_9194_regions_simple.csv";
    std::string out_dir = "path/to/17_balances/hmp_pca/";

    LOG_INFO << "Reading data";
    auto dfr = DataframeReader<double>();
    dfr.csv_reader().separator_chars( "," );
    auto meta = dfr.read( from_file( in_meta ));
    LOG_INFO << "meta: " << meta.rows() << " x " << meta.cols();
    auto samples = JplaceReader().read( from_files( utils::dir_list_files( in_jplace, ".*\\.jplace" )));
    LOG_INFO << "samples: " << samples.size();

    LOG_INFO << "Converting to Mass Trees";
    auto trees = convert_sample_set_to_mass_trees( samples, false ).first;

    for( size_t i = 0; i < trees.size(); ++i ) {
        // LOG_DBG << "tree " << i << " " << trees.name_at(i);
        // make_rooted( trees[i] );
        ladderize( trees[i] );
    }

    // Copy the meta data in the correct sample order.
    assert( meta.rows() == trees.size() );
    auto locations = std::vector<double>( meta.rows() );
    auto regions = std::vector<double>( meta.rows() );
    for( size_t i = 0; i < trees.size(); ++i ) {
        locations[i] = meta[ "Location" ].as<double>()[ trees.name_at(i) ];
        regions[i]   = meta[ "Region" ].as<double>()[ trees.name_at(i) ];
    }

    LOG_INFO << "Computing Balance data";
    BalanceSettings settings;
    settings.pseudo_count_summand_all = 0.65;
    settings.tendency = BalanceSettings::WeightTendency::kNone;
    settings.norm = BalanceSettings::WeightNorm::kNone;
    auto const bal_data = mass_balance_data( trees, settings );

    // Will need for all trees.
    auto lm = LayoutParameters();
    lm.ladderize = false;
    lm.stroke.width = 10;

    // Write edge masses and taxon weights
    utils::MatrixWriter<double>().to_file(
        bal_data.edge_masses, out_dir + "edge_masses.csv", trees.names()
    );
    assert( bal_data.edge_masses.cols() == bal_data.tree.edge_count() );
    auto emtotal = std::vector<double>( bal_data.edge_masses.cols(), 0.0 );
    for( size_t r = 0; r < bal_data.edge_masses.rows(); ++r ) {
        for( size_t c = 0; c < bal_data.edge_masses.cols(); ++c ) {
            emtotal[c] += bal_data.edge_masses(r,c);
        }
    }
    assert( emtotal.size() == bal_data.tree.edge_count() );
    write_color_tree_to_svg_file(
        bal_data.tree, lm, emtotal,
         ColorMap( color_list_viridis() ), ColorNormalizationLinear( emtotal ),
        out_dir + "edge_masses.svg"
    );
    // LOG_DBG << "taxon weights " << utils::join( bal_data.taxon_weights );
    // if( bal_data.taxon_weights.size() > 0 ) {
    //     write_color_tree_to_svg_file(
    //         bal_data.tree, lm, bal_data.taxon_weights, ColorMap( color_list_viridis() ), ColorNormalizationLinear(bal_data.taxon_weights),
    //         out_dir + "taxon_weights.svg"
    //     );
    // }

    // Balances computation
    LOG_INFO << "Computing Edge Balances";
    auto const bals = edge_balances( bal_data, true );
    LOG_INFO << "bals: " << bals.rows() << " x " << bals.cols();
    assert( bals.rows() == trees.size() );
    assert( meta.rows() == bals.rows() );
    assert( bals.cols() == trees[0].edge_count() );
    utils::MatrixWriter<double>().to_file(
        bals, out_dir + "balances_raw.csv", trees.names()
    );

    // Remove constant columns.
    auto bals_non_const = bals;
    auto non_const_edges = filter_constant_columns( bals_non_const );
    utils::MatrixWriter<double>().to_file(
        bals_non_const, out_dir + "balances_nc.csv", trees.names()
    );
    LOG_INFO << "bals_non_const: " << bals_non_const.rows() << " x " << bals_non_const.cols();

    size_t const num_comp = 2;
    LOG_INFO << "Balance pca";
    auto const pca = utils::principal_component_analysis(
        bals_non_const, num_comp, utils::PcaStandardization::kCovariance
    );
    assert( pca.eigenvalues.size()  == num_comp );
    // assert( pca.eigenvectors.rows() == trees[0].edge_count() );
    assert( pca.eigenvectors.cols() == num_comp );
    assert( pca.projection.rows()   == trees.size() );
    assert( pca.projection.cols()   == num_comp );

    LOG_INFO << "Writing out balance pca results";

    utils::MatrixWriter<double>().to_file(
        pca.eigenvectors, out_dir + "eigenvectors.csv"
    );
    utils::MatrixWriter<double>().to_file(
        pca.projection, out_dir + "projection.csv", trees.names()
    );

    auto proj_meta = utils::Matrix<double>( trees.size(), num_comp + 2, 0.0 );
    for( size_t i = 0; i < num_comp; ++i ) {
        proj_meta.col(i) = pca.projection.col(i);
    }
    proj_meta.col( num_comp     ) = locations;
    proj_meta.col( num_comp + 1 ) = regions;
    utils::MatrixWriter<double>().to_file(
        proj_meta, out_dir + "proj_meta.csv", trees.names()
    );

    for( size_t i = 0; i < num_comp; ++i ) {
        auto edge_vals = std::vector<double>( bal_data.tree.edge_count(), 0.0 );
        for( size_t j = 0; j < non_const_edges.size(); ++j ) {
            edge_vals[non_const_edges[j]] = pca.eigenvectors(j,i);
        }

        auto cm = ColorMap( color_list_spectral() );
        cm.mask_color( utils::color_from_hex( "#dfdfdf" ));
        auto cn = ColorNormalizationDiverging( edge_vals );
        cn.make_centric();

        write_color_tree_to_svg_file(
            bal_data.tree, lm, edge_vals, cm, cn,
            out_dir + "eigenvector_" + std::to_string(i) + ".svg"
        );
        //
        // auto cm = ColorMap( color_list_spectral() );
        // cm.mask_color( utils::color_from_hex( "#dfdfdf" ));
        // auto cn = ColorNormalizationDiverging( pca.eigenvectors.col(i) );
        // cn.make_centric();
        //
        // write_color_tree_to_svg_file(
        //     bal_data.tree, lm, pca.eigenvectors.col(i), cm, cn,
        //     out_dir + "eigenvector_" + std::to_string(i) + ".svg"
        // );
    }

    LOG_INFO << "Pairwise distances";
    auto const pd_mat = euclidean_distance_matrix( bals );
    LOG_INFO << "pd_mat: " << pd_mat.rows() << " x " << pd_mat.cols();
    utils::MatrixWriter<double>().to_file(
        pd_mat, out_dir + "balance_distances.csv"
    );

    LOG_INFO << "Finished";

    return 0;
}
