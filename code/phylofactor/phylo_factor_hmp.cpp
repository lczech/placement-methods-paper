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
    std::string in_jplace = "path/to/06_samples/hmp_Bact_Unconstr_backbone/samples_jplace/";
    std::string in_meta = "path/to/17_balances/hmp_phylofactor/meta_9194_regions.csv";
    std::string out_dir = "path/to/17_balances/hmp_phylofactor/";

    LOG_INFO << "Reading data";
    auto dfr = DataframeReader<std::string>();
    dfr.csv_reader().separator_chars( "," );
    dfr.col_names_from_first_row( false );
    auto meta = dfr.read( from_file( in_meta ));
    LOG_INFO << "meta: " << meta.rows() << " x " << meta.cols();
    auto const jfiles = utils::dir_list_files( in_jplace, true, ".*\\.jplace" );
    LOG_INFO << "found " << jfiles.size() << " jplace files";

    // {
    //     auto const fit = std::find(
    //         jfiles.begin(), jfiles.end(),
    //         std::string( "path/to/06_samples/hmp_Bact_Unconstr_backbone/samples_jplace/srp_SRS098352.jplace" )
    //     );
    //     auto pos = std::distance( jfiles.begin(), fit );
    //     LOG_DBG << "pos " << pos;
    //
    //     LOG_DBG << "file_basename " << file_basename( jfiles[pos] );
    //     LOG_DBG << "is_file( jfiles[pos] ) " << ( is_file( jfiles[pos] ) ? "si" : "no" );
    //     LOG_DBG << "file_exists( jfiles[pos] ) " << ( file_exists( jfiles[pos] ) ? "si" : "no" );
    //     LOG_DBG << "file_size " << file_size( jfiles[pos] );
    //     auto const cnt = file_read( jfiles[pos] );
    //     LOG_DBG << "length " << cnt.size();
    //
    //     auto const testsmp = JplaceReader().read( from_file( jfiles[pos] ));
    //     LOG_DBG << "sample size " << testsmp.size();
    // }

    auto samples = JplaceReader().read( from_files( jfiles ));
    LOG_INFO << "samples: " << samples.size();

    LOG_INFO << "Converting to Mass Trees";
    auto trees = convert_sample_set_to_mass_trees( samples, false ).first;

    for( size_t i = 0; i < trees.size(); ++i ) {
        // LOG_DBG << "tree " << i << " " << trees.name_at(i);
        // make_rooted( trees[i] );
        ladderize( trees[i] );
    }

    // Copy the meta data in the correct sample order.
    LOG_INFO << "Get meta data";
    auto const& meta_2_col = meta[2].as<std::string>();
    auto const factor = glm_factor( meta_2_col.begin(), meta_2_col.end(), { "stool", "mouth" }, {"NA"} );
    // auto const factor = glm_factor( meta[2].begin(), meta[2].end(), {}, {"NA"} );

    auto const factor_smry = glm_factor_summary( factor );
    for( size_t i = 0; i < factor_smry.size(); ++i ) {
        LOG_DBG1 << factor.levels[i] << ": " << factor_smry[i];
    }
    auto const indicator = glm_indicator_variables( factor, meta.row_names() );
    LOG_DBG << "indicator.rows() " << indicator.rows() << " indicator.cols() " << indicator.cols() << " trees.size() " << trees.size();
    Matrix<double> meta_vals( indicator.rows(), indicator.cols(), 0.0 );
    assert( meta.rows() == trees.size() );
    for( size_t i = 0; i <  indicator.rows(); ++i ) {
        for( size_t j = 0; j <  indicator.cols(); ++j ) {
            meta_vals(i,j) = indicator[ j ].as<double>()[ trees.name_at(i) ];
        }
    }
    utils::MatrixWriter<double>().to_file(
        meta_vals, out_dir + "meta_vals.csv", trees.names()
    );

    LOG_INFO << "Computing Balance data";
    BalanceSettings settings;
    // settings.tendency = BalanceSettings::WeightTendency::kNone;
    // settings.norm = BalanceSettings::WeightNorm::kNone;
    settings.pseudo_count_summand_all = 0.65;
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
    // write_color_tree_to_svg_file(
    //     bal_data.tree, lm, bal_data.taxon_weights,
    //      ColorMap( color_list_viridis() ), ColorNormalizationLinear( bal_data.taxon_weights ),
    //     out_dir + "taxon_weights.svg"
    // );

    LOG_INFO << "Running Phylo Factorization";

    size_t const num_factors = 5;
    auto const factors = phylogenetic_factorization(
        bal_data,
        [&]( std::vector<double> const& balances ){
            auto const fit = glm_fit( meta_vals, balances, glm_family_gaussian() );
            if( !fit.converged ) {
                LOG_DBG << "did not converge";
            }

            if( ! std::isfinite(fit.null_deviance) ) {
                LOG_DBG << "fit.null_deviance " << fit.null_deviance;
            }
            if( ! std::isfinite(fit.deviance) ) {
                LOG_DBG << "fit.deviance " << fit.deviance;
            }
            if( ! std::isfinite(fit.null_deviance) || ! std::isfinite(fit.deviance) ) {
                LOG_DBG << "meta_vals " << join( meta_vals );
                LOG_DBG << "balances " << join( balances );
            }

            return fit.null_deviance - fit.deviance;
        },
        num_factors,
        []( size_t iteration, size_t max_iterations ){
            LOG_DBG1 << "iteration " << iteration << " of " << max_iterations;
        }
    );
    auto const found_factors = factors.size();

    auto balances = utils::Matrix<double>( meta_vals.rows(), found_factors );
    // balances.col(found_factors) = meta_vals;

    for( size_t i = 0; i < found_factors; ++i ) {
        auto const& factor = factors[i];

        LOG_INFO << "edge " << factor.edge_index << " ( " << factor.edge_indices_primary.size() << " | " << factor.edge_indices_secondary.size() << " ), ov " << factor.objective_value;

        // Make a tree with the edges of that factor.
        auto edge_cols = phylo_factor_single_factor_colors( bal_data.tree, factors, i );
        write_color_tree_to_svg_file(
            bal_data.tree, lm, edge_cols,
            out_dir + "factor_edges_" + std::to_string(i) + ".svg"
        );

        // LOG_DBG << "factor.balances.size() " << factor.balances.size();
        // LOG_DBG << "balances.rows() " << balances.rows();
        assert( factor.balances.size() == balances.rows() );
        balances.col(i) = factor.balances;

        auto cm = ColorMap( color_list_viridis() );
        cm.mask_color( Color( 0.8, 0.8, 0.8 ));
        write_color_tree_to_svg_file(
            bal_data.tree, lm, factor.all_objective_values,
             cm, ColorNormalizationLinear( factor.all_objective_values ),
            out_dir + "ovs_" + std::to_string(i) + ".svg"
        );
    }

    // Write balances of the factors.
    utils::MatrixWriter<double>().to_file(
        balances, out_dir + "factor_balances.csv", trees.names()
    );

    // // special test for correlation: the second factor splits away a cherry edge.
    // // let's see how much mass there is there.
    // auto const& sf = factors[1];
    // auto const ei1 = bal_data.tree.edge_at(sf.edge_index).secondary_link().next().edge().index();
    // auto const ei2 = bal_data.tree.edge_at(sf.edge_index).secondary_link().next().next().edge().index();
    // auto sfec = std::vector<Color>( bal_data.tree.edge_count(), Color( 0.8,0.8,0.8 ));
    // sfec[ei1] = sfec[ei2] = Color( 0.8, 0.0, 0.0 );
    // write_color_tree_to_svg_file(
    //     bal_data.tree, lm, sfec, out_dir + "sfec.svg"
    // );
    // assert( bal_data.edge_masses.rows() == meta_vals.size() );
    // auto sfmat = Matrix<double>( meta_vals.size(), 5 );
    // sfmat.col(0) = bal_data.edge_masses.col(ei1);
    // sfmat.col(1) = bal_data.edge_masses.col(ei2);
    // for( size_t i = 0; i < sfmat.rows(); ++i ) {
    //     sfmat(i,2) = bal_data.edge_masses(i, ei1) + bal_data.edge_masses(i, ei2);
    // }
    // sfmat.col(3) = sf.balances;
    // sfmat.col(4) = meta_vals;
    // utils::MatrixWriter<double>().to_file(
    //     sfmat, out_dir + "sfem.csv", trees.names()
    // );

    LOG_INFO << "Finished";

    return 0;
}
