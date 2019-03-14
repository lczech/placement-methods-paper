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

Tree translate_labels( Tree const& tree, std::unordered_map<std::string, std::string> label_map )
{
    auto result = tree;
    for( auto& node : result.nodes() ) {
        auto& name = node.data<CommonNodeData>().name;
        // if( ! name.empty() ) {
        //     if( label_map.count( name ) > 0 ) {
        //         LOG_DBG1 << name << " << " << label_map[name];
        //     } else {
        //         LOG_DBG1 << name;
        //     }
        // }
        if( label_map.count( name ) > 0 ) {
            name = label_map[name];
        }
    }
    return result;
}

void write_edge_mass_matrix(
    Matrix<double> const& mat,
    std::vector<double> const& meta_vals,
    Tree const& tree,
    std::string const& out_dir
) {
    auto write_files = [](
        Matrix<double> const& mat,
        std::string const& out_name
    ){
        utils::MatrixWriter<double>().to_file(
            mat, out_name + ".csv"
        );

        auto const mm = matrix_minmax(mat);
        auto map = ColorMap( color_list_ylorrd() );
        auto norm = ColorNormalizationLinear( mm.min, mm.max );

        auto mat_col = Matrix<Color>( mat.rows(), mat.cols() );
        for( size_t r = 0; r < mat.rows(); ++r ) {
            for( size_t c = 0; c < mat.cols(); ++c ) {
                mat_col( r, c ) = map( norm, mat(r,c) );
            }
        }

        BmpWriter().to_file( mat_col, out_name + ".bmp" );
    };

    // Write first mat
    // write_files( mat, out_dir + "edge_masses_raw" );

    // Sort matrix by nugent score
    auto mat_b = Matrix<double>( mat.rows(), mat.cols() );
    auto const si = sort_indices( meta_vals.begin(), meta_vals.end() );
    for( size_t i = 0; i < si.size(); ++i ) {
        auto const sii = si[i];
        mat_b.row(i) = mat.row(sii).to_vector();
    }
    // write_files( mat_b, out_dir + "edge_masses_rows" );

    // Sort matrix by tree node order
    auto mat_c = Matrix<double>( mat.rows(), mat.cols() );
    size_t c = 0;
    for( auto it : postorder(tree) ) {
        if( it.is_last_iteration() ) {
            continue;
        }
        mat_c.col(c) = mat_b.col( it.edge().index() ).to_vector();
        // LOG_DBG << "c " << c << " index " << it.edge().index() << ": " << join( mat_c.col(c), ", " );
        ++c;
    }
    assert( c == mat_c.cols() );

    // center log ratio transform per sample.
    for( size_t r = 0; r < mat_c.rows(); ++r ) {
        auto log_vec = mat_c.row(r).to_vector();
        for( auto& e : log_vec ) {
            e = std::log(e);
        }
        auto const avg = arithmetic_mean(log_vec);
        for( auto& e : log_vec ) {
            e -= avg;
        }
        mat_c.row(r) = log_vec;
    }

    write_files( mat_c, out_dir + "edge_masses_sorted" );
}

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
    // std::string in_jplace = "path/to/03_bv/03_epa/orig_queries_jplace_clean_15";
    std::string in_jplace = "path/to/03_bv/03_epa/orig_queries_jplace_clean";
    std::string in_meta = "path/to/data/bv/meta/meta_simple.csv";
    std::string label_file = "path/to/data/bv/meta/sequence_names.csv";
    std::string out_dir = "path/to/philr_bv/factors/";

    // ----------------------------------------------------------------
    //     Reading data
    // ----------------------------------------------------------------

    LOG_INFO << "Reading data";
    auto dfr = DataframeReader<double>();
    dfr.csv_reader().separator_chars( "\t" );
    auto meta = dfr.read( from_file( in_meta ));
    LOG_INFO << "meta: " << meta.rows() << " x " << meta.cols();
    auto samples = JplaceReader().read( from_files( utils::dir_list_files( in_jplace, ".*\\.jplace" )));
    LOG_INFO << "samples: " << samples.size();

    LOG_INFO << "Converting to Mass Trees";
    auto trees = convert_sample_set_to_mass_trees( samples, false ).first;
    assert( trees.size() == meta.rows() );

    for( size_t i = 0; i < trees.size(); ++i ) {
        // LOG_DBG << "tree " << i << " " << trees.name_at(i);
        // make_rooted( trees[i] );
        ladderize( trees[i] );
    }

    // Copy the meta data in the correct sample order.
    LOG_INFO << "Get meta vector";
    std::vector<double> meta_vals( trees.size(), 0.0 );
    assert( meta.rows() == trees.size() );
    for( size_t i = 0; i < trees.size(); ++i ) {
        meta_vals[i] = meta[ "nugent" ].as<double>()[ trees.name_at(i) ];
    }

    Matrix<double> meta_vals_mat( trees.size(), meta.cols() );
    for( size_t i = 0; i < trees.size(); ++i ) {
        for( size_t c = 0; c < meta.cols(); ++c ) {
            meta_vals_mat( i, c ) = meta[ c ].as<double>()[ trees.name_at(i) ];
        }
    }

    // ----------------------------------------------------------------
    //     Computing Balance data
    // ----------------------------------------------------------------

    LOG_INFO << "Computing Balance data";
    BalanceSettings settings;
    settings.tendency = BalanceSettings::WeightTendency::kNone;
    settings.norm = BalanceSettings::WeightNorm::kNone;
    settings.pseudo_count_summand_all = 0.65;
    auto const bal_data = mass_balance_data( trees, settings );

    // tree label translation
    auto csv_reader = CsvReader();
    csv_reader.separator_chars("\t");
    auto labels = csv_reader.read( from_file( label_file ));
    std::unordered_map<std::string, std::string> label_map;
    for( auto const& line : labels ) {
        assert( line.size() == 2 );
        assert( label_map.count( line[0] ) == 0 );
        label_map[ line[0] ] = line[1];
    }
    auto const label_tree = translate_labels( bal_data.tree, label_map );

    // Will need for all trees.
    auto lm = LayoutParameters();
    lm.ladderize = false;
    lm.stroke.width = 10;

    // ----------------------------------------------------------------
    //     Write edge masses and taxon weights
    // ----------------------------------------------------------------

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
        label_tree, lm, emtotal,
        ColorMap( color_list_viridis() ), ColorNormalizationLinear( emtotal ),
        out_dir + "edge_masses.svg"
    );
    // write_color_tree_to_svg_file(
    //     bal_data.tree, lm, bal_data.taxon_weights,
    //      ColorMap( color_list_viridis() ), ColorNormalizationLinear( bal_data.taxon_weights ),
    //     out_dir + "taxon_weights.svg"
    // );

    // ----------------------------------------------------------------
    //     Write heatmap stuff
    // ----------------------------------------------------------------

    auto layout = RectangularLayout( label_tree, LayoutType::kPhylogram, false );
    auto const spreading = LayoutSpreading::kAllNodesButRoot;
    layout.inner_node_spreading( spreading );

    // Set colourful edges.
    // std::vector<utils::SvgStroke> strokes;
    // for( size_t i = 0; i < tree.edge_count(); ++i ) {
    //     auto const& cv = scheme[ i % scheme.size() ];
    //     strokes.push_back( utils::SvgStroke() );
    //     strokes.back().color = utils::color_from_name_web( cv );
    //
    //     // auto& node = tree.edge_at(i).secondary_link().node();
    //     // if( is_leaf(node) ) {
    //     //     node.data<CommonNodeData>().name = cv;
    //     // }
    // }
    // layout.set_edge_strokes( strokes );

    // Set label alignment.
    layout.align_labels( true );
    auto spacer_stroke = utils::SvgStroke( utils::Color( 0.8, 0.8, 0.8 ));
    spacer_stroke.dash_array = std::vector<double>({ 2.0, 0.5 });
    spacer_stroke.dash_offset = 2.0;
    layout.set_label_spacer_strokes( spacer_stroke, spreading );

    // Write to file.
    std::ofstream spread_tree_os;
    utils::file_output_stream( out_dir + "spread_tree.svg", spread_tree_os );
    layout.to_svg_document().write( spread_tree_os );
    spread_tree_os.close();

    write_edge_mass_matrix( bal_data.edge_masses, meta_vals, label_tree, out_dir );

    // ----------------------------------------------------------------
    //     Run!
    // ----------------------------------------------------------------

    LOG_INFO << "Running Phylo Factorization";

    auto obj_fct_lin_reg = [&]( std::vector<double> const& balances ){
        // return std::abs( spearmans_rank_correlation_coefficient( meta_vals, balances ));

        // Null deviance.
        auto const ms = mean_stddev( balances );
        auto const var = ms.stddev * ms.stddev;

        // Deviance of a linear regression model.
        auto const lf = simple_linear_regression(
            meta_vals.begin(), meta_vals.end(),
            balances.begin(), balances.end()
        );
        auto const err = mean_squared_error(
            meta_vals.begin(), meta_vals.end(),
            balances.begin(), balances.end(),
            lf
        );
        // auto const fvu = fraction_of_variance_unexplained(
        //     meta_vals.begin(), meta_vals.end(),
        //     balances.begin(), balances.end(),
        //     lf
        // );
        //
        // // LOG_DBG2 << "f " << fvu.first << " s " << fvu.second;
        // return 1.0 - fvu;

        // LOG_DBG2 << "err " << err << " var " << var << " e-v " << (err-var) << " v-e " << (var-err);

        // Objective is the difference between the two.
        // return var;
        return var - err;
        // return var / err;
        // return err / var;
        // return err;
    };

    auto obj_fct_glm = [&]( std::vector<double> const& balances ){
        auto const fit = glm_fit( meta_vals_mat, balances, glm_family_gaussian() );
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
            LOG_DBG << "meta_vals_mat " << join( meta_vals_mat );
            LOG_DBG << "balances " << join( balances );
        }

        return fit.null_deviance - fit.deviance;
    };

    size_t const num_factors = 10;
    auto const factors = phylogenetic_factorization(
        bal_data,
        obj_fct_glm,
        num_factors,
        []( size_t iteration, size_t max_iterations ){
            LOG_DBG1 << "iteration " << iteration << " of " << max_iterations;
        }
    );
    auto const found_factors = factors.size();

    assert( meta_vals.size() == trees.size() );
    auto balances = utils::Matrix<double>( meta_vals.size(), found_factors + 1 );
    balances.col(found_factors) = meta_vals;

    std::ofstream factor_taxa_of;
    file_output_stream( out_dir + "factor_taxa.txt", factor_taxa_of );

    for( size_t i = 0; i < found_factors; ++i ) {
        auto const& factor = factors[i];

        LOG_INFO << "edge " << factor.edge_index << " ( " << factor.edge_indices_primary.size() << " | " << factor.edge_indices_secondary.size() << " ), ov " << factor.objective_value;

        // Make a tree with the edges of that factor.
        auto edge_cols = phylo_factor_single_factor_colors( bal_data.tree, factors, i );
        write_color_tree_to_svg_file(
            label_tree, lm, edge_cols,
            out_dir + "factor_edges_" + std::to_string(i) + ".svg"
        );

        assert( factor.balances.size() == balances.rows() );
        balances.col(i) = factor.balances;

        // write objective value trees
        auto cm = ColorMap( color_list_viridis() );
        cm.mask_color( Color( 0.8, 0.8, 0.8 ));
        write_color_tree_to_svg_file(
            label_tree, lm, factor.all_objective_values,
            cm, ColorNormalizationLinear( factor.all_objective_values ),
            out_dir + "ovs_" + std::to_string(i) + ".svg"
        );

        // write out which taxa are on the secondary part of the factor.
        factor_taxa_of << "Factor " << i << ":\n";
        std::unordered_set<std::string> edge_names;
        for( auto const ei : factor.edge_indices_secondary ) {
            auto const& ed = label_tree.edge_at(ei).secondary_link().node().data<CommonNodeData>();
            if( ! ed.name.empty() ) {
                edge_names.insert(ed.name);
            }
        }
        for( auto const& en : edge_names ) {
            factor_taxa_of << en << "\n";
        }
        factor_taxa_of << "\n";
    }

    // Write balances of the factors.
    utils::MatrixWriter<double>().to_file(
        balances, out_dir + "factor_balances.csv", trees.names()
    );

    // Make a tree with the edges of that factor.
    PhyloFactorCladeColors clade_cols;
    std::vector<utils::Color> clade_colors = {{
        { 0.301961, 0.686275, 0.290196 }, // green
        { 0.215686, 0.494118, 0.721569 }, // blue
        { 0.894118, 0.101961, 0.109804 }, // red
        { 0.596078, 0.305882, 0.639216 }, // purple
        { 1.000000, 0.498039, 0.000000 }, // orange
        // { 1.000000, 1.000000, 0.200000 }, // yellow
        { 0.650980, 0.337255, 0.156863 }, // brown
        { 0.968627, 0.505882, 0.749020 }, // pink
        // { 0.600000, 0.600000, 0.600000 }  // grey
        { 0.400000, 0.760784, 0.647059 }, // teal
        // { 0.988235, 0.552941, 0.384314 }, // orange
        { 0.552941, 0.627451, 0.796078 }, // lilac
        { 0.905882, 0.541176, 0.764706 }, // magenta
    }};
    clade_cols.clade_colors = clade_colors;
    auto all_edge_cols = phylo_factor_clade_colors( bal_data.tree, factors, 0, clade_cols );
    write_color_tree_to_svg_file(
        label_tree, lm, all_edge_cols,
        out_dir + "factors_tree.svg"
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