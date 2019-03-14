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
#include <cassert>
#include <cmath>
#include <string>
#include <unordered_map>
#include <unordered_set>

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

// =================================================================================================
//     translate_tree_labels
// =================================================================================================

Tree translate_tree_labels( Tree const& tree, std::string const& in_label_file )
{
    if( in_label_file.empty() ) {
        // result.label_tree = result.bal_data.tree;
        return tree;
    }

    auto csv_reader = CsvReader();
    csv_reader.separator_chars("\t");
    auto labels = csv_reader.read( from_file( in_label_file ));

    std::unordered_map<std::string, std::string> label_map;
    for( auto const& line : labels ) {
        assert( line.size() == 2 );
        assert( label_map.count( line[0] ) == 0 );
        label_map[ line[0] ] = line[1];
    }

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

// =================================================================================================
//     read_data
// =================================================================================================

struct PFData
{
    // std::vector<double> meta_vals;
    Matrix<double> meta_vals_mat;
    std::vector<size_t> col_sort_indices;
    std::vector<std::string> col_names;
    int main_col_idx = -1;

    SampleSet samples;
    TreeSet trees;
    Tree label_tree;
    BalanceData bal_data;
    LayoutParameters lm;
};

PFData read_data(
    std::string const& in_jplace,
    std::string const& in_meta,
    BalanceSettings settings,
    std::string const& in_label_file = "",
    std::string        main_col_name = "",
    double line_width = 10.0,
    std::string meta_sep_char = "\t",
    bool use_only_main_col = false
) {
    PFData data;

    LOG_INFO << "Reading meta data";
    // auto dfr = DataframeReader<double>();
    auto dfr = DataframeReader<std::string>();
    dfr.csv_reader().separator_chars( meta_sep_char );
    auto meta_str = dfr.read( from_file( in_meta ));
    LOG_INFO << "meta_str: " << meta_str.rows() << " x " << meta_str.cols();

    if( use_only_main_col ) {
        if( main_col_name.empty() ) {
            throw std::runtime_error( "Cannot use main_col_name if not given." );
        }

        Dataframe tmp;
        for( auto const& name : meta_str.row_names() ) {
            tmp.add_row( name );
        }
        tmp.add_col( main_col_name, meta_str[ main_col_name ].as<std::string>().to_vector() );
        meta_str = tmp;
    }

    std::string report;
    auto meta = glm_prepare_dataframe( meta_str, report );
    LOG_INFO << report;
    LOG_INFO << "meta: " << meta.rows() << " x " << meta.cols();

    LOG_INFO << "Reading jplace files";
    data.samples = JplaceReader().read( from_files( utils::dir_list_files( in_jplace, ".*\\.jplace" )));
    LOG_INFO << "samples: " << data.samples.size();

    LOG_INFO << "Converting to Mass Trees";
    data.trees = convert_sample_set_to_mass_trees( data.samples, false ).first;

    for( size_t i = 0; i < data.trees.size(); ++i ) {
        // LOG_DBG << "tree " << i << " " << data.trees.name_at(i);
        // make_rooted( data.trees[i] );
        ladderize( data.trees[i] );
    }

    // Copy the meta data in the correct sample order.
    LOG_INFO << "Get meta vector";
    if( main_col_name.empty() && meta.cols() == 1 ) {
        main_col_name = meta[0].name();
    }
    if( ! main_col_name.empty() ) {
        auto meta_vals = std::vector<double>( data.trees.size(), 0.0 );
        assert( meta.rows() == data.trees.size() );
        for( size_t i = 0; i < data.trees.size(); ++i ) {
            meta_vals[i] = meta[ main_col_name ].as<double>()[ data.trees.name_at(i) ];
        }
        data.col_sort_indices = sort_indices( meta_vals.begin(), meta_vals.end() );
        data.main_col_idx = meta[ main_col_name ].index();
    } else {
        data.col_sort_indices = std::vector<size_t>( data.trees.size() );
        std::iota( data.col_sort_indices.begin(), data.col_sort_indices.end(), 0 );
    }

    data.meta_vals_mat = Matrix<double>( data.trees.size(), meta.cols() );
    data.col_names = std::vector<std::string>( meta.cols() );
    assert( meta.rows() == data.trees.size() );
    for( size_t i = 0; i < data.trees.size(); ++i ) {
        for( size_t c = 0; c < meta.cols(); ++c ) {
            data.meta_vals_mat( i, c ) = meta[ c ].as<double>()[ data.trees.name_at(i) ];
        }
    }
    for( size_t c = 0; c < meta.cols(); ++c ) {
        data.col_names[c] = meta[ c ].name();
    }

    LOG_INFO << "Computing Balance data";
    data.bal_data = mass_balance_data( data.trees, settings );
    data.label_tree = translate_tree_labels( data.bal_data.tree, in_label_file );

    // Will need for all trees.
    data.lm.ladderize = false;
    // data.lm.stroke.width = 4;
    // data.lm.stroke.width = 10;
    data.lm.stroke.width = line_width;

    return data;
}

// =================================================================================================
//     calc_and_write_balances_mat_and_pca
// =================================================================================================

void calc_and_write_balances_mat_and_pca( PFData const& data, std::string const& out_dir )
{
    LOG_DBG << "calc_and_write_balances_mat_and_pca";

    // Test for bv dataset
    // {
    //     auto edge_colors = std::vector<Color>( data.bal_data.tree.edge_count(), Color( 0.8,0.8,0.8 ) );
    //     edge_colors[1015] = Color( 1,0,0 );
    //     edge_colors[362] = Color( 0,1,0 );
    //     write_color_tree_to_svg_file(
    //         data.bal_data.tree, data.lm, edge_colors,
    //         out_dir + "interesting_edges.svg"
    //     );
    // }

    // Preprae data and write out matrices
    auto mat = edge_balances( data.bal_data );
    utils::MatrixWriter<double>().to_file(
        mat, out_dir + "edge_balances_all.csv", data.trees.names()
    );

    // correlation per edge, if we have a main meta data feature to correlate with.
    if( data.main_col_idx > -1 ) {
        LOG_DBG1 << "correlation";

        assert( mat.cols() == data.label_tree.edge_count() );
        assert( mat.rows() == data.meta_vals_mat.rows() );

        // Test for bv dataset
        // {
        //     std::ofstream ie_of;
        //     utils::file_output_stream( out_dir + "interesting_edges.txt", ie_of );
        //     ie_of << "edge_362\tedge_1015\tnugent\n";
        //     for( size_t r = 0; r < mat.rows(); ++r ) {
        //         ie_of << mat( r, 362 ) << "\t" << mat( r, 1015 ) << "\t" << data.meta_vals_mat( r, data.main_col_idx ) << "\n";
        //     }
        //
        //     auto fit_edge = [](
        //         Matrix<double> const&      x_predictors,
        //         std::vector<double> const& y_response
        //     ){
        //         auto const fit = glm_fit( x_predictors, y_response );
        //         if( !fit.converged ) {
        //             LOG_DBG << "did not converge";
        //         }
        //
        //         if( ! std::isfinite(fit.null_deviance) ) {
        //             LOG_DBG << "fit.null_deviance " << fit.null_deviance;
        //         }
        //         if( ! std::isfinite(fit.deviance) ) {
        //             LOG_DBG << "fit.deviance " << fit.deviance;
        //         }
        //         if( ! std::isfinite(fit.null_deviance) || ! std::isfinite(fit.deviance) ) {
        //             LOG_DBG << "data.meta_vals_mat " << join( x_predictors );
        //             LOG_DBG << "balances " << join( y_response );
        //         }
        //
        //         LOG_DBG << "fit.null_deviance " << fit.null_deviance << " fit.deviance " << fit.deviance;
        //         LOG_DBG << "fit nd - d " << (fit.null_deviance - fit.deviance);
        //     };
        //
        //     auto mat_ie = Matrix<double>( data.meta_vals_mat.rows(), 1 );
        //     mat_ie.col(0) = data.meta_vals_mat.col( data.main_col_idx ).to_vector();
        //     auto bal_362  = mat.col( 362 );
        //     auto bal_1015 = mat.col( 1015 );
        //
        //     LOG_DBG << "fit 362";
        //     fit_edge( mat_ie, bal_362 );
        //     LOG_DBG << "fit 1015";
        //     fit_edge( mat_ie, bal_1015 );
        // }

        // Write out the edges that have a strong correlation,
        // for testing purposes only.
        // std::ofstream corr_of;
        // file_output_stream( out_dir + "edge_balances_correlation_strong.txt", corr_of );
        // corr_of << "meta:\n";
        // corr_of << join( data.meta_vals_mat.col(data.main_col_idx).to_vector() ) << "\n\n";

        auto cm = ColorMap( color_list_spectral() );
        cm.mask_color( utils::color_from_hex( "#dfdfdf" ));
        auto cn = ColorNormalizationDiverging( -1.0, 0.0, 1.0 );

        auto corr = std::vector<double>( mat.cols(), 0.0 );
        for( size_t i = 0; i < mat.cols(); ++i ) {
            corr[i] = spearmans_rank_correlation_coefficient(
                mat.col(i).begin(), mat.col(i).end(),
                data.meta_vals_mat.col(data.main_col_idx).begin(),
                data.meta_vals_mat.col(data.main_col_idx).end()
            );

            assert(( -1.0 <= corr[i] && corr[i] <= 1.0 ) || ( ! std::isfinite(corr[i]) ) );
            // if( corr[i] > 0.5 || corr[i] < -0.5 ) {
            //     corr_of << "At edge " << i << " correlation " << corr[i] << "\n";
            //     corr_of << join( mat.col(i).to_vector() ) << "\n\n";
            // }
        }
        // corr_of << "All corrs: " << join(corr) << "\n";

        write_color_tree_to_svg_file(
            data.label_tree, data.lm, corr, cm, cn,
            out_dir + "edge_balances_correlation.svg"
        );
    }

    // filter out const cols
    auto const non_const_edges = filter_constant_columns( mat );
    utils::MatrixWriter<double>().to_file(
        mat, out_dir + "edge_balances_filtered.csv"
    );

    // Run PCA
    LOG_DBG1 << "pca";
    size_t components = 10;
    auto pca = utils::principal_component_analysis(
        mat, components, utils::PcaStandardization::kCovariance
    );

    LOG_DBG1 << "results";

    // Write result matrices
    utils::MatrixWriter<double>().to_file(
        pca.eigenvectors, out_dir + "edge_balances_pca_eigenvectors.csv"
    );
    utils::MatrixWriter<double>().to_file(
        pca.projection, out_dir + "edge_balances_pca_projection.csv"
    );
    auto proj_meta = utils::Matrix<double>( data.trees.size(), components + data.meta_vals_mat.cols(), 0.0 );
    assert( proj_meta.rows() == data.meta_vals_mat.rows() );
    for( size_t i = 0; i < components; ++i ) {
        proj_meta.col(i) = pca.projection.col(i).to_vector();
    }
    for( size_t i = 0; i < data.meta_vals_mat.cols(); ++i ) {
        proj_meta.col( components + i ) = data.meta_vals_mat.col(i).to_vector();
    }
    utils::MatrixWriter<double>().to_file(
        proj_meta, out_dir + "edge_balances_pca_proj_meta.csv", data.trees.names()
    );

    // Write trees
    for( size_t i = 0; i < components; ++i ) {
        auto edge_vals = std::vector<double>( data.bal_data.tree.edge_count(), 0.0 );
        for( size_t j = 0; j < non_const_edges.size(); ++j ) {
            edge_vals[non_const_edges[j]] = pca.eigenvectors(j,i);
        }
        for( auto const& edge : data.bal_data.tree.edges() ) {
            // edge_vals[edge.index()] = pca.eigenvectors( edge.index(),i );
            if( is_leaf(edge) ) {
                edge_vals[edge.index()] = std::numeric_limits<double>::quiet_NaN();
            }
        }

        auto cm = ColorMap( color_list_spectral() );
        cm.mask_color( utils::color_from_hex( "#dfdfdf" ));
        auto cn = ColorNormalizationDiverging( edge_vals );
        cn.make_centric();

        write_color_tree_to_svg_file(
            data.bal_data.tree, data.lm, edge_vals, cm, cn,
            out_dir + "edge_balances_pca_eigenvector_" + std::to_string(i) + ".svg"
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
}

// =================================================================================================
//     Write edge masses and taxon weights
// =================================================================================================

void write_edge_masses_taxon_weights( PFData const& data, std::string const& out_dir )
{
    LOG_DBG << "write_edge_masses_taxon_weights";

    utils::MatrixWriter<double>().to_file(
        data.bal_data.edge_masses, out_dir + "edge_masses.csv", data.trees.names()
    );
    assert( data.bal_data.edge_masses.cols() == data.bal_data.tree.edge_count() );
    auto emtotal = std::vector<double>( data.bal_data.edge_masses.cols(), 0.0 );
    auto emrtotal = std::vector<double>( data.bal_data.edge_masses.cols(), 0.0 );
    for( size_t r = 0; r < data.bal_data.edge_masses.rows(); ++r ) {
        for( size_t c = 0; c < data.bal_data.edge_masses.cols(); ++c ) {
            emtotal[c] += data.bal_data.edge_masses(r,c);
            emrtotal[c] += data.bal_data.raw_edge_masses(r,c);
        }
    }
    assert( emtotal.size() == data.bal_data.tree.edge_count() );
    write_color_tree_to_svg_file(
        data.label_tree, data.lm, emtotal,
        ColorMap( color_list_viridis() ), ColorNormalizationLinear( emtotal ),
        out_dir + "edge_masses.svg"
    );
    write_color_tree_to_svg_file(
        data.label_tree, data.lm, emrtotal,
        ColorMap( color_list_viridis() ), ColorNormalizationLinear( emrtotal ),
        out_dir + "edge_masses_raw.svg"
    );

    // Only write weights if we have them.
    if( ! data.bal_data.taxon_weights.empty() ) {
        auto const mm = minimum_maximum(
            data.bal_data.taxon_weights.begin(), data.bal_data.taxon_weights.end()
        );
        if( mm.min < mm.max ) {
            write_color_tree_to_svg_file(
                data.label_tree, data.lm, data.bal_data.taxon_weights,
                ColorMap( color_list_viridis() ), ColorNormalizationLinear( data.bal_data.taxon_weights ),
                out_dir + "taxon_weights.svg"
            );
        }
    }
}

// =================================================================================================
//     write_edge_mass_matrix
// =================================================================================================

void write_edge_mass_matrix(
    Matrix<double> const& mat,
    // std::vector<double> const& meta_vals,
    std::vector<size_t> col_sort_indices,
    Tree const& tree,
    std::string const& out_dir,
    std::string const& suffix = ""
) {
    LOG_DBG << "write_edge_mass_matrix";

    auto write_files = [](
        Matrix<double> const& mat,
        std::string const& out_name
    ){
        utils::MatrixWriter<double>().to_file(
            mat, out_name + ".csv"
        );

        auto const mm = matrix_minmax(mat);
        auto map = ColorMap( color_list_ylorrd() );
        // map.reverse(true);
        auto norm = ColorNormalizationLinear( mm.min, mm.max );

        LOG_DBG << out_name << " min: " << mm.min << " max: " << mm.max;

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
    // auto const si = sort_indices( meta_vals.begin(), meta_vals.end() );
    auto const& si = col_sort_indices;
    for( size_t i = 0; i < si.size(); ++i ) {
        auto const sii = si[i];
        mat_b.row(i) = mat.row(sii).to_vector();
    }
    // write_files( mat_b, out_dir + "edge_masses_rows" );

    // Sort matrix by tree node order
    auto mat_c = Matrix<double>( mat.rows(), mat.cols() );
    // size_t c = 0;
    // for( auto it : postorder(tree) ) {
    //     if( it.is_last_iteration() ) {
    //         continue;
    //     }
    //     mat_c.col(c) = mat_b.col( it.edge().index() ).to_vector();
    //     // LOG_DBG << "c " << c << " index " << it.edge().index() << ": " << join( mat_c.col(c), ", " );
    //     ++c;
    // }
    // assert( c == mat_c.cols() );

    // size_t c = 0;
    // for( auto it : postorder(tree) ) {
    //     if( it.is_last_iteration() ) {
    //         continue;
    //     }
    //     mat_c.col(c) = mat_b.col( it.edge().index() ).to_vector();
    //     // LOG_DBG << "c " << c << " index " << it.edge().index() << ": " << join( mat_c.col(c), ", " );
    //     ++c;
    // }
    // assert( c == mat_c.cols() );
    size_t node_counter = 0;
    auto visits = std::vector<size_t>( tree.node_count(), 0 );
    for( auto it : eulertour( tree )) {
        auto const node_index = it.node().index();

        // Count the how many-th visit this is. As we have a bifurcating tree,
        // it can never surpass 3 visits.
        ++visits[ node_index ];
        assert( visits[ node_index ] <= 3 );

        if( is_root( it.node() )) {
            continue;
        } else if( is_leaf( it.node() ) || visits[ node_index ] == 2 ) {
            mat_c.col(node_counter) = mat_b.col( it.edge().index() ).to_vector();
            ++node_counter;
        }
    }

    // center log ratio transform per sample.
    for( size_t r = 0; r < mat_c.rows(); ++r ) {
        auto log_vec = mat_c.row(r).to_vector();
        for( auto& e : log_vec ) {
            e = std::log(e);
        }
        auto const avg = arithmetic_mean(log_vec);
        for( size_t c = 0; c < log_vec.size(); ++c ) {
            auto& e = log_vec[c];
            e -= avg;

            // Correction for taxon weight pseudo cound additions...
            // if( mat_c(r,c) == 0.65 ) {
            //     e = 0.0;
            // }
        }
        mat_c.row(r) = log_vec;
    }

    write_files( mat_c, out_dir + "edge_masses_sorted" + suffix );
}

// =================================================================================================
//     write_heat_tree
// =================================================================================================

void write_heat_tree(
    Matrix<double> const& raw_edge_masses,
    // std::vector<double> const& meta_vals,
    std::vector<size_t> col_sort_indices,
    Tree const& tree,
    std::string const& out_dir
) {
    LOG_DBG << "write_heat_tree";

    // assert( raw_edge_masses.rows() == meta_vals.size() );
    assert( raw_edge_masses.rows() == col_sort_indices.size() );
    assert( raw_edge_masses.cols() == tree.edge_count() );

    // Transpose and sort the matrix by edges/nodes.
    std::vector<size_t> tmp;
    for( auto const& e : tree.edges() ) {
        tmp.push_back( e.secondary_link().node().index() );
    }
    auto const sorting = utils::sort_indices( tmp.begin(), tmp.end() );
    assert( sorting.size() == tree.node_count() - 1 );
    auto matrix_t = Matrix<double>( raw_edge_masses.cols(), raw_edge_masses.rows() );
    for( size_t i = 0; i < sorting.size(); ++i ) {
        assert( sorting[i] < matrix_t.rows() );
        matrix_t.row( sorting[i] ) = raw_edge_masses.col(i).to_vector();
    }

    // Sort columns by nugent score.
    auto matrix = Matrix<double>( matrix_t.rows(), matrix_t.cols() );
    // auto const si = sort_indices( meta_vals.begin(), meta_vals.end() );
    auto const& si = col_sort_indices;
    assert( si.size() == matrix_t.cols() );
    for( size_t i = 0; i < si.size(); ++i ) {
        auto const sii = si[i];
        matrix.col(i) = matrix_t.col(sii).to_vector();
    }

    // center log ratio transform per sample.
    for( size_t c = 0; c < matrix.cols(); ++c ) {
        auto log_vec = matrix.col(c).to_vector();
        for( auto& e : log_vec ) {
            e = std::log(e);
        }
        auto const avg = arithmetic_mean(log_vec);
        for( size_t r = 0; r < log_vec.size(); ++r ) {
            auto& e = log_vec[r];
            e -= avg;
        }
        matrix.col(c) = log_vec;
        // auto const mm = minimum_maximum( log_vec );
        // LOG_DBG1 << "heat mat col " << c << " min " << mm.min << " max " << mm.max << " avg " << avg;
    }

    // Calc tree masses
    auto const edge_masses = matrix_col_sums( raw_edge_masses );
    assert( edge_masses.size() == tree.edge_count() );

    // Make color maps and norms.
    auto const mm    = matrix_minmax(matrix);
    auto matrix_norm = ColorNormalizationLinear( mm.min, mm.max );
    auto tree_norm   = ColorNormalizationLogarithmic( edge_masses.begin(), edge_masses.end() );

    auto write_heattree_ = [&](
        ColorMap const& map_mat, ColorMap const& map_tree, std::string const& mname
    ){

        // Make color matrix
        auto mat_col = Matrix<Color>( matrix.rows(), matrix.cols() );
        for( size_t r = 0; r < matrix.rows(); ++r ) {
            for( size_t c = 0; c < matrix.cols(); ++c ) {
                mat_col( r, c ) = map_mat( matrix_norm, matrix(r,c) );
            }
        }

        // Prepare heat tree
        HeatTreeParameters params;
        params.tree = tree;
        params.ladderize = false;
        params.type = LayoutType::kPhylogram;
        params.color_per_branch = map_tree( tree_norm, edge_masses );
        params.matrix = mat_col;

        // params.column_labels =
        auto const doc = heat_tree( params, map_mat , matrix_norm, map_tree, tree_norm );
        std::ostringstream out;
        doc.write( out );
        utils::file_write( out.str(), out_dir + "heat_tree_" + mname + ".svg" );
    };

    write_heattree_( ColorMap( color_list_viridis() ), ColorMap( color_list_viridis() ), "viridis" );
    write_heattree_( ColorMap( color_list_heat() ), ColorMap( color_list_heat() ), "heat" );
    write_heattree_( ColorMap( color_list_heat() ), ColorMap( color_list_viridis() ), "viridis_heat" );
    write_heattree_( ColorMap( color_list_heat() ), ColorMap( color_list_bupubk() ), "bupubk_heat" );
    // write_heattree_( ColorMap( color_list_orrd() ), "orrd" );
    // write_heattree_( ColorMap( color_list_ylorrd() ), "ylorrd" );
}

// =================================================================================================
//     write_heat_map_stuff
// =================================================================================================

void write_heat_map_stuff( PFData const& data, std::string const& out_dir )
{
    LOG_DBG << "write_heat_map_stuff";

    auto layout = RectangularLayout( data.label_tree, LayoutType::kPhylogram, false );
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

    write_heat_tree(
        data.bal_data.raw_edge_masses, data.col_sort_indices, data.label_tree, out_dir
    );
}

// =================================================================================================
//     write_factor_tree
// =================================================================================================

void write_factor_tree(
    PFData const& data, std::vector<PhyloFactor> const& factors, std::string const& out_dir
) {
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

    // dummy fill so that we have enough colors for experiments that want more than 10 factors.
    for( size_t i = clade_cols.clade_colors.size(); i < factors.size(); ++i ) {
        clade_cols.clade_colors.push_back( clade_colors[ i % clade_colors.size() ] );
    }

    auto all_edge_cols = phylo_factor_clade_colors( data.bal_data.tree, factors, 0, clade_cols );
    write_color_tree_to_svg_file(
        data.label_tree, data.lm, all_edge_cols,
        out_dir + "factors_tree.svg"
    );
}

// =================================================================================================
//     write_factor_edges
// =================================================================================================

void write_factor_edges(
    PFData const& data, std::vector<PhyloFactor> const& factors, std::string const& out_dir
) {
    auto const found_factors = factors.size();
    assert( data.col_sort_indices.size() == data.trees.size() );
    auto balances = utils::Matrix<double>(
        data.trees.size(), found_factors + data.meta_vals_mat.cols()
    );
    auto col_names = std::vector<std::string>( data.meta_vals_mat.cols() + found_factors );

    // Add meta data cols as well
    // balances.col(found_factors) = data.meta_vals;
    for( size_t c = 0; c < data.meta_vals_mat.cols(); ++c ) {
        balances.col( found_factors + c ) = data.meta_vals_mat.col(c).to_vector();
        col_names[ found_factors + c ] = data.col_names[c];
    }

    std::ofstream factor_taxa_of;
    file_output_stream( out_dir + "factor_taxa.txt", factor_taxa_of );

    for( size_t i = 0; i < found_factors; ++i ) {
        auto const& factor = factors[i];

        LOG_INFO << "edge " << factor.edge_index << " ( " << factor.edge_indices_primary.size() << " | " << factor.edge_indices_secondary.size() << " ), ov " << factor.objective_value;

        // Make a tree with the edges of that factor.
        auto edge_cols = phylo_factor_single_factor_colors( data.bal_data.tree, factors, i );
        write_color_tree_to_svg_file(
            data.label_tree, data.lm, edge_cols,
            out_dir + "factor_edges_" + std::to_string(i) + ".svg"
        );

        assert( factor.balances.size() == balances.rows() );
        balances.col(i) = factor.balances;
        col_names[ i ] = "factor_" + std::to_string(i);

        // write objective value trees
        auto cm = ColorMap( color_list_viridis() );
        cm.mask_color( Color( 0.8, 0.8, 0.8 ));
        write_color_tree_to_svg_file(
            data.label_tree, data.lm, factor.all_objective_values,
            cm, ColorNormalizationLinear( factor.all_objective_values ),
            out_dir + "ovs_" + std::to_string(i) + ".svg"
        );

        // write out which taxa are on the secondary part of the factor.
        factor_taxa_of << "Factor " << i << ":\n";
        std::unordered_set<std::string> edge_names;
        for( auto const ei : factor.edge_indices_secondary ) {
            auto const& ed = data.label_tree.edge_at(ei).secondary_link().node().data<CommonNodeData>();
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
        balances, out_dir + "factor_balances.csv", data.trees.names(), col_names, "sample"
    );
}

// =================================================================================================
//     Objective Functions
// =================================================================================================

double obj_fct_corr( PFData const& data, std::vector<double> const& balances )
{
    if( data.main_col_idx < 0 ) {
        throw std::runtime_error( "No main col given for correlation objective function." );
    }

    auto const& meta_vals = data.meta_vals_mat.col( data.main_col_idx );
    return std::abs( spearmans_rank_correlation_coefficient( meta_vals, balances ));
}

double obj_fct_lin_reg( PFData const& data, std::vector<double> const& balances )
{
    // return std::abs( spearmans_rank_correlation_coefficient( meta_vals, balances ));

    if( data.main_col_idx < 0 ) {
        throw std::runtime_error( "No main col given for simple linear regression." );
    }

    // Null deviance.
    auto const ms = mean_stddev( balances );
    auto const var = ms.stddev * ms.stddev;

    auto const& meta_vals = data.meta_vals_mat.col( data.main_col_idx );

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
}

double obj_fct_glm( PFData const& data, std::vector<double> const& balances )
{
    GlmExtras extras;
    // extras.mean_deviance = true;
    auto const fit = glm_fit( data.meta_vals_mat, balances, glm_family_gaussian(), extras );
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
        LOG_DBG << "data.meta_vals_mat " << join( data.meta_vals_mat );
        LOG_DBG << "balances " << join( balances );
    }

    return fit.null_deviance - fit.deviance;
}

double obj_fct_phyca( PFData const& data, std::vector<double> const& balances )
{
    assert( balances.size() == data.bal_data.edge_masses.rows() );
    auto var_vec = std::vector<double>( balances.size(), 0.0 );
    for( size_t i = 0; i < balances.size(); ++i ) {
        for( size_t j = 0; j < data.bal_data.edge_masses.cols(); ++j ) {
            var_vec[i] += balances[i] * data.bal_data.edge_masses( i, j );
        }
    }
    auto const ms = mean_stddev( var_vec );
    return ms.stddev * ms.stddev;
}

// =================================================================================================
//     main
// =================================================================================================

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.date = true;
    utils::Logging::details.time = true;
    utils::Options::get().number_of_threads( 40 );
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    // ----------------------------------------------------------------
    //     Input Files Tara Unconstr Euk Tree
    // ----------------------------------------------------------------

    // std::string in_jplace = "path/to/06_samples/tara_Euks_Unconstr_backbone/samples_jplace/";
    // std::string in_meta = "path/to/17_balances/tara_place/data_df.csv";
    // std::string in_label_file = "";
    // std::string main_col = "";
    // std::string out_dir = "path/to/17_balances/tara_place/";
    //
    // double line_width = 10.0;
    // std::string meta_sep_char = ",";
    // bool use_only_main_col = true;

    // ----------------------------------------------------------------
    //     Input Files Tara Tree Berneyse
    // ----------------------------------------------------------------

    // std::string in_jplace = "path/to/18_tara/05_samples/samples/";
    // std::string in_meta = "path/to/18_tara/06_pf/data_df.csv";
    // std::string in_label_file = "";
    // std::string main_col = "";
    // std::string out_dir = "path/to/18_tara/06_pf/";
    //
    // double line_width = 10.0;
    // std::string meta_sep_char = ",";

    // ----------------------------------------------------------------
    //     Input Files HMP
    // ----------------------------------------------------------------

    // std::string in_jplace = "path/to/06_samples/hmp_Bact_Unconstr_backbone/samples_jplace_of/";
    // std::string in_meta = "path/to/17_balances/hmp_pf_of/meta_9194_regions_of_clean.csv";
    // std::string in_label_file = "";
    // std::string main_col = "";
    // std::string out_dir = "path/to/17_balances/hmp_pf_of/";
    //
    // double line_width = 10.0;
    // std::string meta_sep_char = ",";

    // ----------------------------------------------------------------
    //     Input Files BV
    // ----------------------------------------------------------------

    // BV files.
    std::string in_label_file = "path/to/data/bv/meta/sequence_names.csv";
    std::string in_meta = "path/to/data/bv/meta/meta_simple.csv";
    // std::string in_meta = "path/to/data/bv/meta/meta_simple_nugent.csv";
    // std::string main_col = "nugent";
    std::string main_col = "";

    // Settings BV normal.
    std::string in_jplace = "path/to/03_bv/03_epa/orig_queries_jplace_clean";
    std::string out_dir = "path/to/philr_bv/factors/";

    double line_width = 10.0;
    std::string meta_sep_char = "\t";
    // bool use_only_main_col = true;
    bool use_only_main_col = false;

    // Settings BV small.
    // std::string in_jplace = "path/to/17_balances/bv_place_small/03_place/queries";
    // std::string out_dir = "path/to/17_balances/bv_place_small/04_pf/";
    // double line_width = 4.0;

    // ----------------------------------------------------------------
    //     Settings
    // ----------------------------------------------------------------

    BalanceSettings settings;
    settings.pseudo_count_summand_all = 0.65;

    if( argc < 2 || argc > 3 ) {
        throw std::runtime_error( "expecting 1 or 2 arguments" );
    }

    std::string weight_prefix;
    if( std::string( argv[1] ) == "no" ) {
        LOG_INFO << "no taxon weighing";
        settings.tendency = BalanceSettings::WeightTendency::kNone;
        settings.norm = BalanceSettings::WeightNorm::kNone;
        weight_prefix = "no_tw";
    } else if( std::string( argv[1] ) == "yes" ) {
        LOG_INFO << "with taxon weighing";
        weight_prefix = "tw";
    } else {
        throw std::runtime_error( "taxon weights yes or no" );
    }
    if( argc == 3 ) {
        main_col = std::string( argv[2] );
    }

    out_dir += ( main_col.empty() ? "" : main_col + "_" ) + weight_prefix + "/";
    LOG_INFO << "target dir: " << out_dir;
    dir_create( out_dir );

    if( main_col.empty() ) {
        use_only_main_col = false;
    }

    // if( main_col.empty() ) {
    //     LOG_INFO << "!!! running GLM";
    // } else {
    //     LOG_INFO << "!!! running SLR";
    // }

    size_t const num_factors = 10;

    // ----------------------------------------------------------------
    //     Run!
    // ----------------------------------------------------------------

    auto const data = read_data(
        in_jplace, in_meta, settings, in_label_file, main_col,
        line_width, meta_sep_char, use_only_main_col
    );

    calc_and_write_balances_mat_and_pca( data, out_dir );
    write_edge_masses_taxon_weights( data, out_dir );

    // write_edge_mass_matrix( data.bal_data.edge_masses, data.col_sort_indices, data.label_tree, out_dir );
    // write_edge_mass_matrix( data.bal_data.raw_edge_masses, data.col_sort_indices, data.label_tree, out_dir, "_raw" );
    write_heat_map_stuff( data, out_dir );

    LOG_INFO << "Running Phylo Factorization";
    auto const factors = phylogenetic_factorization(
        data.bal_data,
        [&]( std::vector<double> const& balances ){
            // return obj_fct_corr( data, balances );
            // return obj_fct_lin_reg( data, balances );
            // return obj_fct_phyca( data, balances );
            return obj_fct_glm( data, balances );

            // if( main_col.empty() ) {
            //     // LOG_INFO << "!!! running GLM";
            //     return obj_fct_glm( data, balances );
            // } else {
            //     // LOG_INFO << "!!! running SLR";
            //     return obj_fct_lin_reg( data, balances );
            // }
        },
        num_factors,
        []( size_t iteration, size_t max_iterations ){
            LOG_DBG1 << "iteration " << iteration << " of " << max_iterations;
        }
    );

    write_factor_tree(  data, factors, out_dir );
    write_factor_edges( data, factors, out_dir );

    LOG_INFO << "Finished";

    return 0;
}
