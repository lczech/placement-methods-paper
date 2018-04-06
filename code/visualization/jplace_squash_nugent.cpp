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

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

// =================================================================================================
//     Write Squash Tree To SVG
// =================================================================================================

/**
 * @brief Write an SVG file containing a tree with colored branches.
 */
void write_squash_tree_to_svg(
    placement::PlacementTree const&  tree,
    std::vector<utils::Color> const& colors_per_node,
    std::string const&               svg_filename
) {
    auto copy = tree;
    ladderize(copy);
    // LOG_DBG << "Ladderize tree: " << ( validate_topology(copy) ? "valid!" : "INVALID!" );

    // Make a layout tree.
    // auto layout = tree::RectangularLayout( tree );
    auto layout = tree::RectangularLayout( copy, tree::LayoutType::kPhylogram );
    // layout.scaler_x( 300 );

    // Set colourful node shapes.
    if( colors_per_node.size() != tree.node_count() ) {
        throw std::runtime_error( "Node colors vec of wrong size" );
    }

    std::vector<utils::SvgGroup> node_shapes;
    node_shapes.resize( tree.node_count() );
    for( size_t i = 0; i < tree.node_count(); ++i ) {
        if( colors_per_node[i] == Color(1,0,1) ) {
            continue;
        }
        // Add a shape according to the sample color.
        node_shapes[i].add( utils::SvgCircle(
            utils::SvgPoint( 0, 0 ),
            6.0,
            utils::SvgStroke( SvgStroke::Type::kNone ),
            utils::SvgFill( colors_per_node[i] )
        ));
    }
    layout.set_node_shapes( node_shapes );

    // Write to file
    std::ostringstream out;
    layout.to_svg_document().write( out );
    utils::file_write( out.str(), svg_filename + ".svg" );
}

void make_squash_tree(
    std::unordered_map<std::string, double>  const& nugent,
    DefaultTree const& sc_tree,
    std::vector<std::string> const& sample_names,
    std::string const& outfile
)
{
    // auto const& ref_assignments = massignments;
    // auto const& ref_assignments = iassignments;

    if( nugent.size() != leaf_node_count(sc_tree) ) {
        LOG_WARN << nugent.size() << " == nugent.size() != leaf_node_count(sc_tree) == " << leaf_node_count(sc_tree);
        return;
    }

    if( nugent.size() != sample_names.size() ) {
        LOG_WARN << nugent.size() << " == nugent.size() != sample_names.size() == " << sample_names.size();
        return;
    }

    auto norm = utils::ColorNormalizationLinear( 0.0, 10.0 );
    auto map = ColorMap({ Color( 0, 0, 1 ), Color( 1, 0, 0 ) });
    // auto colors_per_node = map( norm, nugent );

    auto colors_per_node = std::vector<utils::Color>( sc_tree.node_count(), Color(1,0,1) );
    size_t tmp_sum = 0;
    for( size_t i = 0; i < sc_tree.node_count(); ++i ) {
        auto const& name = sc_tree.node_at(i).data<DefaultNodeData>().name;
        if( name.empty() ) {
            continue;
        }

        // Fake test to see if it is a num
        bool is_num = false;
        try{
            // LOG_DBG << stoi(name);
            tmp_sum += stoi(name);
            is_num = true;
        } catch(...){}
        if( is_num ) {
            continue;
        }

        // find index in sample name list
        // auto match = std::find(sample_names.begin(), sample_names.end(), name);
        // size_t index;
        // if(match != sample_names.end()) {
        //     index = match - sample_names.begin();
        // } else {
        //     LOG_WARN << "cannot find sample name " << name;
        //     continue;
        // }
        //
        // colors_per_node[ i ] = map( norm, nugent[index] );

        colors_per_node[ i ] = map( norm, nugent.at(name) );
    }

    //  do not write node names
    auto cpy = sc_tree;
    for( auto& node : cpy.nodes() ) {
        node->data<DefaultNodeData>().name = "";
    }

    // write_color_tree_to_nexus( cpy, {}, outdir + "/cluster_kmeans.nexus" );
    write_squash_tree_to_svg( cpy, colors_per_node, outfile );
}

std::unordered_map<std::string, double> nugent_map( std::string const& csv_file )
{
    std::unordered_map<std::string, double> result;

    auto csv_reader = CsvReader();
    auto table = csv_reader.from_file( csv_file );
    for( auto const& line : table ) {
        if( line.size() != 4 ) {
            throw std::runtime_error( "line.size() != 4" );
        }

        result[ line[0] ] = stod( line[3] );
    }
    return result;
}

// =================================================================================================
//      Main
// =================================================================================================

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;

    utils::Options::get().number_of_threads( 4 );
    // LOG_BOLD << utils::Options::get().info();
    // LOG_BOLD;

    LOG_INFO << "Started";

    // Check if the command line contains the right number of arguments.
    if (argc != 4) {
        throw std::runtime_error(
            "Need to provide three arguments.\n"
        );
    }

    // In out dirs.
    auto epadir = utils::dir_normalize_path( std::string( argv[1] ));
    auto nugentfile = std::string( argv[2] );
    auto outdir = utils::dir_normalize_path( std::string( argv[3] ));
    utils::dir_create(outdir);

    utils::Options::get().command_line( argc, argv );
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;


    // -------------------------------------------------------------------------
    //     Load meta data
    // -------------------------------------------------------------------------

    auto nugents = nugent_map( nugentfile );

    // for( auto const& n : nugents ) {
    //     LOG_DBG << n.first << " " << n.second;
    // }

    // -------------------------------------------------------------------------
    //     Find and read sample jplace files.
    // -------------------------------------------------------------------------

    // Get all jplace files, either in the epa dir, or in its subdirs
    auto jplace_filenames = utils::dir_list_files( epadir, ".*\\.jplace" );
    std::sort( jplace_filenames.begin(), jplace_filenames.end() );

    // Output list of samples in the order that we use them for the matrix
    // LOG_INFO << "Writing jplace sample order list";
    // std::ofstream jplace_files_of( outdir + "jplace_order.list" );
    // for( auto const& elem : jplace_filenames ) {
    //     jplace_files_of << elem << "\n";
    // }
    // jplace_files_of.close();
    // LOG_INFO << "done";
    // LOG_INFO;

    SampleSet sset;
    JplaceReader jplace_reader;

    // Process all jplace files.
    LOG_INFO << "reading " << jplace_filenames.size() << " jplace sample files";
    sset = jplace_reader.from_files( jplace_filenames );

    // size_t file_count = 0;
    // for( auto const& jplace_filename : jplace_filenames ) {
    //
    //     // Progress output.
    //     // LOG_INFO << "at jplace file " << file_count << " (" << jplace_filename << ")";
    //     ++file_count;
    //
    //     sset.add( jplace_reader.from_file( epadir + jplace_filename ), jplace_filename );
    // }

    // Final output for jplace reading
    assert( sset.size() == jplace_filenames.size() );
    LOG_INFO << "finished reading sample jplace files";
    LOG_INFO;

    // -------------------------------------------------------------------------
    //     Squash
    // -------------------------------------------------------------------------

    // Options::get().random_seed( 1623390399 );

    std::vector<std::string> sample_names;
    for( auto const& smp : sset ) {
        sample_names.push_back( smp.name );
        // LOG_DBG << smp.name;
    }

    LOG_INFO << "Converting Trees";
    auto mass_trees = convert_sample_set_to_mass_trees( sset );
    // sset.clear();

    LOG_INFO << "Starting squash clustering";
    auto sc = tree::SquashClustering();
    sc.run( std::move( mass_trees.first ) );
    LOG_INFO << "Finished squash clustering";

    // LOG_DBG << "Clusters";
    // for( size_t i = 0; i < sc.clusters.size(); ++i ) {
    //     LOG_DBG1 << i << ": " << ( sc.clusters[i].active ? "act " : "dea " )
    //              << sc.clusters[i].count << " "
    //              << sc.clusters[i].distances.size() << " "
    //              << sc.clusters[i].tree.node_count();
    // }
    //
    // LOG_DBG << "Mergers";
    // for( size_t i = 0; i < sc.mergers.size(); ++i ) {
    //     LOG_DBG1 << i << ": "
    //              << sc.mergers[i].index_a << " " << sc.mergers[i].distance_a << " | "
    //              << sc.mergers[i].index_b << " " << sc.mergers[i].distance_b;
    // }

    // Remove jplace part from file names
    // for( auto& fn : jplace_filenames ) {
    //     fn = utils::replace_all( fn, ".jplace", "" );
    // }

    LOG_INFO << "Writing cluster tree";
    std::ofstream file_clust_tree;
    utils::file_output_stream( outdir + "/cluster.newick",  file_clust_tree );
    auto const sct_str = sc.tree_string( sample_names );
    file_clust_tree << sct_str;

    // LOG_INFO << "Writing fat trees";
    // for( size_t i = 0; i < sc.clusters.size(); ++i ) {
    //     auto const& cc = sc.clusters[i];
    //
    //     auto const cv = mass_tree_mass_per_edge( cc.tree );
    //     auto const colors = counts_to_colors(cv, false);
    //
    //     write_color_tree_to_nexus( avg_tree, colors, outdir + "/tree_" + std::to_string(i) + ".nexus" );
    //     write_color_tree_to_svg( avg_tree, colors, outdir + "/tree_" + std::to_string(i) );
    // }

    // -------------------------------------------------------------------------
    //     Squash Tree with kmeans node colors
    // -------------------------------------------------------------------------

    auto sc_tree = DefaultTreeNewickReader().from_string( sct_str );

    // auto const& ref_assignments = massignments;
    // auto const& ref_assignments = iassignments;

    make_squash_tree( nugents, sc_tree, sample_names, outdir + "/cluster_nugent" );

    LOG_INFO << "Finished ";
    return 0;
}
