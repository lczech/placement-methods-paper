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
#include <string>
#include <unordered_map>

using namespace genesis;
using namespace genesis::placement;

/**
 * Some tests of EMD on the BV data.
 */
int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    LOG_INFO << "Started " << utils::current_time();
    size_t count = 0;

    // In out dirs.
    std::string indir  = "bacterial_vaginosis/03_epa_magny/orig_queries/";
    std::string outdir = "bacterial_vaginosis/04_dists_b/";

    // Reading and writing fasta.
    auto reader = JplaceReader();
    SampleSet sample_set;
    tree::TreeSet tset;

    // Process all fasta files, filter them, chunk them.
    auto sample_names = utils::split( utils::file_read( indir + "../../00_data/sample_names" ), "\n" );
    for( auto const& sample_name : sample_names ) {
        if( count > 0 && count % 10 == 0 ) {
            LOG_DBG << "reading at " << count;
        }
        ++count;

        // Read and add to set.
        auto filename = indir + "sample_" + sample_name + "/RAxML_portableTree.orig_queries.jplace";
        auto sample = reader.from_file( filename );
        sample_set.add( sample, sample_name );
        tset.add( sample_name, sample.tree() );
    }
    auto avg_length_tree = tree::average_branch_length_tree( tset );
    LOG_INFO << "reading finished at " << count;
    count = 0;

    // Remap all samples to the avg branch length tree.
    for( size_t i = 0; i < sample_set.size(); ++i ) {
        auto new_smp = Sample( avg_length_tree );
        for( auto const& query : sample_set[i].sample.pqueries() ) {
            new_smp.add(query);
        }
        sample_set[i].sample = new_smp;
    }

    LOG_INFO << "Using " << utils::Options::get().number_of_threads() << " threads.";
    LOG_INFO << "Matrix calculation started at " << utils::current_time();
    // auto emd_matrix = utils::Matrix<double>( sample_set.size(), sample_set.size(), 0.0 );
    // auto nhd_matrix = utils::Matrix<double>( sample_set.size(), sample_set.size(), 0.0 );
    // auto nhd_matrix = node_histogram_distance( sample_set, 50 );

    auto emd_matrix = earth_movers_distance( sample_set );

    // for( size_t i = 0; i < sample_set.size(); ++i ) {
    //
    //     LOG_DBG << "EMD at " << count;
    //     // if( count > 0 && count % 10 == 0 ) {
    //     // }
    //     ++count;
    //
    //     // The result is symmetric - we only calculate the upper triangle.
    //     for( size_t j = i + 1; j < sample_set.size(); ++j ) {
    //         emd_matrix(i, j) = earth_movers_distance(
    //             sample_set[i].sample,
    //             sample_set[j].sample
    //         );
    //         emd_matrix(j, i) = emd_matrix(i, j);
    //
    //         // nhd_matrix(i, j) = node_histogram_distance(
    //         //     sample_set[i].sample,
    //         //     sample_set[j].sample,
    //         //     50
    //         // );
    //         // nhd_matrix(j, i) = nhd_matrix(i, j);
    //     }
    // }

    LOG_INFO << "Matrix calculation finished at " << utils::current_time();
    utils::file_write( utils::to_string( emd_matrix ), outdir + "emd.mat" );
    // utils::file_write( utils::to_string( nhd_matrix ), outdir + "nhd.mat" );

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
