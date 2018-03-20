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

#include "genesis/utils/core/logging.hpp"

#include "genesis/taxonomy/formats/taxopath_generator.hpp"
#include "genesis/taxonomy/formats/taxopath_parser.hpp"
#include "genesis/taxonomy/functions/taxonomy.hpp"
#include "genesis/taxonomy/functions/taxopath.hpp"
#include "genesis/taxonomy/iterator/preorder.hpp"
#include "genesis/taxonomy/iterator/postorder.hpp"
#include "genesis/taxonomy/taxon.hpp"
#include "genesis/taxonomy/taxonomy.hpp"
#include "genesis/taxonomy/taxopath.hpp"

#include "genesis/utils/core/fs.hpp"
#include "genesis/utils/core/std.hpp"
#include "genesis/utils/formats/csv/reader.hpp"
#include "genesis/utils/io/input_stream.hpp"
#include "genesis/utils/text/string.hpp"

#include "genesis/utils/math/histogram.hpp"
#include "genesis/utils/math/histogram/accumulator.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace genesis;
using namespace genesis::taxonomy;

/**
 * @brief This was a short test to count some sizes in the silva taxonomy.
 */
int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();

    auto tax = Taxonomy();
    auto map = std::unordered_map< std::string, size_t >();
    auto sum_map = std::unordered_map< std::string, size_t >();

    auto csv_reader = utils::CsvReader();
    csv_reader.separator_chars( "\t" );

    auto taxopath_parser = TaxopathParser();

    std::string infile = "path/to/silva/taxmap_slv_ssu_ref_nr_123.1.txt";
    std::string outfile = "tax_vec.txt";

    std::ifstream ifs (infile);
    // auto it = utils::CountingIstream(ifs);
    utils::InputStream it( utils::make_unique< utils::StreamInputSource >( ifs ));

    LOG_TIME << "start reading taxonomy";
    while( it ) {
        // if( it.line() % 10000 == 0 ) {
        //     LOG_DBG << "At line " << it.line();
        // }

        auto fields = csv_reader.parse_line( it );
        if( fields.size() < 3 ) {
            LOG_WARN << "Malformed taxon at line " << it.line() - 1;
            continue;
        }

        auto tax_name = fields[3];
        if( tax_name.empty() || tax_name.back() != ';' ) {
            LOG_WARN << "Malformed taxon at line " << it.line() - 1 << ": " << tax_name;
            continue;
        }

        // tax_name += fields[4];

        add_from_taxopath(
            tax,
            taxopath_parser.from_string( tax_name )
        );

        ++map[ tax_name ];
    }

    LOG_TIME << "finished reading taxonomy";
    sort_by_name( tax );
    // std::cout << tax;

    // LOG_DBG << "validate: " << ( validate( tax ) ? "true" : "false" );

    auto make_sum_map = [&] ( Taxon& t ) {
        auto gen = TaxopathGenerator();
        auto name = gen(t);

        sum_map[ name ] = map[ name + ";" ];
        for( auto const & c : t ) {
            auto cname = gen(c);

            if( sum_map.count( cname ) > 0 ) {
                sum_map[ name ] += sum_map[ cname ];
            } else {
                LOG_WARN << "postorder broken";
                LOG_DBG << "at " << name;
                LOG_DBG << "with c " << cname;

                return;
            }
        }
    };
    postorder_for_each( tax, make_sum_map );

    LOG_TIME << "write tax_count_file";
    std::ofstream tax_count_file;
    tax_count_file.open( "tax_counts" );

    auto print_with_count = [&] ( Taxon& t ) {
        auto gen = TaxopathGenerator();
        auto name = gen(t);
        if( map.count( name + ";" ) == 0 ) {
            map[ name + ";" ] = 0;
            // tax_count_file << "not in map: " << name << ";" << "\n";
            // tax_count_file << name << "\t0\n";
            // return;
        }

        tax_count_file << name;
        tax_count_file << "\t" << taxon_level(t);
        tax_count_file << "\t" << map[ name + ";" ];
        tax_count_file << "\t" << sum_map[ name ];
        tax_count_file << "\n";

        // size_t child_count = total_taxa_count(t);
        // tax_count_file << "\t" << child_count  << "\n";
    };

    preorder_for_each( tax, print_with_count );
    tax_count_file.close();
    LOG_TIME << "finished tax_count_file";

    auto vec = std::vector< std::pair<std::string, int> >(
        map.begin(),
        map.end()
    );

    auto sort_by_val = []( std::pair<std::string, int> const& lhs, std::pair<std::string, int> const& rhs ) {
        return lhs.second > rhs.second;
    };

    std::sort( vec.begin(), vec.end(), sort_by_val );

    LOG_DBG << "map has " << map.size() << " entries.";
    LOG_DBG << "vec has " << vec.size() << " entries.";

    // Output sorted by freq
    LOG_TIME << "write tax_vec_file";
    std::ofstream tax_vec_file;
    tax_vec_file.open( outfile );

    auto ha = utils::HistogramAccumulator();
    for( auto const& p : vec ) {
        tax_vec_file << p.first << "  -->  " << p.second << "\n";
        ha.accumulate( p.second, 1.0 );
    }
    // for( auto const& e : ha ) {
    //     tax_vec_file << e.first << "\t" << e.second << "\n";
    // }
    tax_vec_file.close();

    LOG_TIME << "finished tax_vec_file";
    LOG_INFO << "Finished.";
    return 0;
}
