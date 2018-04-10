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
using namespace genesis::sequence;

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    LOG_INFO << "Started " << utils::current_time();

    // Check if the command line contains the right number of arguments.
    if (argc != 2) {
        throw std::runtime_error(
            "Need to provide a fsa/fasta file.\n"
        );
    }

    // Input file
    auto infile = std::string( argv[1] );
    auto name = utils::file_filename( infile );;

    // Read file
    auto set = FastaReader().from_file( infile );

    // Write output
    std::ofstream out;
    out.open( name + ".tsv" );

    size_t count = 0;
    // std::unordered_map<std::string, size_t> counts;
    for( auto const& seq : set ) {
        auto const s = seq.label().find( "description=\"" );
        auto const e = seq.label().find( " ", s );
        if( e == std::string::npos ) {
            throw std::runtime_error( "Wrong label" );
        }

        auto const label = seq.label().substr( s + 13, e - s - 13 );
        // LOG_DBG << seq.label() << ": =" << label << "=";
        // break;

        auto sites = seq.sites();
        sites.erase(
            std::remove_if(
                sites.begin(), sites.end(),
                []( char chr ){
                    return chr != 'A' && chr != 'C' && chr != 'G' && chr != 'T';
                }
            ),
            sites.end()
        );

        // out << label << "\t" << utils::to_string_leading_zeros( counts[label], 4 ) << "\t";
        out << label << "\t" << utils::to_string_leading_zeros( count, 4 ) << "\t";

        out << sites << "\n";

        // ++counts[label];
        ++count;
    }

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
