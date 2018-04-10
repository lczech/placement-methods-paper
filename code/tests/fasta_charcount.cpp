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
#include <map>

using namespace genesis;
using namespace genesis::sequence;

int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;

    utils::Options::get().number_of_threads( 4 );
    // LOG_BOLD << utils::Options::get().info();
    // LOG_BOLD;

    LOG_INFO << "Started";

    std::string infile = "data/silva/SILVA_123.1_SSURef_Nr99_tax_silva_full_align_trunc.fasta";
    utils::InputStream instr( utils::make_unique< utils::FileInputSource >( infile ));
    auto it = FastaInputIterator( instr );

    std::map<char, size_t> result;
    while( it ) {
        for( auto const& s : *it ) {
            ++result[s];
        }
        ++it;
    }

    for( auto e : result ) {
        LOG_DBG << "'" << e.first << "': " << e.second;
    }

    LOG_INFO << "Finished";
    return 0;
}
