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
#include <vector>

using namespace genesis;
using namespace genesis::utils;

int main( int argc, char** argv )
{
    if (argc != 2) {
        throw std::runtime_error( "Need jplace file" );
    }
    utils::Logging::log_to_stdout();
    LOG_INFO << "Started " << utils::current_time();

    auto const file = std::string( argv[1] );
    auto cont = file_read( file );

    std::string const f_ckcomma = "\n,\n  ],\n";
    std::string const goodcomma = "\n  ],\n";

    auto const cnt = count_substring_occurrences( cont, f_ckcomma );
    LOG_INFO << "Found f_ckcomma " << cnt << " times.";

    if( cnt == 0 ) {
        // nothing. already replaced it
    } else if( cnt == 1 ) {
        cont = replace_all( cont, f_ckcomma, goodcomma );
        Options::get().allow_file_overwriting(true);
        file_write( cont, file );
    } else {
        LOG_WARN << "f_ckcomma occurs too often! baaaad!";
    }

    LOG_INFO << "Finished " << utils::current_time();
    return 0;
}
