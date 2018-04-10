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

#include "genesis/sequence.hpp"
#include "genesis/utils.hpp"

using namespace genesis;
using namespace genesis::sequence;

// =================================================================================================
//     Main
// =================================================================================================

int main( int argc, char** argv )
{
    // Stupidly simple CLI
    if (argc != 2) {
        throw std::runtime_error( "Need to provide the path to a fasta file." );
    }

    // Input
    auto const in_file = std::string( argv[1] );
    auto fasta_it = FastaInputIterator().from_file( in_file );

    // Process
    while( fasta_it ) {
        std::cout << fasta_it->label() << "\n";
        ++fasta_it;
    }

    return 0;
}
