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
using namespace genesis::tree;

/**
 * The placement tree produced by sativa (i i remember correctly) contains a prefix "r_"
 * for every reference taxon. remove them.
 */
int main( int argc, char** argv )
{
    (void) argc;
    (void) argv;

    // Activate logging.
    utils::Logging::log_to_stdout();

    Tree t = DefaultTreeNewickReader().from_file( "tax_pruning/00_ref/Gen_tax_constr_mfu.tre" );

    // LOG_BOLD << PrinterCompact().print(t);

    for( auto& n  : t.nodes() ) {
        auto& d = n->data<DefaultNodeData>();
        if( utils::starts_with( d.name, "r_" ) ) {
            d.name = d.name.substr( 2 );
        } else if( ! d.name.empty() ) {
            LOG_WARN << "Node name: " << d.name;
        }
    }

    auto newick_writer = DefaultTreeNewickWriter();
    newick_writer.enable_names(true);
    newick_writer.enable_branch_lengths(true);
    newick_writer.branch_length_precision(10);

    newick_writer.to_file( t, "tax_pruning/00_ref/constraint.newick" );

    return 0;
}
