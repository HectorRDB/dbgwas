/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

// We include the header file for the tool
#include "build_dbg.hpp"
#include "map_reads.hpp"
#include "statistical_test.h"
#include "generate_output.h"
#include "global.h"
#include "Utils.h"
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
/********************************************************************************/

int main (int argc, char* argv[])
{
    //get the template path
    pathToExecParent = argv[0];
    pathToExecParent.replace(pathToExecParent.rfind("DBGWAS"), 6, "");

    //TODO: we should remove this and deal better with the relative/absolute paths...
    //TODO: we should remove this and deal better with the relative/absolute paths...
    if (pathToExecParent!="./")
        fatalError("Sorry, for the moment, to execute this tool you have to cd to the bin folder and execute as ./DBGWAS <parameters>.");
    //TODO: we should remove this and deal better with the relative/absolute paths...
    //TODO: we should remove this and deal better with the relative/absolute paths...

    try
    {
        //Build DBG
        cerr << "Step 1. Building DBG and mapping strains on the DBG..." << endl;
        build_dbg().run(argc, argv); //this call will set up graph and nodeIdToUnitigId
        map_reads().run(argc, argv);
        cerr << "Done!" << endl;

        //Run the statistical test
        cerr << "Step 2. Running statistical test (bugwas + gemma)..." << endl;
        statistical_test().run(argc, argv);
        cerr << "Done!" << endl;

        //Find the neighbourhood around significant unitigs...
        cerr << "Step 3. Building visualisation around significant unitigs..." << endl;
        generate_output().run(argc, argv);
        cerr << "Done!" << endl;
    }
    catch (Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

