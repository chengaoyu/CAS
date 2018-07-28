/*
This file is part of CASToR.

    CASToR is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.

    CASToR is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along with
    CASToR (in file GNU_GPL.TXT). If not, see <http://www.gnu.org/licenses/>.

Copyright 2017-2018 all CASToR contributors listed below:

    --> current contributors: Thibaut MERLIN, Simon STUTE, Didier BENOIT, Claude COMTAT, Marina FILIPOVIC, Mael MILLARDET
    --> past contributors: Valentin VIELZEUF

This is CASToR version 2.0.1.
*/

/*!
  \file
  \ingroup utils
  \brief This code allows the conversion of a GATE macro file defining a cylindricalPET or ecat geometry into a CASToR ASCII .geom file. \n
         This script currently only supports 'cylindricalPET' and 'ecat' systems. \n
         Two arguments are required : \n
         1. The gate .mac file : enter file.mac. We suppose that there is all the geometry and the digitizer information in this file. \n
         2. The output name : it will be a .geom file.
*/

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <math.h>

#include "gDataConversionUtilities.hh"

using namespace std;

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      ShowHelp()
  \param   a_returnCode
  \brief   Display main command line options for castor-GATEMacToGeom
*/
void showHelp(int a_returnCode)
{
  cout << endl;
  cout << " Usage:  castor-GATEMacToGeom -m input_mac_file" << endl;
  cout << "                              -o output_geom_file" << endl;
  cout << "                              -v verbose" << endl;
  cout << endl;
  cout << "[Input settings]:" << endl;
  cout << "  -m input_mac_file     : give the path to a GATE macro file describing a cylindricalPET, ecat or SPECThead system" << endl;
  cout << "  -o output_geom_file   : give the alias of the output geom file which will be generated in the scanner repository (default location: config/scanner/)" << endl;
  //cout << "  -v verbose            : verbosity level" << endl;
  cout << endl;
  #ifdef CASTOR_VERSION
  cout << "  This program is part of the CASToR release version " << CASTOR_VERSION << "." << endl;
  cout << endl;
  #endif
  Exit(a_returnCode);
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
// Basic checking and processing the input

int main(int argc, char *argv[])
{
  // Path to mac file
  string path_to_mac = "";
  // Path to the output geometric file
  string path_to_geom = "";
  
  if (argc < 3)
  {
    cerr << "Invalid number of arguments. Please read the help section." << endl;
    showHelp(0);
  }
  else
  {
    for (int i=1; i<argc; i++)
    {
      string option = (string)argv[i];

      // .mac file
      if (option=="-m")
      {
        if (argv[i+1] == NULL)
        {
          Cerr("***** castor-GATEMacToGeom :: argument missing for option: " << option << endl);
          Exit(EXIT_FAILURE);
        }
        else
          path_to_mac = argv[i+1];

        // Check if the file exists
        ifstream fcheck(path_to_mac.c_str());
        if(!fcheck.good())
        {
          Cerr("***** castor-GATEMacToGeom :: cannot read provided mac file : " << path_to_mac << endl);
          Exit(EXIT_FAILURE);
        }
        fcheck.close();
        i++;
      }
      
      // output file
      else if (option=="-o")
      {
        if (argv[i+1] == NULL)
        {
          Cerr("***** castor-GATEMacToGeom :: argument missing for option: " << option << endl);
          Exit(EXIT_FAILURE);
        }
        else
          path_to_geom = argv[i+1];
        i++;
      }
      
      else
      {
        Cerr(endl << "***** castor-GATERootToCastor :: Unknown option '" << option << "' !" << endl);
        showHelp(0);
        Exit(EXIT_FAILURE);
      }
    }
  }
  

  // ============================================================================================================
  // SOutputManager object initialisation:
  // ============================================================================================================
    
  sOutputManager* p_outputManager = sOutputManager::GetInstance();  

  // Set path to the config directory
  if (p_outputManager->CheckConfigDir(""))
  {
    Cerr("***** castor-GATEMacToGeom :: A problem occured while checking for the config directory path !" << endl);
    Exit(EXIT_FAILURE);
  }
  
  string scanner_repository = sOutputManager::GetInstance()->GetPathToConfigDir() + "scanner" + OS_SEP;
  path_to_geom = scanner_repository + path_to_geom + ".geom";
  
    cout << "path mac " << path_to_mac << endl;
  // ============================================================================================================
  // Get & check system type 
  // ============================================================================================================
  switch ( GetGATESystemType(path_to_mac) ) 
  {
    case GATE_SYS_CYLINDRICAL:
      cout << endl << " --- CylindricalPET system detected. Proceeding to conversion... --- " << endl << endl;
      
      if(CreateGeomWithCylindrical(path_to_mac , path_to_geom) )
      {
        Cerr("***** castor-GATEMacToGeom :: An error occured while trying to process mac file for cylindrical system: " << path_to_mac << endl);
        Exit(EXIT_FAILURE);
      }
      break;
      
    case GATE_SYS_ECAT:
      cout << endl << " --- ECAT system detected. Proceeding to conversion... --- " << endl;
      if(CreateGeomWithECAT(path_to_mac , path_to_geom) )
      {
        Cerr("***** castor-GATEMacToGeom :: An error occured while trying to process mac file for ecat system: " << path_to_mac  << endl);
        Exit(EXIT_FAILURE);
      }
      break;

    case GATE_SYS_SPECT:
      cout << endl << " --- SPECThead system detected. Proceeding to conversion... --- " << endl;
      if(CreateGeomWithSPECT(path_to_mac , path_to_geom) )
      {
        Cerr("***** castor-GATEMacToGeom :: An error occured while trying to process mac file for SPECT system: " << path_to_mac  << endl);
        Exit(EXIT_FAILURE);
      }
      break;
      
    default: // Unknown system
      Cerr("***** castor-GATEMacToGeom :: System type not supported : " << endl);
      Cerr("                   This script only supports conversion for cylindricalPET ecat and SPECThead systems"  << endl);
      Cerr("                   The system type is recovered from the lines '/gate/systems/...'"  << endl);
      Exit(EXIT_FAILURE);  
      break;
  }
  cout << endl << " --- Conversion completed --- " << endl << endl;
  return 0;
}
