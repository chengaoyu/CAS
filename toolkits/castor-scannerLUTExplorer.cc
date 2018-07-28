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
  \brief   This program can be used to visualize a scanner Look-Up-Table (geometric position and orientation of the system crystals/pixels), element by element.
  \details The program read a CASToR Look-Up-Table file (.hscan/.lut), or generate the LUT from a generic geometry file (.geom), and allow to explore it, element by element. \n
           For CT scanners or SPECT cameras described by a .geom file, a datafile must also be provided in order to get information regarding the projection angles.
*/

#include "gVariables.hh"
#include "gOptions.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "iDataFilePET.hh"
#include "iDataFileSPECT.hh"
#include "iDataFileCT.hh"
#include "sOutputManager.hh"
#include "sScannerManager.hh"
#include "sRandomNumberGenerator.hh"
#include "iScannerPET.hh"
#include "sAddonManager.hh"


/**
 * @defgroup INPUT_FILE_TYPE input file type
 *
 *    \brief Input file type for castor-scannerLUTExplorer (lut or geom) \n
 *           Defined in castor-scannerLUTExplorer.cc
 * @{
 */
/** Constant corresponding to an undefined type (default) (=-1) */
#define IF_TYPE_UNKNOWN -1
/** Constant corresponding to a Look-Up-Table file (.hscan/.lut) (=0) */
#define IF_TYPE_LUT 0
/** Constant corresponding to a generic geometry file (.geom )file (=1) */
#define IF_TYPE_GEOM 1
/** @} */


// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================
//                                                        H E L P     F U N C T I O N S
// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================


/*!
  \fn      ShowHelp()
  \brief   Display main command line options for castor-recon
*/
void ShowHelp()
{
  // Show help
  cout << endl;
  cout << "Usage:  castor-scannerLUTExplorer  -df datafile.cdh  [settings]" << endl;
  cout << endl;
  cout << "This program can be used to explore a datafile and get some info about it. By default, it simply prints general information recovered" << endl;
  cout << "from the reader. If the '-e' option is supply, an element-by-element exploration will be performed: information about each element will be" << endl;
  cout << "displayed." << endl;
  cout << "For CT scanners or SPECT cameras described by a .geom file, a datafile must also be provided in order to get information regarding the projection angles." << endl;
  cout << endl;
  cout << "[Mandatory parameters]:" << endl;
  cout << "  -sf   scan_file      : Give the path to the header of a single lut file, or to a geom file." << endl;
  cout << endl;
  cout << "[Options]:" << endl;
  cout << "  -e                   : Flag for scanner element by element exploration (will list the content of each element)." << endl;
  cout << "  -i                   : Flag for interactive one by one scanner element exploration when -e option is supplied." << endl;
  cout << "  -g                   : Flag for global exploration (geometric characteristics of all elements will be displayed)." << endl;
  cout << "  -o   path_out_file   : Path to an output file in which the data will be printed." << endl;
  cout << "  -df  data_file       : For CT or SPECT system, using a .geom file as input. Give the path to the header of a datafile in order to get projection informations." << endl;
  cout << "  -vb value            : Give the verbosity level, from 0 (no verbose) to 2 (default: 1)" << endl;
  cout << "  --help,-h,-help      : Print out this help page." << endl; // managed by main
  cout << endl;
  #ifdef BUILD_DATE
  cout << "  Build date: " << BUILD_DATE << endl;
  cout << endl;
  #endif
  #ifdef CASTOR_VERSION
  cout << "  This program is part of the CASToR release version " << CASTOR_VERSION << "." << endl;
  cout << endl;
  #endif
}




/*!
  \fn      GetInputFileType(string& ap_pathToScanFilename, string& ap_pathToBinary, int& ap_modality, int& ap_inputFileType)
  \param   ap_pathToScanFilename : path to the provided scanner file
  \param   ap_pathToBinary : path to a binary file
  \param   a_modality : system modality flag
  \param   ap_inputFileType : type of the scanner file (lut/geom)
  \brief   Identify type of the provided scanner file
*/
void GetInputFileType(string& ap_pathToScanFilename, string& ap_pathToBinary, int a_modality, int& ap_inputFileType)
{
  // Check if we have a LUT or GEOM file
  
  // checking for LUT
  if(ap_pathToBinary.find(".hscan") != string::npos
  || ap_pathToBinary.find(".lut")   != string::npos)
  {
    ap_inputFileType = IF_TYPE_LUT;
    if( ap_pathToBinary.find_last_of(".hscan") != string::npos ) // header has been provided
      ap_pathToBinary = ap_pathToBinary.substr(ap_pathToBinary.find_last_of(".hscan")) + ".lut";
    else // lut has been provided
      ap_pathToScanFilename = ap_pathToBinary.substr(ap_pathToBinary.find_last_of(".lut")) + ".hscan";
    
    Cout("  A scanner LUT binary file has been provided " << endl);
  }
  // Checking for GEOM (look for a mandatory key)
  else
  {
    // Try to recover a mandatory geom file key, to check whether we have a geom scanner file
    FLTNB key = 0.;
    string key_str = "";
    
    // Get a key specific to the modality
    if (a_modality == SCANNER_PET)
      key_str = "number of rsectors";
    else if (a_modality == SCANNER_SPECT_CONVERGENT)
      key_str = "number of detector heads";
    else if (a_modality == SCANNER_CT)
      key_str = "detector radius";
    else
    {
      Cerr("***** castor-scannerLUTExplorer() -> Error, this modality is not implemented yet !" << endl);
      Exit(EXIT_FAILURE);
    }
    
    int rvalue = ReadDataASCIIFile(ap_pathToScanFilename, key_str, &key, 1, KEYWORD_MANDATORY);
    
    if (rvalue > 1)
    {
      Cerr("***** castor-scannerLUTExplorer() -> An error occurred while trying to read a mandatory parameter for a .geom file in the scanner header file !" << endl);
      Exit(EXIT_FAILURE);
    }
    else if (rvalue == 1)
    {
      Cerr("***** castor-scannerLUTExplorer() -> The provided file has not been identified as a scanner binary file (no .hscan or .lut extension), nor a .geom file (a mandatory key is missing) !" << endl);
      Exit(EXIT_FAILURE);
    }
    else
    {
      ap_inputFileType = IF_TYPE_GEOM;
      Cout("  A generic geometry file has been provided " << endl);
    }
  }
}




// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================
//                                                        M A I N     P R O G R A M
// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================

int main(int argc, char** argv)
{
  // ============================================================================================================
  // MPI stuff (we make all instances but the first one returning 0 directly)
  // ============================================================================================================
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  int mpi_size = 1;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  if (mpi_rank!=0) return 0;
  #endif

  // No argument, then show help
  if (argc==1)
  {
    ShowHelp();
    Exit(EXIT_SUCCESS);
  }

  // ============================================================================================================
  // The few parameters
  // ============================================================================================================

  // String of a single path to a scanner file name
  string path_to_scan_filename = "";
  
  // String of a single path to a data file name (required for SPECT/CT with geom)
  string path_to_data_filename = "";
  
  // Modality
  int modality = SCANNER_UNKNOWN;

  // General verbose level
  int verbose = 1;
  // Flag for element-by-element exploration
  bool elt_by_elt_flag = false;
  // Flag for global exploration
  bool global_flag = false;
  // Flag for input file type
  int input_file_type = IF_TYPE_UNKNOWN;
  // Output file (log) string
  string path_fout = "";


  // ============================================================================================================
  // Read command-line parameters
  // ============================================================================================================

  // Must manually increment the option index when an argument is needed after an option
  for (int i=1; i<argc; i++)
  {
    // Get the option as a string
    string option = (string)argv[i];
    // Show help
    if (option=="-h" || option=="--help" || option=="-help")
    {
      ShowHelp();
      Exit(EXIT_SUCCESS);
    }
    // General verbosity level
    else if (option=="-vb")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-scannerLUTExplorer() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose))
      {
        Cerr("***** castor-scannerLUTExplorer() -> Exception when trying to read provided verbosity level '" << verbose << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // ScannerFile
    else if (option=="-sf") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-scannerLUTExplorer() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_to_scan_filename = (string)argv[i+1];
      i++;
    }
    // DataFile
    else if (option=="-df") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-scannerLUTExplorer() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_to_data_filename = (string)argv[i+1];
      i++;
    }
    // Element by element exploration
    else if (option=="-e")
    {
      elt_by_elt_flag = true;
    }
    // Global exploration
    else if (option=="-g")
    {
      global_flag = true;
    }
    else if (option=="-o")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-scannerLUTExplorer()) -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_fout = argv[i+1];
      i++;
    }
    
    // Unknown option!
    else
    {
      Cerr("***** castor-scannerLUTExplorer() -> Unknown option '" << option << "' !" << endl);
      Exit(EXIT_FAILURE);
    }
  }

  // ============================================================================================================
  // Some checks
  // ============================================================================================================

  // Check for a datafile name
  if (path_to_scan_filename=="")
  {
    Cerr("***** castor-scannerLUTExplorer() -> Please provide a datafile name !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Verbose
  if (verbose<=0)
  {
    Cerr("***** castor-scannerLUTExplorer() -> Verbose less than 1 has no sense for a program used to solely print out information on screen !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Restrict verbose level to 2 to avoid having to much useless information during initialization
  if (verbose>2) verbose = 2;
  // Cannot have element-by-element and global both enables
  if (global_flag && elt_by_elt_flag)
  {
    Cerr("***** castor-scannerLUTExplorer() -> Cannot use global and element-by-element mode at the same time !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ----------------------------------------------------------------------------------------
  // Check input file type if LUT or geom
  // ----------------------------------------------------------------------------------------
  
  
  // Check if we have two files for the scanner basename (.hscan/.lut) -> we have a lut
  // Get scanner basename 
  string path_to_binary = path_to_scan_filename;
  
  
  // ============================================================================================================
  // Initializations
  // ============================================================================================================

  // Verbose (we know it is at least 1 so we don't check it)
  Cout("==============================================================" << endl);
  Cout("castor-scannerLUTExplorer() -> Initialization starts" << endl);

  // Get user endianness (interfile I/O)
  GetUserEndianness();
  
  // ----------------------------------------------------------------------------------------
  // Create sScannerManager
  // ----------------------------------------------------------------------------------------

  // Create the scanner manager
  sScannerManager* p_ScannerManager = sScannerManager::GetInstance();  
  p_ScannerManager->SetVerbose(verbose);
  
  // First, recover the modality
  string modality_str = "";
  
  if( ReadDataASCIIFile(path_to_scan_filename, "modality", &modality_str, 1, KEYWORD_MANDATORY) )
  {
    Cerr("***** castor-scannerLUTExplorer() -> Error when trying to read the 'modality' key in the scanner file !" << endl);
    Exit(EXIT_FAILURE);
  }
  
  modality = p_ScannerManager->GetModalityFromString(modality_str);
  
  // Identify input file type and system 
  GetInputFileType(path_to_scan_filename, path_to_binary, modality, input_file_type);
  
  // Initialized output Manager, if output log is enabled
  if(path_fout != "") 
  {
    sOutputManager* p_outputManager; 
    p_outputManager = sOutputManager::GetInstance();
      
    p_outputManager->SetVerbose(verbose);
  
    //p_outputManager->SetDataFileName(vpath_to_initial_img);
    
    // Set path to the config directory
    if (p_outputManager->CheckConfigDir(""))
    {
      Cerr("***** castor-proj() -> A problem occured while checking for the config directory path !" << endl);
      Exit(EXIT_FAILURE);
    }
    
    // Initialize output directory and base name
    if (p_outputManager->InitOutputDirectory(path_fout, ""))
    {
      Cerr("***** castor-proj() -> A problem occured while initializing output directory !" << endl);
      Exit(EXIT_FAILURE);
    }
    // Log command line
    if (p_outputManager->LogCommandLine(argc,argv))
    {
      Cerr("***** castor-proj() -> A problem occured while logging command line arguments !" << endl);
      Exit(EXIT_FAILURE);
    }
  }
  



  // Get system name from the dataFile
  string scanner_name = "";
  if (ReadDataASCIIFile(path_to_scan_filename, "scanner name", &scanner_name, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-scannerLUTExplorer() -> A problem occured while trying to find the system name in the datafile header !" << endl);
    Exit(EXIT_FAILURE);
  } 
  
  
  if (p_ScannerManager->InitScannerWithFile(path_to_scan_filename, scanner_name, input_file_type) )
  {
    Cerr("***** castor-scannerLUTExplorer() -> A problem occurred during scanner object initialization ! !" << endl);
    Exit(EXIT_FAILURE);
  }
  
  if (p_ScannerManager->BuildScannerObject() )
  {
    Cerr("***** castor-scannerLUTExplorer() -> A problem occurred during scanner object construction ! !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->InstantiateScanner() )
  {
    Cerr("***** castor-scannerLUTExplorer() -> A problem occurred while creating Scanner object !" << endl);
    Exit(EXIT_FAILURE);
  }
  
  // For SPECT / CT with geom, check that a datafile has been provide (mandatory)
  if(input_file_type == IF_TYPE_GEOM && p_ScannerManager->GetScannerType() != SCANNER_PET )
  {
    if(path_to_data_filename == "")
    {
      Cerr("***** castor-scannerLUTExplorer() -> Error : For SPECT / CT systems, a datafile is required to compute the look-up-table from a geom file (i.e projection angles) !" << endl);
      Exit(EXIT_FAILURE);
    }
    else if (p_ScannerManager->GetGeometricInfoFromDataFile(path_to_data_filename))
    {
      Cerr("***** castor-scannerLUTExplorer() -> A problem occurred while retrieving scanner fields from the datafile header !" << endl);
      Exit(EXIT_FAILURE);
    }
  }
  if (p_ScannerManager->BuildLUT() )
  {
    Cerr("***** castor-scannerLUTExplorer() -> A problem occurred while generating/reading the LUT !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->CheckParameters())
  {
    Cerr("***** castor-scannerLUTExplorer() -> A problem occured while checking scanner manager parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->Initialize())
  {
    Cerr("***** castor-scannerLUTExplorer() -> A problem occured while initializing scanner !" << endl);
    Exit(EXIT_FAILURE);
  }

  // Verbose
  Cout("castor-scannerLUTExplorer() -> End of initialization" << endl);
  Cout("==============================================================" << endl);

  // ============================================================================================================
  // Actions
  // ============================================================================================================

  // -------------------------------------------------------------------
  // Ask the scanner to describe itself
  // -------------------------------------------------------------------
  p_ScannerManager->Describe();
  Cout("==============================================================" << endl);

  // -------------------------------------------------------------------
  // Then if the element-by-element option is on, describe each element
  // -------------------------------------------------------------------

  if (elt_by_elt_flag)
  {
    Cout("castor-scannerLUTExplorer() -> Start exploration of all scanner elements" << endl);
    // Get index start and stop
    int64_t index_start = 0; int64_t index_stop = p_ScannerManager->GetSystemNbElts();

    // Launch the loop on all elements
    int64_t index = index_start;
    while (index>=index_start && index<index_stop)
    {
      // Verbose
      Cout("-------------------------  Element index " << index << "  -------------------------" << endl);
      
      // Get the current element
      FLTNB position1[3], position2[3];
      FLTNB orientation1[3], orientation2[3];
                                             
      p_ScannerManager->GetScannerObject()->GetPositionsAndOrientations(0, index, position1, position2, orientation1, orientation2);
      
      // Describe the element
      Cout("Scanner element center location (x,y,z): "<< position2[0] << " ; " << position2[1] << " ; " << position2[2] << " ." ;);
      Cout("Orientation (x,y,z): "<< orientation2[0] << " ; " << orientation2[1] << " ; " << orientation2[2] << endl;);
      
      // The user may provide a specific element index or simply press enter to get to the next
      cout << "--------> Give an element index or simply press enter for next element: " << flush;
      // Read the answer
      string answer = ""; getline(cin,answer);
      cout << endl;
      // If empty answer, then go to the next element
      if (answer=="") index++;
      // Otherwise
      else
      {
        // Convert the answer to int64_t
        int64_t next_index = stoll(answer);
        if (next_index<index_start || next_index>=index_stop)
        {
          Cerr("***** castor-scannerLUTExplorer() -> The provided element index (" << next_index << ") is out of datafile range"
            << " [" << index_start << ":" << index_stop << "[ !" << endl);
          break;
        }
        else index = next_index;
      }

    }
    Cout("==============================================================" << endl);
  }

  // -------------------------------------------------------------------
  // Or perform a global exploration
  // -------------------------------------------------------------------

  else if (global_flag)
  {
    // Verbose
    Cout("castor-scannerLUTExplorer() -> Start global exploration of all elements" << endl);

    // Get index start and stop
    int64_t index_start = 0; int64_t index_stop = p_ScannerManager->GetSystemNbElts();

    // Launch the loop on all elements
    for (int64_t index=index_start ; index<index_stop ; index++)
    {
      // Get the current element
      FLTNB position1[3], position2[3];
      FLTNB orientation1[3], orientation2[3];
                                             
      p_ScannerManager->GetScannerObject()->GetPositionsAndOrientations(0, index, position1, position2, orientation1, orientation2);

      Cout("Scanner element center location (x,y,z): "<< position2[0] << " ; " << position2[1] << " ; " << position2[2] << " ." ;);
      Cout("Orientation (x,y,z): "<< orientation2[0] << " ; " << orientation2[1] << " ; " << orientation2[2] << endl;);
    }

    Cout("==============================================================" << endl);
    
    // -------------------------------------------------------------------
    // Ask the scanner to describe itself
    // -------------------------------------------------------------------
    p_ScannerManager->Describe();
    Cout("==============================================================" << endl);
  }

  // ============================================================================================================
  // End
  // ============================================================================================================
    
  // Ending
  #ifdef CASTOR_MPI
  MPI_Finalize();
  #endif
  return EXIT_SUCCESS;
}
