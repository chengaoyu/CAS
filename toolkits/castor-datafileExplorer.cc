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
  \brief  This program can be used to explore a datafile, event by event.
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
  cout << "Usage:  castor-datafileExplorer  -df datafile.cdh  [settings]" << endl;
  cout << endl;
  cout << "This program can be used to explore a datafile and get some info about it. By default, it simply prints general information recovered" << endl;
  cout << "from the reader. If the '-e' option is supply, an event-by-event exploration will be performed: information about each event will be" << endl;
  cout << "displayed." << endl;
  cout << endl;
  cout << "[Mandatory parameters]:" << endl;
  cout << "  -df   datafile.cdh   : Give the path to a single datafile." << endl;
  cout << endl;
  cout << "[Options]:" << endl;
  cout << "  -e                   : Flag for event by event exploration (will list the content of each event)." << endl;
  cout << "  -i                   : Flag for interactive one by one event exploration when -e option is supplied." << endl;
  cout << "  -g                   : Flag for global exploration (will not list the content of each event but will supply a summary)." << endl;
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

  // String of a single path to a data file name
  string path_to_data_filename;
  // General verbose level
  int verbose = 1;
  // Flag for event-by-event exploring
  bool event_by_event_flag = false;
  // Flag for interactive mode for event-by-event
  bool interactive_flag = false;
  // Flag for global exploration
  bool global_flag = false;

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
        Cerr("***** castor-datafileExplorer() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose))
      {
        Cerr("***** castor-datafileExplorer() -> Exception when trying to read provided verbosity level '" << verbose << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // DataFile
    else if (option=="-df") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-datafileExplorer() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_to_data_filename = (string)argv[i+1];
      i++;
    }
    // Event by event exploration
    else if (option=="-e")
    {
      event_by_event_flag = true;
    }
    // Interactive mode for event-by-event exploration
    else if (option=="-i")
    {
      interactive_flag = true;
    }
    // Global exploration
    else if (option=="-g")
    {
      global_flag = true;
    }
    // Unknown option!
    else
    {
      Cerr("***** castor-datafileExplorer() -> Unknown option '" << option << "' !" << endl);
      Exit(EXIT_FAILURE);
    }
  }

  // ============================================================================================================
  // Some checks
  // ============================================================================================================

  // Check for a datafile name
  if (path_to_data_filename=="")
  {
    Cerr("***** castor-datafileExplorer() -> Please provide a datafile name !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Verbose
  if (verbose<=0)
  {
    Cerr("***** castor-datafileExplorer() -> Verbose less than 1 has no sense for a program used to solely print out information on screen !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Restrict verbose level to 2 to avoid having to much useless information during initialization
  if (verbose>2) verbose = 2;
  // Cannot have event-by-event and global both enables
  if (global_flag && event_by_event_flag)
  {
    Cerr("***** castor-datafileExplorer() -> Cannot use global and event-by-event mode at the same time !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ============================================================================================================
  // Initializations
  // ============================================================================================================

  // Verbose (we know it is at least 1 so we don't check it)
  Cout("==============================================================" << endl);
  Cout("castor-datafileExplorer() -> Initialization starts" << endl);

  // ----------------------------------------------------------------------------------------
  // Create sScannerManager
  // ----------------------------------------------------------------------------------------

  // Get user endianness (interfile I/O)
  GetUserEndianness();

  // Create the scanner manager
  sScannerManager* p_ScannerManager = sScannerManager::GetInstance();  
  p_ScannerManager->SetVerbose(verbose);

  // Get system name from the dataFile
  string scanner_name = "";
  if (ReadDataASCIIFile(path_to_data_filename, "Scanner name", &scanner_name, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-datafileExplorer() -> A problem occured while trying to find the system name in the datafile header !" << endl);
    Exit(EXIT_FAILURE);
  } 
  if (p_ScannerManager->FindScannerSystem(scanner_name) )
  {
    Cerr("***** castor-datafileExplorer() -> A problem occurred while searching for scanner system !" << endl);
    Exit(EXIT_FAILURE);
  } 
  if (p_ScannerManager->BuildScannerObject() )
  {
    Cerr("***** castor-datafileExplorer() -> A problem occurred during scanner object construction ! !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->InstantiateScanner() )
  {
    Cerr("***** castor-datafileExplorer() -> A problem occurred while creating Scanner object !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->GetGeometricInfoFromDataFile(path_to_data_filename))
  {
    Cerr("***** castor-datafileExplorer() -> A problem occurred while retrieving scanner fields from the datafile header !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->BuildLUT() )
  {
    Cerr("***** castor-datafileExplorer() -> A problem occurred while generating/reading the LUT !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->CheckParameters())
  {
    Cerr("***** castor-datafileExplorer() -> A problem occured while checking scanner manager parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->Initialize())
  {
    Cerr("***** castor-datafileExplorer() -> A problem occured while initializing scanner !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ----------------------------------------------------------------------------------------
  // Create vDataFile
  // ----------------------------------------------------------------------------------------

  // Create a default image dimensions and quantification object
  oImageDimensionsAndQuantification* p_ID = new oImageDimensionsAndQuantification();
  p_ID->SetDefault();
  p_ID->SetVerbose(verbose);
  // Virtual datafile
  vDataFile* p_DataFile = NULL;
  // Create specific data file
  if (p_ScannerManager->GetScannerType() == SCANNER_PET) p_DataFile = new iDataFilePET();
  else if (p_ScannerManager->GetScannerType() == SCANNER_SPECT_CONVERGENT) p_DataFile = new iDataFileSPECT(); 
  else if (p_ScannerManager->GetScannerType() == SCANNER_CT) p_DataFile = new iDataFileCT(); 
  // Unknown scanner
  else
  {
    Cerr("***** castor-datafileExplorer() -> Unknown scanner type (" << p_ScannerManager->GetScannerType() << ") for datafile construction ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  // Set data file name
  p_DataFile->SetHeaderDataFileName(path_to_data_filename);
  // Few parameters
  p_DataFile->SetImageDimensionsAndQuantification(p_ID);
  p_DataFile->SetBedIndex(0);
  p_DataFile->SetVerbose(verbose);
  // Read information from the header, without affecting the quantification (as no ImageDimensionsAndQuantification object is used here)
  bool do_not_affect_quantification = false;
  if (p_DataFile->ReadInfoInHeader(do_not_affect_quantification))
  {
    Cerr("***** castor-datafileExplorer() -> A problem occurred during datafile header reading ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_DataFile->CheckParameters())
  {
    Cerr("***** castor-datafileExplorer() -> A problem occurred while checking datafile parameters ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_DataFile->ComputeSizeEvent())
  {
    Cerr("***** castor-datafileExplorer() -> A problem occurred in datafile initialization ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_DataFile->InitializeMappedFile())
  {
    Cerr("***** castor-datafileExplorer() -> A problem occurred in datafile initialization ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_DataFile->PrepareDataFile())
  {
    Cerr("***** castor-datafileExplorer() -> A problem occured in datafile preparation ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }

  // Verbose
  Cout("castor-datafileExplorer() -> End of initialization" << endl);
  Cout("==============================================================" << endl);

  // ============================================================================================================
  // Actions
  // ============================================================================================================

  // -------------------------------------------------------------------
  // Ask the datafile to describe itself
  // -------------------------------------------------------------------
  p_DataFile->Describe();
  Cout("==============================================================" << endl);

  // -------------------------------------------------------------------
  // Then if the event-by-event option is on, describe each event
  // -------------------------------------------------------------------

  if (event_by_event_flag)
  {
    Cout("castor-datafileExplorer() -> Start exploration of all events" << endl);
    // Get index start and stop
    int64_t index_start = 0; int64_t index_stop = 0;
    p_DataFile->GetEventIndexStartAndStop(&index_start, &index_stop);
    // Launch the loop on all events
    int64_t index = index_start;
    while (index>=index_start && index<index_stop)
    {
      // Verbose
      Cout("-------------------------  Event index " << index << "  -------------------------" << endl);
      // Get the current event
      vEvent* p_event = p_DataFile->GetEvent(index);
      if (p_event==NULL)
      {
        Cerr("***** castor-datafileExplorer() -> An error occured while getting the event from index " << index << " !" << endl);
        Exit(EXIT_FAILURE);
      }
      // Describe the event
      p_event->Describe();
      // Interactive mode
      if (interactive_flag)
      {
        // The user may provide a specific event index or simply press enter to get to the next
        cout << "--------> Give an event index or simply press enter for next event: " << flush;
        // Read the answer
        string answer = ""; getline(cin,answer);
        cout << endl;
        // If empty answer, then go to the next event
        if (answer=="") index++;
        // Otherwise
        else
        {
          // Convert the answer to int64_t
          int64_t next_index = stoll(answer);
          if (next_index<index_start || next_index>=index_stop)
          {
            Cerr("***** castor-datafileExplorer() -> The provided event index (" << next_index << ") is out of datafile range"
              << " [" << index_start << ":" << index_stop << "[ !" << endl);
            break;
          }
          else index = next_index;
        }
      }
      else index++;
    }
    Cout("==============================================================" << endl);
  }

  // -------------------------------------------------------------------
  // Or perform a global exploration
  // -------------------------------------------------------------------

  else if (global_flag)
  {
    // Verbose
    Cout("castor-datafileExplorer() -> Start global exploration of all events" << endl);
    // Progression (increments of 2%)
    cout << "0   10%  20%  30%  40%  50%  60%  70%  80%  90%  100%" << endl;
    cout << "|----|----|----|----|----|----|----|----|----|----|" << endl;
    cout << "|" << flush;
    int progression_percentage_old = 0;
    int progression_nb_bars = 0;
    uint64_t progression_printing_index = 0;
    // Get index start and stop
    int64_t index_start = 0; int64_t index_stop = 0;
    p_DataFile->GetEventIndexStartAndStop(&index_start, &index_stop);
    // Cumulative values
    HPFLTNB total_nb_data = 0.;
    // Launch the loop on all events
    for (int64_t index=index_start ; index<index_stop ; index++)
    {
      // Print progression
      if (progression_printing_index%1000==0)
      {
        int progression_percentage_new = ((int)( (((float)(index-index_start+1))/((float)(index_stop-index_start)) ) * 100.));
        if (progression_percentage_new>=progression_percentage_old+2) // Increments of 2%
        {
          int nb_steps = (progression_percentage_new-progression_percentage_old)/2;
          for (int i=0; i<nb_steps; i++)
          {
            cout << "-" << flush;
            progression_nb_bars++;
          }
          progression_percentage_old += nb_steps*2;
        }
      }
      progression_printing_index++;
      // Get the current event
      vEvent* p_event = p_DataFile->GetEvent(index);
      if (p_event==NULL)
      {
        Cerr("***** castor-datafileExplorer() -> An error occured while getting the event from index " << index << " !" << endl);
        Exit(EXIT_FAILURE);
      }
      // Number of data
      for (INTNB b=0; b<p_event->GetNbValueBins(); b++) total_nb_data += ((HPFLTNB)(p_event->GetEventValue(b)));
    }
    // End of progression printing (do not log out with Cout here)
    int progression_total_bars = 49;
    for (int i=0; i<progression_total_bars-progression_nb_bars; i++) cout << "-";
    cout << "|" << endl;
    Cout("  --> Total number of data: " << total_nb_data << endl);
    Cout("==============================================================" << endl);
  }

  // ============================================================================================================
  // End
  // ============================================================================================================

  // Delete objects in the inverse order in which they were created
  delete p_DataFile;
    
  // Ending
  #ifdef CASTOR_MPI
  MPI_Finalize();
  #endif
  return EXIT_SUCCESS;
}
