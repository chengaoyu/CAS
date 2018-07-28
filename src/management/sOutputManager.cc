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
  \ingroup management

  \brief Implementation of class sOutputManager
*/

#include "sOutputManager.hh"
#include "oImageSpace.hh"
#include "oImageDimensionsAndQuantification.hh"
#include <iomanip> 

#ifdef CASToR_USE_CMAKE
  #include "oCASToRConfig.hh"
#endif

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

// Implementation of the exit function that must be used instead of the standard exit
void Exit(int code)
{
  if (code!=0) Cerr("***** Exit function called. Abort with code " << code << "." << endl);
  #ifdef CASTOR_MPI
  MPI_Finalize();
  #endif
  exit(code);
}

// Singleton : set pointer to object to NULL
sOutputManager *sOutputManager::mp_Instance = NULL;

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \brief sOutputManager constructor.
  \details It is private at this class is singleton. 
           It should be instanciated using the GetInstance() function 
           Initialize the member variables to their default values.
*/
sOutputManager::sOutputManager()
{
  // Initialize all members as default
  m_verbose = 0;
  m_mpiRank = 0;
  m_baseName = "";
  m_pathName = "";
  m_pathToConfigDir = "";
  m_mergeOutputDynImgFlag = false;
  mp_dataFileName = {""};
  // Set scientific numbers
  cout << std::scientific;
  cerr << std::setprecision(std::numeric_limits<FLTNB>::digits10+1);
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \brief sOutputManager destructor. 
*/
sOutputManager::~sOutputManager()
{
  // Just have to close the log file
  if (m_mpiRank==0 && m_logFile) m_logFile.close();
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn CheckConfigDir
  \param a_path
  \brief Set the path to the CASTOR config directory 
  \details Set the path to the CASTOR config directory from the given path if it is provided
           Otherwise, try to get it from the environment variable CASTOR_CONFIG
  \return 0 if success, positive value otherwise
*/
int sOutputManager::CheckConfigDir(const string& a_path)
{
  // Case 1: a path is provided
  if (a_path!="")
  {
    if (m_verbose>=3) Cout("sOutputManager::CheckConfigDir() -> Directory selected from option as '" << a_path << "'" << endl); 
    m_pathToConfigDir = a_path + OS_SEP;
  }
  // Case 2: no path provided so we look after the CASTOR_CONFIG environment variable
  else
  {
    #ifdef CASToR_USE_CMAKE
    string tmp_path = CASTOR_CONFIG;
    if (tmp_path.empty()) // throw error if empty
    #elif defined(CASTOR_USE_MINGW)
    #ifdef CASTOR_CONFIG
    // This macro CASTOR_CONFIG_STRING is declared in sOutputManager.hh and automatically convert the value from the environment
    // CASTOR_CONFIG into a string (i.e. including the double quotes so that it can be used here in the affectation)
    string tmp_path = ConvertAllSlashOccurrencesToBackSlash(CASTOR_CONFIG_STRING);
    if (tmp_path.empty())
    #else
    // Here the CASTOR_CONFIG variable must have been defined before cross-compilation, so we write a message that will make the
    // compiler crash and display this fake line! Please let it as is.
    When cross-compiling, you must define the CASTOR_CONFIG environment variable before that. This compilation error you are seeing is normal!
    #endif
    #else
    char* tmp_path = getenv("CASTOR_CONFIG"); 
    if (tmp_path==NULL) // throw error if empty
    #endif
    {
      Cerr("***** sOutputManager::CheckConfigDir() -> No path nor CASTOR_CONFIG variable provided !" << endl);
      return 1;
    }
    if (m_verbose>=3) Cout("sOutputManager::CheckConfigDir() -> Directory selected from environment variable as '" << ((string)tmp_path) << "'" << endl);
    m_pathToConfigDir = ((string)tmp_path) + OS_SEP;
  }
  // End
  return 0;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetPathToConfigDir
  \brief Return the path to the CASTOR config directory 
  \details Just return the path if it has already been initialized
           Otherwise, the function recovers the path from environnement variables
           If any error, the working directory is returned instead
  \return a string containing path to the CASToR configuration directory
*/
const string& sOutputManager::GetPathToConfigDir()
{
  // If the config directory has already been initialized, then simply return its current value
  if (m_pathToConfigDir!="") return m_pathToConfigDir;
  // Otherwise, this means that this function is probably called from before the singleton initialization.
  else
  {
    if (m_verbose>=3) Cout("sOutputManager::GetPathToConfigDir ..."<< endl); 
    // We get the directory from the environment variable
    #if defined(CASToR_USE_CMAKE)
    string tmp_path = CASTOR_CONFIG;
    if (tmp_path.empty())
    #elif defined(CASTOR_USE_MINGW)
    #ifdef CASTOR_CONFIG
    // This macro CASTOR_CONFIG_STRING is declared in sOutputManager.hh and automatically convert the value from the environment
    // CASTOR_CONFIG into a string (i.e. including the double quotes so that it can be used here in the affectation)
    string tmp_path = ConvertAllSlashOccurrencesToBackSlash(CASTOR_CONFIG_STRING);
    if (tmp_path.empty())
    #else
    // Here the CASTOR_CONFIG variable must have been defined before cross-compilation, so we write a message that will make the
    // compiler crash and display this fake line! Please let it as is.
    When cross-compiling, you must define the CASTOR_CONFIG environment variable before that. This compilation error you are seeing is normal!
    #endif
    #else
    char* tmp_path = getenv("CASTOR_CONFIG");
    if (tmp_path==NULL)
    #endif
    {
      Cerr("***** sOutputManager::CheckConfigDir() -> No path nor CASTOR_CONFIG variable provided ! Try working directory instead." << endl);
      m_pathToConfigDir = ".";
    }
    else
      #if defined(CASToR_USE_CMAKE) || defined(CASTOR_USE_MINGW)
      m_pathToConfigDir = tmp_path + OS_SEP;
      #else
      m_pathToConfigDir = ((string)tmp_path) + OS_SEP;
      #endif
  
    return m_pathToConfigDir;
  }
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InitOutputDirectory
  \param a_pathFout : path to an output file as provided by the user
  \param a_pathDout : path to an output directory as provided by the user
  \brief Create the output directory if any, extract the base name and create the log file.
  \return 0 if success, and positive value otherwise.
*/
int sOutputManager::InitOutputDirectory(const string& a_pathFout, const string& a_pathDout)
{
  #ifdef CASTOR_VERBOSE
  if (m_verbose>=4) Cout("+++++ sOutputManager::InitOutputDirectory() -> Enter"<< endl); 
  #endif

  // Check unicity of the path
  if (a_pathFout!="" && a_pathDout!="")
  {
    Cerr("***** sOutputManager::InitOutputDirectory() -> Either a file path (-fout) or a directory path (dout) should be provided.  cannot be both provided, make your choice !" << endl);
    return 1;
  }

  // Check if the provided path ends with '.' or '/' or '\', then alert and crash
  string the_path = a_pathFout;
  if (a_pathDout!="") the_path = a_pathDout;
  string last_char = the_path.substr(the_path.length()-1);
  if (last_char=="." || last_char==OS_SEP)
  {
    Cerr("***** sOutputManager::InitOutputDirectory() -> Please provide a path not finishing by '.', '"<< OS_SEP <<"' character !" << endl);
    return 1;
  }

  // Verbose
  if (m_verbose>=1) Cout("sOutputManager::InitOutputDirectory() -> Output path is '" << a_pathFout << a_pathDout << "'" << endl);

  // -------------------------------------------------
  // First case: a file path is provided
  // -------------------------------------------------

  if (a_pathFout!="")
  {
    // Get the last slash position
    size_t last_slash_pos = a_pathFout.find_last_of(OS_SEP);
    // No slash
    if (last_slash_pos==string::npos) 
    {
      m_pathName = "";
      m_baseName = a_pathFout;
    }
    // Some slashes
    else
    {
      // Everything before the slash becomes the path
      m_pathName = a_pathFout.substr(0,last_slash_pos+1);
      // Everything after the slash becomes the base name
      m_baseName = a_pathFout.substr(last_slash_pos+1);
    }
  }

  // -------------------------------------------------
  // Second case: a directory path is provided
  // -------------------------------------------------

  else if (a_pathDout!="")
  {
    // Get the last slash position
    size_t last_slash_pos = a_pathDout.find_last_of(OS_SEP);
    // No slash
    if (last_slash_pos==string::npos) 
    {
      m_pathName = a_pathDout + OS_SEP;
      m_baseName = a_pathDout;
    }
    // Some slashes
    else
    {
      // The whole becomes the path
      m_pathName = a_pathDout + OS_SEP;
      // Everything after the slash becomes the base name
      m_baseName = a_pathDout.substr(last_slash_pos+1);
    }
    // Create directory only for first MPI instance
    if (m_mpiRank==0)
    {
      #ifdef _WIN32
      string instruction = "if not exist " + m_pathName + " mkdir " + m_pathName;
      #else
      string instruction = "mkdir -p " + m_pathName;
      #endif
      int error = system(instruction.c_str());
      if (error)
      {
        Cerr("***** sOutputManager::InitOutputDirectory() -> Failed to create output directory with name '" << m_pathName << "' !" << endl);
        return 1;
      }
    }
    // Verbose
    if (m_verbose>=3) Cout("  --> Output files will be written inside directory '" << m_pathName << "'" << endl);
  }

  // -------------------------------------------------
  // Third: create the log file
  // -------------------------------------------------

  // Only first MPI instance deals with this
  if (m_mpiRank==0)
  {
    // Create file name
    string log_file_name = m_pathName + m_baseName + ".log";
    // Open and check
    m_logFile.open(log_file_name.c_str());
    if (!m_logFile)
    {
      Cerr("***** sOutputManager::InitOutputDirectory() -> Failed to create output log file as '" << log_file_name << "' ! Are you sure the provided output path exists ?" << endl);
      return 1;
    }
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int sOutputManager::LogCommandLine(int argc, char** argv)
{
  if (m_verbose>=3) Cout("sOutputManager::LogCommandLine() -> Write command line into the log file"<< endl); 
  
  // Exit function is MPI rank is anything other than 0
  if (m_mpiRank!=0) return 0;

  if (m_logFile)
  {
    m_logFile << "==================================================================================================" << endl;
    m_logFile << "                                      COMMAND LINE CONTEXT" << endl;
    m_logFile << "==================================================================================================" << endl;
    // Print command line
    m_logFile << "Command line: ";
    for (int i=0; i<argc; i++) m_logFile << argv[i] << " ";
    m_logFile << endl;
    #ifdef _WIN32
    char pwd[MAX_PATH];
    GetCurrentDirectory(MAX_PATH,pwd);
    m_logFile << "Working directory: " << pwd << endl;
    #else
    m_logFile << "Working directory: " << getenv("PWD") << endl;
    #endif
    std::time_t date_of_execution = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    m_logFile << "Data of execution: " << std::ctime(&date_of_execution);
    m_logFile << "Float numbers precision in bytes (for matrices and some computations): " << sizeof(FLTNB) << endl;
    m_logFile << "Float numbers precision in bytes (for sensitive computations that require at least double precision): " << sizeof(HPFLTNB) << endl;
    m_logFile << "Float numbers precision in bytes (for datafile reading/writing): " << sizeof(FLTNBDATA) << endl;
    m_logFile << "Float numbers precision in bytes (for scanner LUT reading/writing): " << sizeof(FLTNBLUT) << endl;
    m_logFile << "CASToR version: " << CASTOR_VERSION << endl;
    m_logFile << "==================================================================================================" << endl << flush;
  }
  else 
    return 1;

  // End
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      inline int  sOutputManager::SetOutNbPrec()
  \param   a_format : string containing output format and precision
  \brief   Set the output format and precision used for numeric display
  \return  0 if success, positive value otherwise 
*/
int sOutputManager::SetOutNbPrec(string a_format)
{
  string format;
  int32_t precision;
  string option = "-onbc";
  
  // parsing   
  size_t pos = a_format.find_first_of(",");      
  if (pos == string::npos)
  {
    Cerr("***** sOutputManager::SetOutNbPrec() -> Error, format must be (format,precision) (e.g = (f,5) or (s,0)!" << endl);
    return 1;
  }
  
  format = a_format.substr(0,pos);
  if (format == "f")
    cout << std::fixed;

  if (ReadStringOption(a_format.substr(pos+1), &precision, 1, ",", option))
  {
    Cerr("***** sOutputManager::SetOutNbPrec() -> Invalid argument " << a_format << " for option " << option << " !" << endl);
    Cerr("***** sOutputManager::SetOutNbPrec() -> Format must be (format,precision) (e.g = (f,5) or (s,0)!" << endl);
    return 1;
  }

  if(precision>0)
    cout << std::setprecision(precision);
  else 
    cout << std::setprecision(std::numeric_limits<FLTNB>::digits10+1);
      
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
