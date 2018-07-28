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
  \brief Declaration of class sOutputManager
*/

#ifndef SOUTPUTMANAGER_HH
#define SOUTPUTMANAGER_HH 1

#include "gVariables.hh"
class oImageDimensionsAndQuantification;

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// Cross-platform directory separator, if required
// Windows might auto-convert '/' to '\' though
#ifdef _WIN32
#define OS_SEP "\\"
#else
#define OS_SEP "/"
#endif

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// Use to convert a macro variable provided by the pre-processor, to a string
#define QUOTES(value) #value
#define TOSTRING(macro) QUOTES(macro)
// The CASTOR_CONFIG environment variable converted to a string
#ifdef CASTOR_CONFIG
#define CASTOR_CONFIG_STRING TOSTRING(CASTOR_CONFIG)
#endif

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

#ifdef CASTOR_MPI
// Macros for logging messages on standard output (with MPI)
#define Cout(MESSAGE)                                           \
  do                                                            \
  {                                                             \
    int mpi_rank; MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);     \
    if (mpi_rank==0)                                            \
    {                                                           \
      std::cout << MESSAGE;                                     \
      sOutputManager* instance = sOutputManager::GetInstance(); \
      if (instance!=NULL)                                       \
      {                                                         \
        ofstream& logMac = instance->GetLogFile();              \
        if (logMac) logMac << MESSAGE;                          \
      }                                                         \
    }                                                           \
  } while(0)
// Macros for logging messages on standard error (with MPI)
#define Cerr(MESSAGE)                                           \
  do                                                            \
  {                                                             \
    int mpi_rank; MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);     \
    std::cerr << "***** The following error message was issued" \
              << " by MPI instance " << mpi_rank << endl;       \
    std::cerr << MESSAGE;                                       \
    sOutputManager* instance = sOutputManager::GetInstance();   \
    if (instance!=NULL)                                         \
    {                                                           \
      ofstream& logMac = instance->GetLogFile();                \
      if (logMac)                                               \
      {                                                         \
        logMac << "***** The following error message was issued"\
               << " by MPI instance " << mpi_rank << endl;      \
        logMac << MESSAGE;                                      \
      }                                                         \
    }                                                           \
  } while(0)
#else
// Macros for logging messages on standard output (without MPI)
#define Cout(MESSAGE)                                         \
  do                                                          \
  {                                                           \
    std::cout << MESSAGE;                                     \
    sOutputManager* instance = sOutputManager::GetInstance(); \
    if (instance!=NULL)                                       \
    {                                                         \
      ofstream& logMac = instance->GetLogFile();              \
      if (logMac) logMac << MESSAGE;                          \
    }                                                         \
  } while(0)
// Macros for logging messages on standard error (without MPI)
#define Cerr(MESSAGE)                                         \
  do                                                          \
  {                                                           \
    std::cerr << MESSAGE;                                     \
    sOutputManager* instance = sOutputManager::GetInstance(); \
    if (instance!=NULL)                                       \
    {                                                         \
      ofstream& logMac = instance->GetLogFile();              \
      if (logMac) logMac << MESSAGE;                          \
    }                                                         \
  } while(0)
#endif

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// Verbose levels definitions
#define VERBOSE_LIGHT  1
#define VERBOSE_NORMAL 2
#define VERBOSE_DETAIL 3
#define VERBOSE_DEBUG_LIGHT  4
#define VERBOSE_DEBUG_NORMAL 5
#define VERBOSE_DEBUG_EVENT  6
#define VERBOSE_DEBUG_MAX    7

// Macro for CASTOR_VERBOSE to print out the entering of each function.
// This may be used for debugging. When this compilation variable is not
// set, this does nothing.

// If CASTOR_VERBOSE is not set, we define this macro as doing nothing ';'
#ifndef CASTOR_VERBOSE

  #define DEBUG_VERBOSE(IGNORED1,IGNORED2) ;

// If CASTOR_VERBOSE is set, we differentiate between windows/unix systems
#else
  // Case for pure windows
  #if defined(_WIN32) && !defined(CASTOR_USE_MINGW)
    #define DEBUG_VERBOSE(CURRENT_LEVEL,APPLICABLE_LEVEL)  \
      do                                                   \
      {                                                    \
        if (CURRENT_LEVEL >= APPLICABLE_LEVEL)             \
        {                                                  \
          cout << "+++++ " << __FUNCSIG__                  \
               << " -> Entering this function" << endl;    \
        }                                                  \
      } while(0);
  // Case for unix and cross-compilation from unix to windows
  #else
    #define DEBUG_VERBOSE(CURRENT_LEVEL,APPLICABLE_LEVEL)  \
      do                                                   \
      {                                                    \
        if (CURRENT_LEVEL >= APPLICABLE_LEVEL)             \
        {                                                  \
          cout << "+++++ " << __PRETTY_FUNCTION__          \
               << " -> Entering this function" << endl;    \
        }                                                  \
      } while(0);
  #endif
#endif

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// Function used to make an clean exit taking MPI finalize into account
/*!
  \fn      void Exit()
  \param   int code
  \brief   This function must be used to exit the program whenever needed
  \details It does a clean exit taking MPI finalize into account if applicable.
*/
void Exit(int code);

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

/*!
  \class   sOutputManager
  \brief   Singleton class that manages output writing on disk (images, sinograms, etc). \n
           It also manages logging and printing on screen.
*/
class sOutputManager
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      static sOutputManager* sOutputManager::GetInstance()
      \brief   Instanciate the singleton object and Initialize member variables if not already done, 
               return a pointer to this object otherwise
      \return  instance of the sOutputManager singleton
    */
    static sOutputManager* GetInstance() 
    { 
      if (mp_Instance == NULL) mp_Instance = new sOutputManager();
      return mp_Instance;
    }


  // -------------------------------------------------------------------
  // Get & Set functions
  public:
    /*!
      \fn      inline void sOutputManager::SetVerbose()
      \param   a_verbose
      \brief   set verbosity
    */
    inline void SetVerbose( int a_verbose )
           {m_verbose = a_verbose;}
    /*!
      \fn      inline void sOutputManager::SetDataFileName
      \return  set the datafile names
    */
    inline void SetDataFileName( vector<string> ap_dataFileName )
           {mp_dataFileName = ap_dataFileName;}
    /*!
      \fn      inline ofstream& sOutputManager::GetLogFile
      \return  the ofstream object for the log file
    */
    inline ofstream& GetLogFile()
           {return m_logFile;}
    /*!
      \fn      inline const string& sOutputManager::GetBaseName()
      \return  the string containing the output image base name
    */
    inline const string& GetBaseName()
           {return m_baseName;}
    /*!
      \fn      inline const string& sOutputManager::GetPathName()
      \return  the string containing the path to the output directory
    */
    inline const string& GetPathName()
           {return m_pathName;}
    /*!
      \fn      inline vector<string> sOutputManager::GetDataFileName
      \return  recover datafile names
    */
    inline vector<string> GetDataFileName()
           {return mp_dataFileName;}
    /*!
      \fn      inline void sOutputManager::SetMPIRank()
      \param   a_mpiRank
      \brief   Initialize the machine index for MPI
    */
    inline void SetMPIRank( int a_mpiRank )
           {m_mpiRank = a_mpiRank;}
    /*!
      \fn      inline void sOutputManager::SetMergeDynImagesFlag()
      \param   a_flag
      \brief   Set to on the flag indicating that a dynamic serie of 3D images should be
               written on disk in one file instead of one file for each 3D image 
               associated with a metaheader, or not
    */
    inline void SetMergeDynImagesFlag( bool a_flag )
           {m_mergeOutputDynImgFlag = a_flag;}
    /*!
      \fn      inline bool sOutputManager::MergeDynImages()
      \brief   Indicate if a dynamic serie of 3D images should be merged in one file (true)
               or written on disk as one file for each 3D image, associated with an Interfile metaheader
      \return  true if images should be merged in one file, false otherwise.
    */
    inline bool MergeDynImages()
           {return m_mergeOutputDynImgFlag;}


    /*!
      \fn      inline int  sOutputManager::SetOutNbPrec()
      \param   a_format : string containing output format and precision
      \brief   Set the output format and precision used for numeric display
      \return  0 if success, positive value otherwise 
    */
    int SetOutNbPrec(string a_format);


  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      int sOutputManager::CheckConfigDir()
      \param   a_path
      \brief   Set the path to the CASTOR config directory from the given path if not empty or through
               the existence of the environment variable CASTOR_CONFIG
      \return  0 if success, positive value otherwise
    */
    int CheckConfigDir( const string& a_path );
    /*!
      \fn      const string& sOutputManager::GetPathToConfigDir()
      \brief   Return the path to the CASTOR config directory 
      \details Just return the path if it has already been initialized
               Otherwise, the function recovers the path from environnement variables
               If any error, the working directory is returned instead
      \return  astring containing path to the CASToR configuration directory
    */
    const string& GetPathToConfigDir();
    /*!
      \fn      int sOutputManager::InitOutputDirectory()
      \param   a_pathFout : path to an output file as provided by the user
      \param   a_pathDout : path to an output directory as provided by the user
      \brief   Create the output directory if any, extract the base name and create the log file.
      \return  0 if success, and positive value otherwise.
    */
    int InitOutputDirectory( const string& a_pathFout, const string& a_pathDout );
    /*!
      \fn      int sOutputManager::LogCommandLine()
      \param   argc
      \param   argv
      \brief   Write log file header with the provided command line options and different informations
      \return  0 if success, positive value otherwise.
    */
    int LogCommandLine( int argc, char** argv );
    
    
  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \brief   sOutputManager constructor.
      \details It is private at this class is singleton. 
               It should be instanciated using the GetInstance() function 
               Initialize the member variables to their default values.
    */
    sOutputManager();
    /*!
      \brief   sOutputManager destructor. 
    */
    ~sOutputManager();
    // Prevent the compiler to generate methods to copy the object
    sOutputManager(sOutputManager const&){};     
    void operator=(sOutputManager const&){};


  // -------------------------------------------------------------------
  // Data members
  private:
    static sOutputManager* mp_Instance; /*!< Pointer to this singleton object */
    int m_verbose;                      /*!< Verbosity */
    int m_mpiRank;                      /*!< Machine index for MPI */
    string m_baseName;                  /*!< String containing the image base name */
    string m_pathName;                  /*!< String containing the path to the output directory */
    string m_pathToConfigDir;           /*!< String containing the path to the CASToR configuration directory */
    vector<string> mp_dataFileName;     /*!< String containing the datafile name(s) */
    ofstream m_logFile;                 /*!< File object for the output log file */
    bool m_mergeOutputDynImgFlag;       /*!< Flag indicating if dynamic image should be written on disk as one file (true),
                                             or a serie of 3D image associated with a metaheader (false). Default: false */
};

#endif
