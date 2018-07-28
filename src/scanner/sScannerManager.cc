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
  \ingroup  scanner
  \brief    Implementation of class sScannerManager
*/

#include "sScannerManager.hh"
#include "sAddonManager.hh"

// Singleton : set pointer to object to NULL
sScannerManager *sScannerManager::mp_Instance = NULL;

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \brief sScannerManager constructor.
  \details It is private at this class is singleton. 
           This class should be instanciated using the GetInstance() function 
           Initialize the member variables to their default values.
*/
sScannerManager::sScannerManager()
{
  mp_Scanner = NULL;
  mp_ID = NULL;
  m_verbose = -1;
  m_pathToScannerFile = "";
  m_scannerName = "";
  m_hasUserScannerFile = false;
  m_hasGenericScannerFile = false;
  m_allParametersChecked = false;
  m_saveLUTFlag = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \brief sScannerManager destructor. 
*/
sScannerManager::~sScannerManager()
{
  if (mp_Scanner != NULL) delete mp_Scanner;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn CheckParameters
  \brief Check if all parameters have been correctly initialized, and call the CheckParameters function of the scanner object
  \return 0 if success. Positive value otherwise
*/
int sScannerManager::CheckParameters()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // Check mandatory parameters
  if (mp_Scanner == NULL)
  {
    Cerr("***** sScannerManager::CheckParameters() -> Scanner object not initialized !" << endl);
    return 1;
  }
  if (!m_hasUserScannerFile && !m_hasGenericScannerFile)
  {
    Cerr("***** sScannerManager::CheckParameters() -> Scanner file type (generic or user provided) not initialized !" << endl);
    return 1;
  }
  if (m_scannerName.empty())
  {
    Cerr("***** sScannerManager::CheckParameters() -> Scanner name not initialized" << endl);
    return 1;
  }
  if (m_pathToScannerFile.empty())
  {
    Cerr("***** sScannerManager::CheckParameters() -> Path to scanner file not initialized !" << endl);
    return 1;
  }
  if (m_verbose<0)
  {
    Cerr("***** sScannerManager::CheckParameters() -> Verbosity level not initialized !" << endl);
    return 1;
  }
  // Check parameters of the scanner object
  if (mp_Scanner->CheckParameters())
  {
    Cerr("***** sScannerManager::CheckParameters() -> A problem occurred while checking Scanner Object parameters !" << endl);
    return 1;
  }
  // All parameters are checked
  m_allParametersChecked = true;
  // End
  return 0;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      sScannerManager::Describe()
  \brief   Call the eponym function from the Scanner object (if initialized)
*/
void sScannerManager::Describe()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  
  if (mp_Scanner == NULL)
  {
    Cerr("***** WARNING: sScannerManager::Describe() -> Trying to describe the scanner object while it is not initialized !" << endl);
  }
  else
  {
    Cout(endl << "---------------------------------------------------------------------------" << endl);
    Cout("sScannerManager::Describe() -> Here is some generic content of the scanner:" << endl);
    Cout("  --> Scanner name: " << m_scannerName << endl);
    Cout("  --> Scanner file: " << m_pathToScannerFile << endl);
    mp_Scanner->Describe();
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn Initialize
  \brief Initialization :
         - check if all parameters of the manager have been checked
         - call the initialization function of the scanner object 
  \return 0 if success. Positive value otherwise
*/
int sScannerManager::Initialize()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // Parameters checked ?
  if (!m_allParametersChecked)
  {
     Cerr("***** sScannerManager::Initialize() -> Parameters have not been checked !" << endl);
     return 1;
  }
  // Verbose
  if (m_verbose>=VERBOSE_LIGHT) Cout("sScannerManager::Initialize() -> From scanner " << m_scannerName << endl);
  // Initialize the scanner object
  if (mp_Scanner->Initialize())
  {
    Cerr("***** sScannerManager::Initialize() -> A problem occurred while initializing Scanner Object !" << endl);
    return 1;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ShowScannersDescription
  \brief Get the description associated to the different scanners and print all on screen.
         Walk through the scanner repository and look for the keyword "description" in .geom file and .hscan file.
  \return 0 if success, positive value otherwise
  \todo Check everything output correctly for all implemented scanners
*/
int sScannerManager::ShowScannersDescription()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return 0;
  #endif

  // Verbose
  if (m_verbose>=VERBOSE_DETAIL) Cout("sScannerManager::ShowScannersDescription ..."<< endl);

  // Gather all the available scanner from the repository
  vector<string> list_scanner_names;
  if (GetAvailableScanners(&list_scanner_names))
  {
    Cerr("***** sScannerManager::ShowScannersDescription() -> A problem occurred while recovering scanner names from the scanner repository !" << endl);
    return 1;
  }

  // Print-out descriptions
  cout << endl << "Here is the list of all available scanners in the repository along with their options:" << endl << endl;
  for (unsigned int iter = 0; iter!=list_scanner_names.size(); iter++)
  {
    // Get the actual scanner in temporary local variable
    string scanner_file = sOutputManager::GetInstance()->GetPathToConfigDir() + "scanner" + OS_SEP;
    // Print out the name of this scanner
    cout << "------------------------------------------------------------------" << endl;
    cout << "-----  \"" << list_scanner_names[iter] << "\"" << endl;
    cout << "------------------------------------------------------------------" << endl;
    // Proceed
    scanner_file.append(list_scanner_names[iter]);
    // Look for the "description" keyword in the .geom or .hscan file
    string description;
    ReadDataASCIIFile(scanner_file, "description", &description, 1, 1);
    // Print out description
    cout << description << endl << endl;
  }
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn FindScannerSystem
  \param a_scannerName : string containing name of the required scanner
  \brief Look for a file matching with the scanner name in parameter inside the scanner repository 
  \return 0 if success (scanner found). Positive value otherwise
*/
int sScannerManager::FindScannerSystem(string a_scannerName)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Check emptiness
  if (a_scannerName.empty())
  {
    Cerr("***** sScannerManager::FindScannerSystem() -> scanner name is empty !" << endl);
    return 1;
  }

  // Verbose
  if (m_verbose>=VERBOSE_DETAIL) Cout("sScannerManager::FindScannerSystem() -> Search for scanner " << a_scannerName << " in the configuration directory" << endl);

  // Get the list of scanner from the repository
  vector<string> repository_scanner_names;
  if (GetAvailableScanners(&repository_scanner_names))
  {
    Cerr("***** sScannerManager::FindScannerSystem() -> A problem occurred while recovering scanner names from the scanner repository !" << endl);
    return 1;
  }

  // Set the scanner name
  m_scannerName = a_scannerName;
  m_pathToScannerFile = sOutputManager::GetInstance()->GetPathToConfigDir() + "scanner" + OS_SEP;

  // String to recover generic/user-made scanner file
  string generic_scanner_file, user_scanner_file;

  // Loop over the scanner available in the repository, and check if the scanner name match with any of these
  for (unsigned int index_name=0 ; index_name<repository_scanner_names.size() ; index_name++)
  {
    string gfile = a_scannerName;
    string ufile = a_scannerName;
    // Check if the scanner match a geom file (generic_scanner_file)
    if (gfile.append(".geom") == repository_scanner_names[index_name])
    {
      generic_scanner_file = a_scannerName;
      generic_scanner_file.append(".geom");
      m_hasGenericScannerFile = true;
      if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Matched generic file for the scanner system: " << repository_scanner_names[index_name] << endl);
    }
    // Check if the scanner match an user file (user provided LUT)
    if (ufile.append(".hscan") == repository_scanner_names[index_name])
    {
      user_scanner_file = a_scannerName;
      user_scanner_file.append(".hscan");
      m_hasUserScannerFile = true;
      if (m_verbose>=VERBOSE_DETAIL) Cout("  --> matched custom LUT for the scanner system: " << repository_scanner_names[index_name] << endl);
    }
  }

  // When both a generic file and user-provided LUT exist for this system, the generic file is selected by default
  if (m_hasGenericScannerFile && m_hasUserScannerFile)
  {
    if (m_verbose>=VERBOSE_DETAIL)
    {
      Cerr("***** WARNING sScannerManager::FindScannerSystem() -> Both a generic file and user-provided LUT have been detected for the " << a_scannerName << " system." << endl);
      Cerr("                                                      The generic file (*.geom) will be selected by default." << endl);
    }
    m_hasUserScannerFile = false;
  }

  // Initialize m_pathToScannerFile member variable
  if (m_hasGenericScannerFile) m_pathToScannerFile.append(generic_scanner_file.c_str());
  else if (m_hasUserScannerFile) m_pathToScannerFile.append(user_scanner_file.c_str());
  // Unknown scanner, output scanner description and throw error
  else
  {
    Cerr("***** sScannerManager::FindScannerSystem() -> Scanner '"<< a_scannerName << 
         "' is not known in the scanner repository. Please provide a LUT/generic file for this scanner in the scanner repository in: " << 
         sOutputManager::GetInstance()->GetPathToConfigDir() << "scanner" << OS_SEP << endl << endl;);
    ShowScannersDescription();
    return 1;
  }
  // End  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      int sScannerManager::InitScannerWithFile()
  \param   a_pathScanFile : string containing the path to the scanner file
  \param   a_scannerName  : string containing the name of the required scanner
  \param   a_fileTypeFlag : string containing the type of the scanner file (0=lut, 1=geom)
  \brief   Initialize member variables (file path, file type, and scanner name) with the provided arguments
  \return  0 if success (scanner found). Positive value otherwise
*/
int sScannerManager::InitScannerWithFile(string a_pathScanFile, string a_scannerName, int a_fileTypeFlag)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("sScannerManager::InitScannerWithFile()" << endl);

  // Set the scanner name
  m_scannerName = a_scannerName;
  m_pathToScannerFile = a_pathScanFile;
  
  if( a_fileTypeFlag == 0 )
    m_hasUserScannerFile = true;
  else if ( a_fileTypeFlag == 1 )
    m_hasGenericScannerFile = true;
  else
  {
    Cerr("***** sScannerManager::InitScannerWithFile() -> Error : File type provided in argument must be == 0 (lut) or == 1 (geom). Current value is  '" << a_fileTypeFlag << " !" << endl);
    return 1;
  }
  
  // End  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn BuildScannerObject()
  \brief Instantiate the specific scanner object related to the modality, and set verbosity of scanner object
  \todo delete the check on modality once all scanner classes will be completely implemented ? 
  \return 0 if success. Positive value otherwise
*/
int sScannerManager::BuildScannerObject()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("sScannerManager::BuildScannerObject() -> Start building"<< endl); 
    
  // Check scanner modality
  string scanner_type;
  typedef vScanner *(*maker_scanner) ();

  // Get the system type from the value of the modality field inside the scanner header
  if (ReadDataASCIIFile(m_pathToScannerFile, "modality", &scanner_type, 1, KEYWORD_MANDATORY) == 1)
  {
    Cerr("***** sScannerManager::BuildScannerObject() -> 'Modality' field not found in the header of the scanner configuration file at :"<< endl);
    Cerr("*****                                          " << m_pathToScannerFile << endl);
    return 1;
  }

  // Get scanner's list from addon manager
  std::map <string,maker_scanner> list = sAddonManager::GetInstance()->mp_listOfScannerTypes;

  // Instanciate the scanner class corresponding to the modality, throw error if no match
  if (list[scanner_type]) mp_Scanner = list[scanner_type]();
  else
  {
    Cerr("***** sScannerManager::BuildScannerObject() -> Modality '" << scanner_type << "' is unknown !" << endl);
    sAddonManager::GetInstance()->ShowHelpScanner();
    return 1;
  }

  // Set scanner verbosity
  mp_Scanner->SetVerbose(m_verbose);

  // Set scanner image dimensions and quantification pointer
  mp_Scanner->SetImageDimensionsAndQuantification(mp_ID);

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetGeometricInfoFromDataFile
  \param a_path : string containing the path to datafile header
  \brief Call the specialized function of the scanner object in order to 
         get geometric informations from the datafile header
  \return 0 if success. Positive value otherwise
*/
int sScannerManager::GetGeometricInfoFromDataFile(string a_path)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("sScannerManager::GetGeometricInfoFromDataFile() -> Look for acquisition specific settings into datafile header"<< endl); 

  // Get information
  if (mp_Scanner->GetGeometricInfoFromDataFile(a_path))
  {
    Cerr("***** sScannerManager::GetGeometricInfoFromDataFile() -> An error occurred while getting information from the datafile." << endl);          
    return 1;
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InstantiateScanner
  \brief Instantiate scanner using the related function in the scanner classes
  \todo delete the check on scanner type once all scanner classes will be completely implemented. 
  \return 0 if success. Positive value otherwise
*/
int sScannerManager::InstantiateScanner()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("sScannerManager::InstantiateScanner() -> Instantiate the scanner geometry structure"<< endl); 
  
  // Check if the scanner object is known (TODO : delete this check once all scanner classes will be completely implemented)
  if (mp_Scanner->GetScannerType()<0)
  {
    Cerr("***** sScannerManager::BuildScannerObject() -> Unknow scanner type !" << endl);
    return 1;
  }

  // Instantiate geometry
  if (mp_Scanner->Instantiate(m_hasUserScannerFile))
  {
    Cerr("***** sScannerManager::InstantiateScanner() -> An error occurred while instanciating the Scanner object !" << endl);
    return 1;
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn BuildLUT()
  \brief Call the eponym function of the scanner class
  \return 0 if success. Positive value otherwise
*/
int sScannerManager::BuildLUT()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("sScannerManager::BuildLUT() -> Generate the geometric Look-Up Table"<< endl); 
    
  // Generate LUT
  if (mp_Scanner->BuildLUT(m_hasUserScannerFile))
  {
    Cerr("***** sScannerManager::BuildLUT() -> A problem occurred while generating/reading the LUT !" << endl);          
    return 1;
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetAvailableScanners
  \param ap_scannerNames : vector list of string to recover the available scanner names
  \brief Gather all the names of the header files (.geom & .hscan) 
         in the repository folder in the vector<string> passed in parameter
  \return 0 if sucess, positive value otherwise
*/
int sScannerManager::GetAvailableScanners(vector<string> *ap_scannerNames)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_DETAIL) Cout("sScannerManager::GetAvailableScanners() -> Scan the scanner configuration directory"<< endl); 

  // Initialize directory
  DIR *respository_dir;
  struct dirent *ent;
  string str_geom(".geom"), str_hscan(".hscan");
  string scanner_repository = sOutputManager::GetInstance()->GetPathToConfigDir() + "scanner" + OS_SEP;

  // Open directory
  if ((respository_dir = opendir(scanner_repository.c_str())) != NULL) 
  {
    // Print all the files and directories within the repository directory
    while ((ent = readdir (respository_dir)) != NULL) 
    {
      string scanner_name = ent->d_name;
      // Get rid of backup files in linux
      if(scanner_name.at(scanner_name.size()-1) != '~')
      // .geom or .hscan file found
      if(scanner_name.find(str_geom)!=string::npos || scanner_name.find(str_hscan)!=string::npos )
      {
        ap_scannerNames->push_back(scanner_name);
      }
    }
    // Close directory
    closedir (respository_dir);
  } 
  else 
  {
    Cerr("***** sScannerManager::GetAvailableScanners() -> Could not open the repository directory at: " << scanner_repository << endl);
    return 1;
  }

  // End
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      int sScannerManager::GetModalityFromString()
  \param   a_systemStr : String corresponding to a system (PET, CT, etc..)
  \brief   A simple utility function which returns the integer \n
           corresponding to the system string  passed in parameter
  \return  The integer corresponding to the scaner, as defined in the SCANNER_TYPE macro 
*/
int sScannerManager::GetModalityFromString(string a_systemStr)
{
  if( a_systemStr == "PET" ) 
    return SCANNER_PET;
  else if( a_systemStr == "SPECT_PINHOLE" )
    return SCANNER_SPECT_PINHOLE;
  else if( a_systemStr == "SPECT_CONVERGENT" )
    return SCANNER_SPECT_CONVERGENT;
  else if( a_systemStr == "CT" )
    return SCANNER_CT;
  else if( a_systemStr == "SINOGRAM" )
    return SCANNER_SINOGRAM;
  else
    return SCANNER_UNKNOWN;
  
  return SCANNER_UNKNOWN;
}
    
    
    
    
// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetScannerLayerNbRings
  \param a_layer : layer index
  \brief Ask the number of rings to the scanner object for a specific layer. 
         Returns an error if this information is not available for the scanner type of the object (eg : SPECT systems) 
  \return The number of rings in the system if success. NEGATIVE value otherwise

int sScannerManager::GetScannerLayerNbRings(int a_layer)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_DETAIL) Cout("sScannerManager::GetScannerLayerNbRings() -> For layer " << a_layer << endl); 

  // Works only for PET scanners
  if (mp_Scanner->GetScannerType() == SCANNER_PET) 
  {
    // Get the number of rings
    int nb_rings = mp_Scanner->GetScannerLayerNbRings(a_layer);
    // Check for error
    if (nb_rings <= 0)
    {
      Cerr("***** sScannerManager::GetScannerLayerNbRings() -> Error when trying to get the number of rings in the system for the crystal layer " << a_layer+1 << " !" << endl);
      return -1;
    }
    else return nb_rings;
  }
  else
  {
    Cerr("***** sScannerManager::GetScannerLayerNbRings() -> This function is only available for PET scanners !" << endl);
    return -1;
  }
}
*/

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int sScannerManager::GetSPECTSpecificParameters(uint16_t* ap_nbOfProjections, 
                                                uint16_t* ap_nbHeads,
                                                FLTNB* ap_acquisitionZoom,
                                                uint16_t* ap_nbOfBins, 
                                                FLTNB* ap_pixSizeXY,
                                                FLTNB*& ap_angles, 
                                                FLTNB*& ap_CORtoDetectorDistance,
                                                int* ap_headRotDirection)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_DETAIL) Cout("sScannerManager::GetSPECTSpecificParameters() -> Get acquisition dependent parameters"<< endl); 
    
  // Check modality
  if (mp_Scanner->GetScannerType() != SCANNER_SPECT_CONVERGENT && mp_Scanner->GetScannerType() != SCANNER_SPECT_PINHOLE)
  {
    Cerr("***** sScannerManager::GetSPECTSpecificParameters()-> The scanner object is not of SPECT type !" << endl);
    return 1;
  }
  else
  {
    // TODO: SS: remove this function from the vScanner, let it be specific to the SPECT scanner and then do a dynamic_cast here
    //           and do the same for all functions declared in vScanner as virtual but that are entirely specific !
    if (mp_Scanner->GetSPECTSpecificParameters(ap_nbOfProjections, ap_nbHeads, ap_acquisitionZoom, ap_nbOfBins, ap_pixSizeXY, ap_angles, ap_CORtoDetectorDistance, ap_headRotDirection) )
    {
      Cerr("***** sScannerManager::GetSPECTSpecificParameters()-> A problem occurred while retrieving SPECT parameters from the scanner object !" << endl);
      return 1;
    }
  }
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn GetCTSpecificParameters (Reconstruction)
  \param ap_nbOfProjections : number of views/projections
  \param ap_angles : an array containing angles for each view
  \param ap_headRotDirection : head rotation direction
  \brief Transfer geometric information recovered from the datafile to the scanner object, by copy of pointers
  \return 0 if success, positive value otherwise
*/
int sScannerManager::GetCTSpecificParameters(uint16_t* ap_nbOfProjections, 
                                             FLTNB*& ap_angles, 
                                             int* ap_headRotDirection)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // Verbose
  if (m_verbose>=VERBOSE_DETAIL) Cout("sScannerManager::GetCTSpecificParameters() -> Get acquisition dependent parameters"<< endl); 
    
  // Check modality
  if (mp_Scanner->GetScannerType() != SCANNER_CT)
  {
    Cerr("***** sScannerManager::GetCTSpecificParameters()-> The scanner object is not of CT type !" << endl);
    return 1;
  }
  else
  {
    // TODO: SS: remove this function from the vScanner, let it be specific to the SPECT scanner and then do a dynamic_cast here
    //           and do the same for all functions declared in vScanner as virtual but that are entirely specific !
    if (mp_Scanner->GetCTSpecificParameters(ap_nbOfProjections, ap_angles, ap_headRotDirection) )
    {
      Cerr("***** sScannerManager::GetCTSpecificParameters()-> A problem occurred while retrieving CT parameters from the scanner object !" << endl);
      return 1;
    }
  }
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_GetModalityStopValueMainLoop
  \brief Get the stop value for the main loop of analytic projection depending on the modality 
  \return the required stop value if success, NEGATIVE value otherwise
  \todo Cases for CT and sinogram scanner types
*/
int64_t sScannerManager::PROJ_GetModalityStopValueMainLoop()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  if (m_verbose>=VERBOSE_DETAIL) Cout("sScannerManager::PROJ_GetModalityStopValueMainLoop ..."<< endl); 
  
  if (mp_Scanner->GetScannerType() == SCANNER_PET)
  {
    return mp_Scanner->GetSystemNbElts();
  }
  else if(mp_Scanner->GetScannerType() == SCANNER_SPECT_PINHOLE ||
          mp_Scanner->GetScannerType() == SCANNER_SPECT_CONVERGENT )
  {
    return (int64_t)mp_Scanner->PROJ_GetSPECTNbProjections(); // cast from uint16_t to int
  }
  else if(mp_Scanner->GetScannerType() == SCANNER_CT)
  {
    // TODO
    Cerr("sScannerManager::PROJ_GetModalityStopValueMainLoop()-> Not implemented for CT yet!" << endl);
    return -1;
  }
  else if(mp_Scanner->GetScannerType() == SCANNER_SINOGRAM)
  {
    // TODO + should call a GetType() like function for scanner defined as a sinogram
    Cerr("sScannerManager::PROJ_GetModalityStopValueMainLoop()-> Not implemented for Sinogram scanner yet!" << endl);
    return -1;
  }
  else
  {
    Cerr("sScannerManager::PROJ_GetModalityStopValueMainLoop()-> Unknown scanner type!" << endl);
    return -1;
  }
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_GetModalityStartValueInnerLoop
  \param a_elt1 : Current nb of processed crystals (PET), projections (SPECT)
  \brief Get the start value for the inner loop of analytic projection depending on the modality 
  \return the required stop value if success, NEGATIVE value otherwise
  \todo Cases for CT and sinogram scanner types
  \todo Precise with SS on which index SPECT inner loop should be done
*/
int64_t sScannerManager::PROJ_GetModalityStartValueInnerLoop(int64_t a_elt1)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  if (mp_Scanner->GetScannerType() == SCANNER_PET)
  {
    return a_elt1+1;
  }
  else if(mp_Scanner->GetScannerType() == SCANNER_SPECT_PINHOLE ||
          mp_Scanner->GetScannerType() == SCANNER_SPECT_CONVERGENT )
  {
    return 0; // (first crystal)
    // TODO SS: here the inner loop should be done on the number of views
  }
  else if(mp_Scanner->GetScannerType() == SCANNER_CT)
  {
    // TODO
    Cerr("sScannerManager::PROJ_GetModalityStopValueMainLoop()-> Not implemented for CT yet!" << endl);
    return -1;
  }
  else if(mp_Scanner->GetScannerType() == SCANNER_SINOGRAM)
  {
    // TODO + should call a GetType() like function for scanner defined as a sinogram
    Cerr("sScannerManager::PROJ_GetModalityStopValueMainLoop()-> Not implemented for Sinogram scanner yet!" << endl);
    return -1;
  }
  else
  {
    Cerr("sScannerManager::PROJ_GetModalityStopValueMainLoop()-> Unknown scanner type!" << endl);
    return -1;
  }
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_GetCurrentProgression
  \param a_elt1 : Current nb of processed #1 crystals (PET), projections (SPECT)
  \param a_elt2 : Current nb of processed #2 crystals (PET), crystals (SPECT)
  \param ap_nbEltsArray : Total number of elements processed for each #1 crystals (PET/CT systems) 
  \param a_nbRGates
  \param a_nbCGates
  \param a_fr
  \param a_rg
  \param a_cg
  \brief Get numerator value according to the modality to compute percent progression during the analytical projection process
  \return the required progression value if success, negative value otherwise
  \todo Cases for CT and sinogram scanner types
  \todo Optimize this, for now it's quite a lot of operation for each couple of elements
  \todo Check everything is ok for 3D/4D PET and SPECT 
*/
int64_t sScannerManager::PROJ_GetCurrentProgression(int64_t a_elt1, int64_t a_elt2, int64_t* ap_nbEltsArray, 
                                                     int a_nbRGates, int a_nbCGates,
                                                     int a_fr, int a_rg, int a_cg)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  // TODO : optimization, maybe too much operations for each couple of elements
  if (mp_Scanner->GetScannerType() == SCANNER_PET)
  {
    int64_t nb_total_elts = mp_Scanner->GetSystemNbElts();

    return ap_nbEltsArray[a_elt1]+a_elt2+1 + 
           (int64_t)(a_fr*a_nbRGates*a_nbCGates + a_rg*a_nbCGates + a_cg)* nb_total_elts*nb_total_elts/2;
  }
  else if(mp_Scanner->GetScannerType() == SCANNER_SPECT_PINHOLE ||
          mp_Scanner->GetScannerType() == SCANNER_SPECT_CONVERGENT )
  {
    return a_elt2 + a_elt1*mp_Scanner->GetSystemNbElts() + 
           a_fr * a_nbRGates * a_nbCGates * mp_Scanner->GetSystemNbElts() +
           a_rg * a_nbCGates * mp_Scanner->GetSystemNbElts() +
           a_cg * mp_Scanner->GetSystemNbElts();
  }
  else if(mp_Scanner->GetScannerType() == SCANNER_CT)
  {
    //TODO
    Cerr("sScannerManager::PROJ_GetCurrentProgression()-> Not implemented for CT yet!" << endl);
    return -1;
  }
  else if(mp_Scanner->GetScannerType() == SCANNER_SINOGRAM)
  {
    //TODO + should call a GetType() like function for scanner defined as a sinogram
    Cerr("sScannerManager::PROJ_GetCurrentProgression()-> Not implemented for Sinogram scanner yet!" << endl);
    return -1;
  }
  else
  {
    Cerr("sScannerManager::PROJ_GetCurrentProgression()-> Unknown scanner type!" << endl);
    return -1;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_GetProgressionFinalValue
  \brief Get numerator value according to the modality to compute percent progression during the projection process
  \return the required progression value if success, negative value otherwise
  \todo Cases for CT and sinogram scanner types
*/
int64_t sScannerManager::PROJ_GetProgressionFinalValue()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  if (mp_Scanner->GetScannerType() == SCANNER_PET)
  {
    //return mp_Scanner->GetSystemNbElts();
    int64_t nb_total_elts = mp_Scanner->GetSystemNbElts();
    return nb_total_elts*nb_total_elts/2;
  }
  else if(mp_Scanner->GetScannerType() == SCANNER_SPECT_PINHOLE ||
          mp_Scanner->GetScannerType() == SCANNER_SPECT_CONVERGENT )
  {
    return (int64_t)mp_Scanner->PROJ_GetSPECTNbProjections()*mp_Scanner->GetSystemNbElts(); // cast from uint16_t to int
  }
  else if(mp_Scanner->GetScannerType() == SCANNER_CT)
  {
    //TODO
    Cerr("sScannerManager::PROJ_GetProgressionFinalValue()-> Not implemented for CT yet!" << endl);
    return -1;
  }
  else if(mp_Scanner->GetScannerType() == SCANNER_SINOGRAM)
  {
    //TODO + should call a GetType() like function for scanner defined as a sinogram
    Cerr("sScannerManager::PROJ_GetProgressionFinalValue()-> Not implemented for Sinogram scanner yet!" << endl);
    return -1;
  }
  else
  {
    Cerr("sScannerManager::PROJ_GetProgressionFinalValue()-> Unknown scanner type!" << endl);
    return -1;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_SetPETSpecificParameters (Analytic Projection)
  \param a_maxAxialDiffmm : max axial difference in mm between 2 crystals forming a lor
  \brief Deliver to the PET scanner object all informations provided from the datafile header
  \return 0 if success, positive value otherwise
  \todo How to handle systems with several layer of rings ?
*/
int sScannerManager::PROJ_SetPETSpecificParameters(FLTNB a_maxAxialDiffmm)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  if (m_verbose>=VERBOSE_DETAIL) Cout("sScannerManager::PROJ_SetPETSpecificParameters ..."<< endl); 
  
  if (mp_Scanner->GetScannerType() != SCANNER_PET)
  {
    Cerr("***** sScannerManager::SetPETSpecificParameters() -> The scanner object is not of PET type !" << endl);
    return 1;
  }
  else
  {
    mp_Scanner->SetPETMaxAxialDiffmm(a_maxAxialDiffmm);
  }
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_GetPETSpecificParameters (Analytic Projection)
  \param ap_maxRingDiff : maximal axial difference in mm between 2 crystals forming a lor
  \brief Transfer addresses to each geometric parameter of the PET scanner objets to the corresponding pointer of the datafile passed as argument
  \return 0 if success, positive value otherwise
  \todo How to handle systems with several layer of rings ?
*/
int sScannerManager::PROJ_GetPETSpecificParameters(FLTNB* ap_maxRingDiff)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  if (m_verbose>=VERBOSE_DETAIL) Cout("sScannerManager::PROJ_GetPETSpecificParameters ..."<< endl); 
    
  if (mp_Scanner->GetScannerType() != SCANNER_PET )
  {
    Cerr("***** sScannerManager::PROJ_GetPETSpecificParameters()-> The scanner object is not of PET type !" << endl);
    return 1;
  }
  else
  {
    // Negative value means no ring diff restriction
    FLTNB max_ring_diff = -1.;
    
    if (mp_Scanner->PROJ_GetPETSpecificParameters(&max_ring_diff) )
    {
      Cerr("***** sScannerManager::PROJ_GetPETSpecificParameters()-> A problem occurred while retrieving PET parameters from the scanner object !" << endl);
      return 1;
    }
    
    *ap_maxRingDiff = max_ring_diff;
  }
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_SetSPECTSpecificParameters (Analytic Projection)
  \param ap_nbOfBins : 2 elements array containing transaxial number of bins
  \param a_nbOfProjections : number of views/projections
  \param a_firstAngle : angle of the first view
  \param a_lastAngle : angle of the last view
  \param ap_projectionAngles : an array containing angles for each view
  \param a_CORtoDetectorDistance : a distance between the center of rotation and the detector
  \param a_RotDirection : Rotation direction of the head (clockwise/counter-clockwise)
  \brief Deliver to the SPECT scanner object all informations provided from the acquisition parameters
  \details For analytical projection, this data is provided from the command-line options 
  \return 0 if success, positive value otherwise
*/
int sScannerManager::PROJ_SetSPECTSpecificParameters(uint16_t* ap_nbOfBins, 
                                                     uint32_t a_nbOfProjections, 
                                                     FLTNB a_firstAngle, 
                                                     FLTNB a_stepAngle, 
                                                     FLTNB* ap_projectionAngles, 
                                                     FLTNB a_CORtoDetectorDistance,
                                                     string a_rotDirection)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)

  if (m_verbose>=VERBOSE_DETAIL) Cout("sScannerManager::PROJ_SetSPECTSpecificParameters ..."<< endl); 
    
  // Return error if wrong modality
  if (mp_Scanner->GetScannerType() != SCANNER_SPECT_CONVERGENT && mp_Scanner->GetScannerType() != SCANNER_SPECT_PINHOLE)
  {
    Cerr("***** sScannerManager::SetPETSpecificParameters()-> The scanner object is not of PET type !" << endl);
    return 1;
  }
  else
  {
    // Nb bins
    mp_Scanner->PROJ_SetSPECTNbBins(ap_nbOfBins);
    
    // Nb projections
    if(a_nbOfProjections == 0)
    {
      Cerr("***** sScannerManager::SetSPECTSpecificParameters()-> Error, SPECT analytic projection requires an user-specified number of projections !" << endl);
      return 1;
    }
    
    mp_Scanner->PROJ_SetSPECTNbProjections(a_nbOfProjections);


    // Projection angles initialized with either each angles, or calculated from first/last angle
    // By default, we chose the user-provided custom initialization for SPECT projection angles
    // If not initialized, we computed them from first/last angle
    if(ap_projectionAngles == NULL)
    {
      // If no custom initialization have been provided, then we check if we have first and last angles values for automatic initialization of the projection angles
      if(a_firstAngle<0 || a_stepAngle<0)
      {
        Cerr("***** sScannerManager::SetSPECTSpecificParameters()-> Error, SPECT projection requires to set the projection angles using either the '-SPECT_ang' or '-SPECT_c_ang' options !" << endl);
        return 1;
      }
      else
      {
        // Fill the SPECT_projection_angles array
        ap_projectionAngles = new FLTNB[a_nbOfProjections];
        ap_projectionAngles[0] = a_firstAngle;
        
        for(uint32_t a=1 ; a<a_nbOfProjections ; a++)
          ap_projectionAngles[a] = ap_projectionAngles[a-1] + a_stepAngle;
      }
    }
    
    if(mp_Scanner->PROJ_SetSPECTAngles(ap_projectionAngles) )
    {
      Cerr("***** sScannerManager::SetSPECTAngles() -> An error occured while trying to initialize SPECT projection angles !" << endl);
      return 1;
    }

    // CORtoDetectorDistance
    if(mp_Scanner->PROJ_SetSPECTCORtoDetectorDistance(a_CORtoDetectorDistance) )
    {
      Cerr("***** sScannerManager::SetPETSpecificParameters() -> An error occured while trying to initialize SPECT distance between center of rotation to detectors !" << endl);
      return 1;
    }
    
    // Set rotation direction
    mp_Scanner->SetRotDirection(a_rotDirection);
  }
  
  return 0;
}
