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
  \brief    Declaration of class sScannerManager
*/

#ifndef SSCANNERMANAGER_HH
#define SSCANNERMANAGER_HH 1

#include "gVariables.hh"
#include "gOptions.hh"
#include "sOutputManager.hh"
#ifdef _WIN32
#include "oDirentWin32.hh"
#else
#include "dirent.h"
#endif
#include "vScanner.hh"

/**
 * @defgroup SCANNER_TYPE system type
 *
 *    \brief System type (PET, CT, SPECT, ...) \n
 *           Defined in sScannerManager.hh
 * @{
 */
/** Constant corresponding to an unknown system (=-1) */
#define SCANNER_UNKNOWN -1
/** Constant corresponding to a PET system (=0) */
#define SCANNER_PET 0
/** Constant corresponding to a SPECT system with pinhole collimator (=1) */
#define SCANNER_SPECT_PINHOLE 1
/** Constant corresponding to a SPECT system with a convergent/parallel collimator (=2) */
#define SCANNER_SPECT_CONVERGENT 2
/** Constant corresponding to a CT system (=3) */
#define SCANNER_CT 3
/** Constant corresponding to any system whose geometry is defined by sinograms (=4) */
#define SCANNER_SINOGRAM 4
/** @} */


/*!
  \class   sScannerManager
  \brief   Singleton class that Instantiate and initialize the scanner object.
  \details This class Instantiate and initialize the scanner object depending on its initilization file (generic .geom file or user LUT). \n
           It holds several informations on the system and can be accessed from every class.
*/
class sScannerManager
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      static sScannerManager* sScannerManager::GetInstance()
      \brief   Instanciate the singleton object and Initialize member variables if not already done, 
               return a pointer to this object otherwise
      \return  instance of the sScannerManager singleton
    */
    static sScannerManager* GetInstance() 
    {
      if (mp_Instance == NULL)
        mp_Instance = new sScannerManager();
      return mp_Instance;
    }

  // -----------------------------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      sScannerManager::Describe()
      \brief   Call the eponym function from the Scanner object (if initialized)
    */
    void Describe();
    /*!
      \fn      int sScannerManager::CheckParameters()
      \brief   Check if all parameters have been correctly initialized, and call the CheckParameters function of the scanner object
      \return  0 if success. Positive value otherwise
    */
    int CheckParameters();
    /*!
      \fn      int sScannerManager::Initialize()
      \brief   Initialization : \n
               - check if all parameters of the manager have been checked \n
               - call the initialization function of the scanner object 
      \return  0 if success. Positive value otherwise
    */
    int Initialize();
    /*!
      \fn      int sScannerManager::ShowScannersDescription()
      \brief   Get the description associated to the different scanners and print all on screen. \n
               Walk through the scanner repository and look for the keyword "description" in .geom file and .hscan file.
      \return  0 if success, positive value otherwise
      \todo    Check everything output correctly for all implemented scanners
    */
    int ShowScannersDescription();
    /*!
      \fn      int sScannerManager::FindScannerSystem()
      \param   a_scannerName : string containing name of the required scanner
      \brief   Look for a file matching with the scanner name in parameter inside the scanner repository 
      \return  0 if success (scanner found). Positive value otherwise
    */
    int FindScannerSystem(string a_scannerName);
    /*!
      \fn      int sScannerManager::InitScannerWithFile()
      \param   a_pathScanFile : string containing the path to the scanner file
      \param   a_scannerName  : string containing the name of the required scanner
      \param   a_fileTypeFlag : string containing the type of the scanner file (0=lut, 1=geom)
      \brief   Initialize member variables (file path, file type, and scanner name) with the provided arguments
      \return  0 if success (scanner found). Positive value otherwise
    */
    int InitScannerWithFile(string a_pathScanFile, string a_scannerName, int a_fileTypeFlag);
    /*!
      \fn      int sScannerManager::BuildScannerObject()
      \brief   Instantiate the specific scanner object related to the modality, and set verbosity of scanner object
      \todo    delete the check on scanner type once all scanner classes will be completely implemented. 
      \return  0 if success. Positive value otherwise
    */
    int BuildScannerObject();
    /*!
      \fn      int sScannerManager::GetGeometricInfoFromDataFile()
      \param   a_path : string containing the path to datafile header
      \brief   Call the specialized function of the scanner object in order to 
               get geometric informations from the datafile header
      \return  0 if success. Positive value otherwise
    */
    int GetGeometricInfoFromDataFile(string a_pathToDataFilename);
    /*!
      \fn      int sScannerManager::InstantiateScanner()
      \brief   Instantiate scanner using the related function in the scanner classes
      \todo    delete the check on scanner type once all scanner classes will be completely implemented. 
      \return  0 if success. Positive value otherwise
    */
    int InstantiateScanner();
    /*!
      \fn      int sScannerManager::BuildLUT()
      \brief   Call the eponym function of the scanner class
      \return  0 if success. Positive value otherwise
    */
    int BuildLUT();
    /*!
      \fn      int sScannerManager::GetScannerLayerNbRings()
      \param   a_layer: layer index
      \brief   DEPRECATED
               Ask the number of rings to the scanner object for a specific layer. \n
               Returns an error if this information is not available for the scanner type of the object (eg : SPECT systems) 
      \return  The number of rings in the system if success. NEGATIVE value otherwise
    
    int GetScannerLayerNbRings(int a_layer);*/
    /*!
      \fn      int sScannerManager::GetModalityFromString()
      \param   a_systemStr : String corresponding to a system (PET, CT, etc..)
      \brief   A simple utility function which returns the integer \n
               corresponding to the system string  passed in parameter
      \return  The integer corresponding to the scaner, as defined in the SCANNER_TYPE macro 
    */
    int GetModalityFromString(string a_systemStr);


  // -------------------------------------------------------------------
  // Get & Set functions
  public:
    /*!
      \fn      inline bool sScannerManager::HasUserScannerFile()
      \return  true if the scanner file is an user LUT, false otherwise (generic file is used)
    */
    inline bool HasUserScannerFile()
           {return m_hasUserScannerFile;}
    /*!
      \fn      inline string sScannerManager::GetPathToScannerFile()
      \return  the path to the scanner file
    */
    inline string GetPathToScannerFile()
           {return m_pathToScannerFile;}
    /*!
      \fn      inline string sScannerManager::GetScannerName()
      \return  the scanner name
    */
    inline string GetScannerName()
           {return m_scannerName;}
    /*!
      \fn      inline int sScannerManager::GetScannerType()
      \return  the scanner type as returned by the scanner object,
               or SCANNER_UNKNOWN if the scanner object has not been initialized
    */
    inline int GetScannerType()
           {return (mp_Scanner!=NULL) ? mp_Scanner->GetScannerType() : SCANNER_UNKNOWN ;}
    /*!
      \fn      inline vScanner* sScannerManager::GetScannerObject()
      \return  the vScanner object
    */
    inline vScanner* GetScannerObject()
           {return mp_Scanner;}
    /*!
      \fn      inline int sScannerManager::GetSystemNbElts()
      \return  the number of elements in the system
    */
    inline int GetSystemNbElts()
           {return mp_Scanner->GetSystemNbElts();}
    /*!
      \fn      inline void sScannerManager::SetVerbose()
      \param   a_verboseLevel
      \brief   set verbosity
    */
    inline void SetVerbose(int a_verboseLevel)
           {m_verbose = a_verboseLevel;}
    /*!
      \fn      inline void sScannerManager::SetImageDimensionsAndQuantification()
      \param   oImageDimensionsAndQuantification* ap_ID
      \brief   Set the pointer to the image dimensions and quantification object
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ID)
           {mp_ID = ap_ID;}
    /*!
      \fn      inline void sScannerManager::SetSaveLUTFlag()
      \param   a_flag
      \brief   Set to on the flag indicating a LUT generated by a geom file
               should be written on disk or not.
    */
    inline void SetSaveLUTFlag(bool a_flag)
           {m_saveLUTFlag = a_flag;}
    /*!
      \fn      inline bool sScannerManager::SaveLUTFlag()
      \brief   Get flag indicating a LUT generated by a geom file
               should be written on disk or not.
      \return  true if scanner LUT has to be saved, false otherwise
    */
    inline bool SaveLUTFlag()
           {return m_saveLUTFlag;}
    
    
  // -------------------------------------------------------------------
  // Analytical projection functions
  public:
    /*!
      \fn      int64_t sScannerManager::PROJ_GetModalityStopValueMainLoop()
      \brief   Get the stop value for the main loop of analytic projection depending on the modality 
      \return  the required stop value if success, NEGATIVE value otherwise
      \todo    Cases for CT and sinogram scanner types
    */
    int64_t PROJ_GetModalityStopValueMainLoop();
    /*!
      \fn      int64_t sScannerManager::PROJ_GetModalityStartValueInnerLoop()
      \param   a_elt1 : Current nb of processed crystals (PET), projections (SPECT)
      \brief   Get the start value for the inner loop of analytic projection depending on the modality 
      \return  the required stop value if success, NEGATIVE value otherwise
      \todo    Cases for CT and sinogram scanner types
    */
    int64_t PROJ_GetModalityStartValueInnerLoop(int64_t a_elt1);
    /*!
      \fn      int64_t sScannerManager::PROJ_GetCurrentProgression()
      \param   a_elt1 : Current nb of processed #1 crystals (PET), projections (SPECT)
      \param   a_elt2 : Current nb of processed #2 crystals (PET), crystals (SPECT)
      \param   ap_nbEltsArray : Total number of elements processed for each #1 crystals (PET/CT systems) 
      \param   a_nbRGates
      \param   a_nbCGates
      \param   a_fr
      \param   a_rg
      \param   a_cg
      \brief   Get numerator value according to the modality to compute percent progression during the analytical projection process
      \return  the required progression value if success, negative value otherwise
      \todo    Cases for CT and sinogram scanner types
      \todo    Optimize this, for now it's quite a lot of operation for each couple of elements
      \todo    Check everything is ok for 3D/4D PET and SPECT 
    */
    int64_t PROJ_GetCurrentProgression(int64_t a_elt1, int64_t a_elt2, int64_t* ap_nbEltsArray, int a_nbRGates, int a_nbCGates, int a_fr, int a_rg, int a_cg);
    /*!
      \fn      int64_t sScannerManager::PROJ_GetProgressionFinalValue()
      \brief   Get numerator value according to the modality to compute percent progression during the projection process
      \return  the required progression value if success, negative value otherwise
      \todo    Cases for CT and sinogram scanner types
    */
    int64_t PROJ_GetProgressionFinalValue();


  // -------------------------------------------------------------------
  // Functions dedicated to transfer geometric information for Projection, from the main class to the scanner class (Set)
  //                                                                       from the scanner class to the datafile (Get, for datafile header writing)
  public:
    // RECONSTRUCTION
    /*!
      \fn      int sScannerManager::GetSPECTSpecificParameters()
      \param   ap_nbOfProjections : number of views/projections
      \param   ap_nbHeads : number of heads in the SPECT system
      \param   ap_acquisitionZoom : zoom during acquisition for monolithic detectors
      \param   ap_nbOfBins : 2 elements array containing transaxial number of bins
      \param   ap_pixSizeXY : 2 elements array containing transaxial/axial pixel sizes
      \param   ap_angles : an array containing angles for each view
      \param   ap_CORtoDetectorDistance : a distance between the center of rotation and the detector  
      \param   ap_headRotDirection : head rotation direction
      \brief   Transfer geometric information recovered from the datafile to the scanner object
      \return  0 if success, positive value otherwise
    */
    int GetSPECTSpecificParameters(uint16_t* ap_nbOfProjections, 
                                   uint16_t* ap_nbHeads,
                                   FLTNB* ap_acquisitionZoom,
                                   uint16_t* ap_nbOfBins, 
                                   FLTNB*  ap_pixSizeXY,
                                   FLTNB*& ap_angles, 
                                   FLTNB*& ap_CORtoDetectorDistance,
                                   int* ap_headRotDirection);
    /*!
      \fn      int sScannerManager::GetCTSpecificParameters()
      \param   ap_nbOfProjections : number of views/projections
      \param   ap_angles : an array containing angles for each view
      \param   ap_headRotDirection : head rotation direction
      \brief   Transfer geometric information recovered from the datafile to the scanner object
      \return  0 if success, positive value otherwise
    */
    int GetCTSpecificParameters(uint16_t* ap_nbOfProjections, 
                                FLTNB*& ap_angles, 
                                int* ap_headRotDirection);
    // ANALYTICAL PROJECTION
    /*!
      \fn      int sScannerManager::PROJ_GetPETSpecificParameters()
      \param   ap_maxAxialDiffmm : max axial difference in mm between 2 crystals forming a lor
      \brief   Transfer addresses to each geometric parameter of the PET scanner objets to the corresponding pointer of the datafile passed as argument
      \return  0 if success, positive value otherwise
      \todo    How to handle systems with several layer of rings ?
    */
    int PROJ_GetPETSpecificParameters(FLTNB* ap_maxAxialDiffmm);
    /*!
      \fn      int sScannerManager::PROJ_SetPETSpecificParameters()
      \param   a_maxAxialDiffmm : max axial difference in mm between 2 crystals forming a lor
      \brief   Deliver to the PET scanner object all informations provided from the datafile header
      \return  0 if success, positive value otherwise
      \todo    How to handle systems with several layer of rings ?
    */
    int PROJ_SetPETSpecificParameters(FLTNB a_maxAxialDiffmm);
    /*!
      \fn      int sScannerManager::PROJ_SetSPECTSpecificParameters()
      \param   ap_nbOfBins : 2 elements array containing transaxial number of bins
      \param   a_nbOfProjections : number of views/projections
      \param   a_firstAngle : angle of the first view
      \param   a_lastAngle : angle of the last view
      \param   ap_projectionAngles : an array containing angles for each view
      \param   a_CORtoDetectorDistance : a distance between the center of rotation and the detector
      \param   a_RotDirection : Rotation direction of the head (clockwise/counter-clockwise)
      \brief   Deliver to the SPECT scanner object all informations provided from the acquisition parameters
      \details For analytical projection, this data is provided from the command-line options 
      \return  0 if success, positive value otherwise
    */
    int PROJ_SetSPECTSpecificParameters(uint16_t* ap_nbOfBins,
                                        uint32_t a_nbOfProjections, 
                                        FLTNB a_firstAngle, 
                                        FLTNB a_stepAngle, 
                                        FLTNB* ap_projectionAngles, 
                                        FLTNB a_CORtoDetectorDistance,
                                        string a_rotDirection);


  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \brief   sScannerManager constructor.
      \details It is private at this class is singleton. \n
               It should be instanciated using the GetInstance() function \n 
               Initialize the member variables to their default values.
    */
    sScannerManager();
    /*!
      \brief   sScannerManager destructor. 
    */
    ~sScannerManager();
    // Prevent the compiler to generate methods to copy the object :
    sScannerManager(sScannerManager const&){};
    void operator=(sScannerManager const&){};
    /*!
      \fn      int sScannerManager::GetAvailableScanners()
      \param   ap_scannerNames : vector list of string to recover the available scanner names
      \brief   Gather all the names of the header files (.geom & .hscan) 
               in the repository folder in the vector<string> passed in parameter
      \return  0 if sucess, positive value otherwise
    */
    int GetAvailableScanners(vector<string> *ap_scannerNames);


  // -------------------------------------------------------------------
  // Data members
  private:
    static sScannerManager* mp_Instance;       /*!< Pointer to the instance of this class */
    vScanner* mp_Scanner;                      /*!< Pointer to the Scanner object */
    oImageDimensionsAndQuantification* mp_ID;  /*!< Pointer to the image dimensions and quantification object */
    int m_verbose;                             /*!< Verbosity */
    string m_pathToScannerFile;                /*!< String containing the path to the scanner file */
    string m_scannerName;                      /*!< String containing the scanner name */
    bool m_hasUserScannerFile;                 /*!< Boolean indicating if the scanner geometry if defined by a user LUT */
    bool m_hasGenericScannerFile;              /*!< Boolean indicating if the scanner geometry if defined by a generic file */
    bool m_allParametersChecked;               /*!< Boolean indicating if all variables of the class have been checked */
    bool m_saveLUTFlag;                        /*!<Flag indicating a LUT generated by a geom file should be written on disk or not.*/
};

#endif
