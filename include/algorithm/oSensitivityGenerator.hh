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
  \ingroup  algorithm
  \brief    Declaration of class oSensitivityGenerator
*/

#ifndef OSENSITIVITYGENERATOR_HH
#define OSENSITIVITYGENERATOR_HH 1

#include "gVariables.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "oProjectorManager.hh"
#include "oImageConvolverManager.hh"
#include "oDeformationManager.hh"
#include "oDynamicModelManager.hh"
#include "oImageSpace.hh"
#include "vDataFile.hh"
#include "vDeformation.hh"
#include "sOutputManager.hh"
#include "vScanner.hh"



/*!
  \class   oSensitivityGenerator
  \brief   This class is designed to manage the computation of the sensitivity image
  \details The sensitivity computation can be done in two ways: \n
           (i) using a normalization data file if provided, where the loop for computation will be done on all events included in the datafile, with included normalization factors; \n
           (ii) if no normalization file is provided, using a loop on all scanner elements. \n
           A mumap in cm-1 can also be provided in order to include the attenuation correction in the sensitivity image.
*/
class oSensitivityGenerator
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public oSensitivityGenerator::oSensitivityGenerator()
      \brief   The constructor of oSensitivityGenerator
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    oSensitivityGenerator();
    /*!
      \fn      public oSensitivityGenerator::~oSensitivityGenerator()
      \brief   The destructor of oSensitivityGenerator
      \details This is the default and unique destructor. It does not take any parameter and 
               its role is only to free or delete all structures that were built by this class.
    */
    ~oSensitivityGenerator();


  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      public int oSensitivityGenerator::CheckParameters()
      \brief   A public function used to check the parameters settings
      \details This function does not take any parameter and is used to check that all mandatory
               members were correctly parameterized.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckParameters();
    /*!
      \fn      public int oSensitivityGenerator::Initialize()
      \brief   A public function used to initialize the sensitivity generator
      \details This function does not take any parameter and is used to initialize everything that
               should be initialized.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */
    int Initialize();
    /*!
      \fn      public int oSensitivityGenerator::Launch()
      \brief   A public function used to launch the sensitivity generator (compute the sensitivity image)
      \details This function does not take any parameter and is used to launch the computation of the
               sensitivity image. In fact, it simply switches between the CPU and GPU versions of the
               function (which are private).
      \return  An integer reflecting the computation status; 0 if no problem, another value otherwise.
    */
    int Launch();


  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      private int oSensitivityGenerator::InitializeAttenuationFiles()
      \brief   Initialize the attenuation images provided for sensitivity computation
      \details This function is called by the Initialize() function. It checks and reads the provided
               attenuation images to be taken into account for sensitivity computation.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */
    int InitializeAttenuationFiles();
    /*!
      \fn      private int oSensitivityGenerator::InitializeNormalizationFiles()
      \brief   Initialize the normalization datafiles provided for sensitivity computation
      \details This function is called by the Initialize() function. It checks and reads the provided
               normalization datafiles to be taken into account for sensitivity computation.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */
    int InitializeNormalizationFiles();
    /*!
      \fn      private int oSensitivityGenerator::LaunchCPU()
      \brief   Launch the computation of the sensitivity image (CPU version)
      \details This function calls either the ComputeSensitivityFromNormalizationFile() or the
               ComputeSensitivityFromScanner() function based on normalization data files provided or not. \n
               It does the loop over all bed positions.
      \return  An integer reflecting the computation status; 0 if no problem, another value otherwise.
    */
    int LaunchCPU();
    #ifdef CASTOR_GPU
    /*!
      \fn      private int oSensitivityGenerator::LaunchGPU()
      \brief   Launch the computation of the sensitivity image (GPU version)
      \details This function is not implemented
      \return  An integer reflecting the computation status; 0 if no problem, another value otherwise.
    */
    int LaunchGPU();
    #endif
    /*!
      \fn      private int oSensitivityGenerator::ComputeSensitivityFromHistogramDataFile()
      \param   int a_bed
      \brief   Launch the computation of the sensitivity image for this bed, based on the input histogram data files
      \details The loop is done over the events of the histogram data files, which represent the
               collection of all data channels for the different dynamic frames.
      \return  An integer reflecting the computation status; 0 if no problem, another value otherwise.
    */
    int ComputeSensitivityFromHistogramDataFile(int a_bed);
    /*!
      \fn      private int oSensitivityGenerator::ComputeSensitivityFromNormalizationFile()
      \param   int a_bed
      \brief   Launch the computation of the sensitivity image for this bed, based on normalization data files
      \details The loop is done over the events of the normalization data files, which represent the
               collection of all data channels for the different dynamic frames.
      \return  An integer reflecting the computation status; 0 if no problem, another value otherwise.
    */
    int ComputeSensitivityFromNormalizationFile(int a_bed);
    /*!
      \fn      private int oSensitivityGenerator::ComputeSensitivityFromScanner()
      \param   int a_bed
      \brief   Launch the computation of the sensitivity image for this bed, based on a loop over all
               scanner elements
      \return  An integer reflecting the computation status; 0 if no problem, another value otherwise.
    */
    int ComputeSensitivityFromScanner(int a_bed);
    /*!
      \fn      private int oSensitivityGenerator::ProcessThisLine()
      \param   oProjectionLine* ap_Line
      \param   vEvent* ap_Event
      \param   int a_bed
      \param   int a_frame
      \param   int a_respGate
      \param   int a_cardGate
      \param   int a_thread
      \brief   This function manages the computation of the sensitivity contribution of the projection line provided as a parameter.
      \details For the given line, the given event, the current bed, frame, respiratory and cardiac gates, the function performs the
               forward projectio into the attenuation images (if any), compute the sensitivity contribution, and back-projects this
               contribution into the sensitivity image.
      \return  An integer reflecting the projection status; 0 if no problem, another value otherwise.
    */
    int ProcessThisLine(oProjectionLine* ap_Line, vEvent* ap_Event, int a_bed, int a_frame, int a_respGate, int a_cardGate, int a_thread);


  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      public string oSensitivityGenerator::GetPathToSensitivityImage()
      \brief   This function return the path to the sensitivity image.
      \return  m_pathToSensitivityImage
    */
    string GetPathToSensitivityImage();
    /*!
      \fn      public inline void oSensitivityGenerator::SetImageDimensionsAndQuantification()
      \param   oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification
      \brief   This function is used to set the pointer to the oImageDimensionsAndQuantification in use.
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
           {mp_ImageDimensionsAndQuantification = ap_ImageDimensionsAndQuantification;};
    /*!
      \fn      public inline void oSensitivityGenerator::SetImageSpace()
      \param   oImageSpace* ap_ImageSpace
      \brief   This function is used to set the pointer to the oImageSpace in use.
    */
    inline void SetImageSpace(oImageSpace* ap_ImageSpace)
           {mp_ImageSpace = ap_ImageSpace;};
    /*!
      \fn      public inline void oSensitivityGenerator::SetScanner()
      \param   vScanner* ap_Scanner
      \brief   This function is used to set the pointer to the vScanner in use.
    */
    inline void SetScanner(vScanner* ap_Scanner)
           {mp_Scanner = ap_Scanner;}
    /*!
      \fn      public inline void oSensitivityGenerator::SetProjectorManager()
      \param   oProjectorManager* ap_ProjectorManager
      \brief   This function is used to set the pointer to the oProjectorManager in use.
    */
    inline void SetProjectorManager(oProjectorManager* ap_ProjectorManager)
           {mp_ProjectorManager = ap_ProjectorManager;};
    /*!
      \fn      public inline void oSensitivityGenerator::SetImageConvolverManager()
      \param   oImageConvolverManager* ap_ImageConvolverManager
      \brief   This function is used to set the pointer to the oImageConvolverManager in use.
    */
    inline void SetImageConvolverManager(oImageConvolverManager* ap_ImageConvolverManager)
           {mp_ImageConvolverManager = ap_ImageConvolverManager;}
    /*!
      \fn      public inline void oSensitivityGenerator::SetDeformationManager()
      \param   oDeformationManager* ap_DeformationManager
      \brief   This function is used to set the pointer to the oDeformationManager in use.
    */
    inline void SetDeformationManager(oDeformationManager* ap_DeformationManager)
           {mp_DeformationManager = ap_DeformationManager;};
    /*!
      \fn      public inline void oSensitivityGenerator::SetDataFile()
      \param   vDataFile** a2p_DataFile
      \brief   This function is used to set the pointer to the vDataFile array in use.
    */
    inline void SetDataFile(vDataFile** a2p_DataFile)
           {m2p_DataFile = a2p_DataFile;};
    /*!
      \fn      public inline void oSensitivityGenerator::SetComputeFromHistogramFlag()
      \param   bool a_computeFromHistogramFlag
      \brief   This function is used to set the m_computeFromHistogramFlag
    */
    inline void SetComputeFromHistogramFlag(bool a_computeFromHistogramFlag)
           {m_computeFromHistogramFlag = a_computeFromHistogramFlag;}
    /*!
      \fn      public inline void oSensitivityGenerator::SetPathToNormalizationFileName()
      \param   vector<string> ap_pathToNormalizationFileName
      \param   bool a_inverseDataFileOrderFlag
      \brief   This function is used to set the path to the normalization file names, and the flag saying if their order must be reversed.
    */
    inline void SetPathToNormalizationFileName(vector<string> ap_pathToNormalizationFileName, bool a_inverseDataFileOrderFlag)
           { mp_pathToNormalizationFileName = ap_pathToNormalizationFileName;
             m_inverseDataFileOrderFlag = a_inverseDataFileOrderFlag; }
    /*!
      \fn      public inline void oSensitivityGenerator::SetPathToAttenuationImage()
      \param   string a_pathToAttenuationImage
      \brief   This function is used to set the path to the attenuation image.
    */
    inline void SetPathToAttenuationImage(string a_pathToAttenuationImage)
           {m_pathToAttenuationImage = a_pathToAttenuationImage;};
    /*!
      \fn inline public oSensitivityGenerator::SetPathToMaskImage
      \param a_pathToMaskImage
      \brief Set path to a mask image
    */
    inline void SetPathToMaskImage(string a_pathToMaskImage) {m_pathToMaskImg = a_pathToMaskImage;};
    /*!
      \fn      public inline void oSensitivityGenerator::SetAtnGateImageDimensionsAndQuantification()
      \param   int a_nbAtnRespGateImages
      \param   int a_nbAtnCardGateImages
      \brief   This function is used to set the number of attenuation images for respiratory and cardiac gates.
    */
    inline void SetNumberOfAtnGateImages(int a_nbAtnRespGateImages, int a_nbAtnCardGateImages)
           {m_nbAtnRespGateImages = a_nbAtnRespGateImages; m_nbAtnCardGateImages = a_nbAtnCardGateImages;};
    /*!
      \fn      public inline void oSensitivityGenerator::SetGPUflag()
      \param   bool a_flagGPU
      \brief   This function is used to set the GPU flag; do we use GPU or not.
    */
    inline void SetGPUflag(bool a_flagGPU)
           {m_flagGPU = a_flagGPU;};
    /*!
      \fn      public inline void vProjector::SetVerbose()
      \param   int a_verboseLevel
      \brief   Set the verbose level.
    */
    inline void SetVerbose(int a_verboseLevel)
           {m_verbose = a_verboseLevel;};


  // -------------------------------------------------------------------
  // Data members
  private:
    string m_pathToSensitivityImage;                  /*!< The actual path to the created sensitivity image */
    oImageDimensionsAndQuantification* 
                 mp_ImageDimensionsAndQuantification; /*!< Pointer to the oImageDimensionsAndQuantification object */
    oImageSpace* mp_ImageSpace;                       /*!< Pointer to the Image Space object */
    vScanner* mp_Scanner;                             /*!< Pointer to the Scanner object */
    oProjectorManager* mp_ProjectorManager;           /*!< Pointer to the Projector Manager object */
    oImageConvolverManager* mp_ImageConvolverManager; /*!< Pointer to the Image Convolver Manager object */
    oDeformationManager* mp_DeformationManager;       /*!< Pointer to the Deformation Manager object */
    vDataFile** m2p_DataFile;                         /*!< Pointer to the array of vDataFile objects */
    bool m_computeFromHistogramFlag;                  /*!< A boolean to specify if the datafile is histogram and that the sensitivity must be computed from it */
    vector<string> mp_pathToNormalizationFileName;    /*!< Vector of normalization file names */
    bool m_inverseDataFileOrderFlag;                  /*!< Flag to say that the provided normalization data file names order must be reversed */
    vDataFile*** m3p_NormalizationDataFile;           /*!< Pointer to the double array of vDataFile objects for normalization */
    bool m_oneNormalizationFileForAllBeds;            /*!< A boolean to specify if only one normalization file has been provided for all beds */
    bool m_oneNormalizationFileForAllFrames;          /*!< A boolean to specify if only one normalization file has been provided for all frames */
    string m_pathToAttenuationImage;                  /*!< String of the path to the attenuation images */
    string m_pathToMaskImg;                           /*!< String of the path to the mask image */
    bool m_mumapAttenuationFlag;                      /*!< Flag to say if an attenuation image has to been provided */
    bool m_forwardProjectAttenuation;                 /*!< Flag to say that we need to forward project the provided attenuation map (many conditions for that) */
    int m_nbAtnRespGateImages;                        /*!< Number of attenuation images corresponding to the respiratory gates */
    int m_nbAtnCardGateImages;                        /*!< Number of attenuation images corresponding to the cardiac gates */
    uint64_t* mp_lineCounter;                         /*!< Counters of the number of effectively projected lines, 1 per thread */
    bool m_flagGPU;                                   /*!< Do we use GPU or not (default=false) */ 
    int m_verbose;                                    /*!< Verbosity (default=-1)*/
    bool m_checked;                                   /*!< Boolean that says if the parameters were checked or not */
    bool m_initialized;                               /*!< Boolean that says if the projector was initialized or not */
};

#endif













