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
  \ingroup  projector
  \brief    Declaration of class oProjectorManager
*/

#ifndef OPROJECTORMANAGER_HH
#define OPROJECTORMANAGER_HH 1

#include "gVariables.hh"
#include "vScanner.hh"
#include "vProjector.hh"
#include "oSystemMatrix.hh"
#include "oProjectionLine.hh"
#include "vDataFile.hh"
#include "vEvent.hh"

/*!
  \class   oProjectorManager
  \brief   This class is designed to manage the projection part of the reconstruction
  \details As each manager class, it is created in the main program, all parameters are
           then set, checked, and the manager is initialized. The manager is then used by the
           algorithm itself, where the function ComputeProjectionLine() is called to compute
           the system matrix elements for the provided event, which are stored in a
           oProjectionLine. For multi-threading implementation, each thread have its own
           oProjectionLine. The manager will make use of a vProjector or a oSystemMatrix
           to compute the system matrix elements for forward and backward projections.
*/
class oProjectorManager
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public oProjectorManager::oProjectorManager()
      \brief   The constructor of oProjectorManager
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    oProjectorManager();
    /*!
      \fn      public oProjectorManager::~oProjectorManager()
      \brief   The destructor of oProjectorManager
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built by this class.
    */
    ~oProjectorManager();


  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      public int oProjectorManager::CheckParameters()
      \brief   A function used to check the parameters settings
      \details This function does not take any parameter and is used to check that all mandatory
               members were correctly parameterized.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckParameters();
    /*!
      \fn      public int oProjectorManager::CheckSPECTAttenuationCompatibility()
      \param   const string& a_pathToAttenuationImage
      \brief   A function used to check specific compatibility with SPECT and attenuation correction
      \details This function checks that the projector is compatible with SPECT and attenuation correction
               and that the computation strategy is not image-based; in case the provided path to an attenuation
               image is not empty.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckSPECTAttenuationCompatibility(const string& a_pathToAttenuationImage);
    /*!
      \fn      public int oProjectorManager::Initialize()
      \brief   A function used to initialize the manager and the projectors or system matrices it manages
      \details This function does not take any parameter and is used to initialize everything that
               should be initialized. In a few words, it parses the options, then creates and
               initializes the projectors or system matrices based on the provided options.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    int Initialize();
    /*!
      \fn      public void oProjectorManager::ApplyBedOffset()
      \param   int a_bed
      \brief   Compute the bed offset from the provided bed index and apply it to all projection lines
    */
    void ApplyBedOffset(int a_bed);
    /*!
      \fn      public void oProjectorManager::SetSensitivityModeOn()
      \brief   Say that the projector will be used to compute the global sensitivity
      \details It change the TOF application, as well as POI, and set the sensitivity mode on for the vProjector
    */
    void SetSensitivityModeOn();
    /*!
      \fn      public void oProjectorManager::SetSensitivityModeOff()
      \brief   Say that the projector will no longer be used to compute the global sensitivity
      \details It restores the TOF application, as well as POI, and set the sensitivity mode off for the vProjector
    */
    void SetSensitivityModeOff();
    /*!
      \fn      public oProjectionLine* oProjectorManager::ComputeProjectionLine()
      \param   vEvent* ap_Event
      \param   int a_th
      \brief   This function is used to compute system matrix elements from the associated projector or pre-computed system matrix
               given a vEvent, a bed position and a thread index. As this function returns a pointer to a oProjectionLine (filled
               from its own oProjectionLine buffers), if a problem occurs, it returns the NULL pointer.
      \return  A pointer to the oProjectionLine containing the system matrix elements.
    */
    oProjectionLine* ComputeProjectionLine(vEvent* ap_Event, int a_th);


  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      public inline void oProjectorManager::SetVerbose()
      \param   int a_verboseLevel
      \brief   Set the verbose level
    */
    inline void SetVerbose(int a_verboseLevel)
           {m_verbose = a_verboseLevel;}
    /*!
      \fn      public inline void oProjectorManager::SetScanner()
      \param   vScanner* ap_Scanner
      \brief   Set the scanner in use
    */
    inline void SetScanner(vScanner* ap_Scanner)
           {mp_Scanner = ap_Scanner;}
    /*!
      \fn      public inline void oProjectorManager::SetImageDimensionsAndQuantification()
      \param   oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification
      \brief   Set the image dimensions in use
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
           {mp_ImageDimensionsAndQuantification = ap_ImageDimensionsAndQuantification;}
    /*!
      \fn      public inline void oProjectorManager::SetDataFile()
      \param   vDataFile* ap_DataFile
      \brief   Set a data file in use to later recover some information from it
    */
    inline void SetDataFile(vDataFile* ap_DataFile)
           {mp_DataFile = ap_DataFile;}
    /*!
      \fn      public inline void oProjectorManager::SetComputationStrategy()
      \param   int a_computationStrategy
      \brief   Set the computation strategy for the system matrix elements storage
    */
    inline void SetComputationStrategy(int a_computationStrategy)
           {m_computationStrategy = a_computationStrategy;}
    /*!
      \fn      public inline void oProjectorManager::SetOptionsForward()
      \param   const string& a_optionsForward
      \brief   Set the forward projection options contained in the provided string
    */
    inline void SetOptionsForward(const string& a_optionsForward)
           {m_optionsForward = a_optionsForward;}
    /*!
      \fn      public inline void oProjectorManager::SetOptionsBackward()
      \param   const string& a_optionsBackward
      \brief   Set the backward projection options contained in the provided string
    */
    inline void SetOptionsBackward(const string& a_optionsBackward)
           {m_optionsBackward = a_optionsBackward;}
    /*!
      \fn      public inline void oProjectorManager::SetOptionsCommon()
      \param   const string& a_optionsCommon
      \brief   Set the common projection options contained in the provided string
    */
    inline void SetOptionsCommon(const string& a_optionsCommon)
           {m_optionsCommon = a_optionsCommon;}
    /*!
      \fn      public inline int oProjectorManager::GetNbTOFBins()
      \brief   Get the number of TOF bins associated to the projector
      \return  The number of TOF bins m_nbTOFBins
    */
    inline int GetNbTOFBins()
           {return m_nbTOFBins;}
    /*!
      \fn      public inline int oProjectorManager::GetComputationStrategy()
      \brief   Get the computation strategy for the storage in the projection line
      \return  The computation strategy m_computationStrategy
    */
    inline int GetComputationStrategy()
           {return m_computationStrategy;}
    /*!
      \fn      public bool oProjectorManager::IsForwardOperatorCompatibleWithSPECTAttenuationCorrection()
      \return  Return the compatibility with SPECT attenuation correction of the forward operator
    */
    bool IsForwardOperatorCompatibleWithSPECTAttenuationCorrection();
    /*!
      \fn      public bool oProjectorManager::IsBackwardOperatorCompatibleWithSPECTAttenuationCorrection()
      \return  Return the compatibility with SPECT attenuation correction of the backward operator
    */
    bool IsBackwardOperatorCompatibleWithSPECTAttenuationCorrection();
    /*!
      \fn      public inline void oProjectorManager::ProcessAndSetMask()
      \param   FLTNB* ap_maskImage
      \brief   Process and set the provided mask image for projector masking
    */
    int ProcessAndSetMask(FLTNB* ap_maskImage);

  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      private int oProjectorManager::ParseOptionsAndInitializeProjectors()
      \brief   Parse forward and backward projection options contained in the previously provided
               strings. This function is called inside the Initialize() function.
      \details Parse forward and backward projection options contained in the m_optionsForward and
               m_optionsBackward strings. Specific pure virtual functions of the vProjector are used to read
               parameters and initialize them. If a oSystemMatrix is used, specific functions are also called
               to read and initialize it. This function is called inside the Initialize() function. The syntax
               for the declaration of the projector and associated options is described inside the main program.
      \return  An integer reflecting the parsing and initialization status; 0 if no problem,
               another value otherwise.
    */
    int ParseOptionsAndInitializeProjectors();


  // -------------------------------------------------------------------
  // Data members
  private:
    // Scanner and image dimensions
    vScanner* mp_Scanner;                   /*!< Pointer to the vScanner object in use */
    oImageDimensionsAndQuantification* 
      mp_ImageDimensionsAndQuantification;  /*!< Pointer to the oImageDimensionsAndQuantification object in use */
    // DataFile
    vDataFile* mp_DataFile;                 /*!< Pointer to a vDataFile object in use */
    // TOF and POI options
    int m_applyTOF;                        /*!< Integer tagging the use of type of TOF use */
    bool m_applyPOI;                        /*!< Boolean that says if we apply POI or not */
    int m_nbTOFBins;                        /*!< The number of TOF bins in use */
    // Computation strategy for projection lines
    int m_computationStrategy;              /*!< The integer describing the computation strategy for system matrix elements storage */
    // Forward and backward options for projectors
    string m_optionsForward;                /*!< The string containing options for the forward projections */
    string m_optionsBackward;               /*!< The string containing options for the backward projections */
    // Common options for projectors
    string m_optionsCommon;                 /*!< The string containing common options for the projectors */
    // Forward and backward projectors
    string m_forwardProjectorName;          /*!< The name of the forward projector provided in the options */
    string m_backwardProjectorName;         /*!< The name of the backward projector provided in the options */
    oSystemMatrix* mp_SystemMatrixForward;  /*!< The pointer to the oSystemMatrix used for forward projections */
    oSystemMatrix* mp_SystemMatrixBackward; /*!< The pointer to the oSystemMatrix used for backward projections */
    vProjector* mp_ProjectorForward;        /*!< The pointer to the vProjector used for forward projections */
    vProjector* mp_ProjectorBackward;       /*!< The pointer to the vProjector used for backward projections */
    bool m_useSystemMatrixForward;          /*!< Boolean that says if a oSystemMatrix is used for forward projections */
    bool m_useSystemMatrixBackward;         /*!< Boolean that says if a oSystemMatrix is used for backward projections */
    bool m_useProjectorForward;             /*!< Boolean that says if a vProjector is used for forward projections */
    bool m_useProjectorBackward;            /*!< Boolean that says if a vProjector is used for backward projections */
    bool m_useMatchedProjectors;            /*!< Boolean that says if matched projectors are used */
    // Forward and backward projection lines (as many as threads)
    oProjectionLine** m2p_ProjectionLines;  /*!< The table of pointers to the oProjectionLines, one per thread */
    // Verbose level
    int m_verbose;                          /*!< The verbose level */
    // Has been checked ?
    bool m_checked;                         /*!< Boolean that says if the parameters were checked or not */
    // Has been initialized ?
    bool m_initialized;                     /*!< Boolean that says if the manager was initialized or not */
    bool m_applyMask;                       /*!< Apply the voxel mask to projectors */
    bool* mp_mask;                          /*!< Mask for voxels: only true voxels will be taken into account for projector coefficients computation (currently 3D - XYZ)*/
};

#endif
