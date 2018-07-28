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
  \ingroup  optimizer
  \brief    Declaration of class oOptimizerManager
*/

#ifndef OOPTIMIZERMANAGER_HH
#define OOPTIMIZERMANAGER_HH 1

#include "gVariables.hh"
#include "vOptimizer.hh"
#include "vPenalty.hh"
#include "oImageSpace.hh"
#include "vDataFile.hh"

/*!
  \class   oOptimizerManager
  \brief   This class is designed to manage the optimization part of an iterative reconstruction
  \details As each manager class, it is created in the main program, all parameters are then set,
           checked, and the manager is initialized. The manager is then used by the algorithm
           itself, where the function DataUpdateStep() function is called for each event to apply
           the forward projection, perform operations in the data space and apply the back-proj,
           based on a oProjectionLine, a vEvent and a vOptimizer. The ImageUpdateStep() function
           is called after the loop on events to apply the update operations in the image space,
           based on the back-projected correction images.
*/
class oOptimizerManager
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public oOptimizerManager::oOptimizerManager()
      \brief   The constructor of oOptimizerManager
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    oOptimizerManager();
    /*!
      \fn      public oOptimizerManager::~oOptimizerManager()
      \brief   The destructor of oOptimizerManager
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built by this class.
    */
    ~oOptimizerManager();


  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      public int oOptimizerManager::CheckParameters()
      \brief   A function used to check the parameters settings
      \details This function does not take any parameter and is used to check that all mandatory
               members were correctly parameterized.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckParameters();
    /*!
      \fn      public int oOptimizerManager::Initialize()
      \brief   A function used to initialize the manager and the optimizer it manages
      \details This function does not take any parameter and is used to initialize everything that
               should be initialized. In a few words, it parses the options, then creates and
               initializes the optimizer based on the provided options.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    int Initialize();
    /*!
      \fn      public int oOptimizerManager::PreDataUpdateStep()
      \param   int a_iteration
      \param   int a_nbIterations
      \param   int a_subset
      \param   int* ap_nbSubsets
      \brief   A function that simply calls the eponym function from the vOptimizer
      \return  An integer reflecting the execution status; 0 if no problem, another value otherwise.
    */
    int PreDataUpdateStep(int a_iteration, int a_nbIterations, int a_subset, int* ap_nbSubsets);
    /*!
      \fn      public int oOptimizerManager::PostDataUpdateStep()
      \param   int a_iteration
      \param   int a_nbIterations
      \param   int a_subset
      \param   int* ap_nbSubsets
      \brief   A function that simply calls the eponym function from the vOptimizer
      \return  An integer reflecting the execution status; 0 if no problem, another value otherwise.
    */
    int PostDataUpdateStep(int a_iteration, int a_nbIterations, int a_subset, int* ap_nbSubsets);
    /*!
      \fn      public int oOptimizerManager::DataUpdateStep()
      \param   oProjectionLine* ap_Line
      \param   vEvent* ap_Event
      \param   int a_bed
      \param   int a_timeFrame
      \param   int a_respGate
      \param   int a_cardGate
      \param   int a_iteration
      \param   int a_thread
      \brief   A function dedicated to the update step in the data space (for each event inside the loop)
      \details This function will call many functions of the vOptimizer that split up the update step in the data space
               in smaller steps: 1. forward projection, 2. optional step, 3. backproj sensitivity for histogram mode,
               4. optional step, 5. compute corrections, 6. optional step, 7. backproj corrections, 8. compute FOMs.
               This function is called by the iterative algorithm, inside the loop on all events.
      \return  An integer reflecting the execution status; 0 if no problem, another value otherwise.
    */
    int DataUpdateStep( oProjectionLine* ap_Line, vEvent* ap_Event,
                        int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                        int a_iteration, int a_thread );
    /*!
      \fn      public int oOptimizerManager::ImageUpdateStep()
      \param   int a_iteration
      \param   int a_nbSubsets
      \brief   A function dedicated to the update step in the image space (performed after the loop on events)
      \details This function will update the visited voxels first (see in vOptimizer for details), manage the call
               for penalty computation and call the image update step function specific to the optimizer.
      \return  An integer reflecting the execution status; 0 if no problem, another value otherwise.
    */
    int ImageUpdateStep( int a_iteration, int a_nbSubsets );


  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      private int oOptimizerManager::ParseOptionsAndInitializeProjectors()
      \brief   Parse optimizer and penalty options contained in the previously provided
               strings. This function is called inside the Initialize() function.
      \details Parse optimizer and penalty options contained in the m_optionsOptimizer and
               m_optionsPenalty strings. Specific pure virtual functions of the vOptimizer
               and vPenalty are used to read parameters and initialize them. This function
               is called inside the Initialize() function. The syntax for the declaration of
               the optimizer, penalty and associated options is described inside the main program.
      \return  An integer reflecting the parsing and initialization status; 0 if no problem,
               another value otherwise.
    */
    int ParseOptionsAndInitializeOptimizerAndPenalty();


  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      public inline void oOptimizerManager::SetVerbose()
      \param   int a_verboseLevel
      \brief   Set the verbose level
    */
    inline void SetVerbose(int a_verboseLevel)
           {m_verbose = a_verboseLevel;}
    /*!
      \fn      public inline void oOptimizerManager::SetOptionsOptimizer()
      \param   const string& a_optionsOptimizer
      \brief   Set the optimizer projection options contained in the provided string
    */
    inline void SetOptionsOptimizer(const string& a_optionsOptimizer)
           {m_optionsOptimizer = a_optionsOptimizer;}
    /*!
      \fn      public inline void oOptimizerManager::SetOptionsPenalty()
      \param   const string& a_optionsPenalty
      \brief   Set the penalty projection options contained in the provided string
    */
    inline void SetOptionsPenalty(const string& a_optionsPenalty)
           {m_optionsPenalty = a_optionsPenalty;}
    /*!
      \fn      public inline void oOptimizerManager::SetImageDimensionsAndQuantification()
      \param   oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification
      \brief   Set the image dimensions in use
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
           {mp_ImageDimensionsAndQuantification = ap_ImageDimensionsAndQuantification;}
    /*!
      \fn      public inline void oOptimizerManager::SetImageSpace()
      \param   oImageSpace* ap_ImageSpace
      \brief   Set the image space in use
    */
    inline void SetImageSpace(oImageSpace* ap_ImageSpace)
           {mp_ImageSpace = ap_ImageSpace;}
    /*!
      \fn      public inline void oOptimizerManager::SetNbTOFBins()
      \param   int a_nbTOFBins
      \brief   Set the number of TOF bins in use
    */
    inline void SetNbTOFBins(int a_nbTOFBins)
           {m_nbTOFBins = a_nbTOFBins;}
    /*!
      \fn      public inline void oOptimizerManager::SetDataMode()
      \param   int a_dataMode
      \brief   Set the mode of the data (histogram, list-mode)
    */
    inline void SetDataMode(int a_dataMode)
           {m_dataMode = a_dataMode;}
    /*!
      \fn      public inline void oOptimizerManager::SetDataType()
      \param   int a_dataType
      \brief   Set the type of the data (pet, spect, etc)
    */
    inline void SetDataType(int a_dataType)
           {m_dataType = a_dataType;}
    /*!
      \fn      public inline void oOptimizerManager::SetDataSpec()
      \param   int a_dataSpec
      \brief   Set the specificity of the data: EMISSION or TRANSMISSION
    */
    inline void SetDataSpec(int a_dataSpec)
           {m_dataSpec = a_dataSpec;}
    /*!
      \fn      public inline void oOptimizerManager::SetOptimizerFOMFlag()
      \param   bool a_optimizerFOMFlag
      \brief   Set the optimizer FOM flag that specifies if some figures-of-merit (FOM) will be computed in the data space
    */
    inline void SetOptimizerFOMFlag(bool a_optimizerFOMFlag)
           {m_optimizerFOMFlag = a_optimizerFOMFlag;}
    /*!
      \fn      public inline void oOptimizerManager::SetOptimizerImageStatFlag()
      \param   bool a_optimizerImageStatFlag
      \brief   Set the optimizer image stat flag that specifies if some basic statistics about image update is computed
    */
    inline void SetOptimizerImageStatFlag(bool a_optimizerImageStatFlag)
           {m_optimizerImageStatFlag = a_optimizerImageStatFlag;}
    /*!
      \fn      public inline int oOptimizerManager::GetNbBackwardImages()
      \brief   Return the number of backward images used by the optimizer, explaining why the eponym function of vOptimizer is called
    */
    inline int GetNbBackwardImages()
           {return mp_Optimizer->GetNbBackwardImages();}
    /*!
      \fn      public inline FLTNB oOptimizerManager::GetInitialValue()
      \brief   Return the initial image value used by the optimizer, explaining why the eponym function of vOptimizer is called
    */
    inline FLTNB GetInitialValue()
           {return mp_Optimizer->GetInitialValue();}


  // -------------------------------------------------------------------
  // Data members
  private:
    oImageDimensionsAndQuantification* 
      mp_ImageDimensionsAndQuantification; /*!< Pointer to the oImageDimensionsAndQuantification object in use */
    oImageSpace* mp_ImageSpace;                 /*!< Pointer to the oImageSpace object in use */
    int m_dataMode;                        /*!< Flag indicating if the data is List (=0) or Histogram (=1) mode*/
    int m_dataType;                        /*!< Flag indicating if the data is PET (=0),SPECT (=1) or TRANSMISSION type (=2) */
    int m_dataSpec;                        /*!< Flag indicating if the data is EMISSION or TRANSMISSION */
    int m_nbTOFBins;                       /*!< The number of TOF bins in use */
    string m_optionsOptimizer;             /*!< The string containing options for the optimizer projections */
    string m_optionsPenalty;               /*!< The string containing options for the penalty projections */
    bool m_optimizerFOMFlag;               /*!< Flag that says if some figures-of-merit will be computed in the data space */
    bool m_optimizerImageStatFlag;         /*!< Flag that says if some basic statistics about the image update will be computed */
    vOptimizer* mp_Optimizer;              /*!< The actual optimizer in use */
    vPenalty* mp_Penalty;                  /*!< The actual penalty in use (optional) */
    int m_verbose;                         /*!< The verbose level */
};

#endif
