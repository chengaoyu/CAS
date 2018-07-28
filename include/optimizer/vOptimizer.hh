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
  \brief    Declaration of class vOptimizer
*/

#ifndef VOPTIMIZER_HH
#define VOPTIMIZER_HH 1

#include "gVariables.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "oImageSpace.hh"
#include "vDataFile.hh"
#include "oProjectionLine.hh"

/*!
  \class   vOptimizer
  \brief   This class is designed to generically described any iterative optimizer
  \details This class is an abstract one, in the sense that it cannot be used on its own
           because several pure virtual functions belong to it. Its children are
           implementations of actual optimizers. Everywhere in the code, this parent class
           should be used instead of any of its children. The optimizer is created by the
           oOptimizerManager is used during the iterative process of the algorithm through the
           call of multiple functions. All functions prefixed by DataStepX are called within the
           DataUpdateStep() function of the oOptimizerManager. The other main function is the
           ImageUpdateStep() function called by the eponym function from the manager. Anyway,
           all these aforementioned functions are virtual, but are designed to be flexible enaugh
           for any kind of optimizer. Finally, a particular optimizer (a child of this vOptimizer)
           will be characterized by the implementation of some pure virtual functions, some for
           options and initialization management, and three particular ones that really describe
           what this particular optimizer does: SensitivitySpecificOperations(),
           DataSpaceSpecificOperations() and ImageSpaceSpecificOperations(). Read the description
           of these functions below to get more details.
*/
class vOptimizer
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public vOptimizer::vOptimizer()
      \brief   The constructor of vOptimizer
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    vOptimizer();
    /*!
      \fn      virtual public vOptimizer::~vOptimizer()
      \brief   The destructor of vOptimizer
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built by this class.
               It is virtual, so that it is automatically called when a child object is deleted.
    */
    virtual ~vOptimizer();


  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      public int vOptimizer::ShowHelp()
      \brief   A function used to show help about the optimizer
      \details This function simply calls the ShowHelpSpecific() function implemented by children.
    */
    void ShowHelp();
    /*!
      \fn      public int vOptimizer::CheckParameters()
      \brief   A public function used to check the parameters settings
      \details This function does not take any parameter and is used to check that all mandatory
               members were correctly parameterized. At the end, it calls the pure virtual
               CheckSpecificParameters() function implemented by children.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckParameters();
    /*!
      \fn      public int vOptimizer::Initialize()
      \brief   A public function used to initialize the optimizer
      \details This function does not take any parameter and is used to initialize everything that
               should be initialized. At the end, it calls the pure virtual InitializeSpecific()
               function implemented by children.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */
    int Initialize();
    /*!
      \fn      public int vOptimizer::UpdateVisitedVoxels()
      \brief   A public function used to update the 'visited' voxels after each subset
      \details This function is called at the end of each subset. Based on the sensitivity image, it labels
               which voxel has been visited or not during each subset. At the end of the iteration, all
               voxels that have not been visited in any subset are set to 0. There is a dedicated matrix in
               the oImageSpace to keep track of visited voxels. Practically, this avoids having those never
               visited voxels trapped in their initial value if the latter is different from 0.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    int UpdateVisitedVoxels();
    /*!
      \fn      public int vOptimizer::PreDataUpdateStep()
      \param   int a_iteration
      \param   int a_nbIterations
      \param   int a_subset
      \param   int* ap_nbSubsets
      \brief   A public function used to do stuff that need to be done at the beginning of a subset (before
               the data update step; i.e. the loop over all events)
      \details It does some reseting for the FOM computation. At the end, it calls the PreDataUpdateSpecificStep()
               function, which does nothing but being virtual, so it can be overloaded to perform specific stuff
               at this step by specific optimizers that would need to do so.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    int PreDataUpdateStep(int a_iteration, int a_nbIterations, int a_subset, int* ap_nbSubsets);
    /*!
      \fn      public int vOptimizer::PostDataUpdateStep()
      \param   int a_iteration
      \param   int a_nbIterations
      \param   int a_subset
      \param   int* ap_nbSubsets
      \brief   A public function used to do stuff that need to be done at the beginning of a subset (before
               the data update step; i.e. the loop over all events)
      \details It does some reseting for the FOM computation. At the end, it calls the PostDataUpdateSpecificStep()
               function, which does nothing but being virtual, so it can be overloaded to perform specific stuff
               at this step by specific optimizers that would need to do so.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    int PostDataUpdateStep(int a_iteration, int a_nbIterations, int a_subset, int* ap_nbSubsets);


  // -------------------------------------------------------------------
  // Virtual public member functions that may be overloaded by specific child optimizers
  public:
    /*!
      \fn      public virtual int vOptimizer::DataStep1ForwardProjectModel()
      \param   oProjectionLine* ap_Line
      \param   vEvent* ap_Event
      \param   int a_bed
      \param   int a_timeFrame
      \param   int a_respGate
      \param   int a_cardGate
      \param   int a_thread
      \brief   A public function used to compute the model: forward projection of the provided event
      \details Inside the DataUpdateStep() of the oOptimizerManager, this is the first function called. This function performs
               the forward projection of the provided event, following the provided projection line. It takes all dynamic
               dimensions into account, including the intrinsic basis functions. Note that this function uses the ForwardProject()
               function that automatically deals with attenuation for SPECT and includes all multiplicative terms theoretically
               included in the system matrix. Then, all additive corrections are added to the forward projection. So the unit
               of the resulting model is the same as the data. The result of this forward projection is put in the
               m2p_forwardValues table (for the given thread). The other dimension of this table is for the TOF bins. 
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int DataStep1ForwardProjectModel( oProjectionLine* ap_Line, vEvent* ap_Event,
                                              int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                              int a_thread );
    /*!
      \fn      public virtual int vOptimizer::DataStep2Optional()
      \param   oProjectionLine* ap_Line
      \param   vEvent* ap_Event
      \param   int a_bed
      \param   int a_timeFrame
      \param   int a_respGate
      \param   int a_cardGate
      \param   int a_iteration
      \param   int a_thread
      \brief   A public function which does nothing but being virtual.
      \details Inside the DataUpdateStep() of the oOptimizerManager, this is the second function called. This function can be
               overloaded by specific optimizers if needed.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int DataStep2Optional( oProjectionLine* ap_Line, vEvent* ap_Event,
                                   int a_bed, int a_timeFrame, int a_respGate, int a_cardGate, int a_iteration,
                                   int a_thread );
    /*!
      \fn      public virtual int vOptimizer::DataStep3BackwardProjectSensitivity()
      \param   oProjectionLine* ap_Line
      \param   vEvent* ap_Event
      \param   int a_bed
      \param   int a_timeFrame
      \param   int a_respGate
      \param   int a_cardGate
      \param   int a_thread
      \brief   A public function used to back-project the sensitivity terms for the provided event
      \details Inside the DataUpdateStep() of the oOptimizerManager, this is the third function called. This function is only
               called when using histogram data. Before performing the back-projection into the sensitivity image, it calls
               the pure virtual SensitivitySpecificOperations() function whose role is simply to provide the weight associated
               to the projection line.  Note that this function uses the BackwardProject() function that automatically deals
               with attenuation for SPECT and includes all multiplicative terms theoretically included in the system matrix.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int DataStep3BackwardProjectSensitivity( oProjectionLine* ap_Line, vEvent* ap_Event,
                                                     int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                                     int a_thread );
    /*!
      \fn      public virtual int vOptimizer::DataStep4Optional()
      \param   oProjectionLine* ap_Line
      \param   vEvent* ap_Event
      \param   int a_bed
      \param   int a_timeFrame
      \param   int a_respGate
      \param   int a_cardGate
      \param   int a_iteration
      \param   int a_thread
      \brief   A public function which does nothing but being virtual.
      \details Inside the DataUpdateStep() of the oOptimizerManager, this is the fourth function called. This function can be
               overloaded by specific optimizers if needed.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int DataStep4Optional( oProjectionLine* ap_Line, vEvent* ap_Event,
                                   int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                   int a_iteration, int a_thread );
    /*!
      \fn      public virtual int vOptimizer::DataStep5ComputeCorrections()
      \param   oProjectionLine* ap_Line
      \param   vEvent* ap_Event
      \param   int a_bed
      \param   int a_timeFrame
      \param   int a_respGate
      \param   int a_cardGate
      \param   int a_thread
      \brief   A public function used to compute the correction terms in the data space, for the provided event
      \details Inside the DataUpdateStep() of the oOptimizerManager, this is the fifth function called. Its role is to compute
               the correction terms in the data space, based on the forward model and the data. In order to be specific to each
               optimizer, it calls the pure virtual function DataSpaceSpecificOperations(), where the computation is done. The
               correction terms are put in the m3p_backwardValues (a dimension for threads, one for the number of backward images
               and the last for TOF bins).
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int DataStep5ComputeCorrections( oProjectionLine* ap_Line, vEvent* ap_Event,
                                             int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                             int a_thread );
    /*!
      \fn      public virtual int vOptimizer::DataStep6Optional()
      \param   oProjectionLine* ap_Line
      \param   vEvent* ap_Event
      \param   int a_bed
      \param   int a_timeFrame
      \param   int a_respGate
      \param   int a_cardGate
      \param   int a_iteration
      \param   int a_thread
      \brief   A public function which does nothing but being virtual.
      \details Inside the DataUpdateStep() of the oOptimizerManager, this is the sixth function called. This function can be
               overloaded by specific optimizers if needed.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int DataStep6Optional( oProjectionLine* ap_Line, vEvent* ap_Event,
                                   int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                   int a_iteration, int a_thread );
    /*!
      \fn      public virtual int vOptimizer::DataStep7BackwardProjectCorrections()
      \param   oProjectionLine* ap_Line
      \param   vEvent* ap_Event
      \param   int a_bed
      \param   int a_timeFrame
      \param   int a_respGate
      \param   int a_cardGate
      \param   int a_thread
      \brief   A public function used to back-project the correction terms into the backward correction image
      \details Inside the DataUpdateStep() of the oOptimizerManager, this is the seventh function called. This function performs
               the backward projection of the correction terms held into the m3p_backwardValues, for the provided event and
               following the provided projection line. It takes all dynamic dimensions into account, including the intrinsic
               basis functions. Note that this function uses the BackwardProject() function that automatically deals with
               attenuation for SPECT and automatically includes all multiplicative terms theoretically included in the system
               matrix.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int DataStep7BackwardProjectCorrections( oProjectionLine* ap_Line, vEvent* ap_Event,
                                                     int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                                     int a_thread );
    /*!
      \fn      public virtual int vOptimizer::DataStep8ComputeFOM()
      \param   oProjectionLine* ap_Line
      \param   vEvent* ap_Event
      \param   int a_timeFrame
      \param   int a_respGate
      \param   int a_cardGate
      \param   int a_thread
      \brief   A public function used to update the computation of figures-of-merit in the data space
      \details Inside the DataUpdateStep() of the oOptimizerManager, this is the eighth function called. This function updates
               the computation of figures-of-merit in the data space for the provided event. It uses the m4p_FOMXXX structures.
               The computation is done only if the FOM flag is set. For the moment, log-likelihood and RMSE are computed.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int DataStep8ComputeFOM( oProjectionLine* ap_Line, vEvent* ap_Event,
                                     int a_timeFrame, int a_respGate, int a_cardGate,
                                     int a_thread );
    /*!
      \fn      public virtual int vOptimizer::ImageUpdateStep()
      \param   int a_iteration
      \param   int a_nbSubsets
      \brief   A public function used to perform the image update step of the optimizer
      \details This function is called by the eponym function from the oOptimizerManager. It will manage the dynamic loops,
               compute the sensitivity using the private ComputeSensitivity() function, and update each voxel according to
               the specific optimizer by calling the pure virtual ImageSpaceSpecificOperations() function.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int ImageUpdateStep( int a_iteration, int a_nbSubsets );


  // -----------------------------------------------------------------------------------------
  // Virtual private member functions that may be overloaded by specific child optimizers
  private:
    /*!
      \fn      private virtual int vOptimizer::PreDataUpdateSpecificStep()
      \param   int a_iteration
      \param   int a_nbIterations
      \param   int a_subset
      \param   int* ap_nbSubsets
      \brief   A private function used to perform any step required by the child optimizer, before the loop on event inside the
               subset loop.
      \details The vOptimizer implementation just does nothing. It is called by the PreDataUpdateStep() function. It is virtual
               so it can be overloaded by the child optimizer if needed.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int PreDataUpdateSpecificStep(int a_iteration, int a_nbIterations, int a_subset, int* ap_nbSubsets);
    /*!
      \fn      private virtual int vOptimizer::PostDataUpdateSpecificStep()
      \param   int a_iteration
      \param   int a_nbIterations
      \param   int a_subset
      \param   int* ap_nbSubsets
      \brief   A private function used to perform any step required by the child optimizer, after the loop on event inside the
               subset loop.
      \details The vOptimizer implementation just does nothing. It is called by the PostDataUpdateStep() function. It is virtual
               so it can be overloaded by the child optimizer if needed.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int PostDataUpdateSpecificStep(int a_iteration, int a_nbIterations, int a_subset, int* ap_nbSubsets);

  // -----------------------------------------------------------------------------------------
  // Pure virtual public member functions that need to be implemented by child optimizers
  public:
    /*!
      \fn      public virtual int vOptimizer::ReadConfigurationFile() = 0
      \param   const string& a_configurationFile
      \brief   A function used to read options from a configuration file
      \details This function implements the reading of all options associated to a child optimizer, from
               a configuration file. It is pure virtual so is implemented by children. It checks the reading
               status but not the options values that will be checked by the CheckSpecificParameters()
               function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    virtual int ReadConfigurationFile( const string& a_configurationFile ) = 0;
    /*!
      \fn      public virtual int vOptimizer::ReadOptionsList() = 0
      \param   const string& a_optionsList
      \brief   A function used to read options from a list of options
      \details This function implements the reading of all options associated to a child optimizer, from
               a list of options. It is pure virtual so is implemented by children. It checks the reading
               status but not the options values that will be checked by the CheckSpecificParameters()
               function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    virtual int ReadOptionsList( const string& a_optionsList ) = 0;


  // -----------------------------------------------------------------------------------------
  // Pure virtual private member functions that need to be implemented by child optimizers
  private:
    /*!
      \fn      private virtual void vOptimizer::ShowHelpSpecific() = 0
      \brief   A function used to show help about the child module
      \details This function must describe what the optimizer does and how to use it. It describes in
               details the different parameters of the optimizer, and how to set them through the use
               of a configuration file or a list of options. It is pure virtual so is implemented by
               children. It is private because called by the public ShowHelp() function.
    */
    virtual void ShowHelpSpecific() = 0;
    /*!
      \fn      private virtual int vOptimizer::CheckSpecificParameters() = 0
      \brief   A private function used to check the parameters settings specific to the child optimizer
      \details This function is used to check that all parameters specific to the optimizer are correctly set
               within allowed values. It is called by the CheckParameters() function. It is pure virtual
               so is implemented by children.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    virtual int CheckSpecificParameters() = 0;
    /*!
      \fn      private virtual int vOptimizer::InitializeSpecific() = 0
      \brief   A private function used to initialize everything specific to the child optimizer
      \details This function is used to initialize everything specific to the optimizer that should be
               initialized. It is called by the Initialize() function. It is pure virtual so is
               implemented by children.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    virtual int InitializeSpecific() = 0;
    /*!
      \fn      private virtual int vOptimizer::SensitivitySpecificOperations() = 0
      \param   FLTNB a_data
      \param   FLTNB a_forwardModel
      \param   FLTNB* ap_weight
      \param   FLTNB a_multiplicativeCorrections
      \param   FLTNB a_additiveCorrections
      \param   FLTNB a_blankValue
      \param   FLTNB a_quantificationFactor
      \param   oProjectionLine* ap_Line
      \brief   A private function used to compute the sensitivity weight associated to the provided data
      \details This function is pure virtual, so must be implemented by specific optimizers. The computed
               sensitivity term must be put at the ap_weight pointer. All potentially useful information
               is provided as parameters. Note that the multiplicate corrections parameter only includes
               such correction specific to the event; in other words, it does not include the quantification
               factor which is given separately.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
                                               FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                               FLTNB a_quantificationFactor, oProjectionLine* ap_Line ) = 0;
    /*!
      \fn      private virtual int vOptimizer::DataSpaceSpecificOperations() = 0
      \param   FLTNB a_data
      \param   FLTNB a_forwardModel
      \param   FLTNB* ap_backwardValues
      \param   FLTNB a_multiplicativeCorrections
      \param   FLTNB a_additiveCorrections
      \param   FLTNB a_blankValue
      \param   FLTNB a_quantificationFactor
      \param   oProjectionLine* ap_Line
      \brief   A private function used to compute the correction term in the data space from the provided data
      \details This function is pure virtual, so must be implemented by specific optimizers. The computed
               correction term must be put at the ap_backwardValues pointer. All potentially useful information
               is provided as parameters. This is the operation that the optimization algorithm performs in the
               data space, in the loop on events. Note that the multiplicate corrections parameter only includes
               such correction specific to the event; in other words, it does not include the quantification
               factor which is given separately.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                             FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                             FLTNB a_quantificationFactor, oProjectionLine* ap_Line ) = 0;
    /*!
      \fn      private virtual int vOptimizer::ImageSpaceSpecificOperations() = 0
      \param   FLTNB a_currentImageValue
      \param   FLTNB* ap_newImageValue
      \param   FLTNB a_sensitivity
      \param   FLTNB* ap_correctionValues
      \param   INTNB a_voxel
      \brief   A private function used to update the image value from the provided data
      \details This function is pure virtual, so must be implemented by specific optimizers. The old values are
               provided, as well as the location of the new that will be calculated. This is the operation that
               the optimization algorithm performs in the image update step, after the loop on events. Note that
               based on the number of backward images used by the specific optimizer, there can be multiple
               correction values in ap_correctionValues, explaining why it is provided as a pointer. Note also
               that this function is only called for non-zero sensitivity, so no need to check it.
               The design of this function will evolve in order to take into account optimizers with penalties and
               alternated optimizers like MLAA which estimate the attenuation at the same time.
      \return  An integer reflecting the process status; 0 if no problem, another value otherwise.
    */
    virtual int ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                              FLTNB a_sensitivity, FLTNB* ap_correctionValues, INTNB a_voxel ) = 0;


  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      public inline void vOptimizer::SetVerbose()
      \param   int a_verboseLevel
      \brief   Set the verbose level.
    */
    inline void SetVerbose(int a_verbose)
           {m_verbose = a_verbose;}
    /*!
      \fn      public inline void vOptimizer::SetImageDimensionsAndQuantification()
      \param   oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification
      \brief   Set the pointer to the image dimensions in use.
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
           {mp_ImageDimensionsAndQuantification = ap_ImageDimensionsAndQuantification;}
    /*!
      \fn      public inline void vOptimizer::SetImageSpace()
      \param   oImageSpace* ap_ImageSpace
      \brief   Set the pointer to the image space in use.
    */
    inline void SetImageSpace(oImageSpace* ap_ImageSpace)
           {mp_ImageSpace = ap_ImageSpace;}
    /*!
      \fn      public inline void vOptimizer::SetNbTOFBins()
      \param   int a_nbTOFBins
      \brief   Set the number of TOF bins in use.
    */
    inline void SetNbTOFBins(int a_nbTOFBins)
           {m_nbTOFBins = a_nbTOFBins;}
    /*!
      \fn      public inline void vOptimizer::SetDataMode()
      \param   int a_dataMode
      \brief   Set the data mode in use.
    */
    inline void SetDataMode(int a_dataMode)
           {m_dataMode = a_dataMode;}
    /*!
      \fn      public inline void vOptimizer::SetDataType()
      \param   int a_dataType
      \brief   Set the data type in use.
    */
    inline void SetDataType(int a_dataType)
           {m_dataType = a_dataType;}
    /*!
      \fn      public inline void vOptimizer::SetDataSpec()
      \param   int a_dataSpec
      \brief   Set the data physical specificity in use.
    */
    inline void SetDataSpec(int a_dataSpec)
           {m_dataSpec = a_dataSpec;}
    /*!
      \fn      public inline void vOptimizer::SetAttenuationImage()
      \param   FLTNB* ap_attenuationImage
      \param   int a_thread
      \brief   Set the attenuation image corresponding to the current thread and current event.
    */
    inline void SetAttenuationImage(FLTNB* ap_attenuationImage, int a_thread)
           {m2p_attenuationImage[a_thread] = ap_attenuationImage;}
    /*!
      \fn      public inline void vOptimizer::SetFOMFlag()
      \param   bool a_optimizerFOMFlag
      \brief   Set the FOM flag specifying if figures-of-merit will be computed or not.
    */
    inline void SetFOMFlag(bool a_optimizerFOMFlag)
           {m_optimizerFOMFlag = a_optimizerFOMFlag;}
    /*!
      \fn      public inline void vOptimizer::SetImageStatFlag()
      \param   bool a_optimizerImageStatFlag
      \brief   Set the image stat flag specifying if basic statistics about image udpate will be computed or not.
    */
    inline void SetImageStatFlag(bool a_optimizerImageStatFlag)
           {m_optimizerImageStatFlag = a_optimizerImageStatFlag;}
    /*!
      \fn      public inline int vOptimizer::GetNbBackwardImages()
      \brief   Get the number of backward images used by the specific optimizer
      \return  m_nbBackwardImages
    */
    inline int GetNbBackwardImages()
           {return m_nbBackwardImages;}
    /*!
      \fn      public inline FLTNB vOptimizer::GetInitialValue()
      \brief   Get the initial image value (for initialization)
      \return  m_initialValue
    */
    inline FLTNB GetInitialValue()
           {return m_initialValue;}
    /*
    inline bool GetAcceptPenalty()
           {return m_penaltyEnergyFunctionDerivativesOrder>0;}
    inline int GetPenaltyEnergyFunctionDerivativesOrder()
           {return m_penaltyEnergyFunctionDerivativesOrder;}
    */


  // -------------------------------------------------------------------
  // Private member functions
  protected:
    /*!
      \fn      protected FLTNB vOptimizer::ComputeSensitivity()
      \param   FLTNB**** a4p_sensitivityImage
      \param   int a_timeBasisFunction
      \param   int a_respBasisFunction
      \param   int a_cardBasisFunction
      \param   int a_voxel
      \param   int a_nbSubsets
      \brief   A function used to compute the sensitivity of a given voxel and a given set of dynamic basis functions
      \details This function computes the sensitivity for a given set of dynamic basis functions and for a given number
               of subsets. It loops over the dynamic frames and gates to compute the sensitivity based on basis functions
               coefficients. In list-mode, as the sensitivity is computed for the whole FOV, the value is divided by the
               current number of subsets of the current iteration, given as a parameter.
      \return  The computed sensitivity
    */
    FLTNB ComputeSensitivity( FLTNB**** a4p_sensitivityImage,
                              int a_timeBasisFunction, int a_respBasisFunction, int a_cardBasisFunction,
                              int a_voxel, int a_nbSubsets );
    /*!
      \fn      protected FLTNB vOptimizer::ForwardProject()
      \param   oProjectionLine* ap_Line
      \param   FLTNB* ap_image = NULL
      \brief   A function used to forward project the provided image (or 1 if NULL), based on the provided oProjectionLine
      \details Based on the data type, the function calls the projection function of the oProjectionLine class to forward
               project the image taking the multiplicative terms into account, and eventually dealing with the SPECT attenuation.
      \return  The computed forward projection
    */
    FLTNB ForwardProject( oProjectionLine* ap_Line, FLTNB* ap_image = NULL );
    /*!
      \fn      protected void vOptimizer::BackwardProject()
      \param   oProjectionLine* ap_Line
      \param   FLTNB* ap_image
      \param   FLTNB a_value
      \brief   A function used to backward project the provided value into the provided image, based on the provided oProjectionLine
      \details Based on the data type, the function calls the projection function of the oProjectionLine class to backward
               project the value into the image taking the multiplicative terms into account, and eventually dealing with the SPECT
               attenuation.
    */
    void BackwardProject( oProjectionLine* ap_Line, FLTNB* ap_image, FLTNB a_value );

  // -------------------------------------------------------------------
  // Data members
  protected:
    int m_verbose;                         /*!< The verbose level */
    int m_nbBackwardImages;                /*!< The number of backward images used by the specific optimizer */
    int m_nbTOFBins;                       /*!< The number of TOF bins in use */
    FLTNB** m2p_forwardValues;             /*!< Buffer for forward projected values, as many as threads and as many as TOF bins */
    FLTNB*** m3p_backwardValues;           /*!< Buffer for values to be back-projected, as many as threads, as many as TOF bins and as many as backward images */
    FLTNB m_initialValue;                  /*!< The initial value of the image for the specific optimizer */
    bool m_listmodeCompatibility;          /*!< A flag saying if the optimizer is compatible with list-mode data */
    bool m_histogramCompatibility;         /*!< A flag saying if the optimizer is compatible with histogram data */
    bool m_emissionCompatibility;          /*!< A flag saying if the optimizer is compatible with emission data */
    bool m_transmissionCompatibility;      /*!< A flag saying if the optimizer is compatible with transmission data */
    oImageDimensionsAndQuantification* 
      mp_ImageDimensionsAndQuantification; /*!< Pointer to the image dimensions and quantification object */
    oImageSpace* mp_ImageSpace;             /*!< Pointer to the image space object in use */
    int m_dataMode;                        /*!< The mode of the data (list-mode or histogram) */
    int m_dataType;                        /*!< The type of the data (PET, SPECT, ...) */
    int m_dataSpec;                        /*!< The physical specificity of the data (EMISSION or TRANSMISSION) */
    FLTNB** m2p_attenuationImage;          /*!< The attenuation image currently in use for the current event, one per thread (used for SPECT) */
    // Optimizer figures-of-merit computation
    bool m_optimizerFOMFlag;               /*!< A flag saying if figures-of-merit are computed in the data space */
    FLTNB**** m4p_FOMLogLikelihood;        /*!< The log-likelihood FOM, as many as frames, respiratory gates, cardiac gates and threads */
    FLTNB**** m4p_FOMRMSE;                 /*!< The RMSE FOM, as many as frames, respiratory gates, cardiac gates and threads */
    uint64_t**** m4p_FOMNbBins;            /*!< The number of bins contributing to the FOM computation, as many as frames, respiratory gates, cardiac gates and threads */
    HPFLTNB**** m4p_FOMNbData;              /*!< The number of data (counts) contributing to the FOM computation, as many as frames, respiratory gates, cardiac gates and threads */
    // Image update statistics
    bool m_optimizerImageStatFlag;         /*!< A flag saying if some basic statistics about image update are computed */
    INTNB* mp_imageStatNbVox;              /*!< The number of voxels contributing to the image statistics, one per thread to be thread safe */
    FLTNB* mp_imageStatMin;                /*!< The minimum image value, one per thread to be thread safe */
    FLTNB* mp_imageStatMax;                /*!< The maximum image value, one per thread to be thread safe */
    HPFLTNB* mp_imageStatMean;              /*!< The mean image value, one per thread to be thread safe */
    HPFLTNB* mp_imageStatVariance;          /*!< The image variance value, one per thread to be thread safe */
    HPFLTNB* mp_correctionStatMean;         /*!< The mean additive update correction value, one per thread to be thread safe */
    HPFLTNB* mp_correctionStatVariance;     /*!< The variance of additive update correction values, one per thread to be thread safe */
    //int m_penaltyEnergyFunctionDerivativesOrder; /*!< The penalty derivatives order accepted by the specific optimizer (0 if none) */
};


// ----------------------------------------------------------------------
// Part of code that manages the auto declaration of children classes
// ----------------------------------------------------------------------

// Macro for the function that creates the object
#define FUNCTION_OPTIMIZER(CLASS) \
  static vOptimizer *make_optimizer() { return new CLASS(); };

// Macro for the class that links the appropriate function to the map of objects
#define CLASS_OPTIMIZER(NAME,CLASS)                                                           \
  class NAME##OptimizerCreator                                                                \
  {                                                                                           \
    public:                                                                                   \
      NAME##OptimizerCreator()                                                                \
        { sAddonManager::GetInstance()->mp_listOfOptimizers[#NAME] = CLASS::make_optimizer; } \
  };                                                                                          \
  static NAME##OptimizerCreator OptimizerCreator##NAME;

#endif
