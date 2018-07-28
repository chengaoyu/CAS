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
  \brief    Declaration of class oIterativeAlgorithm
*/

#ifndef OITERATIVEALGORITHM_HH
#define OITERATIVEALGORITHM_HH 1

#include "gVariables.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "oProjectorManager.hh"
#include "oImageConvolverManager.hh"
#include "oImageProcessingManager.hh"
#include "oOptimizerManager.hh"
#include "oDeformationManager.hh"
#include "oDynamicModelManager.hh"
#include "oImageSpace.hh"
#include "vDataFile.hh"
#include "sChronoManager.hh"



/*!
  \class  oIterativeAlgorithm
  \brief  This is the main class for iterative reconstructions, that manages the iteration loops. \n
          This class manages an iterative reconstruction of any kind, using a vDataFile, and through the use of an oProjector, an oOptimizer, a oConvolver, a oImageSpace. 
*/
class oIterativeAlgorithm
{
  // Constructor & Destructor
  public:
    /*!
      \fn    public oIterativeAlgorithm::oIterativeAlgorithm
      \brief oIterativeAlgorithm constructor. 
             Initialize the member variables to their default values.
    */
    oIterativeAlgorithm();
    /*!
      \fn    public oIterativeAlgorithm::~oIterativeAlgorithm
      \brief oIterativeAlgorithm destructor. 
    */
    virtual ~oIterativeAlgorithm();



  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn public oIterativeAlgorithm::Iterate
      \brief Just call either the IterateCPU or the IterateGPU function as asked for
      \return 0 if success, positive value otherwise
    */
    int Iterate();
    /*!
      \fn public oIterativeAlgorithm::IterateCPU
      \brief Perform the iterative loop of the algorithm, call the different object 
             for optimization, projection, update, etc.
             Function designed to be executed on the CPU only.
      \details Loops over the iterations, subsets, bed position \n
               Call functions related to each steps of iterative reconstruction: \n
                \n
               StepBeforeIterationLoop() \n
               / Loop on iterations \n
               | StepBeforeSubsetLoop(iteration) \n
               |  / Loop on subsets \n
               |  | StepPreProcessInsideSubsetLoop(iteration,subset) \n
               |  | / Loop on bed positions \n
               |  | | StepInnerLoopInsideSubsetLoop(iteration,subset,bed) \n
               |  | StepPostProcessInsideSubsetLoop(iteration,subset) \n
               |  StepAfterSubsetLoop(iteration) \n
               StepAfterIterationLoop() \n
          
      \return 0 if success, positive value otherwise
    */
    int IterateCPU();
    #ifdef CASTOR_GPU
    /*!
      \fn public oIterativeAlgorithm::IterateGPU
      \brief Perform the iterative loop of the algorithm, call the different object 
             for optimization, projection, update, etc. \n
             Function designed to be executed on the GPU only. \n
             This function is NOT yet implemented
      \return 0 if success, positive value otherwise
    */
    int IterateGPU();
    #endif
    /*!
      \fn     inline public void oIterativeAlgorithm::SetSaveSensitivityHistoFlag()
      \param  bool a_saveSensitivityHistoFlag
      \brief  Set the flag that specifies if the sensitivity image in histogram mode has to be saved for each subset/iteration
    */
    inline void SetSaveSensitivityHistoFlag(bool a_saveSensitivityHistoFlag) {m_saveSensitivityHistoFlag = a_saveSensitivityHistoFlag;}
    /*!
      \fn     inline public void oIterativeAlgorithm::SetSaveSubsetImageFlag()
      \param  bool a_saveImageAfterSubsets
      \brief  Set the flag that specifies if the image has to be saved for each subset
    */
    inline void SetSaveSubsetImageFlag(bool a_saveImageAfterSubsets) {m_saveImageAfterSubsets = a_saveImageAfterSubsets;}
    /*!
      \fn inline public oIterativeAlgorithm::SetOptimizerManager
      \param ap_OptimizerManager
      \brief Set the Optimizer Manager Object 
    */
    inline void SetOptimizerManager(oOptimizerManager* ap_OptimizerManager) {mp_OptimizerManager = ap_OptimizerManager;}
    /*!
      \fn inline public oIterativeAlgorithm::SetImageDimensionsAndQuantification
      \param ap_ImageDimensionsAndQuantification
      \brief Set the Image Dimensions and Quantification Object 
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification) {mp_ID = ap_ImageDimensionsAndQuantification;};
    /*!
      \fn inline public oIterativeAlgorithm::SetImageSpace
      \param ap_ImageSpace
      \brief Set the Image Space Object 
    */
    inline void SetImageSpace(oImageSpace* ap_ImageSpace) {mp_ImageSpace = ap_ImageSpace;};
    /*!
      \fn inline public oIterativeAlgorithm::SetProjectorManager
      \param ap_ProjectorManager
      \brief Set the Projector Manager Object 
    */
    inline void SetProjectorManager(oProjectorManager* ap_ProjectorManager) {mp_ProjectorManager = ap_ProjectorManager;};
    /*!
      \fn inline public oIterativeAlgorithm::SetImageConvolverManager
      \param ap_ImageConvolverManager
      \brief Set the Image Convolver Manager Object 
    */
    inline void SetImageConvolverManager(oImageConvolverManager* ap_ImageConvolverManager) {mp_ImageConvolverManager = ap_ImageConvolverManager;}
    /*!
      \fn inline public oIterativeAlgorithm::SetImageProcessingManager
      \param ap_ImageProcessingManager
      \brief Set the Image Processing Manager Object 
    */
    inline void SetImageProcessingManager(oImageProcessingManager* ap_ImageProcessingManager) {mp_ImageProcessingManager = ap_ImageProcessingManager;}
    /*!
      \fn inline public oIterativeAlgorithm::SetDynamicModelManager
      \param ap_DynamicModelManager
      \brief Set the Dynamic Model Manager Object 
    */
    inline void SetDynamicModelManager(oDynamicModelManager* ap_DynamicModelManager) {mp_DynamicModelManager = ap_DynamicModelManager;}
    /*!
      \fn inline public oIterativeAlgorithm::SetDeformationManager
      \param ap_DeformationManager
      \brief Set the Deformation Manager Object 
    */
    inline void SetDeformationManager(oDeformationManager* ap_DeformationManager) {mp_DeformationManager = ap_DeformationManager;}
    /*!
      \fn inline public oIterativeAlgorithm::SetDataFile
      \param a2p_DataFile
      \brief Set the list of DataFile
    */
    inline void SetDataFile(vDataFile** a2p_DataFile) {m2p_DataFile = a2p_DataFile;};
    /*!
      \fn inline public oIterativeAlgorithm::SetGPUflag
      \param a_flagGPU
      \brief Set the GPU flag
    */
    inline void SetGPUflag(bool a_flagGPU) {m_flagGPU = a_flagGPU;};
    /*!
      \fn inline public oIterativeAlgorithm::SetVerbose
      \param a_verboseLevel
      \brief Set Verbosity
    */
    inline void SetVerbose(int a_verboseLevel) {m_verbose = a_verboseLevel;};
    /*!
      \fn inline public oIterativeAlgorithm::SetNbBeds
      \param a_nbBeds
      \brief Set number of beds (bed positions)
    */
    inline void SetNbBeds(int a_nbBeds) {m_nbBeds = a_nbBeds;};
    /*!
      \fn inline public oIterativeAlgorithm::SetPathInitImage
      \param a_pathToInitialImage
      \brief Set path to an initial image
    */
    inline void SetPathInitImage(string a_pathToInitialImage) {m_pathToInitialImg = a_pathToInitialImage;};
    /*!
      \fn      inline public oIterativeAlgorithm::SetPathToAttenuationImage
      \param   string a_pathToAttenuationImage
      \brief   This function is used to set the path to the attenuation image.
    */
    inline void SetPathToAttenuationImage(string a_pathToAttenuationImage) {m_pathToAtnImg = a_pathToAttenuationImage;};
    /*!
      \fn inline public oIterativeAlgorithm::SetPathToSensitivityImage
      \param a_pathToSensitivityImage
      \brief Set path to the sensitivity image
    */
    inline void SetPathToSensitivityImage(string a_pathToSensitivityImage) {m_pathToSensitivityImg = a_pathToSensitivityImage;};
    /*!
      \fn inline public oIterativeAlgorithm::SetPathToMultiModalImage
      \param a_pathToMultiModalImage
      \brief Set path to multimodal images
    */
    inline void SetPathToMultiModalImage(vector<string> a_pathToMultiModalImage) {m_pathToMultiModalImg = a_pathToMultiModalImage;};
    /*!
      \fn inline public oIterativeAlgorithm::SetPathToMaskImage
      \param a_pathToMaskImage
      \brief Set path to a mask image
    */
    inline void SetPathToMaskImage(string a_pathToMaskImage) {m_pathToMaskImg = a_pathToMaskImage;};
    /*!
      \fn public oIterativeAlgorithm::SetNbIterationsAndSubsets
      \param a_nbIterationsSubsets
      \brief Set the number of iterations and subsets
      \details The provided string is a list of iteration:subset pairs 
               separated by commas.
      \return 0 if success, positive value otherwise
    */
    int SetNbIterationsAndSubsets(const string& a_nbIterationsSubsets);
    /*!
      \fn public oIterativeAlgorithm::SetOutputIterations
      \param a_outputIterations
      \brief Set the selected output iterations
      \details The provided string is a list of couple a:b separated by commas. \n
              It means that we save one iteration over a until b is reached. \n
               "b" must be incrementing for each successive couples.
              If the list is empty, we save all iterations by default. \n
              If the list is equal to "-1", then we save only the last iteration.
      \return 0 if success, positive value otherwise
    */
    int SetOutputIterations(const string& a_outputIterations);
    /*!
      \fn public virtual oIterativeAlgorithm::InitSpecificOptions
      \param a_specificOptions
      \brief Set the selected output iterations
      \details Not Implemented yet
      \return 0
    */
    virtual int InitSpecificOptions(string a_specificOptions);



  // -------------------------------------------------------------------
  // Protected member functions
  protected:
    /*!
      \fn oIterativeAlgorithm::StepBeforeIterationLoop
      \brief This function is called at the beginning of the IterateCPU function.
      \details Initialization and memory allocation for the imageSpace and
               some managers
      \return 0 if success, positive value otherwise.
    */
    virtual int StepBeforeIterationLoop();
    /*!
      \fn oIterativeAlgorithm::StepAfterIterationLoop
      \brief This function is called at the end of the IterateCPU function.
      \details Free Memory for the imageSpace and
               some managers
      \return 0 if success, positive value otherwise.
    */
    virtual int StepAfterIterationLoop();
    /*!
      \fn oIterativeAlgorithm::StepBeforeSubsetLoop
      \param a_iteration : iteration index
      \brief This function is called before starting the subset loop.
      \return 0 if success, positive value otherwise.
    */
    virtual int StepBeforeSubsetLoop(int a_iteration);
    /*!
      \fn oIterativeAlgorithm::StepAfterSubsetLoop
      \param a_iteration : iteration index
      \brief This function is called after finishing the subset loop.
      \details Clean the main images from never visited voxels \n
               Apply post-convolution/post-processing if needed \n
               Write output images on disk as requested by the user
      \todo manage output of coeff  (save parametric images using intrinsic basis functions)
      \return 0 if success, positive value otherwise.
    */
    virtual int StepAfterSubsetLoop(int a_iteration);
    /*!
      \fn oIterativeAlgorithm::StepPreProcessInsideSubsetLoop
      \param a_iteration : iteration index
      \param a_subset : subset index
      \brief This function is called right after starting the subset loop.
             Apply any kind of processing on the forward image before projections
      \details Copy current main image into forward image matrix \n
               Reinitialize backward image and 4D gating indices \n
               Apply image processing/convolution on the forward image matrix
               (image to be projeted)
      \todo    mp_DeformationManager->ApplyDeformationsToForwardImage() : Do we keep this function ?
      \return 0 if success, positive value otherwise.
    */
    virtual int StepPreProcessInsideSubsetLoop(int a_iteration, int a_subset);
    /*!
      \fn oIterativeAlgorithm::StepPreProcessInsideSubsetLoop
      \param a_iteration : iteration index
      \param a_subset : subset index
      \brief This function is called right after starting the subset loop. \n
             Apply any kind of image processing on the backward image and
             main image after backprojections
      \details Synchronize parallel results \n
               Apply image deformation/processing/convolution on the backward image \n
               Update the image using the optimizer functions \n
               Apply dynamic model/image processing/convolution on the main image
      \return 0 if success, positive value otherwise.
    */
    virtual int StepPostProcessInsideSubsetLoop(int a_iteration, int a_subset);
    /*!
      \fn oIterativeAlgorithm::StepInnerLoopInsideSubsetLoop
      \param a_iteration : iteration index
      \param a_subset : subset index
      \param a_bed : bed position
      \brief This function is called inside the subset loop and manages the main loop over the events \n
             The loop over the events is multithreaded, and involves a thread lock in the case of image-based deformation
      \details Step 1: Get the current event for that thread index \n
               Step 2: Call the dynamic switch function that updates the current frame and gate numbers, and also detects involuntary patient motion \n
                       Perform image-based deformation if needed (thread lock is required) \n
               Step 3: Compute the projection line \n
               Step 4: Optimize in the data space (forward-proj, update, backward-proj)
      \todo  Check implementation of thread-locked operation on different architectures
      \return 0 if success, positive value otherwise.
    */
    virtual int StepInnerLoopInsideSubsetLoop(int a_iteration, int a_subset, int a_bed);



  // -------------------------------------------------------------------
  // ----- Functions related to OpenMP barrier -----
    /*!
      \fn oIterativeAlgorithm::ThreadBarrierIncrement
      \brief Increment (thread safe) the m_nbThreadsWaiting variable.
    */
    void ThreadBarrierIncrement();
    /*!
      \fn oIterativeAlgorithm::ThreadBarrierFlag
      \brief Check if the m_releaseThreads boolean is enabled or not 
      \return true if the threads can be released, false otherwise
    */
    bool ThreadBarrierFlag();
    /*!
      \fn oIterativeAlgorithm::ThreadBarrierCheck
      \brief Check if all the threads are currently in waiting condition
      \return true if all the threads are waiting, false otherwise
    */
    bool ThreadBarrierCheck();
    /*!
      \fn oIterativeAlgorithm::ThreadBarrierSetOff
      \brief Disable the thread locking variable (thread safe), 
             and reset the number of waiting threads to 0
    */
    void ThreadBarrierSetOff();
    /*!
      \fn oIterativeAlgorithm::ThreadBarrierSetOn
      \brief Enable the thread locking variable (thread safe), 
             and reset the number of waiting threads to 0
    */
    void ThreadBarrierSetOn();



  // -------------------------------------------------------------------
  // Data members
  protected:
    int m_nbIterations;                                 /*!< Number of iterations (default=1)*/ 
    int* mp_nbSubsets;                                  /*!< Number of subsets (default=1)*/  
    bool* mp_outputIterations;                          /*!< A boolean for each iteration saying if we save it or not */  
    int m_verbose;                                      /*!< Verbosity (default=-1)*/  
    bool m_flagGPU;                                     /*!< Do we use GPU or not (default=false) */  
    oImageDimensionsAndQuantification* mp_ID;           /*!< Pointer to the oImageDimensionsAndQuantification object */
    vDataFile** m2p_DataFile;                           /*!< Pointer to the array of vDataFile object */
    oProjectorManager* mp_ProjectorManager;             /*!< Pointer to the Projector Manager object */
    oOptimizerManager* mp_OptimizerManager;             /*!< Pointer to the Optimizer Manager object */
    oDeformationManager* mp_DeformationManager;         /*!< Pointer to the Deformation Manager object */
    oDynamicModelManager* mp_DynamicModelManager;       /*!< Pointer to the Dynamic Model Manager object */
    oImageSpace* mp_ImageSpace;                         /*!< Pointer to the Image Space object */
    oImageConvolverManager* mp_ImageConvolverManager;   /*!< Pointer to the Image Convolver Manager object */
    oImageProcessingManager* mp_ImageProcessingManager; /*!< Pointer to the Image Processing Manager object */
    int m_nbBeds;                                       /*!< number of bed FOVs (1 datafile by bed) (default=-1) */  
    string m_pathToInitialImg;                          /*!< String containing the path to an initialization image */ 
    string m_pathToAtnImg;                              /*!< String of the path to the attenuation images */
    string m_pathToSensitivityImg;                      /*!< String containing the path to a sensitivity image */ 
    vector<string> m_pathToMultiModalImg;               /*!< String vector containing paths to multimodal images */
    string m_pathToMaskImg;                             /*!< String containing the path to a mask image */ 
    bool m_saveSensitivityHistoFlag;                    /*!< Flag specifying that the sensitivity image has to be saved for each subset/iteration in histogram mode */
    bool m_saveImageAfterSubsets;                       /*!< Flag specifying that the image has to be saved after each subset */
    bool m_generalizedImplementation;                   /*!< No purpose yet. For upcoming features */ 
};

#endif













