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
  \ingroup algorithm
  \brief Implementation of class oIterativeAlgorithm
*/

#include "gVariables.hh"
#include "oIterativeAlgorithm.hh"
#include "iEventHistoCT.hh"

bool m_releaseThreads = false; /*!< Boolean indicating that all threads can be released from a lock (dedicated to OMP multithreading management) */ 
int m_nbThreadsWaiting = 0; /*!< Number of theads currently waiting in a lock (dedicated to OMP multithreading management) */ 

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \brief oIterativeAlgorithm constructor. 
         Initialize the member variables to their default values.
*/
oIterativeAlgorithm::oIterativeAlgorithm()
{
  // Default number of iterations and subsets (1 for each)
  m_nbIterations = 1;
  mp_nbSubsets = (int*)malloc(1*sizeof(int));
  mp_nbSubsets[0] = 1;
  m_generalizedImplementation = true;
  mp_outputIterations = NULL; 
  // set all members to default values
  m_verbose = -1;
  m_flagGPU = false;
  mp_ID = NULL;
  m2p_DataFile= NULL;
  mp_ProjectorManager = NULL;
  mp_OptimizerManager = NULL;
  mp_DeformationManager = NULL;
  mp_DynamicModelManager = NULL;
  mp_ImageSpace = NULL;
  mp_ImageConvolverManager = NULL;
  mp_ImageProcessingManager = NULL;
  m_nbBeds = -1;
  m_pathToInitialImg = "";
  m_pathToSensitivityImg = "";
  m_pathToMaskImg = "";
  m_saveSensitivityHistoFlag = false;
  m_saveImageAfterSubsets = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \brief oIterativeAlgorithm destructor. 
*/
oIterativeAlgorithm::~oIterativeAlgorithm()
{
  // Delete some tables
  if (mp_nbSubsets!=NULL) free(mp_nbSubsets);
  if (mp_outputIterations!=NULL) free(mp_outputIterations);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn SetNbIterationsAndSubsets
  \param a_nbIterationsSubsets
  \brief Set the number of iterations and subsets
  \details The provided string is a list of iteration:subset pairs 
           separated by commas.
  \return 0 if success, positive value otherwise
*/
int oIterativeAlgorithm::SetNbIterationsAndSubsets(const string& a_nbIterationsSubsets)
{
  // If the string is empty, then just keep the defaults
  if (a_nbIterationsSubsets=="") return 0;
  // Otherwise, reset the number of iterations to 0
  m_nbIterations = 0;

  // Copy the string
  string buf = a_nbIterationsSubsets;

  // Loop to search for commas
  size_t pos_comma;
  while ((pos_comma=buf.find_first_of(","))!=string::npos)
  {
    // Get the substring before the comma
    string sub_buf = buf.substr(0,pos_comma);
    // Search for columns
    size_t pos_column = sub_buf.find_first_of(":");
    if (pos_column==string::npos || pos_column==0)
    {
      Cerr("***** oIterativeAlgorithm::SetNbIterationsAndSubsets() -> Syntax problem in number of iterations and subsets !" << endl);
      return 1;
    }
    int iter = atoi(  sub_buf.substr(0,pos_column).c_str()  );
    int subs = atoi(  sub_buf.substr(pos_column+1).c_str()  );
    mp_nbSubsets = (int*)realloc(mp_nbSubsets,(m_nbIterations+iter)*sizeof(int));
    for (int it=0; it<iter; it++) mp_nbSubsets[m_nbIterations+it] = subs;
    m_nbIterations += iter;
    buf = buf.substr(pos_comma+1);
  }

  // Last couple of iterations:subsets
  size_t pos_column = buf.find_first_of(":");
  if (pos_column==string::npos || pos_column==0)
  {
    Cerr("***** oIterativeAlgorithm::SetNbIterationsAndSubsets() -> Syntax problem in number of iterations and subsets !" << endl);
    return 1;
  }
  int iter = atoi(  buf.substr(0,pos_column).c_str()  );
  int subs = atoi(  buf.substr(pos_column+1).c_str()  );
  mp_nbSubsets = (int*)realloc(mp_nbSubsets,(m_nbIterations+iter)*sizeof(int));
  for (int it=0; it<iter; it++) mp_nbSubsets[m_nbIterations+it] = subs;
  m_nbIterations += iter;

  if (m_verbose>=3) 
  {
    Cout("oIterativeAlgorithm::SetNbIterationsAndSubsets() ->  Selected numbers of subsets for each iteration:" << endl); 
    Cout("  Iteration / number of subsets "<< endl); 
    for (int it=0 ; it<m_nbIterations ; it++) Cout("    " << it+1 << "  /  " <<  mp_nbSubsets[it] << endl); 
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn SetOutputIterations
  \param a_outputIterations
  \brief Set the selected output iterations
  \details The provided string is a list of couple a:b separated by commas.
          It means that we save one iteration over a until b is reached.
           "b" must be incrementing for each successive couples.
          If the list is empty, we save all iterations by default.
          If the list is equal to "-1", then we save only the last iteration.
  \return 0 if success, positive value otherwise
*/
int oIterativeAlgorithm::SetOutputIterations(const string& a_outputIterations)
{
  if(m_verbose>=2) Cout("oIterativeAlgorithm::SetOutputIterations ..."<< endl); 
    
  // Allocate the output iterations boolean table
  mp_outputIterations = (bool*)malloc(m_nbIterations*sizeof(bool));
  for (int it=0; it<m_nbIterations; it++) mp_outputIterations[it] = false;

  // If the list is empty, we save all iterations by default
  if (a_outputIterations=="")
  {
    for (int it=0; it<m_nbIterations; it++) mp_outputIterations[it] = true;
    return 0;
  }

  // Copy the string
  string buf = a_outputIterations;

  // Loop to search for commas
  size_t pos_comma;
  int current_iteration = -1;
  while ((pos_comma=buf.find_first_of(","))!=string::npos)
  {
    // Get the substring before the comma
    string sub_buf = buf.substr(0,pos_comma);
    // Search for columns
    size_t pos_column = sub_buf.find_first_of(":");
    if (pos_column==string::npos || pos_column==0)
    {
      Cerr("***** oIterativeAlgorithm::SetOutputIterations() -> Syntax problem in output iteration numbers !" << endl);
      return 1;
    }
    int step_iteration = atoi(  sub_buf.substr(0,pos_column).c_str()  );
    if (step_iteration<=0)
    {
      Cerr("***** oIterativeAlgorithm::SetOutputIterations() -> Iteration step must be strictly positive (here it is " << step_iteration << ") !" << endl);
      return 1;
    }
    int stop_iteration = atoi(  sub_buf.substr(pos_column+1).c_str()  );
    if (stop_iteration<=current_iteration)
    {
      Cerr("***** oIterationAlgorithm::SetOutputIterations() -> Output iteration number must be increasing through provided couples !" << endl);
      return 1;
    }
    for (int it=current_iteration+step_iteration; it<stop_iteration; it+=step_iteration) mp_outputIterations[it] = true;
    current_iteration = stop_iteration-1;
    buf = buf.substr(pos_comma+1);
  }

  // Last couple of iterations:subsets
  size_t pos_column = buf.find_first_of(":");
  if (pos_column==string::npos || pos_column==0)
  {
    // Special case if -1 is provided, it means we save the last iteration
    if (buf=="-1")
    {
      mp_outputIterations[m_nbIterations-1] = true;
      // We directly exist here
      return 0;
    }
    else
    {
      Cerr("***** oIterativeAlgorithm::SetOutputIterations() -> Syntax problem in output iteration numbers !" << endl);
      return 1;
    }
  }
  int step_iteration = atoi(  buf.substr(0,pos_column).c_str()  );
  if (step_iteration<=0)
  {
    Cerr("***** oIterativeAlgorithm::SetOutputIterations() -> Iteration step must be strictly positive (here it is " << step_iteration << ") !" << endl);
    return 1;
  }
  int stop_iteration = atoi(  buf.substr(pos_column+1).c_str()  );
  if (stop_iteration<=current_iteration)
  {
    Cerr("***** oIterationAlgorithm::SetOutputIterations() -> Output iteration number must be increasing through provided couples !" << endl);
    return 1;
  }
  for (int it=current_iteration+step_iteration; it<stop_iteration; it+=step_iteration) mp_outputIterations[it] = true;

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn Iterate
  \brief Just call either the IterateCPU or the IterateGPU function as asked for
  \return 0 if success, positive value otherwise
*/
int oIterativeAlgorithm::Iterate()
{
  if(m_verbose>=2) Cout("oIterativeAlgorithm::Iterate ..."<< endl); 
    
  #ifdef CASTOR_GPU
  if (m_flagGPU) 
    return IterateGPU();
  else 
    return IterateCPU();
  #else
  return IterateCPU();
  #endif
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn IterateCPU
  \brief Perform the iterative loop of the algorithm, call the different object 
         for optimization, projection, update, etc.
         Function designed to be executed on the CPU only.
  \details Loops over the iterations, subsets, bed position
           Call functions related to each steps of iterative reconstruction:
           
           StepBeforeIterationLoop()
           / Loop on iterations
           | StepBeforeSubsetLoop(iteration)
           |  / Loop on subsets
           |  | StepPreProcessInsideSubsetLoop(iteration,subset)
           |  | / Loop on bed positions
           |  | | StepInnerLoopInsideSubsetLoop(iteration,subset,bed)
           |  | StepPostProcessInsideSubsetLoop(iteration,subset)
           |  StepAfterSubsetLoop(iteration)
           StepAfterIterationLoop()
           
  \return 0 if success, positive value otherwise
*/
int oIterativeAlgorithm::IterateCPU()
{
  // Synchronize MPI processes
  #ifdef CASTOR_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // Verbose
  if (m_verbose>=1) Cout("oIterativeAlgorithm::IterateCPU() -> Start algorithm for " << m_nbIterations << " iterations" << endl);

  // Initial clock for execution time
  clock_t clock_start_whole = clock();
  time_t time_start_whole = time(NULL);

  // Call before iteration function
  if (StepBeforeIterationLoop())
  {
    Cerr("***** oIterativeAlgorithm::IterateCPU() -> A problem occured while calling StepBeforeIterationLoop() function !" << endl);
    return 1;
  }

  // Loop on iterations
  for (int iteration=0 ; iteration<m_nbIterations ; iteration++)
  {

    // Call before subset function
    if (StepBeforeSubsetLoop(iteration))
    {
      Cerr("***** oIterativeAlgorithm::IterateCPU() -> A problem occured while calling StepBeforeSubsetLoop() function !" << endl);
      return 1;
    }

    // Loop on subsets
    for (int subset=0 ; subset<mp_nbSubsets[iteration] ; subset++)
    {
      // Verbose
      if (m_verbose>=1) 
      {
        Cout("oIterativeAlgorithm::IterateCPU() -> Start iteration " << iteration+1 << "/" << m_nbIterations 
          << " subset " << subset+1 << "/" << mp_nbSubsets[iteration] << endl);
      }
      
      // Clock start for subset execution time
      clock_t clock_start_subset = clock();
      time_t time_start_subset = time(NULL);

      // Call pre-process inside subset loop function
      if (StepPreProcessInsideSubsetLoop(iteration,subset))
      {
        Cerr("***** oIterativeAlgorithm::IterateCPU() -> A problem occured while calling StepPreProcessInsideSubsetLoop() function !" << endl);
        return 1;
      }

      // Synchronize MPI processes
      #ifdef CASTOR_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      #endif

      // Loop on bed positions
      for (int bed=0 ; bed<m_nbBeds ; bed++)
      {
        // Call the inner loop on events inside subset loop function
        if (StepInnerLoopInsideSubsetLoop(iteration,subset,bed))
        {
          Cerr("***** oIterativeAlgorithm::IterateCPU() -> A problem occured while calling StepInnerLoopInsideSubsetLoop() function !" << endl);
          return 1;
        }
      } // End of beds loop

      // Synchronize MPI processes
      #ifdef CASTOR_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      #endif

      // Call post-process inside subset loop function
      if (StepPostProcessInsideSubsetLoop(iteration,subset))
      {
        Cerr("***** oIterativeAlgorithm::IterateCPU() -> A problem occured while calling StepPostProcessInsideSubsetLoop() function !" << endl);
        return 1;
      }

      // Clock stop for subset execution time
      clock_t clock_stop_subset = clock();
      time_t time_stop_subset = time(NULL);
      if (m_verbose>=2) Cout ("oIterativeAlgorithm::IterateCPU() -> Time spent for subset " << subset+1 << " | User: " << time_stop_subset-time_start_subset
                           << " sec | CPU: " << (clock_stop_subset-clock_start_subset)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);


    } // End of subsets loop

              
    // Call after subset function
    if (StepAfterSubsetLoop(iteration))
    {
      Cerr("***** oIterativeAlgorithm::IterateCPU() -> A problem occured while calling StepAfterSubsetLoop() function !" << endl);
      return 1;
    }
  

  } // End of iterations loop

  //  Call after iteration function
  if (StepAfterIterationLoop())
  {
    Cerr("***** oIterativeAlgorithm::IterateCPU() -> A problem occured while calling StepAfterIterationLoop() function !" << endl);
    return 1;
  }

  // Final clock for execution time
  clock_t clock_stop_whole = clock();
  time_t time_stop_whole = time(NULL);
  if (m_verbose>=1) Cout("oIterativeAlgorithm::IterateCPU() -> Total time spent | User: " << time_stop_whole-time_start_whole 
                      << " sec | CPU: " << (clock_stop_whole-clock_start_whole)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);
  

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn StepBeforeIterationLoop
  \brief This function is called at the beginning of the IterateCPU function.
  \details Initialization and memory allocation for the imageSpace and
           some managers
  \return 0 if success, positive value otherwise.
*/
int oIterativeAlgorithm::StepBeforeIterationLoop()
{
  if (m_verbose>=2) Cout("oIterativeAlgorithm::StepBeforeIterationLoop ... " << endl);

  // Check that the image space is already allocated
  if (mp_ImageSpace==NULL || !mp_ImageSpace->Checked())
  {
    Cerr("***** oIterativeAlgorithm::StepBeforeIterationLoop() -> Image space has not been created or was not checked !" << endl);
    return 1;
  }

  // Instantiate all the required images from the oImageSpace
  mp_ImageSpace->InstantiateImage();
  mp_ImageSpace->InstantiateForwardImage();
  mp_ImageSpace->InstantiateBackwardImageFromDynamicBasis(mp_OptimizerManager->GetNbBackwardImages());
  mp_ImageSpace->InstantiateSensitivityImage(m_pathToSensitivityImg);
  mp_ImageSpace->InstantiateOutputImage();
  mp_ImageSpace->InstantiateVisitedVoxelsImage();
  mp_DeformationManager->InstantiateImageForDeformation(mp_ImageSpace);

  // Main image and sensitivity image initialization
  if (mp_ImageSpace->InitImage(m_pathToInitialImg, mp_OptimizerManager->GetInitialValue()) )
  {
    Cerr("***** oIterativeAlgorithm::StepBeforeIterationLoop() -> An error occured while reading the initialization image !" << endl);
    return 1;
  }

  if (mp_ImageSpace->InitAttenuationImage(m_pathToAtnImg) )
  {
    Cerr("***** oIterativeAlgorithm::StepBeforeIterationLoop()-> Error during attenuation image initialization !" << endl);  
    return 1;
  }

  if (mp_ImageSpace->InitSensitivityImage(m_pathToSensitivityImg))
  {
    Cerr("***** oIterativeAlgorithm::StepBeforeIterationLoop() -> An error occured while initializing the sensitivity image !" << endl);
    return 1;
  }

  if (mp_ImageSpace->InitMultiModalImage(m_pathToMultiModalImg))
  {
    Cerr("***** oIterativeAlgorithm::StepBeforeIterationLoop()-> Error during multimodal image initialization !" << endl);
    return 1;
  }

  if (mp_ImageSpace->InitMaskImage(m_pathToMaskImg))
  {
    Cerr("***** oIterativeAlgorithm::StepBeforeIterationLoop()-> Error during mask image initialization !" << endl);
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
  \fn StepBeforeSubsetLoop
  \param a_iteration : iteration index
  \brief This function is called before starting the subset loop.
  \return 0 if success, positive value otherwise.
*/
int oIterativeAlgorithm::StepBeforeSubsetLoop(int a_iteration)
{
  if (m_verbose>=3) Cout("oIterativeAlgorithm::StepBeforeSubsetLoop ... " << endl);
    
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn StepPreProcessInsideSubsetLoop
  \param a_iteration : iteration index
  \param a_subset : subset index
  \brief This function is called right after starting the subset loop.
         Apply any kind of processing on the forward image before projections
  \details Copy current main image into forward image matrix
           Reinitialize backward image and 4D gating indices
           Apply image processing/convolution on the forward image matrix
           (image to be projeted)
  \return 0 if success, positive value otherwise.
*/
int oIterativeAlgorithm::StepPreProcessInsideSubsetLoop(int a_iteration, int a_subset)
{
  if (m_verbose>=3) Cout("oIterativeAlgorithm::StepPreProcessInsideSubsetLoop ... " << endl);

  // Get the chrono manager singleton pointer
  sChronoManager* p_ChronoManager = sChronoManager::GetInstance();

  // Initialize the correction backward image(s)
  mp_ImageSpace->InitBackwardImage();

  // Copy current image in forward-image buffer (apply deformation if needed)
  mp_ImageSpace->PrepareForwardImage() ; 

  // Apply image processing to forward image
  if (mp_ImageProcessingManager->ApplyProcessingForward(mp_ImageSpace))
  {
    Cerr("***** oIterativeAlgorithm::StepPreProcessInsideSubsetLoop() -> A problem occured while applyin image processing to forward image !" << endl);
    return 1;
  }

  // Apply convolver to forward image
  p_ChronoManager->StartConvolution();
  if (mp_ImageConvolverManager->ConvolveForward(mp_ImageSpace))
  {
    Cerr("***** oIterativeAlgorithm::StepPreProcessInsideSubsetLoop() -> A problem occured while applying image convolver to forward image !" << endl);
    return 1;
  }
  p_ChronoManager->StopConvolution();
  
  // Call the pre-event loop function from the optimizer manager
  mp_OptimizerManager->PreDataUpdateStep(a_iteration, m_nbIterations, a_subset, mp_nbSubsets);

  // Initialisation of the backup images for deformation
  // (bu_fward image initialized with current forward image)
  // (bu_bward and bu_sens images set to 0)
  mp_DeformationManager->InitImageForDeformation(mp_ImageSpace);
  
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn StepInnerLoopInsideSubsetLoop
  \param a_iteration : iteration index
  \param a_subset : subset index
  \param a_bed : bed position
  \brief This function is called inside the subset loop and manages the main loop over the events
         The loop over the events is multithreaded, and involves a thread lock in the case of image-based deformation
  \details Step 1: Get the current event for that thread index
           Step 2: Call the dynamic switch function that updates the current frame and gate numbers, and also detects involuntary patient motion
                   Perform image-based deformation if needed (thread lock is required)
           Step 3: Compute the projection line
           Step 4: Optimize in the data space (forward-proj, update, backward-proj)
  \todo  Check the correct implementation of thread-locked operation on different architectures
  \return 0 if success, positive value otherwise.
*/
int oIterativeAlgorithm::StepInnerLoopInsideSubsetLoop(int a_iteration, int a_subset, int a_bed)
{
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    if (m_nbBeds>1) Cout("oIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> Start loop over events for bed " << a_bed+1 << endl << flush);
    else Cout("oIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> Start loop over events" << endl << flush);
  }
    
  // Reinitialize 4D gate indexes
  mp_ID->ResetCurrentDynamicIndices();

  // Apply the bed offset for this bed position
  mp_ProjectorManager->ApplyBedOffset(a_bed);

  // Progression (increments of 2%)
  if (m_verbose>=VERBOSE_NORMAL && mp_ID->GetMPIRank()==0)
  {
    cout << "0   10%  20%  30%  40%  50%  60%  70%  80%  90%  100%" << endl;
    cout << "|----|----|----|----|----|----|----|----|----|----|" << endl;
    cout << "|" << flush;
  }
  int progression_percentage_old = 0;
  int progression_nb_bars = 0;
  uint64_t progression_printing_index = 0;

  // Compute start and stop indices taking MPI into account (the vDataFile does that)
  int64_t index_start = 0;
  int64_t index_stop  = 0;
  m2p_DataFile[a_bed]->GetEventIndexStartAndStop(&index_start, &index_stop, a_subset, mp_nbSubsets[a_iteration]);

  // Multi-threading with OpenMP
  m_releaseThreads = false;
  m_nbThreadsWaiting = 0;

  // Synchronize MPI processes
  #ifdef CASTOR_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // Set the number of threads for projections (right after this loop, we set back the number of threads for image computations)
  #ifdef CASTOR_OMP
  omp_set_num_threads(mp_ID->GetNbThreadsForProjection());
  #endif

  // This boolean is used to report any problem inside the parallel loop
  bool problem = false;

  // Get the chrono manager singleton pointer
  sChronoManager* p_ChronoManager = sChronoManager::GetInstance();

  // Launch the loop with precomputed start and stop and using a step equal to the number of subsets
  int64_t index;
  // Keep the static scheduling with a chunk size at 1, it is important
  #pragma omp parallel for private(index) schedule(static, 1)
  for ( index = index_start  ;  index < index_stop  ;  index += mp_nbSubsets[a_iteration] )
  {              
    // Get the thread index
    int th = 0;
    #ifdef CASTOR_OMP
    th = omp_get_thread_num();
    #endif
    // Print progression (do not log out with Cout here)
    if (m_verbose>=2 && th==0 && mp_ID->GetMPIRank()==0)
    {
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
    }

    // Step 1: Get the current event for that thread index
    p_ChronoManager->StartIterativeDataUpdateStep1(th);
    #ifdef CASTOR_DEBUG
    if (m_verbose>=4) Cout("oIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> Step1: Get current event for that thread index " << endl);
    #endif
    vEvent* event = m2p_DataFile[a_bed]->GetEvent(index, th);
    if (event==NULL)
    {
      Cerr("***** oIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> An error occured while getting the event from index "
        << index << " (thread " << th << ") !" << endl);
      // Specify that there was a problem
      problem = true;
      // We must continue here because we are inside an OpenMP loop
      continue;
    }
    p_ChronoManager->StopIterativeDataUpdateStep1(th);

    // Step 2: Call the dynamic switch function that updates the current frame and gate numbers, and also detects involuntary patient motion
    p_ChronoManager->StartIterativeDataUpdateStep2(th);
    #ifdef CASTOR_DEBUG
    if (m_verbose>=4)
    {
      Cout("oIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> Step2: Check for Dynamic event (frame/gate switch, image-based deformation " << endl);
    }
    #endif
    int dynamic_switch_value = mp_ID->DynamicSwitch(index, event->GetTimeInMs(), a_bed, th);
    // If the DYNAMIC_SWITCH_CONTINUE is returned, then it means that we are not yet at the first frame
    if ( dynamic_switch_value == DYNAMIC_SWITCH_CONTINUE )
    {
      // Then we just skip this event
      continue;
    }
    // Else, if the DYNAMIC_SWITCH_DEFORMATION is returned, then it means that a change of gate or involuntary motion has occured
    // and is associated to a deformation
    else if ( dynamic_switch_value == DYNAMIC_SWITCH_DEFORMATION )
    {
      // With multi-threads, we should wait all threads to be there to perform the deformation
      #ifdef CASTOR_OMP
      // Count the number of threads reaching this point
      ThreadBarrierIncrement();
      // Loop to lock the threads until the m_releaseThreads condition is enabled
      while (ThreadBarrierFlag())
      {
        // Allow thread 0 to enter once all threads are in the loop
        if ( (th==0) && ThreadBarrierCheck() )
        {
          // Perform here any deformation on the forward image.
          mp_DeformationManager->PerformDeformation(mp_ImageSpace);
          // Release the threads and set the thread counter to 0
          ThreadBarrierSetOff();
        }
      }
      // Count the number of threads reaching this point
      ThreadBarrierIncrement();
      // Set the thread lock condition once all threads reached this point
      if (ThreadBarrierCheck()) ThreadBarrierSetOn();
      #else
      // Perform here any resp related deformations on the forward image
      mp_DeformationManager->PerformDeformation(mp_ImageSpace);
      #endif
    }
    p_ChronoManager->StopIterativeDataUpdateStep2(th);

    // Step 3: Compute the projection line
    p_ChronoManager->StartIterativeDataUpdateStep3(th);
    #ifdef CASTOR_DEBUG
    if (m_verbose>=4) Cout("oIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> Step3: Compute the projection line " << endl);
    #endif
    oProjectionLine *line = mp_ProjectorManager->ComputeProjectionLine(event, th);
    if (line==NULL)
    {
      Cerr("***** oIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> A problem occured while computing the projection line !" << endl);
      // Specify that there was a problem
      problem = true;
      // We must continue here because we are inside an OpenMP loop
      continue;
    }
    p_ChronoManager->StopIterativeDataUpdateStep3(th);

    // Step 4: Optimize in the data space (forward-proj, update, backward-proj)
    p_ChronoManager->StartIterativeDataUpdateStep4(th);
    #ifdef CASTOR_DEBUG
    if (m_verbose>=4) Cout("oIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> Step4: Optimize in the data space " << endl);
    #endif
    if (line->NotEmptyLine()) mp_OptimizerManager->DataUpdateStep( line,
                                                                   event,
                                                                   a_bed,
                                                                   mp_ID->GetCurrentTimeFrame(th),
                                                                   mp_ID->GetCurrentRespImage(th),
                                                                   mp_ID->GetCurrentCardImage(th),
                                                                   a_iteration,
                                                                   th);
    p_ChronoManager->StopIterativeDataUpdateStep4(th);
  } // End of indices loop (OpenMP stops here)

  // Synchronize MPI processes
  #ifdef CASTOR_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // Set back the number of threads for image computation
  #ifdef CASTOR_OMP
  omp_set_num_threads(mp_ID->GetNbThreadsForImageComputation());
  #endif

  // End of progression printing (do not log out with Cout here)
  if (m_verbose>=2 && mp_ID->GetMPIRank()==0)
  {
    int progression_total_bars = 49;
    for (int i=0; i<progression_total_bars-progression_nb_bars; i++) cout << "-";
    cout << "|" << endl;
  }

  // If a problem was encountered, then report it here
  if (problem)
  {
    Cerr("***** oIterativeAlgorithm::StepInnerLoopInsideSubsetLoop() -> A problem occured inside the parallel loop over events !" << endl);
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
  \fn StepPreProcessInsideSubsetLoop
  \param a_iteration : iteration index
  \param a_subset : subset index
  \brief This function is called right after starting the subset loop.
         Apply any kind of image processing on the backward image and
         main image after backprojections
  \details Synchronize parallel results
           Apply image deformation/processing/convolution on the backward image
           Update the image using the optimizer functions
           Apply dynamic model/image processing/convolution on the main image
  \return 0 if success, positive value otherwise.
*/
int oIterativeAlgorithm::StepPostProcessInsideSubsetLoop(int a_iteration, int a_subset)
{
  if (m_verbose>=3) Cout("oIterativeAlgorithm::StepPostProcessInsideSubsetLoop ... " << endl);

  // Get the chrono manager singleton pointer
  sChronoManager* p_ChronoManager = sChronoManager::GetInstance();

  // Merge parallel results
  mp_ImageSpace->Reduce();
  
  // If mask provided, apply mask to sensitivity, the CleanNeverVisitedVoxels will take care of the images
  mp_ImageSpace->ApplyMaskToSensitivity();

  // Perform here any image-based deformations related to motion on the backward image.
  mp_DeformationManager->ApplyDeformationsToBackwardImage(mp_ImageSpace);
  
  // Apply convolver to backward images and sensitivity-if-needed
  p_ChronoManager->StartConvolution();
  if (mp_ImageConvolverManager->ConvolveBackward(mp_ImageSpace))
  {
    Cerr("***** oIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occured while applying convolver to backward images !" << endl);
    return 1;
  }
  p_ChronoManager->StopConvolution();
  
  // Call the post data update step function from the optimizer manager
  mp_OptimizerManager->PostDataUpdateStep(a_iteration, m_nbIterations, a_subset, mp_nbSubsets);

  // Optimize in the image space (apply corrections, MAP and sensitivity); pass the number of subsets for list-mode sensitivity scaling
  mp_OptimizerManager->ImageUpdateStep(a_iteration, mp_nbSubsets[a_iteration]);

  // Apply convolver to current estime images
  p_ChronoManager->StartConvolution();
  if (mp_ImageConvolverManager->ConvolveIntra(mp_ImageSpace))
  {
    Cerr("***** oIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occured while applying convolver to current estimate images !" << endl);
    return 1;
  }
  p_ChronoManager->StopConvolution();

  // Post-update Dynamic Modeling step (if enabled). Manage either linear/non linear dynamic model (physiological or not)
  if (mp_DynamicModelManager->ApplyDynamicModel(mp_ImageSpace, a_iteration, mp_nbSubsets[a_iteration]))
  {
    Cerr("***** oIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occured while applying dynamic model to current estimate images !" << endl);
    return 1;
  }
  
  // Apply image processing to current estime images
  if (mp_ImageProcessingManager->ApplyProcessingIntra(mp_ImageSpace))
  {
    Cerr("***** oIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occured while applying image processing to current estimate images !" << endl);
    return 1;
  }

  // Save the sensitivity image in histogram mode, if asked for
  if (mp_ID->GetMPIRank()==0 && mp_outputIterations[a_iteration] && m_saveSensitivityHistoFlag && m2p_DataFile[0]->GetDataMode()==MODE_HISTOGRAM)
  {
    // Get output manager to build the file name
    sOutputManager* p_output_manager = sOutputManager::GetInstance();
    // Build the file name
    string temp_sens = p_output_manager->GetPathName() + p_output_manager->GetBaseName();
    stringstream temp_it; temp_it << a_iteration + 1;
    stringstream temp_ss; temp_ss << a_subset + 1;
    temp_sens.append("_it").append(temp_it.str()).append("_ss").append(temp_ss.str()).append("_sensitivity");
    // Save sensitivity
    mp_ImageSpace->SaveSensitivityImage(temp_sens);
  }

  // Save the current image estimate if asked for
  if (mp_ID->GetMPIRank()==0 && mp_outputIterations[a_iteration] && m_saveImageAfterSubsets)
  {
    // Verbose
    if (m_verbose>=1) Cout("oIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> Save image at iteration " << a_iteration+1 << " for subset " << a_subset+1 << endl);
    // Build image from basis functions
    mp_ImageSpace->ComputeOutputImage();
    // Apply post-convolution if needed
    p_ChronoManager->StartConvolution();
    if (mp_ImageConvolverManager->ConvolvePost(mp_ImageSpace))
    {
      Cerr("***** oIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occured while convolving the output image !" << endl);
      return 1;
    }
    p_ChronoManager->StopConvolution();
    // Apply post-processing if needed
    if (mp_ImageProcessingManager->ApplyProcessingPost(mp_ImageSpace))
    {
      Cerr("***** oIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occured while applying image processing the output image !" << endl);
      return 1;
    }
    // Save Dynamic images
    if (mp_DynamicModelManager->SaveParametricImages(a_iteration, a_subset))
    {
      Cerr("***** oIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occured while saving parametric images related to the dynamic model !" << endl);
      return 1;
    }
    // Apply output transaxial FOV masking
    if (mp_ImageSpace->ApplyOutputFOVMasking())
    {
      Cerr("***** oIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occured while applying output FOV masking !" << endl);
      return 1;
    }
    // Apply output flip
    if (mp_ImageSpace->ApplyOutputFlip())
    {
      Cerr("***** oIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occured while applying output flip !" << endl);
      return 1;
    }
    // Save output image
    if (mp_ImageSpace->SaveOutputImage(a_iteration, a_subset))
    {
      Cerr("***** oIterativeAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occured while saving output image !" << endl);
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
/*
  \fn StepAfterSubsetLoop
  \param a_iteration : iteration index
  \brief This function is called after finishing the subset loop.
  \details Clean the main images from never visited voxels
           Apply post-convolution/post-processing if needed
           Write output images on disk as requested by the user
  \todo : save parametric images using intrinsic basis functions ?
  \return 0 if success, positive value otherwise.
*/
int oIterativeAlgorithm::StepAfterSubsetLoop(int a_iteration)
{
  if (m_verbose>=3) Cout("oIterativeAlgorithm::StepAfterSubsetLoop ... " << endl);

  // Get the chrono manager singleton pointer
  sChronoManager* p_ChronoManager = sChronoManager::GetInstance();

  // Clean never visited voxels
  mp_ImageSpace->CleanNeverVisitedVoxels();
  // Save the main image
  if (mp_ID->GetMPIRank()==0 && mp_outputIterations[a_iteration])
  {
    if (m_verbose>=1) Cout("oIterativeAlgorithm::StepAfterSubsetLoop() -> Save image at iteration " << a_iteration+1 << endl);
    // Build image from basis functions
    mp_ImageSpace->ComputeOutputImage();
    // Apply post-convolution if needed
    p_ChronoManager->StartConvolution();
    if (mp_ImageConvolverManager->ConvolvePost(mp_ImageSpace))
    {
      Cerr("***** oIterativeAlgorithm::StepAfterSubsetLoop() -> A problem occured while convolving the output image !" << endl);
      return 1;
    }
    p_ChronoManager->StopConvolution();
    // Apply post-processing if needed
    if (mp_ImageProcessingManager->ApplyProcessingPost(mp_ImageSpace))
    {
      Cerr("***** oIterativeAlgorithm::StepAfterSubsetLoop() -> A problem occured while applying image processing the output image !" << endl);
      return 1;
    }
    // Save Dynamic images
    if (mp_DynamicModelManager->SaveParametricImages(a_iteration))
    {
      Cerr("***** oIterativeAlgorithm::StepAfterSubsetLoop() -> A problem occured while saving parametric images related to the dynamic model !" << endl);
      return 1;
    }
    // Todo : save parametric images using intrinsic basis functions
    // Apply output transaxial FOV masking
    if (mp_ImageSpace->ApplyOutputFOVMasking())
    {
      Cerr("***** oIterativeAlgorithm::StepAfterSubsetLoop() -> A problem occured while applying output FOV masking !" << endl);
      return 1;
    }
    // Apply output flip
    if (mp_ImageSpace->ApplyOutputFlip())
    {
      Cerr("***** oIterativeAlgorithm::StepAfterSubsetLoop() -> A problem occured while applying output flip !" << endl);
      return 1;
    }
    // Save output image
    if (mp_ImageSpace->SaveOutputImage(a_iteration))
    {
      Cerr("***** oIterativeAlgorithm::StepAfterSubsetLoop() -> A problem occured while saving output image !" << endl);
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
/*
  \fn StepAfterIterationLoop
  \brief This function is called at the end of the IterateCPU function.
  \details Free Memory for the imageSpace and
           some managers
  \return 0 if success, positive value otherwise.
*/
int oIterativeAlgorithm::StepAfterIterationLoop()
{
  if (m_verbose>=2) Cout("oIterativeAlgorithm::StepAfterIterationLoop ... " << endl);
  
  // Deallocate everything
  mp_DeformationManager->DeallocateImageForDeformation(mp_ImageSpace);
  mp_ImageSpace->DeallocateBackwardImageFromDynamicBasis();
  mp_ImageSpace->DeallocateSensitivityImage();
  mp_ImageSpace->DeallocateMultiModalImage();
  mp_ImageSpace->DeallocateMaskImage();
  mp_ImageSpace->DeallocateForwardImage();
  mp_ImageSpace->DeallocateImage();
  mp_ImageSpace->DeallocateOutputImage();
  mp_ImageSpace->DeallocateVisitedVoxelsImage();

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InitSpecificOptions
  \param a_specificOptions
  \brief Set the selected output iterations
  \details Not Implemented yet
  \return 0
*/
int oIterativeAlgorithm::InitSpecificOptions(string a_specificOptions)
{
  if (m_verbose>=2) Cout("oIterativeAlgorithm::InitSpecificOptions ... " << endl);
    
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
// ----- Functions dedicated to barrier for multithreading ----- 
// TODO we should keep in mind that different compiler may deal differently with OpenMP 
// (as the compiler may take liberties and optimize shared variables as register variables, 
// leading to non-concurrent observations of the variable. 
// And this is not cool as this can completely stop the reconstruction process)
// The omp-flush directives should prevent that, but further tests may be required

/*
  \fn ThreadBarrierIncrement
  \brief Increment (thread safe) the m_nbThreadsWaiting variable.
*/
void oIterativeAlgorithm::ThreadBarrierIncrement()
{
  #ifdef CASTOR_DEBUG
  if (m_verbose>=4) Cout("oIterativeAlgorithm::ThreadBarrierIncrement ... " << endl);
  #endif
  
  #pragma omp critical(increment_waiting_threads)
  {
    m_nbThreadsWaiting++;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ThreadBarrierFlag
  \brief Check if the m_releaseThreads boolean is enabled or not 
  \return true if the threads can be released, false otherwise
*/
bool oIterativeAlgorithm::ThreadBarrierFlag()
{
  #ifdef CASTOR_DEBUG
  if (m_verbose>=4) Cout("oIterativeAlgorithm::ThreadBarrierFlag ... " << endl);
  #endif
  
  // Check any change in m_releaseThreads
  #pragma omp flush(m_releaseThreads)
  return !m_releaseThreads;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ThreadBarrierCheck
  \brief Check if all the threads are currently in waiting condition
  \return true if all the threads are waiting, false otherwise
*/
bool oIterativeAlgorithm::ThreadBarrierCheck()
{ 
  #ifdef CASTOR_DEBUG
  if (m_verbose>=4) Cout("oIterativeAlgorithm::ThreadBarrierCheck ... " << endl);
  #endif
  
  #pragma omp flush(m_nbThreadsWaiting) // Check any change in m_nbThreadsWaiting
  return (m_nbThreadsWaiting==mp_ID->GetNbThreadsForProjection());
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ThreadBarrierSetOff
  \brief Disable the thread locking variable (thread safe), and reset the number of waiting threads to 0
*/
void oIterativeAlgorithm::ThreadBarrierSetOff()
{
  #ifdef CASTOR_DEBUG
  if (m_verbose>=4) Cout("oIterativeAlgorithm::ThreadBarrierSetOff ... " << endl);
  #endif
  
  m_releaseThreads = true;
  #pragma omp flush(m_releaseThreads)   // Make sure the value of m_releaseThreads is propagated to other threads
  m_nbThreadsWaiting = 0;
  #pragma omp flush(m_nbThreadsWaiting)   // Same for m_nbThreadsWaiting
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ThreadBarrierSetOn
  \brief Enable the thread locking variable (thread safe), and reset the number of waiting threads to 0
*/
void oIterativeAlgorithm::ThreadBarrierSetOn()
{
  #ifdef CASTOR_DEBUG
  if (m_verbose>=4) Cout("oIterativeAlgorithm::ThreadBarrierSetOn ... " << endl);
  #endif
  
  m_releaseThreads = false;
  #pragma omp flush(m_releaseThreads)   // Make sure the value of m_releaseThreads is propagated to other threads
  m_nbThreadsWaiting = 0;
  #pragma omp flush(m_nbThreadsWaiting)   // Same for m_nbThreadsWaiting
}
