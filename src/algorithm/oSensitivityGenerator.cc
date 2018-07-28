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
  \brief Implementation of class oSensitivityGenerator
*/

#include "gVariables.hh"
#include "oSensitivityGenerator.hh"
#include "iDataFilePET.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oSensitivityGenerator::oSensitivityGenerator()
{
  // Simply default all members to "empty" values
  m_pathToSensitivityImage = "";
  mp_ImageDimensionsAndQuantification = NULL;
  mp_ImageSpace = NULL;
  mp_Scanner = NULL;
  mp_ProjectorManager = NULL;
  mp_ImageConvolverManager = NULL;
  mp_DeformationManager = NULL;
  m2p_DataFile = NULL;
  m_computeFromHistogramFlag = false;
  mp_pathToNormalizationFileName = {};
  m_inverseDataFileOrderFlag = false;
  m3p_NormalizationDataFile = NULL;
  m_oneNormalizationFileForAllBeds = false;
  m_oneNormalizationFileForAllFrames = false;
  m_pathToAttenuationImage = "";
  m_pathToMaskImg = "";
  m_mumapAttenuationFlag = false;
  m_forwardProjectAttenuation = false;
  m_nbAtnRespGateImages = -1;
  m_nbAtnCardGateImages = -1;
  mp_lineCounter = NULL;
  m_flagGPU = false;
  m_verbose = -1;
  m_checked = false;
  m_initialized = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oSensitivityGenerator::~oSensitivityGenerator()
{
  // Delete all normalization data files if any
  if (m3p_NormalizationDataFile)
  {
    for (int bed=0; bed<mp_ImageDimensionsAndQuantification->GetNbBeds(); bed++) if (m3p_NormalizationDataFile[bed])
    {
      for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
        if (m3p_NormalizationDataFile[bed][fr]) delete m3p_NormalizationDataFile[bed][fr];
      free(m3p_NormalizationDataFile[bed]);
    }
    free(m3p_NormalizationDataFile);
  }
  if (mp_lineCounter) free(mp_lineCounter);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oSensitivityGenerator::CheckParameters()
{
  // Check all mandatory parameters
  if (mp_ImageDimensionsAndQuantification==NULL)
  {
    Cerr("***** oSensitivityGenerator::CheckParameters() -> Image dimensions and quantification object is null !" << endl);
    return 1;
  }
  if (mp_ImageSpace==NULL)
  {
    Cerr("***** oSensitivityGenerator::CheckParameters() -> Image space object is null !" << endl);
    return 1;
  }
  if (mp_Scanner==NULL)
  {
    Cerr("***** oSensitivityGenerator::CheckParameters() -> Scanner object is null !" << endl);
    return 1;
  }
  if (mp_ProjectorManager==NULL)
  {
    Cerr("***** oSensitivityGenerator::CheckParameters() -> Projector manager object is null !" << endl);
    return 1;
  }
  if (mp_ImageConvolverManager==NULL)
  {
    Cerr("***** oSensitivityGenerator::CheckParameters() -> Convolver manager object is null !" << endl);
    return 1;
  }
  if (mp_DeformationManager==NULL)
  {
    Cerr("***** oSensitivityGenerator::CheckParameters() -> Deformation manager object is null !" << endl);
    return 1;
  }
  if (m2p_DataFile==NULL)
  {
    Cerr("***** oSensitivityGenerator::CheckParameters() -> Data file array is null !" << endl);
    return 1;
  }
  for (int b=0; b<mp_ImageDimensionsAndQuantification->GetNbBeds(); b++)
  {
    if (m2p_DataFile[b]==NULL)
    {
      Cerr("***** oSensitivityGenerator::CheckParameters() -> Data file object for bed " << b+1 << " is null !" << endl);
      return 1;
    }
  }
  if (m_verbose<0)
  {
    Cerr("***** oSensitivityGenerator::CheckParameters() -> Verbose level is negative !" << endl);
    return 1;
  }
  // If compute the sensitivity from a histogram datafile, then check if the datafile is actually an histogram
  if (m_computeFromHistogramFlag && m2p_DataFile[0]->GetDataMode()!=MODE_HISTOGRAM)
  {
    Cerr("***** oSensitivityGenerator::CheckParameters() -> It was asked to compute global sensitivity from the provided datafile whereas it is not a histogram !" << endl);
    return 1;
  }
  // Sensitivity computation from histogram do not work with dynamism for the moment
  if (m_computeFromHistogramFlag && ( mp_ImageDimensionsAndQuantification->GetNbTimeFrames()>1 ||
                                      mp_ImageDimensionsAndQuantification->GetNbCardGates()>1 ||
                                      mp_ImageDimensionsAndQuantification->GetNbRespGates()>1 ))
  {
    Cerr("***** oSensitivityGenerator::CheckParameters() -> Global sensitivity computation from histogram does not work wth dynamic data yet !" << endl);
    return 1;
  }
  // Now it is checked
  m_checked = true;
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oSensitivityGenerator::Initialize()
{
  // First check that the parameters has been checked !
  if (!m_checked)
  {
    Cerr("***** oSensitivityGenerator::Initialize() -> Must call the CheckParameters() function before initialization !" << endl);
    return 1;
  }
  // Verbose
  if (m_verbose>=2) Cout("oSensitivityGenerator::Initialize() -> Start initialization" << endl);
  // For the moment, SPECT and CT sensitivity are not implemented
  if (m2p_DataFile[0]->GetDataType()==TYPE_SPECT || m2p_DataFile[0]->GetDataType()==TYPE_CT)
  {
    Cerr("***** oSensitivityGenerator::Initialize() -> Sensitivity for SPECT and CT not yet implemented !" << endl);
    return 1;
  }
  
  // Allocate the forward image
  mp_ImageSpace->LMS_InstantiateForwardImage();
  
  // Initialization of normalization files:
  //  --> must be done before attenuation files, because attenuation can be inside the normalization data file
  //  --> is not done if sensitivity generated from the histogram datafile
  if (!m_computeFromHistogramFlag && InitializeNormalizationFiles())
  {
    Cerr("***** oSensitivityGenerator::Initialize() -> A problem occured while initializing the normalization data files !" << endl);
    return 1;
  }
  // Initialization of attenuation related stuff: not done if sensitivity generated from the histogram datafile
  if (!m_computeFromHistogramFlag && InitializeAttenuationFiles())
  {
    Cerr("***** oSensitivityGenerator::Initialize() -> A problem occured while initializing the attenuation files !" << endl);
    return 1;
  }
  // We want to know if we will need to forward project into the attenuation image, so here the full condition
  if (
       m_mumapAttenuationFlag && // 1. We must have an attenuation image that is loaded
       m2p_DataFile[0]->GetDataType()==TYPE_PET && // 2. This is only for PET
       m2p_DataFile[0]->GetDataMode()==MODE_LIST && // 3. This is only for sensitivity generation for list-mode datafile
       ( m3p_NormalizationDataFile==NULL || // 4.1 There is no normalization data file (potentially including ACF)
         !(dynamic_cast<iDataFilePET*>(m3p_NormalizationDataFile[0][0]))->GetAtnCorrectionFlag() ) // 4.2 OR the normalization data file does not contain ACF
     )
  {
    // We use a specific boolean that says we need to forward project into the attenuation image
    m_forwardProjectAttenuation = true;
  }
  else
  {
    // Otherwise, we don't
    m_forwardProjectAttenuation = false;
  }

  // Allocate the backward image for multi-thread
  mp_ImageSpace->InstantiateBackwardImageFromDynamicBins();
  // Allocate the mask image if provided
  if (mp_ImageSpace->InitMaskImage(m_pathToMaskImg))
  {
    Cerr("***** oSensitivityGenerator::Initialize() -> Error during mask image initialization !" << endl);
    return 1;
  }
  // Create sensitivity image file name
  sOutputManager* p_outputManager = sOutputManager::GetInstance();
  m_pathToSensitivityImage = p_outputManager->GetPathName() + p_outputManager->GetBaseName() + "_sensitivity";
  // Allocate line counters
  mp_lineCounter = (uint64_t*)calloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection(),sizeof(uint64_t));
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oSensitivityGenerator::InitializeNormalizationFiles()
{
  // TODO: check the consistency between all normalization files (atnflag, normfactor, etc)

  // If no normalization file name provided, then we will exit, but we do some checks before
  if (mp_pathToNormalizationFileName.size()==0)
  {
    // If the ignore-normalization-correction-flag is on, then we can exit the function
    if (mp_ImageDimensionsAndQuantification->GetIgnoreNormCorrectionFlag()) return 0;
    // In PET, throw an error if normalization correction is in the datafile and is taken into account (only if the datafile is a list-mode)
    if (m2p_DataFile[0]->GetDataType()==TYPE_PET && m2p_DataFile[0]->GetDataMode()==MODE_LIST)
    {
      // Cast the vDataFile into iDataFilePET
      iDataFilePET* p_pet_file = (dynamic_cast<iDataFilePET*>(m2p_DataFile[0]));
      // Check if the normalization data are in the file and we do not ignore it
      if (p_pet_file->GetNormCorrectionFlag())
      {
        Cerr("***** oSensitivityGenerator::InitializeNormalizationFiles() -> Normalization correction is included in the data file while it is not in the sensitivity computation !" << endl);
        return 1;
      }
    }
    // We can exit now
    return 0;
  }

  // Verbose
  if (m_verbose>=2) Cout("oSensitivityGenerator::InitializeNormalizationFiles() -> Initialize normalization files" << endl);

  // Check that the number of normalization files options is the same as the number of beds
  if (mp_pathToNormalizationFileName.size()!=((size_t)(mp_ImageDimensionsAndQuantification->GetNbBeds())))
  {
    // If only one is provided, then we assume that it is the same for all bed positions
    if (mp_pathToNormalizationFileName.size()==1) m_oneNormalizationFileForAllBeds = true;
    // Otherwise, throw an error
    else
    {
      Cerr("***** oSensitivityGenerator::InitializeNormalizationFiles() -> Number of normalization files options should be one or equal to the number of data files !" << endl);
      return 1;
    }
  }

  // Allocate the normalization data files array for the bed positions
  m3p_NormalizationDataFile = (vDataFile***)malloc(mp_ImageDimensionsAndQuantification->GetNbBeds()*sizeof(vDataFile**));

  // Get the scanner manager
  sScannerManager* p_ScannerManager = sScannerManager::GetInstance();

  // Loop on the number of beds
  for (int bed=0; bed<mp_ImageDimensionsAndQuantification->GetNbBeds(); bed++)
  {
    // --------------------------------------------------------------------------
    // Step 0: When only one normalization file names is provided for all beds
    // --------------------------------------------------------------------------
    if (m_oneNormalizationFileForAllBeds)
    {
      // If not the first bed
      if (bed!=0)
      {
        // We copy the pointer of the first bed into this one
        m3p_NormalizationDataFile[bed] = m3p_NormalizationDataFile[0];
        // And we continue the loop
        continue;
      }
    }
    // --------------------------------------------------------------------------
    // Step 1: Get filenames and check consistency with number of beds and frames
    // --------------------------------------------------------------------------
    // We will get the normalization file names for this bed
    vector<string> normalization_files;
    // Compute the bed index associated to the file name; in case of reverse order, we change the index only for datafile management here
    int bed_data_file_name_index = bed;
    if (m_inverseDataFileOrderFlag) bed_data_file_name_index = mp_ImageDimensionsAndQuantification->GetNbBeds() - 1 - bed;
    // Then we search for commas separating the file name associated to each frame
    size_t comma = mp_pathToNormalizationFileName[bed_data_file_name_index].find_first_of(",");
    // Loop to get all file names
    while (comma!=string::npos)
    {
      // Get the file name before the comma
      normalization_files.push_back(mp_pathToNormalizationFileName[bed_data_file_name_index].substr(0,comma));
      // Remove this part from the collection of file names
      mp_pathToNormalizationFileName[bed_data_file_name_index] = mp_pathToNormalizationFileName[bed_data_file_name_index].substr(comma+1);
      // Search for the next comma
      comma = mp_pathToNormalizationFileName[bed_data_file_name_index].find_first_of(",");
    }
    // Get the last file name
    normalization_files.push_back(mp_pathToNormalizationFileName[bed_data_file_name_index]);
    // If we have only one frame, then set the m_oneNormalizationFileForAllFrames to true right now
    if (mp_ImageDimensionsAndQuantification->GetNbTimeFrames()==1) m_oneNormalizationFileForAllFrames = true;
    // Check that the number of file names found is equal to the number of frames
    if (normalization_files.size()!=((size_t)(mp_ImageDimensionsAndQuantification->GetNbTimeFrames())))
    {
      // If only one is provided, then we assume that it is the same for all frames
      if (normalization_files.size()==1) m_oneNormalizationFileForAllFrames = true;
      // Otherwise, throw an error
      else
      {
        Cerr("***** oSensitivityGenerator::InitializeNormalizationFiles() -> Number of normalization files for bed " << bed+1 << " should be one or equal to the number of frames !" << endl);
        return 1;
      }
    }
    // If this is the case, then check for previous beds that it was so
    else if (bed!=0 && !m_oneNormalizationFileForAllFrames) 
    {
      Cerr("***** oSensitivityGenerator::InitializeNormalizationFiles() -> If only one normalization file is provided for all frames, then it must be so for all beds !" << endl);
      return 1;
    }
    // --------------------------------------------------------------------------
    // Step 2: Create the normalization data file objects
    // --------------------------------------------------------------------------
    // Allocate the normalization data file array for the frames
    m3p_NormalizationDataFile[bed] = (vDataFile**)malloc(mp_ImageDimensionsAndQuantification->GetNbTimeFrames()*sizeof(vDataFile*));
    // Loop on the number of frames
    for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
    {
      // If one normalization file for all frames
      if (m_oneNormalizationFileForAllFrames)
      {
        // Then, if not the first frame
        if (fr!=0)
        {
          // We copy the pointer of the first frame into this one
          m3p_NormalizationDataFile[bed][fr] = m3p_NormalizationDataFile[bed][0];
          // And we continue the loop
          continue;
        }
      }
      // Switch on the scanner type and create the specific data files
      if (p_ScannerManager->GetScannerType() == SCANNER_PET) m3p_NormalizationDataFile[bed][fr] = new iDataFilePET();
      else
      {
        Cerr("***** oSensitivityGenerator::InitializeNormalizationFiles() -> Unknown scanner type (" << p_ScannerManager->GetScannerType() << ") or not yet implemented !" << endl);
        return 1;
      }
      m3p_NormalizationDataFile[bed][fr]->SetHeaderDataFileName(normalization_files[fr]);
      m3p_NormalizationDataFile[bed][fr]->SetBedIndex(bed);
      m3p_NormalizationDataFile[bed][fr]->SetVerbose(m2p_DataFile[0]->GetVerbose());
      m3p_NormalizationDataFile[bed][fr]->SetImageDimensionsAndQuantification(mp_ImageDimensionsAndQuantification);
      bool affect_quantification_flag = false;
      if (m3p_NormalizationDataFile[bed][fr]->ReadInfoInHeader(affect_quantification_flag))
      {
        Cerr("***** oSensitivityGenerator::InitializeNormalizationFiles() -> A problem occurred during normalization datafile header reading !" << endl);
        return 1;
      }
      if (m3p_NormalizationDataFile[bed][fr]->CheckParameters())
      {
        Cerr("***** oSensitivityGenerator::InitializeNormalizationFiles() -> A problem occurred while checking normalization datafile parameters !" << endl);
        return 1;
      }
      if (m3p_NormalizationDataFile[bed][fr]->ComputeSizeEvent())
      {
        Cerr("***** oSensitivityGenerator::InitializeNormalizationFiles() -> A problem occurred in normalization datafile event size computation !" << endl);
        return 1;
      }
      if (m3p_NormalizationDataFile[bed][fr]->InitializeMappedFile())
      {
        Cerr("***** oSensitivityGenerator::InitializeNormalizationFiles() -> A problem occurred in normalization datafile initialization !" << endl);
        return 1;
      }
      if (m3p_NormalizationDataFile[bed][fr]->PrepareDataFile())
      {
        Cerr("***** oSensitivityGenerator::InitializeNormalizationFiles() -> A problem occured in normalization datafile preparation !" << endl);
        return 1;
      }
      // Check that the normalization data mode is histogram or normalization
      if (m3p_NormalizationDataFile[bed][fr]->GetDataMode()!=MODE_HISTOGRAM && m3p_NormalizationDataFile[bed][fr]->GetDataMode()!=MODE_NORMALIZATION)
      {
        Cerr("***** oSensitivityGenerator::InitializeNormalizationFiles() -> Normalization data file used for sensitivity computation must be either histogram or normalization mode !" << endl);
        return 1;
      }
    }
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oSensitivityGenerator::InitializeAttenuationFiles()
{
  // If it is asked to ignore attenuation correction, then nothing to do
  if (mp_ImageDimensionsAndQuantification->GetIgnoreAttnCorrectionFlag()) return 0;

  // If attenuation image file name is empty, then we do some checks before
  if (m_pathToAttenuationImage=="")
  {
    // Case for PET: if no normalization data file OR no ACF in the normalization data file
    if (m2p_DataFile[0]->GetDataType()==TYPE_PET && (m3p_NormalizationDataFile==NULL || !(dynamic_cast<iDataFilePET*>(m3p_NormalizationDataFile[0][0]))->GetAtnCorrectionFlag()))
    {
      // Cast the vDataFile into iDataFilePET
      iDataFilePET* p_pet_file = (dynamic_cast<iDataFilePET*>(m2p_DataFile[0]));
      // Throw an error if attenuation correction is in the datafile (at this step we know we do not ignore it)
      if (p_pet_file->GetAtnCorrectionFlag())
      {
        Cerr("***** oSensitivityGenerator::InitializeAttenuationFiles() -> Attenuation correction is included in the data file while it is not in the sensitivity computation !" << endl);
        return 1;
      }
    }
    // We can exit now
    return 0;
  }

  // In PET, if the ACF is already present in the normalization data file, then say it and ignore mumap
  if (m2p_DataFile[0]->GetDataType()==TYPE_PET && m3p_NormalizationDataFile!=NULL && (dynamic_cast<iDataFilePET*>(m3p_NormalizationDataFile[0][0]))->GetAtnCorrectionFlag())
  {
    Cout("oSensitivityGenerator::InitializeAttenuationFiles() -> Ignore provided mumap and consider ACF provided in the normalization data file" << endl);
    return 0;
  }

  // Verbose
  if (m_verbose>=2) Cout("oSensitivityGenerator::InitializeAttenuationFiles() -> Allocate and read attenuation images (assumed to be in cm-1)" << endl);

  // Allocate and read the attenuation image (this puts the image into the m4p_attenuation buffer of oImageSpace)
  if (mp_ImageSpace->InitAttenuationImage(m_pathToAttenuationImage))
  {
    Cerr("***** oSensitivityGenerator::Initialize() -> A problem occured while initializing the attenuation image into the image space !" << endl);  
    return 1;
  }

  // Allocate the foward image, where the attenuation will be copied and deformed if needed
  //mp_ImageSpace->LMS_InstantiateForwardImage();

  // Copy attenuation image into forward-image buffer
  mp_ImageSpace->LMS_CopyAtnToForwardImage(mp_DeformationManager->UseDeformationResp() || mp_DeformationManager->UseDeformationInv(),
                                           mp_DeformationManager->UseDeformationCard() );

  // We do not need the m4p_attenuation buffer of oImageSpace anymore, so free it now
  mp_ImageSpace->LMS_DeallocateAttenuationImage();

  // TODO: Think about adding a field in the attenuation header to specify if the attenuation images have the same spatial resolution as the
  //       scanner in use, in order to apply a smoothing on it. For the moment, we assume it has been done beforehand.

  // TODO: for the moment, we assume that the attenuation image is in cm-1, maybe use a header field to detect that and convert it
  //       automatically.

  // Put attenuation flag on
  m_mumapAttenuationFlag = true;

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oSensitivityGenerator::Launch()
{
  // Verbose
  if (m_verbose>=1) Cout("oSensitivityGenerator::Launch() -> Start the sensitivity computation" << endl);
  // Simply switches between the CPU and GPU versions
  #ifdef CASTOR_GPU
  if (m_flagGPU) return LaunchGPU();
  else return LaunchCPU();
  #else
  return LaunchCPU();
  #endif
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oSensitivityGenerator::LaunchCPU()
{
  // Loop on the bed positions
  for (int bed=0; bed<mp_ImageDimensionsAndQuantification->GetNbBeds(); bed++)
  {
    // Reset line counters
    for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection(); th++) mp_lineCounter[th] = 0;
    // Apply the bed offset for this bed position
    mp_ProjectorManager->ApplyBedOffset(bed);
    // Case where the computation is performed directly from the histogram datafile
    if (m_computeFromHistogramFlag)
    {
      if (ComputeSensitivityFromHistogramDataFile(bed))
      {
        Cerr("***** oSensitivityGenerator::LaunchCPU() -> A problem occured while computing the sensitivity using the histogram data file for bed " << bed+1 << " !" << endl);
        return 1;
      }
    }
    // If a normalization file is provided, then use them to iterate on the elements to be taken into account
    else if (m3p_NormalizationDataFile!=NULL)
    {
      if (ComputeSensitivityFromNormalizationFile(bed))
      {
        Cerr("***** oSensitivityGenerator::LaunchCPU() -> A problem occured while computing the sensitivity using the normalization file for bed " << bed+1 << " !" << endl);
        return 1;
      }
    }
    // Otherwise, use a generic method from the scanner elements
    else
    {
      if (ComputeSensitivityFromScanner(bed))
      {
        Cerr("***** oSensitivityGenerator::LaunchCPU() -> A problem occured while computing the sensitivity from scanner elements for bed " << bed+1 << " !" << endl);
        return 1;
      }
      // Synchronize MPI processes here as the sensitivity computation from scanner elements has not MPI implementation yet
      #ifdef CASTOR_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      #endif
    }
    
    // Sum up the line counter into the first thread
    for (int th=1; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection(); th++) mp_lineCounter[0] += mp_lineCounter[th];
    // Verbose
    if (m_verbose>=2) Cout("  --> Number of effectively projected lines: " << mp_lineCounter[0] << endl);
  }
  
  
  // Set the number of threads for image computations
  #ifdef CASTOR_OMP
  omp_set_num_threads(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation());
  #endif

  // Reduce all results (merging parallel computation)
  for (int fr=0 ; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames() ; fr++)
  {
    //for (int rg=mp_ImageDimensionsAndQuantification->GetPMotionFirstIndexForLMS(fr) ; rg<mp_ImageDimensionsAndQuantification->GetNb1stMotImgsForLMS(fr) ; rg++)
    int frm_nb_1st_mot_imgs = mp_ImageDimensionsAndQuantification->GetNb1stMotImgsForLMS(fr)-mp_ImageDimensionsAndQuantification->GetPMotionFirstIndexForLMS(fr);
    for (int rg=0 ; rg<frm_nb_1st_mot_imgs ; rg++)
      for (int cg=0 ; cg<mp_ImageDimensionsAndQuantification->GetNb2ndMotImgsForLMS() ; cg++)
        {
          // Merge parallel results into the backward image
          int img_0 = 0;
          mp_ImageSpace->ReduceBackwardImage(img_0, fr, rg, cg);
          // Apply mask if provided
          mp_ImageSpace->ApplyMaskToBackwardImage(img_0, fr, rg, cg);
          // Perform any deformations on the backward image
          
          /*
          mp_DeformationManager->ApplyDeformationForSensitivityGeneration(mp_ImageSpace,
                                                                          BACKWARD_DEFORMATION,
                                                                          fr,
                                                                          rg + mp_ImageDimensionsAndQuantification->GetPMotionFirstIndexForLMS(fr),
                                                                          cg);*/
                                                                          
          mp_DeformationManager->ApplyDeformationForSensitivityGeneration(mp_ImageSpace,
                                                                          BACKWARD_DEFORMATION,
                                                                          rg + mp_ImageDimensionsAndQuantification->GetPMotionFirstIndexForLMS(fr),
                                                                          fr,
                                                                          rg,
                                                                          cg);
                                                                          
          //mp_DeformationManager->ApplyDeformationForSensitivityGeneration(mp_ImageSpace, BACKWARD_DEFORMATION, 1, fr, rg, cg);
          //mp_DeformationManager->ApplyDeformationForSensitivityGeneration(mp_ImageSpace, BACKWARD_DEFORMATION, 2, fr, rg, cg);
        }
  }
  
  // Allocate sensitivity image
  mp_ImageSpace->LMS_InstantiateSensitivityImage(); 
  // Copy backward image to sensitivity image
  mp_ImageSpace->LMS_CopyBackwardToSensitivity();
  // Deallocate backward image
  mp_ImageSpace->DeallocateBackwardImageFromDynamicBins();

  // Deallocate the forward image only if attenuation was used
  //if (m_mumapAttenuationFlag) mp_ImageSpace->LMS_DeallocateForwardImage();
  // Deallocate forward image
  mp_ImageSpace->LMS_DeallocateForwardImage();
  // Deallocate mask image
  mp_ImageSpace->DeallocateMaskImage();
  
  // Apply PSF and save only for first instance
  if (mp_ImageDimensionsAndQuantification->GetMPIRank()==0)
  {
    // Apply PSF if needed
    mp_ImageConvolverManager->ConvolveSensitivity(mp_ImageSpace);
    // Save sensitivity image
    mp_ImageSpace->LMS_SaveSensitivityImage(m_pathToSensitivityImage, mp_DeformationManager);
  }
  // Deallocate sensitivity image
  mp_ImageSpace->LMS_DeallocateSensitivityImage();
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oSensitivityGenerator::ComputeSensitivityFromHistogramDataFile(int a_bed)
{
  // For the moment, we restrict the validity of this function to non-gated data
  if (mp_ImageDimensionsAndQuantification->GetNb1stMotImgsForLMS(mp_ImageDimensionsAndQuantification->GetNbTimeFrames()-1)>1 || mp_ImageDimensionsAndQuantification->GetNb2ndMotImgsForLMS()>1)
  {
    Cerr("***** oSensitivityGenerator::ComputeSensitivityFromHistogramDataFile() -> Cannot compute sensitivity from histogram file with gating !" << endl);
    return 1;
  }

  // Verbose
  if (m_verbose>=2)
  {
    if (mp_ImageDimensionsAndQuantification->GetNbBeds()>1) Cout("oSensitivityGenerator::ComputeSensitivityFromHistogramDataFile() -> Start computation for bed " << a_bed+1 << endl);
    else Cout("oSensitivityGenerator::ComputeSensitivityFromHistogramDataFile() -> Start computation" << endl);
  }

  // Initial clock
  clock_t clock_start = clock();
  time_t time_start = time(NULL);

  // In case a problem occurs in the parallel loop
  bool problem = false;
  // Set the number of threads for projections
  #ifdef CASTOR_OMP
  omp_set_num_threads(mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection());
  #endif

  // Get index start and stop
  int64_t index_start = 0;
  int64_t index_stop = 0;
  m2p_DataFile[a_bed]->GetEventIndexStartAndStop(&index_start, &index_stop);

  // Reinitialize 4D gate indexes
  mp_ImageDimensionsAndQuantification->ResetCurrentDynamicIndices();

  // Launch the loop on all events
  int64_t index = 0, printing_index = 0;
  #pragma omp parallel for private(index) schedule(static, 1)
  for ( index = index_start  ;  index < index_stop  ;  index++)
  {
    // Get the thread index
    int th = 0;
    #ifdef CASTOR_OMP
    th = omp_get_thread_num();
    #endif
    // Verbose
    if (th==0 && m_verbose>=2 && mp_ImageDimensionsAndQuantification->GetMPIRank()==0)
    {
      if (printing_index%1000==0)
      {
        int percent = ((int)( ((FLTNB)(index-index_start)) * 100. / ((FLTNB)(index_stop-index_start)) ));
        cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
             << percent << " %                    " << flush;
      }
      printing_index++;
    }
    // Get the current event for that thread index
    vEvent* event = m2p_DataFile[a_bed]->GetEvent(index, th);
    if (event==NULL)
    {
      Cerr("***** oSensitivityGenerator::ComputeSensitivityFromHistogramDataFile() -> An error occured while getting the event from index "
           << index << " (thread " << th << ") !" << endl);
      // A problem was encountered
      problem = true;
      // We must continue here because we are inside an OpenMP loop
      continue;
    }
/*
    // Call the dynamic switch function that updates the current frame and gate numbers, and also detects involuntary patient motion
    int dynamic_switch_value = mp_ImageDimensionsAndQuantification->DynamicSwitch(index, event->GetTimeInMs(), a_bed, th);
    // If the DYNAMIC_SWITCH_CONTINUE is returned, then it means that we are not yet at the first frame
    if (dynamic_switch_value==DYNAMIC_SWITCH_CONTINUE)
    {
      // Then we just skip this event
      continue;
    }
*/
    // Compute the projection line
    oProjectionLine *line = mp_ProjectorManager->ComputeProjectionLine(event, th);
    if (line==NULL)
    {
      Cerr("***** oSensitivityGenerator::ComputeSensitivityFromHistogramDataFile() -> A problem occured while computing the projection line !" << endl);
      // Specify that there was a problem
      problem = true;
      // We must continue here because we are inside an OpenMP loop
      continue;
    }

    // Process this line
    int no_frame = 0;
    int no_resp_gate = 0;
    int no_card_gate = 0;
//    if (line->NotEmptyLine() && ProcessThisLine(line, event, a_bed, mp_ImageDimensionsAndQuantification->GetCurrentTimeFrame(th), no_resp_gate, no_card_gate, th))
    if (line->NotEmptyLine() && ProcessThisLine(line, event, a_bed, no_frame, no_resp_gate, no_card_gate, th))
    {
      Cerr("***** oSensitivityGenerator::ComputeSensitivityFromHistogramDataFile() -> A problem occured while processing a line !" << endl);
    }
  }

  // Synchronize MPI processes
  #ifdef CASTOR_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // End of progression printing (do not log out with Cout here)
  if (m_verbose>=2 && mp_ImageDimensionsAndQuantification->GetMPIRank()==0)
    cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
         << "  --> 100 %                        " << endl;
  // If a problem was encountered, then report it here
  if (problem)
  {
    Cerr("***** oSensitivityGenerator::ComputeSensitivityFromHistogramDataFile() -> A problem occured inside the parallel loop over events !" << endl);
    return 1;
  }

  // Clock total
  clock_t clock_stop = clock();
  time_t time_stop = time(NULL);
  if (mp_ImageDimensionsAndQuantification->GetNbBeds()>1)
    Cout("  --> Time spent for sensitivity generation for bed " << a_bed+1 << " | User: " << time_stop-time_start
         << " sec | CPU: " << (clock_stop-clock_start)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);
  else
    Cout("  --> Time spent for sensitivity generation | User: " << time_stop-time_start
         << " sec | CPU: " << (clock_stop-clock_start)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oSensitivityGenerator::ComputeSensitivityFromNormalizationFile(int a_bed)
{
  // Verbose
  if (m_verbose>=2)
  {
    if (mp_ImageDimensionsAndQuantification->GetNbBeds()>1) Cout("oSensitivityGenerator::ComputeSensitivityFromNormalizationFile() -> Start computation for bed " << a_bed+1 << endl);
    else Cout("oSensitivityGenerator::ComputeSensitivityFromNormalizationFile() -> Start computation" << endl);
  }

  // Initial clock
  clock_t clock_start = clock();
  time_t time_start = time(NULL);

  // Loop on frames
  for (int fr=0 ; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames() ; fr++)
  {
    // Verbose
    if (m_verbose>=2 && mp_ImageDimensionsAndQuantification->GetNbTimeFrames()>1)
      cout << "  --> Processing frame " << fr+1 << endl;
    // Loop on respiratory gates
    for (int rg=0 ; rg<mp_ImageDimensionsAndQuantification->GetNb1stMotImgsForLMS(fr) ; rg++)
    {
      // Verbose
      if (m_verbose>=2 && mp_ImageDimensionsAndQuantification->GetNb1stMotImgsForLMS(fr)>1)
        cout << "  --> Processing respiratory gate " << rg+1 << endl;
      // Loop on cardiac gates
      for (int cg=0 ; cg<mp_ImageDimensionsAndQuantification->GetNb2ndMotImgsForLMS() ; cg++)
      {
        // Verbose
        if (m_verbose>=2 && mp_ImageDimensionsAndQuantification->GetNb2ndMotImgsForLMS()>1)
          cout << "  --> Processing cardiac gate " << cg+1 << endl;
        
        // Perform any forward deformations on the forward image
        //mp_DeformationManager->ApplyDeformationForSensitivityGeneration(mp_ImageSpace, FORWARD_DEFORMATION, fr, rg, cg);
        mp_DeformationManager->ApplyDeformationForSensitivityGeneration(mp_ImageSpace,
                                                                        FORWARD_DEFORMATION,
                                                                        mp_ImageDimensionsAndQuantification->GetPMotionFirstIndexForLMS(fr),
                                                                        fr,
                                                                        rg,
                                                                        cg);

        // In case a problem occurs in the parallel loop
        bool problem = false;
        // Set the number of threads for projections
        #ifdef CASTOR_OMP
        omp_set_num_threads(mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection());
        #endif

        // Get index start and stop
        int64_t index_start = 0;
        int64_t index_stop = 0;
        m3p_NormalizationDataFile[a_bed][fr]->GetEventIndexStartAndStop(&index_start, &index_stop);

        // Index for progression printing
        uint64_t progression_index_total = (index_stop-index_start)
                                         * mp_ImageDimensionsAndQuantification->GetNbTimeFrames()
                                         * mp_ImageDimensionsAndQuantification->GetNb1stMotImgsForLMS(fr)
                                         * mp_ImageDimensionsAndQuantification->GetNb2ndMotImgsForLMS();

        // Launch the loop on all events
        int64_t index = 0, printing_index = 0;
        #pragma omp parallel for private(index) schedule(static, 1)
        for ( index = index_start  ;  index < index_stop  ;  index++)
        {
          // Get the thread index
          int th = 0;
          #ifdef CASTOR_OMP
          th = omp_get_thread_num();
          #endif
          
          // Verbose
          /*
          if (th==0 && m_verbose>=2 && mp_ImageDimensionsAndQuantification->GetMPIRank()==0)
          {
            if (printing_index%1000==0)
            {
              int percent = ((int)( ((FLTNB)(index-index_start)) * 100. / ((FLTNB)(index_stop-index_start)) ));
              cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
                   << percent << " %                    " << flush;
            }
            printing_index++;
          }
          */
          
          if (th==0 && m_verbose>=2 && mp_ImageDimensionsAndQuantification->GetMPIRank()==0)
          {
            if (printing_index%1000==0)
            {
              int64_t progression_chunk = (int64_t)((FLTNB)(index_stop-index_start));
              int64_t progression_index_current = fr * mp_ImageDimensionsAndQuantification->GetNb1stMotImgsForLMS(fr) * mp_ImageDimensionsAndQuantification->GetNb2ndMotImgsForLMS() * progression_chunk
                                                + rg * mp_ImageDimensionsAndQuantification->GetNb2ndMotImgsForLMS() * progression_chunk
                                                + cg * progression_chunk
                                                + index;

              int64_t percent = (int64_t) (( ((FLTNB)progression_index_current)/((FLTNB)progression_index_total) ) * ((FLTNB)100));
              cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
                   << percent << " %                    " << flush;
            }
            printing_index++;
          }

          // Get the current event for that thread index
          vEvent* event = m3p_NormalizationDataFile[a_bed][fr]->GetEvent(index, th);

          if (event==NULL)
          {
            Cerr("***** oSensitivityGenerator::ComputeSensitivityFromNormalizationFile() -> An error occured while getting the event from index "
                 << index << " (thread " << th << ") !" << endl);
            // A problem was encountered
            problem = true;
            // We must continue here because we are inside an OpenMP loop
            continue;
          }
          // Compute the projection line
          oProjectionLine *line = mp_ProjectorManager->ComputeProjectionLine(event, th);
          if (line==NULL)
          {
            Cerr("***** oSensitivityGenerator::ComputeSensitivityFromNormalizationFile() -> A problem occured while computing the projection line !" << endl);
            // Specify that there was a problem
            problem = true;
            // We must continue here because we are inside an OpenMP loop
            continue;
          }
          // Process this line
          if (line->NotEmptyLine() && ProcessThisLine(line, event, a_bed, fr, rg, cg, th))
          {
            Cerr("***** oSensitivityGenerator::ComputeSensitivityFromNormalizationFile() -> A problem occured while processing a line !" << endl);
          }
        }
        // End of progression printing (do not log out with Cout here)
        if (m_verbose>=2 && mp_ImageDimensionsAndQuantification->GetMPIRank()==0)
          cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
               << "  --> 100 %                        " << endl;
        // If a problem was encountered, then report it here
        if (problem)
        {
          Cerr("***** oSensitivityGenerator::ComputeSensitivityFromNormalizationFile() -> A problem occured inside the parallel loop over events !" << endl);
          return 1;
        }
      }
    }
  }

  // Clock total
  clock_t clock_stop = clock();
  time_t time_stop = time(NULL);
  if (mp_ImageDimensionsAndQuantification->GetNbBeds()>1)
    Cout("  --> Time spent for sensitivity generation for bed " << a_bed+1 << " | User: " << time_stop-time_start
         << " sec | CPU: " << (clock_stop-clock_start)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);
  else
    Cout("  --> Time spent for sensitivity generation | User: " << time_stop-time_start
         << " sec | CPU: " << (clock_stop-clock_start)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oSensitivityGenerator::ComputeSensitivityFromScanner(int a_bed)
{
  // No MPI implementation for this function yet, so we exit if not the first instance
  if (mp_ImageDimensionsAndQuantification->GetMPIRank()!=0) return 0;

  // Verbose
  if (m_verbose>=2)
  {
    if (mp_ImageDimensionsAndQuantification->GetNbBeds()>1) Cout("oSensitivityGenerator::ComputeSensitivityFromScanner() -> Start computation for bed " << a_bed+1 << endl);
    else Cout("oSensitivityGenerator::ComputeSensitivityFromScanner() -> Start computation" << endl);
  }

  // Get the scanner manager
  sScannerManager* p_scannerManager = sScannerManager::GetInstance();

  // Total number of elements
  int nb_total_elts = mp_Scanner->GetSystemNbElts();

  // Initial clock
  clock_t clock_start = clock();
  time_t time_start = time(NULL);

  // Initialize main loop start and stop values
  int64_t main_loop_start_index = 0 ;
  int64_t main_loop_stop_index = 0;

  // Check beforehand any issue with the loop start/stop values (not possible in the inner multithreaded loop)
  if (p_scannerManager->PROJ_GetModalityStopValueMainLoop()>0) main_loop_stop_index = p_scannerManager->PROJ_GetModalityStopValueMainLoop();
  else
  {
    Cerr("***** oSensitivityGenerator::ComputeSensitivityFromScanner() -> An error occured when trying to initialize main loop stop index !" << endl); 
    return 1;
  }
  if (p_scannerManager->PROJ_GetModalityStartValueInnerLoop(0)<=0)
  {
    Cerr("***** oSensitivityGenerator::ComputeSensitivityFromScanner() -> An error occured when trying to initialize inner loop start index !" << endl); 
    return 1;
  }

  // Prepare pre-computed sums of events to avoid the exponential evolution of the percentage (for PET)
  // TODO Perhaps replace this with a call to scannerManager for error management
  int64_t* progression_elts_array = new int64_t[nb_total_elts*sizeof(int64_t)];
  progression_elts_array[0] = 0;
  for (int idx_elt=1 ; idx_elt<mp_Scanner->GetSystemNbElts() ; idx_elt++) 
    progression_elts_array[idx_elt] = progression_elts_array[idx_elt-1] + (mp_Scanner->GetSystemNbElts()-idx_elt);

  // Index for progression printing
  uint64_t printing_index = 0;
  uint64_t progression_index_total = p_scannerManager->PROJ_GetProgressionFinalValue()
                                   * mp_ImageDimensionsAndQuantification->GetNb1stMotImgsForLMS(0)
                                   * mp_ImageDimensionsAndQuantification->GetNb2ndMotImgsForLMS();
  
  for (int fr=1 ; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames() ; fr++)
    progression_index_total += p_scannerManager->PROJ_GetProgressionFinalValue()
                            * mp_ImageDimensionsAndQuantification->GetNb1stMotImgsForLMS(fr)
                            * mp_ImageDimensionsAndQuantification->GetNb2ndMotImgsForLMS();
  
  //mp_ImageSpace->LMS_InstantiateForwardImage();
  // Loop on frames
  for (int fr=0 ; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames() ; fr++)
  {
    // Verbose
    if (m_verbose>=2 && mp_ImageDimensionsAndQuantification->GetNbTimeFrames()>1)
      cout << "  --> Processing frame " << fr+1 << endl;
    // Loop on respiratory gates
    //for (int rg=mp_ImageDimensionsAndQuantification->GetPMotionFirstIndexForLMS(fr) ; rg<mp_ImageDimensionsAndQuantification->GetNb1stMotImgsForLMS(fr) ; rg++)
    int frm_nb_1st_mot_imgs = mp_ImageDimensionsAndQuantification->GetNb1stMotImgsForLMS(fr)-mp_ImageDimensionsAndQuantification->GetPMotionFirstIndexForLMS(fr);
    for (int rg=0 ; rg<frm_nb_1st_mot_imgs ; rg++)
    {
      // Verbose
      if (m_verbose>=2 && mp_ImageDimensionsAndQuantification->GetNb1stMotImgsForLMS(fr)>1)
        cout << "  --> Processing respiratory gate " << rg+1 << endl;
      // Loop on cardiac gates
      for (int cg=0 ; cg<mp_ImageDimensionsAndQuantification->GetNb2ndMotImgsForLMS() ; cg++)
      {
        // Verbose
        if (m_verbose>=2 && mp_ImageDimensionsAndQuantification->GetNb2ndMotImgsForLMS()>1)
          cout << "  --> Processing cardiac gate " << cg+1 << endl;
            
        // Loop for involuntary patient motion
        //for (int im=0 ; im<mp_ImageDimensionsAndQuantification->GetNbIPatMotionSubsets() ; im++)
        //{
          // Verbose
        //  if (m_verbose>=2 && mp_ImageDimensionsAndQuantification->GetNbIPatMotionSubsets()>1)
        //    cout << "  --> Processing patient motion data subset " << im+1 << endl;

          // Perform any forward deformations on the forward image (only if attenuation is used)
          //if (m_mumapAttenuationFlag)
          
          /*
          mp_DeformationManager->ApplyDeformationForSensitivityGeneration(mp_ImageSpace,
                                                                          FORWARD_DEFORMATION,
                                                                          fr,
                                                                          rg + mp_ImageDimensionsAndQuantification->GetPMotionFirstIndexForLMS(fr),
                                                                          cg);*/
                                                                          
          mp_DeformationManager->ApplyDeformationForSensitivityGeneration(mp_ImageSpace,
                                                                          FORWARD_DEFORMATION,
                                                                          mp_ImageDimensionsAndQuantification->GetPMotionFirstIndexForLMS(fr),
                                                                          fr,
                                                                          rg,
                                                                          cg);

          // In case a problem occurs in the parallel loop
          bool problem = false;
          // Set the number of threads for projections
          #ifdef CASTOR_OMP
          omp_set_num_threads(mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection());
          #endif
  
          // Start the loop
          int64_t idx_elt1;
  
          #pragma omp parallel for private(idx_elt1) schedule(static, 1)
          for(idx_elt1=main_loop_start_index ; idx_elt1<main_loop_stop_index ; idx_elt1++)
          {
            // Get the thread index
            int th = 0;
            #ifdef CASTOR_OMP
            th = omp_get_thread_num();
            #endif
  
            // Initialize inner loop start and stop values
            int64_t inner_loop_start_index = p_scannerManager->PROJ_GetModalityStartValueInnerLoop(idx_elt1) ;
            int64_t inner_loop_stop_index = mp_Scanner->GetSystemNbElts();
            
            // Inner loop on scanner elements (idx_elt2) which are used to form LOR using the first scanner element (idx_elt1) 
            for (int64_t idx_elt2=inner_loop_start_index ; idx_elt2<inner_loop_stop_index ; idx_elt2++)
            {
              // Print progression (do not log out with Cout here)
              if (th==0)
              {
                if (m_verbose>=2 && printing_index%10000==0)
                {
                 int64_t progression_index_current = p_scannerManager->PROJ_GetCurrentProgression(idx_elt1, idx_elt2, progression_elts_array, 
                                                                                                    mp_ImageDimensionsAndQuantification->GetNb1stMotImgsForLMS(fr), 
                                                                                                    mp_ImageDimensionsAndQuantification->GetNb2ndMotImgsForLMS(), 
                                                                                                                                  fr, rg, cg);
                  int64_t percent = (int64_t) (( ((FLTNB)progression_index_current)/((FLTNB)progression_index_total) ) * ((FLTNB)100));
                  cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
                     << percent << " %                       ";
                }
                printing_index++;
              }
  
              // Allocate an event using the iDataFile
              vEvent* event = m2p_DataFile[a_bed]->PROJ_GenerateEvent(idx_elt1, idx_elt2, th);
  
              // Check from the scanner requirements that this LOR is allowed
              if (!mp_Scanner->IsAvailableLOR(event->GetID1(0), event->GetID2(0))) continue;
  
              // Generate the projection event and compute the projection line 
              oProjectionLine* line = mp_ProjectorManager->ComputeProjectionLine(event, th);
              if (line==NULL)
              {
                Cerr("***** oSensitivityGenerator::ComputeSensitivityFromScanner() -> A problem occured while computing the projection line !" << endl);
                // Specify that there was a problem
                problem = true;
                // We must continue here because we are inside an OpenMP loop
                continue;
              }
  
              // Process this line
              if (line->NotEmptyLine() && ProcessThisLine(line, event, a_bed, fr, rg, cg, th))
              {
                Cerr("***** oSensitivityGenerator::ComputeSensitivityFromScanner() -> A problem occured while processing a line !" << endl);
              }
            }
          }
          // End of progression printing (do not log out with Cout here)
          if (m_verbose>=2 && mp_ImageDimensionsAndQuantification->GetMPIRank()==0)
            cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                 << "  --> 100 %                       " << endl;
          // If a problem was encountered, then report it here
          if (problem)
          {
            Cerr("***** oSensitivityGenerator::ComputeSensitivityFromScanner() -> A problem occured inside the parallel loop over events !" << endl);
            return 1;
          }
        //}
      }
    }
  }

  delete[] progression_elts_array;

  // Clock total
  clock_t clock_stop = clock();
  time_t time_stop = time(NULL);
  if (mp_ImageDimensionsAndQuantification->GetNbBeds()>1)
    Cout("  --> Time spent for sensitivity generation for bed " << a_bed+1 << " | User: " << time_stop-time_start
         << " sec | CPU: " << (clock_stop-clock_start)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);
  else
    Cout("  --> Time spent for sensitivity generation | User: " << time_stop-time_start
         << " sec | CPU: " << (clock_stop-clock_start)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oSensitivityGenerator::ProcessThisLine(oProjectionLine* ap_Line, vEvent* ap_Event, int a_bed, int a_frame, int a_respGate, int a_cardGate, int a_thread)
{
  // ----------------------------------------------------------------------------------------------
  // Part 0: get the multiplicative corrections from the vEvent and the quantification factor
  // ----------------------------------------------------------------------------------------------

  // The multiplicative corrections from the vEvent can be:
  //   1: generated from vDataFile::GenerateEvent(), this will be a default value, so this will be 1.
  //   2: taken from the normalization data file, so this will be the normalization correction factor times the attenuation correction factor
  FLTNB multiplicative_factor = ap_Event->GetMultiplicativeCorrections();
  // If null, then skip this event
  if (multiplicative_factor<=0.) return 0;
  // Get the quantification factor
  FLTNB quantification_factor = mp_ImageDimensionsAndQuantification->GetQuantificationFactor(a_bed, a_frame, a_respGate, a_cardGate);
  // If null, then skip this event
  if (quantification_factor<=0.) return 0;

  // ----------------------------------------------------------------------------------------------
  // Part 1: deal with the attenuation forward projection for PET if any
  // ----------------------------------------------------------------------------------------------

  // The ACF
  FLTNB attenuation_correction_factor = 1.;

  // Do a forward projection only if the attenuation image was provided and the ACF was not included in the normalization data file
  if ( m_forwardProjectAttenuation )
  {
    // Cumulative mu (in cm-1)
    FLTNB cumulative_mu = 0.;
    // Loop on TOF bins
    for (int b=0; b<ap_Line->GetNbTOFBins(); b++)
    {
      // Projection operation (the attenuation image in the m4p_forwardImage buffer of the oImageSpace)
      for (int vl=0 ; vl<ap_Line->GetCurrentNbVoxels(FORWARD, b) ; vl++)
        cumulative_mu += ap_Line->GetVoxelWeights(FORWARD, b, vl) * mp_ImageSpace->m4p_forwardImage[a_frame][a_respGate][a_cardGate][ap_Line->GetVoxelIndex(FORWARD, b,vl)];
    }
    // Update the ACF (converting the cumulative mu in mm-1)
    attenuation_correction_factor = max(((FLTNB)1.),exp(cumulative_mu*((FLTNB)0.1)));
  }

  // ----------------------------------------------------------------------------------------------
  // Part 2: back-projection of the LOR sensitivity
  // ----------------------------------------------------------------------------------------------

  // Compute the global sensitivity (calibration / (ACF * norm))
  FLTNB lor_sensitivity = 1. / (quantification_factor * attenuation_correction_factor * multiplicative_factor);

  // Loop on TOF bins
  for (int b=0; b<ap_Line->GetNbTOFBins(); b++)
  {
    // Backprojection operation into the backward image
    for (int vl=0 ; vl<ap_Line->GetCurrentNbVoxels(BACKWARD, b) ; vl++)
      mp_ImageSpace->m6p_backwardImage[0][a_thread][a_frame][a_respGate][a_cardGate][ap_Line->GetVoxelIndex(BACKWARD, b, vl)] += ap_Line->GetVoxelWeights(BACKWARD, b, vl)
                                                                                                                               * lor_sensitivity;
  }

  // Increment line counters
  mp_lineCounter[a_thread]++;

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

string oSensitivityGenerator::GetPathToSensitivityImage() 
{
  // Get the flag that says if we merge the dynamic files or not
  bool merge_dynamic_files = sOutputManager::GetInstance()->MergeDynImages();
  // If we merge dynamic file, or we are in a static case (no dynamic), we stay with the hdr extension
   if ( ( mp_ImageDimensionsAndQuantification->GetNbCardGates()==1 &&
          mp_ImageDimensionsAndQuantification->GetNbRespGates()==1 &&
          mp_ImageDimensionsAndQuantification->GetNbTimeFrames()==1 ) ||
       merge_dynamic_files )
    return m_pathToSensitivityImage + ".hdr";
  else
    return m_pathToSensitivityImage + ".mhd";
}
