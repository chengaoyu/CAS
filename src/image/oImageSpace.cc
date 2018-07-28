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
  \ingroup  image
  \brief    Implementation of class oImageSpace
*/

#include "oImageSpace.hh"
#include "vOptimizer.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oImageSpace::oImageSpace()
{
  // Initialize the member variables to their default values
  m_loadedSensitivity = false;
  m_loadedMultiModal = false;
  m_loadedMask = false;
  mp_ID = NULL;
  m_verbose = -1;
  m_nbBackwardImages = 1;
  m_nbMiscellaneousImages = 0;
  m2p_miscellaneousImage = NULL;
  // Set all images pointers to NULL
  m4p_image = NULL;
  m4p_forwardImage = NULL;
  m6p_backwardImage = NULL;
  m5p_sensitivity = NULL;
  m4p_attenuation = NULL;
  m5p_refDynBackwardImage = NULL;
  m4p_refDynSensitivityImage = NULL;
  m4p_outputImage = NULL;
  //mp_visitedVoxelsImage = NULL;
  m4p_visitedVoxelsImage = NULL;
  m2p_projectionImage = NULL;
  m5p_multiModalImage = NULL;
  mp_maskImage = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oImageSpace::~oImageSpace()
{
  // Call the deallocation of the miscellaneous images here
  DeallocateMiscellaneousImage();
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageSpace::InstantiateImage()
{
  // Verbose
  if (m_verbose>=3) Cout("oImageSpace::InstantiateImage() -> Initialize to 0."<< endl); 
    
  // For each temporal dimensions, the number of related images is recovered from the number of
  // intrinsic time basis functions stored in mp_ID. If no intrinsic temporal basis functions
  // are initialized (standard reconstruction), this is equivalent to the standard calls to
  // GetNbTimeFrames() / GetNbRespGates() / GetNbCardGates()

  // First pointer is for the number of time basis functions
  m4p_image = (FLTNB****)malloc(mp_ID->GetNbTimeBasisFunctions()*sizeof(FLTNB***));
  for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
  {
    // Second pointer is for the number of respiratory basis functions
    m4p_image[tbf] = (FLTNB***)malloc(mp_ID->GetNbRespBasisFunctions()*sizeof(FLTNB**));
    for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
    {
      // Third pointer is for the number of cardiac basis functions
      m4p_image[tbf][rbf] = (FLTNB**)malloc(mp_ID->GetNbCardBasisFunctions()*sizeof(FLTNB*));
      for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
      {
        // Fourth pointer is for the 3D space
        m4p_image[tbf][rbf][cbf] = (FLTNB*)malloc(mp_ID->GetNbVoxXYZ()*sizeof(FLTNB));
        // Initialize to 0.
        for (INTNB v=0; v<mp_ID->GetNbVoxXYZ(); v++) m4p_image[tbf][rbf][cbf][v] = 0.;
      }
    }
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageSpace::DeallocateImage()
{
  // Check the first pointer
  if (m4p_image)
  {
    // Verbose
    if (m_verbose>=3) Cout("oImageSpace::DeallocateImage() -> Free memory"<< endl); 
    // Loop on time basis functions
    for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++) if (m4p_image[tbf])
    {
      // Loop on respiratory basis functions
      for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++) if (m4p_image[tbf][rbf])
      {
        // Loop on cardiac basis functions
        for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++) if (m4p_image[tbf][rbf][cbf])
        {
          free(m4p_image[tbf][rbf][cbf]);
        }
        free(m4p_image[tbf][rbf]);
      }
      free(m4p_image[tbf]);
    }
    free(m4p_image);
  }
  // Reset the pointer to NULL
  m4p_image = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageSpace::InstantiateForwardImage()
{
  // Verbose
  if (m_verbose>=3) Cout("oImageSpace::InstantiateForwardImage() -> Initialize to 0."<< endl); 
    
  // For each temporal dimensions, the number of related images is recovered from the number of
  // intrinsic time basis functions stored in mp_ID. If no intrinsic temporal basis functions
  // are initialized (standard reconstruction), this is equivalent to the standard calls to
  // GetNbTimeFrames() / GetNbRespGates() / GetNbCardGates()
  
  // First pointer is for the number of time basis functions
  m4p_forwardImage = (FLTNB****)malloc(mp_ID->GetNbTimeBasisFunctions()*sizeof(FLTNB***));
  for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
  {
    // Second pointer is for the number of respiratory basis functions
    m4p_forwardImage[tbf] = (FLTNB***)malloc(mp_ID->GetNbRespBasisFunctions()*sizeof(FLTNB**));
    for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
    {
      // Third pointer is for the number of cardiac basis functions
      m4p_forwardImage[tbf][rbf] = (FLTNB**)malloc(mp_ID->GetNbCardBasisFunctions()*sizeof(FLTNB*));
      for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
      {
        // Fourth pointer is for the 3D space
        m4p_forwardImage[tbf][rbf][cbf] = (FLTNB*)malloc(mp_ID->GetNbVoxXYZ()*sizeof(FLTNB));
        // Initialize to 0.
        for (INTNB v=0; v<mp_ID->GetNbVoxXYZ(); v++) m4p_forwardImage[tbf][rbf][cbf][v] = 0.;
      }
    }
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageSpace::DeallocateForwardImage()
{
  // Check the first pointer
  if (m4p_forwardImage)
  {
    // Verbose
    if (m_verbose>=3) Cout("oImageSpace::DeallocateForwardImage() -> Free memory"<< endl);
    // Loop on time basis functions
    for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++) if (m4p_forwardImage[tbf])
    {
      // Loop on respiratory basis functions
      for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++) if (m4p_forwardImage[tbf][rbf])
      {
        // Loop on cardiac basis functions
        for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++) if (m4p_forwardImage[tbf][rbf][cbf])
        {
          free(m4p_forwardImage[tbf][rbf][cbf]);
        }
        free(m4p_forwardImage[tbf][rbf]);
      }
      free(m4p_forwardImage[tbf]);
    }
    free(m4p_forwardImage);
  }
  
  // Reset the pointer to NULL
  m4p_forwardImage = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageSpace::InstantiateBackwardImageFromDynamicBasis(int a_nbBackwardImages)
{
  // Verbose
  if (m_verbose>=3) Cout("oImageSpace::InstantiateBackwardImageFromDynamicBasis() -> Initialize to 0."<< endl);
    
  // For each temporal dimensions, the number of related images is recovered from the number of
  // intrinsic time basis functions stored in mp_ID. If no intrinsic temporal basis functions
  // are initialized (standard reconstruction), this is equivalent to the standard calls to
  // GetNbTimeFrames() / GetNbRespGates() / GetNbCardGates()

  // First pointer is for the number of backward images (in case the optimizer does need more than 1)
  m_nbBackwardImages = a_nbBackwardImages;
  m6p_backwardImage = (FLTNB******)malloc(m_nbBackwardImages*sizeof(FLTNB*****));
  for (int img=0; img<m_nbBackwardImages; img++)
  {
    // Second pointer is for the number of threads, in order to have thread-safe backward projections
    m6p_backwardImage[img] = (FLTNB*****)malloc(mp_ID->GetNbThreadsForProjection()*sizeof(FLTNB****));
    for (int th=0; th<mp_ID->GetNbThreadsForProjection(); th++)
    {
      // The third pointer is for the number of time basis functions
      m6p_backwardImage[img][th] = (FLTNB****)malloc(mp_ID->GetNbTimeBasisFunctions()*sizeof(FLTNB***));
      for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
      {
        // The fourth pointer is for the number of respiratory basis functions
        m6p_backwardImage[img][th][tbf] = (FLTNB***)malloc(mp_ID->GetNbRespBasisFunctions()*sizeof(FLTNB**));
        for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
        {
          // The fifth pointer is for the number of cardiac basis functions
          m6p_backwardImage[img][th][tbf][rbf] = (FLTNB**)malloc(mp_ID->GetNbCardBasisFunctions()*sizeof(FLTNB*));
          for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
          {
            // The sixth pointer is for the 3D space
            m6p_backwardImage[img][th][tbf][rbf][cbf] = (FLTNB*)malloc(mp_ID->GetNbVoxXYZ()*sizeof(FLTNB));
          }
        }
      }
    }
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageSpace::DeallocateBackwardImageFromDynamicBasis()
{
  // Check the first pointer
  if (m6p_backwardImage)
  {
    // Verbose
    if (m_verbose>=3) Cout("oImageSpace::DeallocateBackwardImageFromDynamicBasis() -> Free memory" << endl);
    // Loop on backward images (in case the optimizer does need more than 1)
    for (int img=0; img<m_nbBackwardImages; img++) if (m6p_backwardImage[img]!=NULL)
    {
      // Loop on threads
      for (int th=0 ; th<mp_ID->GetNbThreadsForProjection() ; th++) if (m6p_backwardImage[img][th])
      {
        // Loop on time basis functions
        for (int tbf=0 ; tbf<mp_ID->GetNbTimeBasisFunctions() ; tbf++) if (m6p_backwardImage[img][th][tbf])
        {
          // Loop on respiratory basis functions
          for (int rbf=0 ; rbf<mp_ID->GetNbRespBasisFunctions() ; rbf++) if (m6p_backwardImage[img][th][tbf][rbf])
          {
            // Loop on cardiac basis functions
            for (int cbf=0 ; cbf<mp_ID->GetNbCardBasisFunctions() ; cbf++) if (m6p_backwardImage[img][th][tbf][rbf][cbf])
            {
              free(m6p_backwardImage[img][th][tbf][rbf][cbf]);
            }
            free(m6p_backwardImage[img][th][tbf][rbf]);
          }
          free(m6p_backwardImage[img][th][tbf]);
        }
        free(m6p_backwardImage[img][th]);
      }
      free(m6p_backwardImage[img]);
    }
    free(m6p_backwardImage);
  }
  // Reset the pointer to NULL
  m6p_backwardImage = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageSpace::InstantiateBackwardImageFromDynamicBins()
{
  // Verbose
  if (m_verbose>=3) Cout("oImageSpace::InstantiateBackwardImageFromDynamicBins() -> Initialize to 0."<< endl);
  // Not related to an optimizer, so we consider only one image
  m6p_backwardImage = (FLTNB******)malloc(1*sizeof(FLTNB*****));
  // Projection threads
  m6p_backwardImage[0] = (FLTNB*****)malloc(mp_ID->GetNbThreadsForProjection()*sizeof(FLTNB****));
  for (int th=0 ; th<mp_ID->GetNbThreadsForProjection() ; th++)
  {
    // Time frames
    m6p_backwardImage[0][th] = (FLTNB****)malloc(mp_ID->GetNbTimeFrames()*sizeof(FLTNB***));
    for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
    {
      // Respiratory gates
      m6p_backwardImage[0][th][fr] = (FLTNB***)malloc(mp_ID->GetNb1stMotImgsForLMS(fr)*sizeof(FLTNB**));
      for (int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr) ; rg++)
      {
        // Cardiac gates
        m6p_backwardImage[0][th][fr][rg] = (FLTNB**)malloc(mp_ID->GetNb2ndMotImgsForLMS()*sizeof(FLTNB*));
        for (int cg=0 ; cg<mp_ID->GetNb2ndMotImgsForLMS() ; cg++)
        {
          // Number of voxels
          m6p_backwardImage[0][th][fr][rg][cg] = (FLTNB*)malloc(mp_ID->GetNbVoxXYZ()*sizeof(FLTNB));
          for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++) m6p_backwardImage[0][th][fr][rg][cg][v] = 0.;
        } 
      }   
    }
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageSpace::DeallocateBackwardImageFromDynamicBins()
{
  // Check the first pointer
  if (m6p_backwardImage)
  {
    // Verbose
    if (m_verbose>=3) Cout("oImageSpace::DeallocateBackwardImageFromDynamicBins() -> Free memory"<< endl);
    // We assumed only one image when considering dynamic bins
    if (m6p_backwardImage[0]!=NULL)
    {
      // Loop on threads
      for (int th=0 ; th<mp_ID->GetNbThreadsForProjection() ; th++) if (m6p_backwardImage[0][th])
      {
        // Loop on time frames
        for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++) if (m6p_backwardImage[0][th][fr])
        {
          // Loop on respiratory gates
          for (int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr) ; rg++) if (m6p_backwardImage[0][th][fr][rg])
          {
            // Loop on cardiac gates
            for (int cg=0 ; cg<mp_ID->GetNb2ndMotImgsForLMS() ; cg++) if (m6p_backwardImage[0][th][fr][rg][cg])
            {
              free(m6p_backwardImage[0][th][fr][rg][cg]);
            }
            free(m6p_backwardImage[0][th][fr][rg]);
          }
          free(m6p_backwardImage[0][th][fr]);
        }
        free(m6p_backwardImage[0][th]);
      }
      free(m6p_backwardImage[0]);
    }
    free(m6p_backwardImage);
  }
  // Reset the pointer to NULL
  m6p_backwardImage = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageSpace::InstantiateSensitivityImage(const string& a_pathToSensitivityImage)
{
  // -----------------------------------------------------------------------------------------------------
  // Case 1: a path to a sensitivity image is provided, meaning that we will not compute the sensitivity
  //         on-the-fly. In that case, we do not take the multi-threading into account. Also note that in
  //         histogram mode, the same sensitivity image will be used in optimization algorithms for all
  //         subsets.
  // -----------------------------------------------------------------------------------------------------

  if (!a_pathToSensitivityImage.empty())
  {
    // Verbose
    if (m_verbose>=3) Cout("oImageSpace::InstantiateSensitivityImage() -> Will be loaded from '" << a_pathToSensitivityImage << "'" << endl);
    // Allocate for all threads first
    m5p_sensitivity = (FLTNB*****)malloc(mp_ID->GetNbThreadsForProjection()*sizeof(FLTNB****));
    // Allocate the rest only for the first thread
    int thread_0 = 0;
    m5p_sensitivity[thread_0] = (FLTNB****)malloc(mp_ID->GetNbTimeFrames()*sizeof(FLTNB***));
    for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
    {
      m5p_sensitivity[thread_0][fr] = (FLTNB***)malloc(mp_ID->GetNbRespGates()*sizeof(FLTNB**));
      for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
      {
        m5p_sensitivity[thread_0][fr][rg] = (FLTNB**)malloc(mp_ID->GetNbCardGates()*sizeof(FLTNB*));
        for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
          m5p_sensitivity[thread_0][fr][rg][cg] = (FLTNB*)malloc(mp_ID->GetNbVoxXYZ()*sizeof(FLTNB));
      }
    }
    // Make all thread pointers pointing to the first
    for (int th=1; th<mp_ID->GetNbThreadsForProjection(); th++)
      m5p_sensitivity[th] = m5p_sensitivity[thread_0];
    // Set the flag for loaded sensitivity to true
    m_loadedSensitivity = true;
  }

  // -----------------------------------------------------------------------------------------------------
  // Case 2: no sensitivity provided, meaning we will commpute it on-the-fly. We thus need it to be
  //         thread-safe, so we allocate for all threads.
  // -----------------------------------------------------------------------------------------------------

  else 
  {
    // Verbose
    if (m_verbose>=3) Cout("oImageSpace::InstantiateSensitivityImage() -> For all threads"<< endl);
    // Allocate for all threads first
    m5p_sensitivity = (FLTNB*****)malloc(mp_ID->GetNbThreadsForProjection()*sizeof(FLTNB****));
    // Allocate the rest only for all threads
    for (int th=0; th<mp_ID->GetNbThreadsForProjection(); th++)
    {
      m5p_sensitivity[th] = (FLTNB****)malloc(mp_ID->GetNbTimeFrames()*sizeof(FLTNB***));
      for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
      {
        m5p_sensitivity[th][fr] = (FLTNB***)malloc(mp_ID->GetNbRespGates()*sizeof(FLTNB**));
        for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
        {
          m5p_sensitivity[th][fr][rg] = (FLTNB**)malloc(mp_ID->GetNbCardGates()*sizeof(FLTNB*));
          for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
          {
            m5p_sensitivity[th][fr][rg][cg] = (FLTNB*)malloc(mp_ID->GetNbVoxXYZ()*sizeof(FLTNB));
            // Initialize to 0.
            for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++) m5p_sensitivity[th][fr][rg][cg][v] = 0.;
          }
        }
      }
    }
    // Set the flag for loaded sensitivity to false
    m_loadedSensitivity = false;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageSpace::DeallocateSensitivityImage()
{
  // Check the first pointer
  if (m5p_sensitivity)
  {
    // Verbose
    if (m_verbose>=3) Cout("oImageSpace::DeallocateSensitivityImage() -> Free memory"<< endl);
    // If loaded sensitivity then fully deallocate only the first thread pointer
    if (m_loadedSensitivity)
    {
      int thread_0 = 0;
      if (m5p_sensitivity[thread_0])
      {
        // Loop on time frames
        for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++) if (m5p_sensitivity[thread_0][fr])
        {
          // Loop on respiratory gates
          for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++) if (m5p_sensitivity[thread_0][fr][rg])
          {
            // Loop on cardiac gates
            for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++) if (m5p_sensitivity[thread_0][fr][rg][cg])
            {
              free(m5p_sensitivity[thread_0][fr][rg][cg]);
            }
            free(m5p_sensitivity[thread_0][fr][rg]);
          }
          free(m5p_sensitivity[thread_0][fr]);
        }
        free(m5p_sensitivity[thread_0]);
      }
    }
    // Otherwise, loop on all threads
    else
    {
      for (int th=0; th<mp_ID->GetNbThreadsForProjection(); th++) if (m5p_sensitivity[th])
      {
        // Loop on time frames
        for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++) if (m5p_sensitivity[th][fr])
        {
          // Loop on respiratory gates
          for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++) if (m5p_sensitivity[th][fr][rg])
          {
            // Loop on cardiac gates
            for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++) if (m5p_sensitivity[th][fr][rg][cg])
            {
              free(m5p_sensitivity[th][fr][rg][cg]);
            }
            free(m5p_sensitivity[th][fr][rg]);
          }
          free(m5p_sensitivity[th][fr]);
        }
        free(m5p_sensitivity[th]);
      }
    }
    // Free the first pointer
    free(m5p_sensitivity);
  }
  // And set the pointer to NULL
  m5p_sensitivity = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB* oImageSpace::AllocateMiscellaneousImage()
{
  // In this function, we want to allocate a new miscellaneous image on the double pointer member
  // m2p_miscellaneousImage, and then we return the pointer to the allocated image.

  // Update the number of miscellaneous images
  m_nbMiscellaneousImages++;
  // Case 1: the current number of image is 0
  if (m_nbMiscellaneousImages==0)
  {
    // Allocate the first pointer
    m2p_miscellaneousImage = (FLTNB**)malloc(m_nbMiscellaneousImages*sizeof(FLTNB*));
  }
  // Case 2: we already have some miscellaneous images
  else
  {
    // Reallocate the first pointer
    m2p_miscellaneousImage = (FLTNB**)realloc(m2p_miscellaneousImage,m_nbMiscellaneousImages*sizeof(FLTNB*));
  }
  // Allocate the image
  m2p_miscellaneousImage[m_nbMiscellaneousImages-1] = (FLTNB*)malloc(mp_ID->GetNbVoxXYZ()*sizeof(FLTNB));
  // Return the pointer
  return m2p_miscellaneousImage[m_nbMiscellaneousImages-1];
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageSpace::DeallocateMiscellaneousImage()
{
  // Test the pointer existence
  if (m2p_miscellaneousImage!=NULL)
  {
    // Verbose
    if (m_verbose>=3) Cout("oImageSpace::DeallocateMiscellaneousImage() -> Free memory"<< endl);
    // Loop over all miscellaneous images
    for (int m=0; m<m_nbMiscellaneousImages; m++)
    {
      // Test the pointer and free it
      if (m2p_miscellaneousImage[m]) free(m2p_miscellaneousImage[m]);
    }
    // Free the first pointer
    free(m2p_miscellaneousImage);
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageSpace::InitMultiModalImage(const vector<string>& a_pathToImage)
{
  // Allocate and read only if a file is provided
  if (!a_pathToImage.empty())
  {
    // Allocate memory (if not already done)
    if (!m_loadedMultiModal)
    {
      m5p_multiModalImage = (FLTNB*****)malloc(mp_ID->GetNbMultiModalImages()*sizeof(FLTNB****));
      for (int nb=0; nb<mp_ID->GetNbMultiModalImages(); nb++)
      {
        // Verbose
        if (m_verbose>=3) Cout("oImageSpace::InitMultiModalImage() -> From file '" << a_pathToImage.at(nb) << "'"<< endl);
        // Open file
        ifstream image_file(a_pathToImage.at(nb).c_str(), ios::in|ios::binary);
        if (!image_file.is_open())
        {
          Cerr("***** oImageSpace::InitMultiModalImage() -> Input file '" << a_pathToImage.at(nb) << "' is missing or corrupted !" << endl);
          // Read failure implies that the multiModal image will never be used, so free all the allocated memory
          DeallocateMultiModalImage();
          return 1;
        }
        m5p_multiModalImage[nb] = (FLTNB****)malloc(mp_ID->GetNbTimeFrames()*sizeof(FLTNB***));
        for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
        {
          m5p_multiModalImage[nb][fr] = (FLTNB***)malloc(mp_ID->GetNbRespGates()*sizeof(FLTNB**));
          for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
          {
            m5p_multiModalImage[nb][fr][rg] = (FLTNB**)malloc(mp_ID->GetNbCardGates()*sizeof(FLTNB*));
            for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
            {
              m5p_multiModalImage[nb][fr][rg][cg] = (FLTNB*)malloc(mp_ID->GetNbVoxXYZ()*sizeof(FLTNB));
            }
          }
        }
        // Interfile image reading (INTF_LERP_ENABLED = interpolation allowed)
        if (IntfReadImage(a_pathToImage.at(nb), m5p_multiModalImage[nb], mp_ID, m_verbose, INTF_LERP_ENABLED))
        {
          Cerr("***** oImageSpace::InitMultiModalImage() -> Error reading interfile image '" << a_pathToImage.at(nb) << "' !" << endl);
          // Read failure implies that the multiModal image will never be used, so free all the allocated memory
          DeallocateMultiModalImage();
          return 1;
        }
        // Close file
        image_file.close();
      }
      // Flag as loaded
      m_loadedMultiModal = true;
    }
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageSpace::DeallocateMultiModalImage()
{
  // Check the first pointer 
  if (m5p_multiModalImage)
  {
    // Verbose
    if (m_verbose>=3) Cout("oImageSpace::DeallocateMultiModalImage() -> Free memory" << endl);
    for (int nb=0; nb<mp_ID->GetNbMultiModalImages(); nb++) if (m5p_multiModalImage[nb])
    {
      for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++) if (m5p_multiModalImage[nb][fr])
      {
        for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++) if (m5p_multiModalImage[nb][fr][rg])
        {
          for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++) if (m5p_multiModalImage[nb][fr][rg][cg])
          {
            free(m5p_multiModalImage[nb][fr][rg][cg]);
          }
          free(m5p_multiModalImage[nb][fr][rg]);
        }
        free(m5p_multiModalImage[nb][fr]);
      }
      free(m5p_multiModalImage[nb]);
    }
    free(m5p_multiModalImage);
  }
  // Reset the pointer to NULL
  m5p_multiModalImage = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageSpace::InitMaskImage(const string& a_pathToImage)
{
  // Allocate and read only if a file is provided
  if (!a_pathToImage.empty())
  {
    // Verbose
    if (m_verbose>=3) Cout("oImageSpace::InitMaskImage() -> From file '" << a_pathToImage << "'"<< endl);
    // Open file
    ifstream image_file(a_pathToImage.c_str(), ios::in|ios::binary);
    if (!image_file.is_open())
    {
      Cerr("***** oImageSpace::InitMaskImage() -> Input file '" << a_pathToImage << "' is missing or corrupted !" << endl);
      return 1;
    }
    // Allocate memory (if not already done)
    if (!m_loadedMask)
    {
      mp_maskImage = (FLTNB*)malloc(mp_ID->GetNbVoxXYZ()*sizeof(FLTNB));
    }
    // Interfile image reading (INTF_LERP_ENABLED = interpolation allowed)
    if (IntfReadImage(a_pathToImage, mp_maskImage, mp_ID, m_verbose, INTF_LERP_ENABLED))
    {
      Cerr("***** oImageSpace::InitMaskImage() -> Error reading interfile image '" << a_pathToImage << "' !" << endl);
      // Read failure implies that the mask image will never be used, so free all the allocated memory
      DeallocateMaskImage();
      return 1;
    }
    // Flag as loaded
    m_loadedMask = true;
    // Close file
    image_file.close();
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageSpace::DeallocateMaskImage()
{
  // Check pointer
  if (mp_maskImage)
  {
    // Verbose
    if (m_verbose>=3) Cout("oImageSpace::DeallocateMaskImage() -> Free memory" << endl);
    free(mp_maskImage);
  }
  // Set the pointer to NULL
  mp_maskImage = NULL;
  // flag as not loaded
  m_loadedMask = false;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InstantiateOutputImage()
  \brief Instanciate Image matrices dedicated to output writing on disk
  \details Additionnal output image matrix is needed in the case the reconstruction uses intrinsic temporal basis functions
           In this case, the image matrices are defined in the temporal image basis functions space, therefore requiring an additional step to recover the images in the regular image-space.
           If no intrinsic temporal basis functions are used, m4p_outputImage just refers to the address of the image matrix containing the regular image.
*/
void oImageSpace::InstantiateOutputImage()
{
  if(m_verbose>=3) Cout("oImageSpace::InstantiateOutputImage ..."<< endl);
    
  // Here, we will allocate the first 3 pointers to the dynamic dimensions
  m4p_outputImage = new FLTNB***[mp_ID->GetNbTimeFrames()];
  for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
  {
    m4p_outputImage[fr] = new FLTNB**[mp_ID->GetNbRespGates()];
    for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
    {
      m4p_outputImage[fr][rg] = new FLTNB*[mp_ID->GetNbCardGates()];
      for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
      {
        // For the last pointer (i.e. the voxels), we allocate them only if we strictly need them, which
        // means in the case where one of the dynamic dimensions makes internal use of basis functions,
        // otherwise, we will just make this output image pointing to the forward image which is useless
        // at the moment of computing the output image.
        if ( mp_ID->GetTimeStaticFlag() &&
             mp_ID->GetRespStaticFlag() &&
             mp_ID->GetCardStaticFlag() )
        {
          // In this case, the time frames and basis functions represent the same thing, same for respiratory/cardiac
          // gates and their equivalent basis functions
          m4p_outputImage[fr][rg][cg] = m4p_forwardImage[fr][rg][cg];
        }
        else
        {
          // We allocate the output image
          m4p_outputImage[fr][rg][cg] = new FLTNB[mp_ID->GetNbVoxXYZ()];
        }
      }
    }
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn DeallocateOutputImage()
  \brief Free memory for the Image matrices dedicated to output writing on disk
*/
void oImageSpace::DeallocateOutputImage()
{
  if(m_verbose>=3) Cout("oImageSpace::DeallocateOutputImage ..."<< endl);
    
  if (m4p_outputImage)
  {
    for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
    {
      if (m4p_outputImage[fr])
      {
        for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
        {
          if (m4p_outputImage[fr][rg])
          {
            for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
            {
              // Distinguish between intrinsic use of dynamic basis functions or not (in the latter,
              // the m4p_outputImage points to the m4p_forwardImage, so we do not touch it)
              if ( (!mp_ID->GetTimeStaticFlag() ||
                    !mp_ID->GetRespStaticFlag() ||
                    !mp_ID->GetCardStaticFlag()) &&
                     m4p_outputImage[fr][rg][cg] )
              {
                delete m4p_outputImage[fr][rg][cg];
              }
            }
            delete[] m4p_outputImage[fr][rg];
          }
        }
        delete[] m4p_outputImage[fr];
      }
    }
    delete[] m4p_outputImage;
  }
  m4p_outputImage = NULL;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InstantiateBwdImageForDeformation()
  \brief Memory allocation for the buffer backward image required for image-based deformation

void oImageSpace::InstantiateRefImagesForDeformation()
{
  if(m_verbose>=3) Cout("oImageSpace::InstantiateBwdImageForDeformation ..."<< endl);
    
  m5p_refDynBackwardImage = new FLTNB****[m_nbBackwardImages];
  
  for (int img=0; img<m_nbBackwardImages; img++)
  {
    m5p_refDynBackwardImage[img] = new FLTNB***[mp_ID->GetNbTimeBasisFunctions()];
    for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
    {
      // The fourth pointer is for the number of respiratory basis functions
      m5p_refDynBackwardImage[img][tbf] = new FLTNB**[mp_ID->GetNbRespBasisFunctions()];
      for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
      {
        // The fifth pointer is for the number of cardiac basis functions
        m5p_refDynBackwardImage[img][tbf][rbf] = new FLTNB*[mp_ID->GetNbCardBasisFunctions()];
        for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
        {
          // The sixth pointer is for the 3D space
          m5p_refDynBackwardImage[img][tbf][rbf][cbf] = new FLTNB[mp_ID->GetNbVoxXYZ()];
        }
      }
    }
  }
  
}
*/

/*
  \fn      oImageSpace::InstantiateRefImagesForDeformation
  \brief   Allocate memory for the buffer sensitivity image required for image-based deformation. This function is called from the Deformation Manager
*/
void oImageSpace::InstantiateRefImagesForDeformation()
{
  if(m_verbose>=3) Cout("oImageSpace::InstantiateRefImagesForDeformation ..."<< endl);


  // Forward image
  m4p_refDynForwardImage = new FLTNB***[mp_ID->GetNbTimeBasisFunctions()];
  for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
  {
    // The fourth pointer is for the number of respiratory basis functions
    m4p_refDynForwardImage[tbf] = new FLTNB**[mp_ID->GetNbRespBasisFunctions()];
    for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
    {
      // The fifth pointer is for the number of cardiac basis functions
      m4p_refDynForwardImage[tbf][rbf] = new FLTNB*[mp_ID->GetNbCardBasisFunctions()];
      for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
      {
        // The sixth pointer is for the 3D space
        m4p_refDynForwardImage[tbf][rbf][cbf] = new FLTNB[mp_ID->GetNbVoxXYZ()];
      }
    }
  }

  
  // Backward image
  m5p_refDynBackwardImage = new FLTNB****[m_nbBackwardImages];
  
  for (int img=0; img<m_nbBackwardImages; img++)
  {
    m5p_refDynBackwardImage[img] = new FLTNB***[mp_ID->GetNbTimeBasisFunctions()];
    for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
    {
      // The fourth pointer is for the number of respiratory basis functions
      m5p_refDynBackwardImage[img][tbf] = new FLTNB**[mp_ID->GetNbRespBasisFunctions()];
      for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
      {
        // The fifth pointer is for the number of cardiac basis functions
        m5p_refDynBackwardImage[img][tbf][rbf] = new FLTNB*[mp_ID->GetNbCardBasisFunctions()];
        for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
        {
          // The sixth pointer is for the 3D space
          m5p_refDynBackwardImage[img][tbf][rbf][cbf] = new FLTNB[mp_ID->GetNbVoxXYZ()];
        }
      }
    }
  }


  // Sensitivity image
  if(m_loadedSensitivity == false)
  {
    m4p_refDynSensitivityImage = new FLTNB***[mp_ID->GetNbTimeFrames()];
    
    for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
    {
      m4p_refDynSensitivityImage[fr] = new FLTNB**[mp_ID->GetNbRespGates()];
      for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
      {
        m4p_refDynSensitivityImage[fr][rg] = new FLTNB*[mp_ID->GetNbCardGates()];
        for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
          m4p_refDynSensitivityImage[fr][rg][cg] = new FLTNB[mp_ID->GetNbVoxXYZ()];
      }
    }
  }

}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn DeallocateBwdImageForDeformation()
  \brief Free memory for the buffer backward image required for image-based deformation

void oImageSpace::DeallocateBwdImageForDeformation()
{
  if(m_verbose>=3) Cout("oImageSpace::DeallocateBwdImageForDeformation ..."<< endl);
    
  for (int img=0; img<m_nbBackwardImages; img++)
  {
    for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
    {
      for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
      {
        for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
        {
          if (m5p_refDynBackwardImage[img][tbf][rbf][cbf]) delete m5p_refDynBackwardImage[img][tbf][rbf][cbf];
        }
        if (m5p_refDynBackwardImage[img][tbf][rbf]) delete m5p_refDynBackwardImage[img][tbf][rbf];
      }
      if (m5p_refDynBackwardImage[img][tbf]) delete[] m5p_refDynBackwardImage[img][tbf];
    }
    if (m5p_refDynBackwardImage[img]) delete[] m5p_refDynBackwardImage[img];
  }
  if (m5p_refDynBackwardImage) delete[] m5p_refDynBackwardImage;
  m5p_refDynBackwardImage = NULL;
}
*/




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      oImageSpace::DeallocateRefImagesForDeformation
  \brief   Free memory for the buffer sensitivity image required for image-based deformation. This function is called from the Deformation Manager
*/
void oImageSpace::DeallocateRefImagesForDeformation()
{
  if(m_verbose>=3) Cout("oImageSpace::DeallocateRefImagesForDeformation ..."<< endl);

  // Forward image
  for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
  {
    for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
    {
      for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
      {
        if (m4p_refDynForwardImage[tbf][rbf][cbf]) delete[] m4p_refDynForwardImage[tbf][rbf][cbf];
      }
      if (m4p_refDynForwardImage[tbf][rbf]) delete[] m4p_refDynForwardImage[tbf][rbf];
    }
    if (m4p_refDynForwardImage[tbf]) delete[] m4p_refDynForwardImage[tbf];
  }
  if (m4p_refDynForwardImage) delete[] m4p_refDynForwardImage;

  m4p_refDynForwardImage = NULL;

  if(m_verbose>=3) Cout("oImageSpace::Deallocatem5p_refDynBackwardImage ..."<< endl);

  // Backward image
  for (int img=0; img<m_nbBackwardImages; img++)
  {
    for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
    {
      for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
      {
        for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
        {
          if(m_verbose>=3) Cout("oImageSpace::Deallocate m5p_refDynBackwardImage[img][tbf][rbf][cbf] ..."<< endl);
          if (m5p_refDynBackwardImage[img][tbf][rbf][cbf]) delete[] m5p_refDynBackwardImage[img][tbf][rbf][cbf];
        }
        if(m_verbose>=3) Cout("oImageSpace::Deallocate m5p_refDynBackwardImage[img][tbf][rbf]..."<< endl);
        if (m5p_refDynBackwardImage[img][tbf][rbf]) delete[] m5p_refDynBackwardImage[img][tbf][rbf];
      }
      if(m_verbose>=3) Cout("oImageSpace::Deallocate m5p_refDynBackwardImage[img][tbf] ..."<< endl);
      if (m5p_refDynBackwardImage[img][tbf]) delete[] m5p_refDynBackwardImage[img][tbf];
    }
    if(m_verbose>=3) Cout("oImageSpace::Deallocate m5p_refDynBackwardImage[img]..."<< endl);
    if (m5p_refDynBackwardImage[img]) delete[] m5p_refDynBackwardImage[img];
  }
  if(m_verbose>=3) Cout("oImageSpace::Deallocate m5p_refDynBackwardImage..."<< endl);
  if (m5p_refDynBackwardImage) delete[] m5p_refDynBackwardImage;
  m5p_refDynBackwardImage = NULL;

  if(m_verbose>=3) Cout("oImageSpace::DeallocateRm4p_refDynSensitivityImage ..."<< endl);
  // Sensitivity image
  if(m_loadedSensitivity == false)
  {
    for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
    {
      for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
      {
        for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
          if (m4p_refDynSensitivityImage[fr][rg][cg]) delete[] m4p_refDynSensitivityImage[fr][rg][cg];
  
        if (m4p_refDynSensitivityImage[fr][rg]) delete[] m4p_refDynSensitivityImage[fr][rg];
      }
      if (m4p_refDynSensitivityImage[fr]) delete[] m4p_refDynSensitivityImage[fr];
    }
    if (m4p_refDynSensitivityImage) delete[] m4p_refDynSensitivityImage;
    m4p_refDynSensitivityImage = NULL;
  }
  
    if(m_verbose>=3) Cout("oImageSpace::DeallocateRm4p_refDynSensitivityImageOK ..."<< endl);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InstantiateSensImageForDeformation()
  \brief Memory allocation for the buffer sensitivity image required for image-based deformation (only for PET HISTOGRAM mode)

void oImageSpace::InstantiateSensImageForDeformation()
{
  if(m_verbose>=3) Cout("oImageSpace::InstantiateSensImageForDeformation ..."<< endl);
    
  m4p_refDynSensitivityImage = new FLTNB***[mp_ID->GetNbTimeFrames()];
  
  for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
  {
    m4p_refDynSensitivityImage[fr] = new FLTNB**[mp_ID->GetNbRespGates()];
    for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
    {
      m4p_refDynSensitivityImage[fr][rg] = new FLTNB*[mp_ID->GetNbCardGates()];
      for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
        m4p_refDynSensitivityImage[fr][rg][cg] = new FLTNB[mp_ID->GetNbVoxXYZ()];
    }
  }
}
*/



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn DeallocateBwdImageForDeformation()
  \brief Free memory for the buffer sensitivity image required for image-based deformation

void oImageSpace::DeallocateSensImageForDeformation()
{
  if(m_verbose>=3) Cout("oImageSpace::DeallocateSensImageForDeformation ..."<< endl);
    
  for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
  {
    for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
    {
      for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
        if (m4p_refDynSensitivityImage[fr][rg][cg]) delete[] m4p_refDynSensitivityImage[fr][rg][cg];

      if (m4p_refDynSensitivityImage[fr][rg]) delete[] m4p_refDynSensitivityImage[fr][rg];
    }
    if (m4p_refDynSensitivityImage[fr]) delete[] m4p_refDynSensitivityImage[fr];
  }
  if (m4p_refDynSensitivityImage) delete[] m4p_refDynSensitivityImage;
  m4p_refDynSensitivityImage = NULL;
}
*/



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InstantiateVisitedVoxelsImage()
  \brief Memory allocation and initialization for the image matrix containing binary information regarding which 3D voxels have been visited during the projection steps.
*/
void oImageSpace::InstantiateVisitedVoxelsImage()
{
  if(m_verbose>=3) Cout("oImageSpace::InstantiateVisitedVoxelsImage ..."<< endl);
    
  /*
  mp_visitedVoxelsImage = new FLTNB[mp_ID->GetNbVoxXYZ()];
  for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++) 
    mp_visitedVoxelsImage[v] = 0.;
  */
  
  m4p_visitedVoxelsImage = new FLTNB***[mp_ID->GetNbTimeFrames()];
  
  for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
  {
    m4p_visitedVoxelsImage[fr] = new FLTNB**[mp_ID->GetNbRespGates()];
    
    for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
    {
      m4p_visitedVoxelsImage[fr][rg] = new FLTNB*[mp_ID->GetNbCardGates()];
      
      for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
      {
        m4p_visitedVoxelsImage[fr][rg][cg] = new FLTNB[mp_ID->GetNbVoxXYZ()];
        
        for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
          m4p_visitedVoxelsImage[fr][rg][cg][v] = 0;
      }
    }
  }
  
  
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn DeallocateVisitedVoxelsImage()
  \brief Free memory for the image matrix containing binary information regarding which 3D voxels have been visited during the projection steps.
*/
void oImageSpace::DeallocateVisitedVoxelsImage()
{
  if(m_verbose>=3) Cout("oImageSpace::InstantiateVisitedVoxelsImage ..."<< endl);
    
  //if (mp_visitedVoxelsImage) delete mp_visitedVoxelsImage;
  //mp_visitedVoxelsImage = NULL;
  

  for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
    for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
      for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
        if(m4p_visitedVoxelsImage[fr][rg][cg]) delete[] m4p_visitedVoxelsImage[fr][rg][cg];
        
  for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
    for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
      if(m4p_visitedVoxelsImage[fr][rg]) delete[] m4p_visitedVoxelsImage[fr][rg];

  for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
    if(m4p_visitedVoxelsImage[fr]) delete[] m4p_visitedVoxelsImage[fr];
      
  if(m4p_visitedVoxelsImage) delete[] m4p_visitedVoxelsImage;
  m4p_visitedVoxelsImage = NULL;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LMS_DeallocateAttenuationImage()
  \brief  Free memory for the Attenuation image matrices (for analytical projection or list-mode sensitivity generation)
*/
void oImageSpace::LMS_DeallocateAttenuationImage()
{
  if(m_verbose>=3) Cout("oImageSpace::LMS_DeallocateAttenuationImage ..."<< endl);
  
  if(m4p_attenuation != NULL)
  {
    for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
    {
      for(int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr); rg++)
      {
        for(int cg=0 ; cg<mp_ID->GetNb2ndMotImgsForLMS(); cg++)
          if(m4p_attenuation[fr][rg][cg] != NULL) delete m4p_attenuation[fr][rg][cg];
          
        if(m4p_attenuation[fr][rg] != NULL) delete m4p_attenuation[fr][rg];
      }
      if(m4p_attenuation[fr] != NULL) delete m4p_attenuation[fr];
    }
    if(m4p_attenuation != NULL) delete[] m4p_attenuation;
    m4p_attenuation = NULL;
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageSpace::InitAttenuationImage(const string& a_pathToAtnImage)
{
  // todo : read the number of dynamic images from Interfile header
  //        Allocate memory and initialize matrices according to the number of fr/rg/cg

  // Allocate only if an image has been provided
  if (a_pathToAtnImage.empty()) return 0;
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL) Cout("oImageSpace::InitAttenuationImage() -> Initialize and read attenuation image from file '" << a_pathToAtnImage << "'"<< endl);
  // Allocate only if not already done
  if (!m4p_attenuation) 
  {
    m4p_attenuation = new FLTNB***[mp_ID->GetNbTimeFrames()];
    for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
    {
      m4p_attenuation[fr] = new FLTNB**[mp_ID->GetNb1stMotImgsForLMS(fr)];
      for(int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr); rg++)
      {
        m4p_attenuation[fr] = new FLTNB**[mp_ID->GetNb1stMotImgsForLMS(fr)];
        for(int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr); rg++)
        {
          m4p_attenuation[fr][rg] = new FLTNB*[mp_ID->GetNb2ndMotImgsForLMS()];
          
          for(int cg=0 ; cg<mp_ID->GetNb2ndMotImgsForLMS(); cg++)
            m4p_attenuation[fr][rg][cg] = new FLTNB[mp_ID->GetNbVoxXYZ()];
        }
      }
    }
  }
  // Open file
  ifstream image_file(a_pathToAtnImage.c_str(), ios::in|ios::binary); // Read the corresponding attenuation image
  // Check opening
  if (!image_file.is_open())
  {
    Cerr("***** oImageSpace::InitAttenuationImage() -> Failed to open attenuation image from file '" << a_pathToAtnImage << "' !" << endl);
    return 1;
  }
  else
  {
    // Interfile image reading (INTF_LERP_ENABLED = interpolation allowed)
    if (IntfReadImage(a_pathToAtnImage, m4p_attenuation, mp_ID, m_verbose, INTF_LERP_ENABLED) )
    {
      Cerr("***** oImageSpace::InitAttenuationImage() -> An error occured while reading from file '" << a_pathToAtnImage << "' !" << endl);  
      return 1;
    }
  }
  // Close file
  image_file.close();
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InitImage
  \param a_pathToInitialImage : path to an existing image
  \param a_value : value to initialize each voxel with, if an input image is not provided
  \brief Initialize the main image, either using:
         - an existing image at path 'a_pathToInitialImage'
         - initialize each voxel with 'a_value'.
  \return 0 if success, positive value otherwise
*/
int oImageSpace::InitImage(const string& a_pathToInitialImage, FLTNB a_value)
{
  if(m_verbose>=3) Cout("oImageSpace::InitImage ..."<< endl);
    
  if (!a_pathToInitialImage.empty()) // Image initiale has been provided
  {
    if (LoadInitialImage(a_pathToInitialImage) )
    {
      Cerr("***** oImageSpace::InitImage()-> Error while trying to read file at " << a_pathToInitialImage << endl);
      return 1;
    }
  }

  else // Uniform initialization
  {
    // We initialize the content of the cylindrical FOV with the provided
    // value, and the exterior with a ratio of the provided value
    FLTNB exterior_fov_value = a_value * 1.; // For the moment we let the same value everywhere
    // Precast half the number of voxels over X and Y minus 1 (for efficiency)
    FLTNB flt_base_x = 0.5*((FLTNB)(mp_ID->GetNbVoxX()-1));
    FLTNB flt_base_y = 0.5*((FLTNB)(mp_ID->GetNbVoxY()-1));
    // Compute FOV elipse radius over X and Y, then squared
    FLTNB squared_radius_x = 0.5 * ((FLTNB)(mp_ID->GetNbVoxX())) * mp_ID->GetVoxSizeX();
    squared_radius_x *= squared_radius_x;
    FLTNB squared_radius_y = 0.5 * ((FLTNB)(mp_ID->GetNbVoxY())) * mp_ID->GetVoxSizeY();
    squared_radius_y *= squared_radius_y;
    // We assume that the computation of the distance from the center for a given
    // voxel and comparing it with the output FOV percentage costs more than performing
    // the loops in an inverse order compared to how the image is stored in memory.
    // Thus we begin the loops over X, then Y, then we test and if test passes, we
    // do the remaining loop over Z and over all dynamic dimensions.
    int x;
    #pragma omp parallel for private(x) schedule(guided)
    for (x=0; x<mp_ID->GetNbVoxX(); x++)
    {
      // Compute X distance from image center, then squared
      FLTNB squared_distance_x = (((FLTNB)x)-flt_base_x) * mp_ID->GetVoxSizeX();
      squared_distance_x *= squared_distance_x;
      // Loop over Y
      for (int y=0; y<mp_ID->GetNbVoxY(); y++)
      {
        // Compute Y distance from image center, then squared
        FLTNB squared_distance_y = (((FLTNB)y)-flt_base_y) * mp_ID->GetVoxSizeY();
        squared_distance_y *= squared_distance_y;
        // The value to affect to the voxel
        FLTNB affected_value = 0.;
        // Test if the voxel is inside the FOV elipse, then we affect the provided value
        if ( squared_distance_x/squared_radius_x + squared_distance_y/squared_radius_y <= 1. ) affected_value = a_value;
        // Else if outside, we affect the exterior value
        else affected_value = exterior_fov_value;
        // Loop over Z
        for (int z=0; z<mp_ID->GetNbVoxZ(); z++)
        {
          // Compute global voxel index
          INTNB index = z*mp_ID->GetNbVoxXY() + y*mp_ID->GetNbVoxX() + x;
          // Loops on dynamic dimensions
          for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
            for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
              for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
                // Affect image value
                m4p_image[tbf][rbf][cbf][index] = affected_value;
        }
      }
    }
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LoadInitialImage()
  \param a_pathToImage : path to an existing image
  \brief Load the initial image provided by the user in the corresponding matrix
  \return 0 if success, positive value otherwise
*/
int oImageSpace::LoadInitialImage(const string& a_pathToImage)
{
  if(m_verbose>=3) Cout("oImageSpace::LoadInitialImage ..."<< endl);
    
  ifstream image_file;
  image_file.open(a_pathToImage.c_str(), ios::in | ios::binary);

  if(!image_file.is_open()) 
  {
    Cerr("***** oImageSpace::LoadInitialImage()-> Error reading file!" << endl);
    return 1;
  }
  else
  {
    // Interfile image reading (INTF_LERP_DISABLED = no interpolation allowed)
    // if(IntfReadImage(a_pathToImage, m4p_image, mp_ID, m_verbose, INTF_LERP_DISABLED))
    if(IntfReadImage(a_pathToImage, m4p_image, mp_ID, m_verbose, INTF_LERP_ENABLED))
    {
      Cerr("***** oImageSpace::LoadInitialImage()-> Error reading Interfile : " << a_pathToImage << " !" << endl);  
      return 1;
    }
  }
  image_file.close();
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InitializationBackwardImage()
  \brief Initialize each voxel of the backward images to 0, also for sensitivity if not loaded (estimated on the fly).
*/
void oImageSpace::InitBackwardImage()
{
  if(m_verbose>=3) Cout("oImageSpace::InitBackwardImage ..." << endl);

  // Reset backward images to 0.
  for (int img=0; img<m_nbBackwardImages; img++)
  {
    int th;
    #pragma omp parallel for private(th) schedule(static, 1)
    for (th=0; th<mp_ID->GetNbThreadsForProjection(); th++)
      for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
        for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
          for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
            for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
            {
              m6p_backwardImage[img][th][tbf][rbf][cbf][v] = 0.;
            }
  }

  // If on-the-fly sensitivity, then reset to 0.
  if (!m_loadedSensitivity)
  {
    int th;
    #pragma omp parallel for private(th) schedule(static, 1)
    for (th=0; th<mp_ID->GetNbThreadsForProjection(); th++)
      for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
        for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
          for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
            for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
              m5p_sensitivity[th][fr][rg][cg][v] = 0.;
  }
}






// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InitSensitivityImage()
  \param a_pathToSensitivityImage : path to the sensitivity image (should be provided only in list-mode reconstruction)
  \brief Initialization for the sensitivity image matrices
  \details Sensitivity image initialization depends on the reconstruction mode :
           - list-mode: Sensitivity image has been computed before reconstruction, it is loaded from the path provided in parameter
           - histogram: Sensitivity image is calculated on the fly during reconstruction.  
           First dimension (thread) is only used in histogram mode, as the on-the-fly sensitivity image computation must be thread safe
  \return  0 if success, positive value otherwise
*/
int oImageSpace::InitSensitivityImage(const string& a_pathToSensitivityImage)
{ 
  if(m_verbose>=3) Cout("oImageSpace::InitSensitivityImage ..."<< endl);
    
  // =====================================================================================================
  // Case 1: a path to a sensitivity image is provided, meaning that we are in listmode and do not compute
  //         the sensitivity on-the-fly. In that case, we do not take the multi-threading into account.
  // =====================================================================================================

  if (!a_pathToSensitivityImage.empty())
  {
    // Open file
    ifstream input_file;
    input_file.open(a_pathToSensitivityImage.c_str(), ios::binary | ios::in);
    // Read file
    if (input_file.is_open())
    { 
      // Interfile image reading (INTF_LERP_DISABLED = no interpolation allowed)
      if(IntfReadImage(a_pathToSensitivityImage, m5p_sensitivity[0], mp_ID, m_verbose, INTF_LERP_DISABLED) )
      {
        Cerr("***** oImageSpace::InitSensitivityImage()-> Error reading Interfile : " << a_pathToSensitivityImage << " !" << endl);  
        return 1;
      }    
      input_file.close();
      m_loadedSensitivity = true;
    }
    else
    {
      Cerr("***** oImageSpace::InitSensitivityImage() -> Input sensitivity file '" << a_pathToSensitivityImage << "' is missing or corrupted !" << endl);
      return 1;
    }
  }

  // =====================================================================================================
  // Case 2: no sensitivity provided, thus we are in histogram mode and will commpute it on-the-fly. We
  //         thus need it to be thread-safe, so we allocated for all threads.
  // =====================================================================================================

  else 
  {
/*
    // Standard initialization
    for (int th=0; th<mp_ID->GetNbThreadsForProjection(); th++)
      for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
        for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
          for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
            for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
              m5p_sensitivity[th][fr][rg][cg][v] = -1.;
* */
    // No need to initialize to any value, as it will be done by the InitBackwardImage function which is called at the beginning
    // of each subset, and deals with the sensitivity initialization to 0 when m_loadedSensitivity is false
    m_loadedSensitivity = false;
  }

  // End
  return 0;
}





// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      oImageSpace::InitRefImagesForDeformation
  \brief   Initialize the references image dedicated to image-based deformation. This function is called from the Deformation Manager
*/
void oImageSpace::InitRefImagesForDeformation()
{
  if(m_verbose>=3) Cout("oDeformationManager::InitBwdImageForDeformation ..." << endl);

  // forward image
  for (int img=0; img<m_nbBackwardImages; img++)
    for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
      for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
        for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
          for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
            m4p_refDynForwardImage[tbf][rbf][cbf][v] = m4p_forwardImage[tbf][rbf][cbf][v];

  // backward image
  for (int img=0; img<m_nbBackwardImages; img++)
    for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
      for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
        for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
          for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
            m5p_refDynBackwardImage[img][tbf][rbf][cbf][v] = 0.;
            
  if(m_loadedSensitivity == false)
  {
    for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
      for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
        for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
          for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
            m4p_refDynSensitivityImage[fr][rg][cg][v] = 0.;
  }
}





// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InitializationSensImageForDeformation()
  \brief Initialize the buffer sensitivity image dedicated to image-based deformation, if required (histogram mode, as sensitivity is not loaded)

void oImageSpace::InitSensImageForDeformation()
{
  if(m_verbose>=3) Cout("oImageSpace::InitSensImageForDeformation ..."<< endl);
  
  if(m_loadedSensitivity == false)
  {
    for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
      for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
        for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
          for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
            m4p_refDynSensitivityImage[fr][rg][cg][v] = 0.;
  }
}
*/




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InstantiateOutputImage()
  \brief Compute output image using the m4p_image matrix and the time/respiratory/cardiac basis functions. Store the result in the m4p_outputImage matrix
  \details If time/respiratory/cardiac basis functions have been initialized, this function has no effect.
*/
void oImageSpace::ComputeOutputImage()
{
  if(m_verbose>=3) Cout("oImageSpace::ComputeOutputImage ..."<< endl);
  
  // First loop on frames
  for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
  {
    // Second loop on respiratory gates
    for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
    {
      // Third loop on cardiac gates
      for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
      {
        // Reset m4p_outputImage to 0
        int v;
        #pragma omp parallel for private(v) schedule(guided)
        for (v=0; v<mp_ID->GetNbVoxXYZ(); v++)
          m4p_outputImage[fr][rg][cg][v] = 0.;              
        // First loop on time basis functions
        for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
        {
          // Get frame/basis coefficient and continue if null
          FLTNB time_basis_coef = mp_ID->GetTimeBasisCoefficient(tbf,fr);
          if (time_basis_coef==0.) continue;
          // Second loop on respiratory basis functions
          for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
          {
            // Get resp_gate/basis coefficient and continue if null
            FLTNB resp_basis_coef = mp_ID->GetRespBasisCoefficient(rbf,rg);
            if (resp_basis_coef==0.) continue;
            // Third loop on cardiac basis functions
            for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
            {
              // Get card_gate_basis coefficient and continue if null
              FLTNB card_basis_coef = mp_ID->GetCardBasisCoefficient(cbf,cg);
              if (card_basis_coef==0.) continue;
              // Compute global coefficient
              FLTNB global_basis_coeff = time_basis_coef * resp_basis_coef * card_basis_coef;
              // Loop on voxel with OpenMP (let the chunk size by default as all values are aligned in memory)
              #pragma omp parallel for private(v) schedule(guided)
              for (v=0; v<mp_ID->GetNbVoxXYZ(); v++)
              {
                // Add contribution from these basis functions
                m4p_outputImage[fr][rg][cg][v] += m4p_image[tbf][rbf][cbf][v] * global_basis_coeff;
              }
            }
          }
        }
      }
    }
  }
  
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageSpace::ApplyOutputFlip()
{
  // If no flip, then return
  if ( !mp_ID->GetFlipOutX() &&
       !mp_ID->GetFlipOutY() &&
       !mp_ID->GetFlipOutZ() ) return 0;
  // Verbose
  if (m_verbose>=2)
  {
    Cout("oImageSpace::ApplyOutputFlip() -> Flip image" << endl);
    if (mp_ID->GetFlipOutX()) Cout("  --> Over X" << endl);
    if (mp_ID->GetFlipOutY()) Cout("  --> Over Y" << endl);
    if (mp_ID->GetFlipOutZ()) Cout("  --> Over Z" << endl);
  }
  // First loop on frames
  for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
  {
    // Second loop on respiratory gates
    for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
    {
      // Third loop on cardiac gates
      for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
      {
        // Flip over Z
        if (mp_ID->GetFlipOutZ())
        {
          // Half loop over Z
          for (INTNB z_1=0; z_1<mp_ID->GetNbVoxZ()/2; z_1++)
          {
            // Compute opposite Z
            INTNB z_2 = mp_ID->GetNbVoxZ() - 1 - z_1;
            // For efficiency
            INTNB base_z1 = z_1 * mp_ID->GetNbVoxXY();
            INTNB base_z2 = z_2 * mp_ID->GetNbVoxXY();
            // Loop over Y
            for (INTNB y=0; y<mp_ID->GetNbVoxY(); y++)
            {
              // For efficiency
              INTNB base_y = y * mp_ID->GetNbVoxX();
              // Loop over X
              for (INTNB x=0; x<mp_ID->GetNbVoxX(); x++)
              {
                // Compute both indices
                INTNB indice_1 = base_z1 + base_y + x;
                INTNB indice_2 = base_z2 + base_y + x;
                // Switch voxels
                FLTNB buffer = m4p_outputImage[fr][rg][cg][indice_1];
                m4p_outputImage[fr][rg][cg][indice_1] = m4p_outputImage[fr][rg][cg][indice_2];
                m4p_outputImage[fr][rg][cg][indice_2] = buffer;
              }
            }
          }
        }
        // Flip over Y
        if (mp_ID->GetFlipOutY())
        {
          // Loop over Z
          for (INTNB z=0; z<mp_ID->GetNbVoxZ(); z++)
          {
            // For efficiency
            INTNB base_z = z * mp_ID->GetNbVoxXY();
            // Half loop over Y
            for (INTNB y_1=0; y_1<mp_ID->GetNbVoxY()/2; y_1++)
            {
              // Compute opposite Y
              INTNB y_2 = mp_ID->GetNbVoxY() - 1 - y_1;
              // For efficiency
              INTNB base_y1 = y_1 * mp_ID->GetNbVoxX();
              INTNB base_y2 = y_2 * mp_ID->GetNbVoxX();
              // Loop over X
              for (INTNB x=0; x<mp_ID->GetNbVoxX(); x++)
              {
                // Compute both indices
                INTNB indice_1 = base_z + base_y1 + x;
                INTNB indice_2 = base_z + base_y2 + x;
                // Switch voxels
                FLTNB buffer = m4p_outputImage[fr][rg][cg][indice_1];
                m4p_outputImage[fr][rg][cg][indice_1] = m4p_outputImage[fr][rg][cg][indice_2];
                m4p_outputImage[fr][rg][cg][indice_2] = buffer;
              }
            }
          }
        }
        // Flip over X
        if (mp_ID->GetFlipOutX())
        {
          // Loop over Z
          for (INTNB z=0; z<mp_ID->GetNbVoxZ(); z++)
          {
            // For efficiency
            INTNB base_z = z * mp_ID->GetNbVoxXY();
            // Loop over Y
            for (INTNB y=0; y<mp_ID->GetNbVoxY(); y++)
            {
              // For efficiency
              INTNB base_y = y * mp_ID->GetNbVoxX();
              // Half loop over X
              for (INTNB x_1=0; x_1<mp_ID->GetNbVoxX()/2; x_1++)
              {
                // Compute opposite X
                INTNB x_2 = mp_ID->GetNbVoxX() - 1 - x_1;
                // Compute both indices
                INTNB indice_1 = base_z + base_y + x_1;
                INTNB indice_2 = base_z + base_y + x_2;
                // Switch voxels
                FLTNB buffer = m4p_outputImage[fr][rg][cg][indice_1];
                m4p_outputImage[fr][rg][cg][indice_1] = m4p_outputImage[fr][rg][cg][indice_2];
                m4p_outputImage[fr][rg][cg][indice_2] = buffer;
              }
            }
          }
        }
      }
    }
  }

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageSpace::ApplyOutputFOVMasking()
{
  // If the output FOV percent is under 0 (the default value) and number of axial slices to be removed is 0, we do not mask anything
  if ( mp_ID->GetFOVOutPercent()<=0. &&
       mp_ID->GetNbSliceOutMask()==0 ) return 0;
  // Verbose
  if (m_verbose>=2)
  {
    Cout("oImageSpace::ApplyOutputFOVMasking() -> Mask output image" << endl);
    if (mp_ID->GetFOVOutPercent()>0.) Cout("  --> Mask transaxial FOV outside " << mp_ID->GetFOVOutPercent() << " %" << endl);
    if (mp_ID->GetNbSliceOutMask()>0) Cout("  --> Mask " << mp_ID->GetNbSliceOutMask() << " slices from both axial FOV limits" << endl);
  }
  // -----------------------------------------------
  // Transaxial FOV masking
  // -----------------------------------------------
  if (mp_ID->GetFOVOutPercent()>0.)
  {
    // Precast half the number of voxels over X and Y minus 1 (for efficiency)
    FLTNB flt_base_x = 0.5*((FLTNB)(mp_ID->GetNbVoxX()-1));
    FLTNB flt_base_y = 0.5*((FLTNB)(mp_ID->GetNbVoxY()-1));
    // Compute FOV elipse radius over X and Y, then squared
    FLTNB squared_radius_x = 0.5 * ((FLTNB)(mp_ID->GetNbVoxX())) * mp_ID->GetVoxSizeX()
                           * mp_ID->GetFOVOutPercent() / 100.;
    squared_radius_x *= squared_radius_x;
    FLTNB squared_radius_y = 0.5 * ((FLTNB)(mp_ID->GetNbVoxY())) * mp_ID->GetVoxSizeY()
                           * mp_ID->GetFOVOutPercent() / 100.;
    squared_radius_y *= squared_radius_y;
    // We assume that the computation of the distance from the center for a given
    // voxel and comparing it with the output FOV percentage costs more than performing
    // the loops in an inverse order compared to how the image is stored in memory.
    // Thus we begin the loops over X, then Y, then we test and if test passes, we
    // do the remaining loop over Z and over all dynamic dimensions.
    int x;
    #pragma omp parallel for private(x) schedule(guided)
    for (x=0; x<mp_ID->GetNbVoxX(); x++)
    {
      // Compute X distance from image center, then squared
      FLTNB squared_distance_x = (((FLTNB)x)-flt_base_x) * mp_ID->GetVoxSizeX();
      squared_distance_x *= squared_distance_x;
      // Loop over Y
      for (int y=0; y<mp_ID->GetNbVoxY(); y++)
      {
        // Compute Y distance from image center, then squared
        FLTNB squared_distance_y = (((FLTNB)y)-flt_base_y) * mp_ID->GetVoxSizeY();
        squared_distance_y *= squared_distance_y;
        // Test if the voxel is inside the FOV elipse, then we skip this voxel
        if ( squared_distance_x/squared_radius_x + squared_distance_y/squared_radius_y <= 1. ) continue;
        // Loop over Z
        for (int z=0; z<mp_ID->GetNbVoxZ(); z++)
        {
          // Compute global voxel index
          INTNB index = z*mp_ID->GetNbVoxXY() + y*mp_ID->GetNbVoxX() + x;
          // First loop on frames
          for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
            // Second loop on respiratory gates
            for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
              // Third loop on cardiac gates
              for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
                // Put 0 into this voxel
                m4p_outputImage[fr][rg][cg][index] = 0.;
        }
      }
    }
  }
  // -----------------------------------------------
  // Axial FOV masking
  // -----------------------------------------------
  if (mp_ID->GetNbSliceOutMask()>0)
  {
    INTNB removed_slices = mp_ID->GetNbSliceOutMask();
    // First loop on frames
    for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
    {
      // Second loop on respiratory gates
      for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
      {
        // Third loop on cardiac gates
        for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
        {
          // Mask slices
          for (int z=0; z<removed_slices; z++)
          {
            // First slices
            INTNB base_z_first = z*mp_ID->GetNbVoxXY();
            // Loop over Y and X
            for (int i=0; i<mp_ID->GetNbVoxXY(); i++)
            {
              INTNB index = base_z_first + i;
              m4p_outputImage[fr][rg][cg][index] = 0.;
            }
            // Last slices
            INTNB base_z_last = (mp_ID->GetNbVoxZ()-1-z)*mp_ID->GetNbVoxXY();
            // Loop over Y and X
            for (int i=0; i<mp_ID->GetNbVoxXY(); i++)
            {
              INTNB index = base_z_last + i;
              m4p_outputImage[fr][rg][cg][index] = 0.;
            }
          }
        }
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

int oImageSpace::ApplyOutputMaskImage()
{
  // If no mask image return
  if ( !IsLoadedMask() ) return 0;
  // Verbose
  if (m_verbose>=2)
  {
    Cout("oImageSpace::ApplyOutputMaskImage() -> Mask output image with the provided input mask image" << endl);
  }
  // -----------------------------------------------
  // Masking
  // -----------------------------------------------
  int v;
  #pragma omp parallel for private(v) schedule(guided)
  for (v=0; v<mp_ID->GetNbVoxXYZ(); v++)
  {
    // First loop on frames
    for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
      // Second loop on respiratory gates
      for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
        // Third loop on cardiac gates
        for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
        {
          if (!(mp_maskImage[v]>0.)) m4p_outputImage[fr][rg][cg][v] = 0.;
        }
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageSpace::ApplyMaskToSensitivity()
{
  // If no mask image return
  if ( !IsLoadedMask() ) return 0;

  // Verbose
  if (m_verbose>=2)
  {
    Cout("oImageSpace::ApplyMaskToSensitivity() -> Mask the sensitivity image with the provided input mask image" << endl);
  }
  // -----------------------------------------------
  // Masking
  // -----------------------------------------------
  int thread_0 = 0;
  int v;
  #pragma omp parallel for private(v) schedule(guided)
  for (v=0; v<mp_ID->GetNbVoxXYZ(); v++)
  {
    // First loop on frames
    for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
      // Second loop on respiratory gates
      for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
        // Third loop on cardiac gates
        for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
        {
          if (!(mp_maskImage[v]>0.)) m5p_sensitivity[thread_0][fr][rg][cg][v] = 0.;
        }
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oImageSpace::ApplyMaskToBackwardImage(int a_imageIndex, int a_timeIndex, int a_respIndex, int a_cardIndex)
{
  // If no mask image return
  if ( !IsLoadedMask() ) return 0;

  // Verbose
  if (m_verbose>=2)
  {
    Cout("oImageSpace::ApplyMaskToBackwardImage() -> Mask the backward image with the provided input mask image" << endl);
  }
  // -----------------------------------------------
  // Masking
  // -----------------------------------------------
  int thread_0 = 0;
  int v;
  #pragma omp parallel for private(v) schedule(guided)
  for (v=0; v<mp_ID->GetNbVoxXYZ(); v++)
  {
    // First loop on frames
    for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
      // Second loop on respiratory gates
      for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
        // Third loop on cardiac gates
        for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
        {
          if (!(mp_maskImage[v]>0.)) m6p_backwardImage[a_imageIndex][thread_0][a_timeIndex][a_respIndex][a_cardIndex][v] = 0.;
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
  \fn SaveOutputImage
  \param a_iteration : current iteration index
  \param   a_subset : current number of subsets (or -1 by default)
  \brief Call the interfile function to write output image on disk
  \return  0 if success, positive value otherwise
*/
int oImageSpace::SaveOutputImage(int a_iteration, int a_subset)
{
  if(m_verbose>=3) Cout("oImageSpace::SaveOutputImage ..."<< endl);
    
  // Get the output manager
  sOutputManager* p_output_manager = sOutputManager::GetInstance();

  // Interfile writing
  // m4p_forwardImage is used here for the special case temporal regularization using intrinsic basis functions is used.
  // In this case, m4p_image contains basis function coefficients and not image values. These are computed in ComputeOutputImage()
  // In any other case, both matrices contain the image values.

  string path_to_img = p_output_manager->GetPathName() + p_output_manager->GetBaseName();

  // Add a suffix for iteration
  if (a_iteration >= 0)
  {
    stringstream ss; ss << a_iteration + 1;
    path_to_img.append("_it").append(ss.str());
  }

  // Add a suffix for subset (if not negative by default), this means that we save a 'subset' image
  if (a_subset >= 0)
  {
    stringstream ss; ss << a_subset + 1;
    path_to_img.append("_ss").append(ss.str());
  }

  // We need one "mode" parameter which indicates if we would like to write 1 image file or several
  // Adjust functions according to that
  if(IntfWriteImgFile(path_to_img, m4p_outputImage, mp_ID, m_verbose) )
  {
    Cerr("***** oImageSpace::SaveOutputImage()-> Error writing Interfile of output image !" << endl);  
    return 1;
  }
  

  // Normal end
  return 0;
}





// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn SaveDebugImage
  \param a_name : output name of the image
  \brief Just a debug function dedicated to write any kind of image on disk in raw format, for debugging purposes
*/
void oImageSpace::SaveDebugImage(const string& a_name)
{
  #ifdef CASTOR_DEBUG
  if(m_verbose>=3) Cout("oImageSpace::SaveDebugImage ..."<< endl);
  #endif
  
  ofstream output_file;
  output_file.open(a_name.c_str(), ios::binary | ios::out);

  //for (int fr=0; fr<ap_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
  //  for (int rg=0; rg<ap_ImageDimensionsAndQuantification->GetNbRespGates(); rg++)
  //    for (int cg=0; cg<ap_ImageDimensionsAndQuantification->GetNbCardGates(); cg++)
  //    {
        // Write file
        output_file.write(reinterpret_cast<char*>(m6p_backwardImage[0][0][0][0][0]), mp_ID->GetNbVoxXYZ()*sizeof(FLTNB));
        // Close file
        output_file.close();
  //    }
  //  }
  //}
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PrepareForwardImage
  \brief Copy current image matrix in the forward-image buffer matrix
*/
void oImageSpace::PrepareForwardImage()
{
  if(m_verbose>=3) Cout("oImageSpace::PrepareForwardImage ..."<< endl);

  for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
    for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
      for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
        for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
          m4p_forwardImage[tbf][rbf][cbf][v] = m4p_image[tbf][rbf][cbf][v];
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn Reduce
  \brief Merge parallel results into the matrix of the backward image matrix of the first thread. Also for MPI.
  \todo Choose computation alternative (do not know yet which one is faster...)
*/
void oImageSpace::Reduce()
{
  if(m_verbose>=3) Cout("oImageSpace::Reduce ..."<< endl);
    
  #if defined(CASTOR_OMP) || defined(CASTOR_MPI)
  if (m_verbose>=3) Cout("oImageSpace::Reduce() -> Merge parallel results" << endl);
  #endif

  // --------------------------------------------------------------------------------
  // Step 1: merge multi-threads results
  // --------------------------------------------------------------------------------

  #ifdef CASTOR_OMP
  if (m_verbose>=3) Cout("  --> Over threads ..." << endl);

  // Special case here where it appears to be always beneficial to use as many threads
  // as the number of threads used for projections. So we set it here, and set it back
  // at the end.
  omp_set_num_threads(mp_ID->GetNbThreadsForProjection());

  // todo : Choose computation alternative (do not know yet which one is faster...)
  int alternative = 2;

  // Alternative 1: standard loops based on the pointers hierarchy but mono-thread
  if (alternative==1 || mp_ID->GetNbThreadsForImageComputation()==1)
  {
    for (int img=0; img<m_nbBackwardImages; img++)
      for (int th=1; th<mp_ID->GetNbThreadsForProjection(); th++) 
        for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
          for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
            for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
              for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
                m6p_backwardImage[img][0][tbf][rbf][cbf][v] += m6p_backwardImage[img][th][tbf][rbf][cbf][v];
  }
  // Alternative 2: multi-thread merging with the voxels loop at first to be thread safe
  // Maybe faster than alternative 1, but the first loop on voxels breaks the memory order... have to try
  else if (alternative==2)
  {
    int v;
    #pragma omp parallel for private(v) schedule(guided)
    for (v=0; v<mp_ID->GetNbVoxXYZ(); v++)
    {
      for (int img=0; img<m_nbBackwardImages; img++)
        for (int th=1; th<mp_ID->GetNbThreadsForProjection(); th++) 
          for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
            for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
              for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
                m6p_backwardImage[img][0][tbf][rbf][cbf][v] += m6p_backwardImage[img][th][tbf][rbf][cbf][v];
    }
  }
  // If on-the-fly sensitivity, then do it also for it.
  if (!m_loadedSensitivity)
  {
    if (alternative==1 || mp_ID->GetNbThreadsForImageComputation()==1)
    {
      for (int th=1; th<mp_ID->GetNbThreadsForProjection(); th++)
        for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
          for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
            for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
              for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
                m5p_sensitivity[0][fr][rg][cg][v] += m5p_sensitivity[th][fr][rg][cg][v];
    }
    else if (alternative==2)
    {
      int v;
      #pragma omp parallel for private(v) schedule(guided)
      for (v=0; v<mp_ID->GetNbVoxXYZ(); v++)
      {
        for (int th=1; th<mp_ID->GetNbThreadsForProjection(); th++)
          for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
            for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
              for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
                m5p_sensitivity[0][fr][rg][cg][v] += m5p_sensitivity[th][fr][rg][cg][v];
      }
    }
  }
  // Set back the number of threads to the one for image computation
  omp_set_num_threads(mp_ID->GetNbThreadsForImageComputation());
  #endif

  // --------------------------------------------------------------------------------
  // Step 2: merge multi-instance results (MPI)
  // --------------------------------------------------------------------------------

  #ifdef CASTOR_MPI
  // We use the MPI_IN_PLACE as the send buffer so that it is the same as the result buffer
  if (m_verbose>=3) Cout("  --> Over instances ..." << endl);

  // Merge backward images
  for (int img=0; img<m_nbBackwardImages; img++)
    for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
      for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
        for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
          MPI_Allreduce( MPI_IN_PLACE, &m6p_backwardImage[img][0][tbf][rbf][cbf][0], mp_ID->GetNbVoxXYZ(), FLTNBMPI, MPI_SUM, MPI_COMM_WORLD );

  // If on-the-fly sensitivity, then do it also for it
  if (!m_loadedSensitivity)
  {
    for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
      for (int rg=0; rg<mp_ID->GetNbRespGates(); rg++)
        for (int cg=0; cg<mp_ID->GetNbCardGates(); cg++)
          MPI_Allreduce( MPI_IN_PLACE, &m5p_sensitivity[0][fr][rg][cg][0], mp_ID->GetNbVoxXYZ(), FLTNBMPI, MPI_SUM, MPI_COMM_WORLD );
  }
  #endif
}





// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn CleanNeverVisitedVoxels
  \brief Based on the visitedVoxelsImage, clean the never visited voxels in the image.
        This function must be called at the end of each iteration
*/
void oImageSpace::CleanNeverVisitedVoxels()
{
  if(m_verbose>=3) Cout("oImageSpace::CleanNeverVisitedVoxels ..."<< endl);
  
  /*
  int v;
  #pragma omp parallel for private(v) schedule(guided)
  for (v=0; v<mp_ID->GetNbVoxXYZ(); v++)
  {
    // Clean this image voxel if the voxel was never visited
    
    if (mp_visitedVoxelsImage[v]==0.)
    {
      for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
        for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
          for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
          {
            m4p_image[tbf][rbf][cbf][v] = 0.;
          }
    }
    // Reset the image
    mp_visitedVoxelsImage[v] = 0.;
  }
  */

  // Clean this image voxel if the voxel was never visited
  for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
    for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
      for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
      {
        int v;
        #pragma omp parallel for private(v) schedule(guided)
        for (v=0; v<mp_ID->GetNbVoxXYZ(); v++)
        {
          if (m4p_visitedVoxelsImage[tbf][rbf][cbf][v]==0.) m4p_image[tbf][rbf][cbf][v] = 0.;
          
          // Reset the image
          m4p_visitedVoxelsImage[tbf][rbf][cbf][v] = 0.;
        }
      }
}









// LIST-MODE SENSITIVITY GENERATION FUNCTIONS
// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LMS_InstantiateImage()
  \brief Allocate memory for the main image matrices (for list-mode sensitivity generation)
  \details This function is dedicated to list-mode sensitivity (LMS) generation.
*/
void oImageSpace::LMS_InstantiateImage()
{
  if(m_verbose>=3) Cout("oImageSpace::LMS_InstantiateImage ..."<< endl);
  
  m4p_image = new FLTNB***[mp_ID->GetNbTimeFrames()];

  for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
  {
    m4p_image[fr] = new FLTNB**[mp_ID->GetNb1stMotImgsForLMS(fr)];
    
    for(int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr) ; rg++)
    {
      m4p_image[fr][rg] = new FLTNB*[mp_ID->GetNb2ndMotImgsForLMS()];
      
      for(int cg=0 ; cg<mp_ID->GetNb2ndMotImgsForLMS() ; cg++)
        m4p_image[fr][rg][cg] = new FLTNB[mp_ID->GetNbVoxXYZ()];
    }
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LMS_DeallocateImage()
  \brief Free memory for the main image matrices (for list-mode sensitivity generation)
  \details This function is dedicated to list-mode sensitivity (LMS) generation.
*/
void oImageSpace::LMS_DeallocateImage()
{
  if(m_verbose>=3) Cout("oImageSpace::LMS_DeallocateImage ..."<< endl);
  
  for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
  {
    for(int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr) ; rg++)
    {
      for(int cg=0 ; cg<mp_ID->GetNb2ndMotImgsForLMS() ; cg++)
        if(m4p_image[fr][rg][cg] != NULL) delete m4p_image[fr][rg][cg];

      if(m4p_image[fr][rg] != NULL) delete[] m4p_image[fr][rg];
    }
    if(m4p_image[fr] != NULL) delete[] m4p_image[fr] ;
  }

  if(m4p_image != NULL) delete[] m4p_image;
  m4p_image = NULL;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LMS_InstantiateForwardImage()
  \brief Allocate memory for the forward image matrices (for list-mode sensitivity generation)
  \details This function is dedicated to list-mode sensitivity (LMS) generation.
*/
void oImageSpace::LMS_InstantiateForwardImage()
{
  if(m_verbose>=3) Cout("oImageSpace::LMS_InstantiateForwardImage ..."<< endl);
  
  m4p_forwardImage = new FLTNB***[mp_ID->GetNbTimeFrames()];
    
  for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
  {
    m4p_forwardImage[fr] = new FLTNB**[mp_ID->GetNb1stMotImgsForLMS(fr)];
    
    for(int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr) ; rg++)
    {
      m4p_forwardImage[fr][rg] = new FLTNB*[mp_ID->GetNb2ndMotImgsForLMS()];
      
      for(int cg=0 ; cg<mp_ID->GetNb2ndMotImgsForLMS() ; cg++)
      {
        m4p_forwardImage[fr][rg][cg] = new FLTNB[mp_ID->GetNbVoxXYZ()];
      }
    }
    
  }
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LMS_DeallocateForwardImage()
  \brief Free memory for the forward image matrices (for list-mode sensitivity generation)
  \details This function is dedicated to list-mode sensitivity (LMS) generation.
*/
void oImageSpace::LMS_DeallocateForwardImage()
{
  if(m_verbose>=3) Cout("oImageSpace::LMS_DeallocateForwardImage ..."<< endl);
  
  for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
  {
    for(int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr) ; rg++)
    {
      for(int cg=0 ; cg<mp_ID->GetNb2ndMotImgsForLMS() ; cg++)
      {
        if(m4p_forwardImage[fr][rg][cg] != NULL) delete m4p_forwardImage[fr][rg][cg];
      }

      if(m4p_forwardImage[fr][rg] != NULL) delete[] m4p_forwardImage[fr][rg];
    }
  
    if(m4p_forwardImage[fr] != NULL) delete[] m4p_forwardImage[fr];
  }

  if(m4p_forwardImage != NULL) delete[] m4p_forwardImage;
  m4p_forwardImage = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LMS_InstantiateSensitivityImage()
  \brief Allocate memory for the sensitivity image matrices (for list-mode sensitivity generation)
  \details This function is dedicated to list-mode sensitivity (LMS) generation.
*/
void oImageSpace::LMS_InstantiateSensitivityImage()
{
  if(m_verbose>=3) Cout("oImageSpace::LMS_InstantiateSensitivityImage ..."<< endl);
  
  m5p_sensitivity = new FLTNB****[1];
  m5p_sensitivity[0] = new FLTNB***[mp_ID->GetNbTimeFrames()];
  for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
  {
    m5p_sensitivity[0][fr] = new FLTNB**[mp_ID->GetNb1stMotImgsForLMS(fr)];
    for (int rg=0; rg<mp_ID->GetNb1stMotImgsForLMS(fr); rg++)
    {
      m5p_sensitivity[0][fr][rg] = new FLTNB*[mp_ID->GetNb2ndMotImgsForLMS()];
      for (int cg=0; cg<mp_ID->GetNb2ndMotImgsForLMS(); cg++)
        m5p_sensitivity[0][fr][rg][cg] = new FLTNB[mp_ID->GetNbVoxXYZ()];
    }
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LMS_DeallocateSensitivityImage()
  \brief Free memory for the sensitivity image matrices (for list-mode sensitivity generation)
  \details This function is dedicated to list-mode sensitivity (LMS) generation.
*/
void oImageSpace::LMS_DeallocateSensitivityImage()
{
  if(m_verbose>=3) Cout("oImageSpace::LMS_DeallocateSensitivityImage ..."<< endl);
  
  for (int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
  {
    for (int rg=0; rg<mp_ID->GetNb1stMotImgsForLMS(fr); rg++)
    {
      for (int cg=0; cg<mp_ID->GetNb2ndMotImgsForLMS(); cg++)
      {
        if(m5p_sensitivity[0][fr][rg][cg] != NULL) delete m5p_sensitivity[0][fr][rg][cg];
      }
      if(m5p_sensitivity[0][fr][rg] != NULL)delete m5p_sensitivity[0][fr][rg];
    }
    if(m5p_sensitivity[0][fr] != NULL)delete m5p_sensitivity[0][fr];
  }
  if(m5p_sensitivity[0] != NULL) delete[] m5p_sensitivity[0];
  
  if(m5p_sensitivity != NULL) delete[] m5p_sensitivity;
  m5p_sensitivity = NULL;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LMS_CopyAtnToImage()
  \brief Copy the attenuation image contained in the 'm2p_attenuation' matrix inside the m2p_image matrix (for list-mode sensitivity generation)
  \details This function is dedicated to list-mode sensitivity (LMS) generation.
*/
void oImageSpace::LMS_CopyAtnToImage()
{
  if(m_verbose>=3) Cout("oImageSpace::LMS_CopyAtnToImage ..."<< endl);
  
  if(m4p_attenuation != NULL)
    for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
      for(int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr) ; rg++)
        for(int cg=0 ; cg<mp_ID->GetNb2ndMotImgsForLMS() ; cg++)
          for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
          {
            m4p_image[fr][rg][cg][v] = m4p_attenuation[fr][rg][cg][v];
          }
  else
    for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
      for(int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr) ; rg++)
        for(int cg=0 ; cg<mp_ID->GetNb2ndMotImgsForLMS() ; cg++)
          for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
          {
            m4p_image[fr][rg][cg][v] = 0;
          }

}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LMS_CopyAtnToForwardImage()
  \param a_use1stMotion
  \param a_use2ndMotion
  \brief Copy the attenuation image contained in the 'm2p_attenuation' matrix inside the m2p_forwardImage matrix (for list-mode sensitivity generation)
  \details This function is dedicated to list-mode sensitivity (LMS) generation.
*/
void oImageSpace::LMS_CopyAtnToForwardImage(bool a_use1stMotion, bool a_use2ndMotion)
{
  if(m_verbose>=3) Cout("oImageSpace::LMS_CopyAtnToForwardImage ..."<< endl);
  
  if(m4p_attenuation != NULL)
    for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
      for(int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr) ; rg++)
        for(int cg=0 ; cg<mp_ID->GetNb2ndMotImgsForLMS() ; cg++)
        {
          // If motion is enabled, copy the attenuation image in all forward image matrices
          int rg_atn = a_use1stMotion ? 0 : rg ;
          int cg_atn = a_use1stMotion ? 0 : cg ;
          
          for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
          {
            m4p_forwardImage[fr][rg][cg][v] = m4p_attenuation[fr][rg_atn][cg_atn][v];
          }
        }
  else
    for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
      for(int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr) ; rg++)
        for(int cg=0 ; cg<mp_ID->GetNb2ndMotImgsForLMS() ; cg++)
          for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
          {
            m4p_forwardImage[fr][rg][cg][v] = 0;
          }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LMS_CopyBackwardToSensitivityImage()
  \brief Copy the backward image containing the result of the sensitivity back-projection, into the sensitivity matrix.
  \details This function is dedicated to list-mode sensitivity (LMS) generation, it is useful to apply convolution and to save the image.
*/
void oImageSpace::LMS_CopyBackwardToSensitivity()
{
  if(m_verbose>=3) Cout("oImageSpace::LMS_CopyBackwardToSensitivity ..."<< endl);
  for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
    for (int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr) ; rg++)
      for (int cg=0 ; cg<mp_ID->GetNb2ndMotImgsForLMS() ; cg++)
        for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
          m5p_sensitivity[0][fr][rg][cg][v] = m6p_backwardImage[0][0][fr][rg][cg][v];
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LMS_PrepareForwardImage()
  \brief Copy current image in forward-image buffer (for list-mode sensitivity generation)
  \details This function is dedicated to list-mode sensitivity (LMS) generation.
*/
void oImageSpace::LMS_PrepareForwardImage()
{
  if(m_verbose>=3) Cout("oImageSpace::LMS_PrepareForwardImage ..."<< endl);
  
  for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
    for (int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr) ; rg++)
      for (int cg=0 ; cg<mp_ID->GetNb2ndMotImgsForLMS() ; cg++)
        for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
          m4p_forwardImage[fr][rg][cg][v] = m4p_image[fr][rg][cg][v];
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void oImageSpace::ReduceBackwardImage(int a_imageIndex, int a_timeIndex, int a_respIndex, int a_cardIndex)
{
  #if defined(CASTOR_OMP) || defined(CASTOR_MPI)
  if (m_verbose>=3) Cout("oImageSpace::ReduceBackwardImage() -> Merge parallel results" << endl);
  #endif

  // --------------------------------------------------------------------------------
  // Step 1: merge multi-threads results
  // --------------------------------------------------------------------------------

  #ifdef CASTOR_OMP
  // Verbose
  if (m_verbose>=3) Cout("  --> Over threads ..." << endl);
  // Do it
  for (int th=1 ; th<mp_ID->GetNbThreadsForProjection() ; th++) 
    for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
      m6p_backwardImage[a_imageIndex][0][a_timeIndex][a_respIndex][a_cardIndex][v] += m6p_backwardImage[a_imageIndex][th][a_timeIndex][a_respIndex][a_cardIndex][v];
  #endif

  // --------------------------------------------------------------------------------
  // Step 2: merge multi-instance results (MPI)
  // --------------------------------------------------------------------------------

  #ifdef CASTOR_MPI
  // We use the MPI_IN_PLACE as the send buffer so that it is the same as the result buffer
  if (m_verbose>=3) Cout("  --> Over instances ..." << endl);
  // Merge backward images
  MPI_Allreduce( MPI_IN_PLACE, &m6p_backwardImage[a_imageIndex][0][a_timeIndex][a_respIndex][a_cardIndex][0], mp_ID->GetNbVoxXYZ(), FLTNBMPI, MPI_SUM, MPI_COMM_WORLD );
  #endif
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LMS_SaveSensitivityImage
  \param a_pathToSensitivityImage : path to the sensitivity image (should be provided only in list-mode reconstruction)
  \param ap_DeformationManager : Pointer to the deformation manager objet (required to retrieve the number of gates in the sensitivity image) 
  \brief Call the interfile function to write the sensitivity image on disk
  \details If image deformation is enabled for respiratory/cardiac gated data, the gated images are summed up into one image and normalize
  \return  0 if success, positive value otherwise
  \todo Interfile management : if the number of sensitivity image to write is different to the actual number of image
                               as it could be the case for dynamic imaging, if one sensitivity image is required for the whole dynamic serie,
                               we will have to give this information to IntfWriteImgFile(), otherwise the Interfile functions will try to write several times the image (leading to segfault or error)
*/
int oImageSpace::LMS_SaveSensitivityImage(const string& a_pathToSensitivityImage, oDeformationManager* ap_DeformationManager)
{
  if(m_verbose>=3) Cout("oImageSpace::LMS_SaveSensitivityImage ..."<< endl);
  
  //sOutputManager* p_output_manager = sOutputManager::GetInstance();
  
  // According to the respiratory motion algorithm, images relatives to resp gates, cardiac gates, or both, have to be summed up and normalize into one image.
  // The resulting image will be summed up into the sensitivity matrix

  int nb_reco_card_images = ap_DeformationManager->GetNbSensImagesCardDeformation(mp_ID->GetNb2ndMotImgsForLMS());    
    
  if (ap_DeformationManager->GetNbSensImagesCardDeformation(mp_ID->GetNb2ndMotImgsForLMS()) == 0) // Average of the card images 
  {
    nb_reco_card_images = 1;
    FLTNB card_normalization = mp_ID->GetNb2ndMotImgsForLMS();
    
    for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
      for(int fr=0 ; fr<mp_ID->GetNbTimeFrames(); fr++)
        for (int rg=0 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr) ; rg++)
        {
          for (int cg=1 ; cg<mp_ID->GetNb2ndMotImgsForLMS() ; cg++)
            m5p_sensitivity[0][fr][rg][0][v] += m5p_sensitivity[0][fr][rg][cg][v];
          m5p_sensitivity[0][fr][rg][0][v] /= card_normalization;
        }
  }
  
  if (ap_DeformationManager->GetNbSensImagesRespDeformation(mp_ID->GetNb1stMotImgsForLMS(0)) == 0) // Average of the resp images 
  {
    for(int fr=0 ; fr<mp_ID->GetNbTimeFrames(); fr++)
    {
      FLTNB resp_normalization = mp_ID->GetNb1stMotImgsForLMS(fr)-mp_ID->GetPMotionFirstIndexForLMS(fr);

      for (int cg=0 ; cg<nb_reco_card_images ; cg++) // Dual gating : Be sure to get the right number of card images
      {          
        for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
        {
          for (int rg=1 ; rg<mp_ID->GetNb1stMotImgsForLMS(fr)-mp_ID->GetPMotionFirstIndexForLMS(fr) ; rg++)
            m5p_sensitivity[0][fr][0][cg][v] += m5p_sensitivity[0][fr][rg][cg][v];
          
          m5p_sensitivity[0][fr][0][cg][v] /= resp_normalization;
        }
      }
    }
  }

  if (SaveSensitivityImage(a_pathToSensitivityImage) )
  {
    Cerr("***** oImageSpace::LMS_SaveSensitivityImage()-> Error writing Sensitivity image !" << endl);  
    return 1;
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn SaveSensitivityImage
  \param a_pathToSensitivityImage : path to the sensitivity image (should be provided only in list-mode reconstruction)
  \brief Call the interfile function to write the sensitivity image on disk
  \return  0 if success, positive value otherwise
*/
int oImageSpace::SaveSensitivityImage(const string& a_pathToSensitivityImage)
{
  if(m_verbose>=3) Cout("oImageSpace::SaveSensitivityImage ..."<< endl);
  
  if (IntfWriteImgFile(a_pathToSensitivityImage, m5p_sensitivity[0], mp_ID, m_verbose) )
  {
    Cerr("***** oImageSpace::SaveSensitivityImage()-> Error writing Interfile of output image !" << endl);  
    return 1;
  }
  
  return 0;
}




// PROJECTION SCRIPT FUNCTIONS
// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_InstantiateProjectionImage()
  \param a_nbProjs : a number of projection slices in the projection
  \param a_nbPixels : a total number of pixels in the projection slices
  \brief Instanciate and initialize projection image for analytical projection
  \details This function is currently only dedicated to SPECT projection, and used by the analytical projection script
  \todo support for PET (requires projection format)
*/
void oImageSpace::PROJ_InstantiateProjectionImage(int a_nbProjs, int a_nbPixels)
{
  if(m_verbose>=3) Cout("oImageSpace::PROJ_InstantiateProjectionImage ..."<< endl);
  
  m2p_projectionImage = new FLTNB*[a_nbProjs];
  
  for(int p=0 ; p<a_nbProjs ; p++)
  {
    m2p_projectionImage[p] = new FLTNB[a_nbPixels];
    for(int px=0 ; px<a_nbPixels ; px++)
      m2p_projectionImage[p][px] = 0;
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_DeallocateProjectionImage()
  \param a_nbProjs : a number of projection slices in the projection
  \brief Free memory for the projection image for analytical projection
  \details This function is currently only dedicated to SPECT projection, and used by the analytical projection script
  \todo support for PET (requires projection format)
*/
void oImageSpace::PROJ_DeallocateProjectionImage(int a_nbProjs)
{
  if(m_verbose>=3) Cout("oImageSpace::PROJ_DeallocateProjectionImage ..."<< endl);
  
  for(int p=0 ; p<a_nbProjs ; p++)
  {
    delete[] m2p_projectionImage[p];
  }
  
  delete[] m2p_projectionImage;
  m2p_projectionImage = NULL;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_InitImage()
  \param a_pathToInitialImage : path to the image to project
  \brief Load the initial image for the analytical projection
  \return  0 if success, positive value otherwise
*/
int oImageSpace::PROJ_InitImage(const string& a_pathToInitialImage)
{
  if(m_verbose>=3) Cout("oImageSpace::PROJ_InitImage ..."<< endl);
  
  if (!a_pathToInitialImage.empty()) // Image initiale
  {
    if(PROJ_LoadInitialImage(a_pathToInitialImage) )
    {
      Cerr("***** oImageSpace::PROJ_InitImage()-> Error while trying to read file at " << a_pathToInitialImage << endl);
      return 1;
    }
  }
  else
  {
    {
      Cerr("***** oImageSpace::PROJ_InitImage()-> No projected image provided ! " << endl);
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
  \fn PROJ_LoadInitialImage()
  \param a_pathToImage : path to the image to project
  \brief Load the initial image for the analytical projection
  \return  0 if success, positive value otherwise
*/
int oImageSpace::PROJ_LoadInitialImage(const string& a_pathToImage)
{
  if(m_verbose>=3) Cout("oImageSpace::PROJ_LoadInitialImage ..."<< endl);

  ifstream image_file;
  image_file.open(a_pathToImage.c_str(), ios::in | ios::binary);

  if(!image_file.is_open()) 
  {
    Cerr("***** oImageSpace::PROJ_LoadInitialImage()-> Error reading file !" << endl);  
    return 1;
  }
  else
  {
    // Interfile image reading (INTF_LERP_DISABLED = no interpolation allowed)
    if(IntfReadImage(a_pathToImage, m4p_image, mp_ID, m_verbose, INTF_LERP_DISABLED))
    {
      Cerr("***** oImageSpace::PROJ_LoadInitialImage()-> Error reading Interfile : " << a_pathToImage << " !" << endl);  
      return 1;
    }
    
  }
  image_file.close();
  return 0;
}





// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PROJ_SaveProjectionImage
  \brief Save an image of the projected data
  \return  0 if success, positive value otherwise
  \todo Support for PET projeted data
        implement interfile file recovery in a scanner class function
*/
int oImageSpace::PROJ_SaveProjectionImage()
{
  if(m_verbose>=3) Cout("oImageSpace::PROJ_SaveProjectionImage ..."<< endl);
  
  string img_name = "_ProjectionImage";

  // Initialize interfile fields structure
  Intf_fields Img_fields;
  IntfKeySetFieldsOutput(&Img_fields , mp_ID);
  
  Img_fields.process_status = "acquired";
  
  // Get specific info about projection
  sScannerManager* p_ScanMgr = sScannerManager::GetInstance();
  
  // Recover the values in the interfile structure.
  // todo : This step for PET projection
  // todo : Implement this step in scanner class functions
  // todo : For spect, deal with systems with several detector heads (currently assume one head)
  if(p_ScanMgr->GetScannerType() == 2)
  {
    uint16_t nb_projs, nb_heads;
    FLTNB zoom;
    uint16_t nb_bins[2];
    FLTNB  pix_size[2];
    FLTNB* p_angles;
    FLTNB* p_radius;
    int dir_rot;
    p_ScanMgr->GetSPECTSpecificParameters(&nb_projs, 
                                          &nb_heads,
                                              &zoom,
                                            nb_bins, 
                                           pix_size,
                                           p_angles, 
                                           p_radius,
                                          &dir_rot);

    Img_fields.nb_projections = nb_projs;
    Img_fields.nb_detector_heads = nb_heads;
    Img_fields.mtx_size[0] = nb_bins[0];
    Img_fields.mtx_size[1] = nb_bins[1];
    Img_fields.vox_size[0] = pix_size[0];
    Img_fields.vox_size[1] = pix_size[1];
    Img_fields.extent_rotation = 360; // 360 by default
    // Check rotation direction from the two first angles
    Img_fields.direction_rotation = (dir_rot == GEO_ROT_CCW)?
                                     "CCW" :
                                     "CW";
    Img_fields.first_angle = p_angles[0];
    
    // Note: In Interfile v3.3, doesn't seem to have a field to provide
    // individually each angle and Center of rotation.
    // Looks like it was planned for ulterior version, at least for the
    // center of rotation
    // For now, we just write each projection angle and radius as an array key
 
    string angles_str = "{";
    string radius_str = "{";
    
    // Flag to check if the radius is similar for each projection angle
    bool has_single_radius = true;
    
    for(uint16_t p=0 ; p<nb_projs ; p++)
    {
      stringstream ss_a, ss_r;
      // projection angles
      ss_a << p_angles[p];
      angles_str.append(ss_a.str());
      (p<nb_projs-1) ? angles_str.append(",") : angles_str.append("}");
      
      // radius
      ss_r << p_radius[p];
      radius_str.append(ss_r.str());
      (p<nb_projs-1) ? radius_str.append(",") : radius_str.append("}");
      if(p_radius[p] != p_radius[0])
        has_single_radius = false;
      
      // Some interfile editors struggle with long line, perhaps because
      // they are stored in char[256] or something like that
      // Hence, break line after a certain number of elements
      if((p+1)%30 == 0) 
      {
        angles_str.append("\n");
        radius_str.append("\n");
      }
    }
    
    // If there is one common radius, write its value in the related string
    if(has_single_radius) 
    {
      stringstream ss;
      ss << p_radius[0];
      radius_str = ss.str();
    }
    
    Img_fields.projection_angles = angles_str;
    Img_fields.radius = radius_str;
  

    // Common code to all modalities
    sOutputManager* p_output_manager = sOutputManager::GetInstance();
    string img_file_name = p_output_manager->GetPathName() + p_output_manager->GetBaseName();
    img_file_name.append("_ProjImage");
  
    if(IntfWriteProjFile(img_file_name, m2p_projectionImage, mp_ID, Img_fields, m_verbose) )
    {
      Cerr("***** oImageSpace::PROJ_SaveProjectionImage()-> Error writing Interfile of output image !" << endl);  
      return 1;
    }
  }
  return 0;
} 


