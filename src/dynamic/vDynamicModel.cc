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
  \ingroup  dynamic
  \brief    Implementation of class vDynamicModel
*/

#include "vDynamicModel.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn vDynamicModel
  \brief Constructor of vDynamicModel. Simply set all data members to default values.
*/
vDynamicModel::vDynamicModel() 
{
  mp_ID = NULL;
  m_nbTimeBF = -1; 
  m_nbModelParam = -1;
  m_verbose = -1;
  m_checked = false;
  m_initialized = false;
  m2p_parametricImages = NULL;
  m2p_modelTACs = NULL;
  m2p_outputParImages = NULL;
  m_saveParImageFlag = true;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn ~vDynamicModel
  \brief Destructor of vDynamicModel.
*/
vDynamicModel::~vDynamicModel() {}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn CheckParameters
  \brief This function is used to check parameters after the latter
         have been all set using Set functions.
  \return 0 if success, positive value otherwise.
*/
int vDynamicModel::CheckParameters()
{
  if(m_verbose>=2) Cout("vDynamicModel::CheckParameters ..."<< endl); 
    
  // Check image dimensions
  if (mp_ID==NULL)
  {
    Cerr("***** vDynamicModel::CheckParameters() -> No image dimensions provided !" << endl);
    return 1;
  }
  
  // Check verbosity
  if (m_verbose<0)
  {
    Cerr("***** vDynamicModel::CheckParameters() -> Wrong verbosity level provided !" << endl);
    return 1;
  }

  // Check number of basis functions
  if (m_nbTimeBF <0)
  {
    Cerr("***** vDynamicModel::CheckParameters() -> Basis functions number has not been initialized !" << endl);
    return 1;
  }

  // Check number of parameters in the model
  if (m_nbModelParam <0)
  {
    Cerr("***** vDynamicModel::CheckParameters() -> Basis functions number has not been initialized !" << endl);
    return 1;
  }
  
  // Check parameters of the child class (if this function is overloaded)
  if (CheckSpecificParameters())
  {
    Cerr("***** vDynamicModel::CheckParameters() -> An error occurred while checking parameters of the child dynamic class !" << endl);
    return 1;
  }
  
  // Normal end
  m_checked = true;
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      vDynamicModel::ComputeOutputParImage
  \brief   Compute output image using the m2p_parametricImages matrix Store the result in the m2p_outputParImages matrix
*/
void vDynamicModel::ComputeOutputParImage()
{
  if(m_verbose >=2) Cout("vDynamicModel::ComputeOutputParImage ..." <<endl);
  
  // If we save parametric image
  if(m_saveParImageFlag)
    for(int p=0 ; p<m_nbModelParam ; p++)
    {
      // If output image matrix is allocated, then copy the current parametric images
      if(m2p_outputParImages    != NULL 
      && m2p_outputParImages[p] != m2p_parametricImages[p])
      {
        for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
          m2p_outputParImages[p][v] = m2p_parametricImages[p][v] ;
      }
      // Point to parametric images otherwise
      else
        m2p_outputParImages = m2p_parametricImages;
    }
    
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      vDynamicModel::ApplyOutputFOVMaskingOnParametricImages
  \brief   Mask the outside of the transaxial FOV based on the m_fovOutPercent
  \details Similar to the eponym function in ImageSpace, but on parametric images
*/
int vDynamicModel::ApplyOutputFOVMaskingOnParametricImages()
{
  // If the output FOV percent is under 0 (the default value) and number of axial slices to be removed is 0, we do not mask anything
  if ( mp_ID->GetFOVOutPercent()<=0. &&
       mp_ID->GetNbSliceOutMask()==0 ) return 0;
       
  // Verbose
  if (m_verbose>=2)
  {
    Cout("vDynamicModel::ApplyOutputFOVMaskingOnParametricImages() -> Mask output image" << endl);
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
          
          for (int b=0 ; b<m_nbTimeBF ; b++)
            m2p_outputParImages[b][index] = 0.;
            
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
        
    // Mask slices
    for (int b=0 ; b<m_nbTimeBF ; b++)
      for (int z=0; z<removed_slices; z++)
      {
        // First slices
        INTNB base_z_first = z*mp_ID->GetNbVoxXY();
        // Loop over Y and X
        for (int i=0; i<mp_ID->GetNbVoxXY(); i++)
        {
          INTNB index = base_z_first + i;
          m2p_outputParImages[b][index] = 0.;
        }
        // Last slices
        INTNB base_z_last = (mp_ID->GetNbVoxZ()-1-z)*mp_ID->GetNbVoxXY();
        // Loop over Y and X
        for (int i=0; i<mp_ID->GetNbVoxXY(); i++)
        {
          INTNB index = base_z_last + i;
          m2p_outputParImages[b][index] = 0.;
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
/*
  \fn      SaveParametricImages
  \param   a_iteration : current iteration index
  \param   a_subset : current number of subsets (or -1 by default)
  \brief   This function is pure virtual so must be implemented by children \n
           Call SaveParametricImages() function of the dynamic model object,
           in order to write on disk any parametric image of the model
  \return  0 if success, positive value otherwise
*/
int vDynamicModel::SaveParametricImages(int a_iteration, int a_subset)
{
  if(m_verbose >=2) Cout("vDynamicModel::SaveParametricImages ..." <<endl);
  
          
  if(m_saveParImageFlag)
  {
    // Get the output manager
    sOutputManager* p_output_manager = sOutputManager::GetInstance();
    
    // Interfile
    string path_to_image = p_output_manager->GetPathName() + p_output_manager->GetBaseName();
    
    // Add a suffix for iteration
    if (a_iteration >= 0)
    {
      stringstream ss; ss << a_iteration + 1;
      path_to_image.append("par_it").append(ss.str());
    }

    // Add a suffix for subset (if not negative by default), this means that we save a 'subset' image
    if (a_subset >= 0)
    {
      stringstream ss; ss << a_subset + 1;
      path_to_image.append("_ss").append(ss.str());
    }
  
    // Write interfile parametric image
    if(IntfWriteImgDynCoeffFile(path_to_image, 
                                m2p_outputParImages, 
                                mp_ID, 
                                m_nbModelParam, 
                                m_verbose) )
    {
      Cerr("***** iPatlakModel::SaveParametricImages()-> Error writing Interfile of output image !" << endl);  
      return 1;
    }
    
  }
  
  return 0;
}
