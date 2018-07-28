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
  \brief    Implementation of class vDeformation
*/

#include "vDeformation.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "oImageSpace.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn vDeformation
  \brief Constructor of vDeformation. Simply set all data members to default values.
*/
vDeformation::vDeformation() 
{
  mp_ID = NULL;
  m_nbTransformations = -1;
  m_verbose = -1;
  m_checked = false;
  m_initialized = false;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ~vDeformation
  \brief Destructor of vDeformation.
*/
vDeformation::~vDeformation() {}




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
int vDeformation::CheckParameters()
{
  if(m_verbose>=2) Cout("vDeformation::CheckParameters ..."<< endl); 
    
  // Check image dimensions
  if (mp_ID==NULL)
  {
    Cerr("***** vDeformation::CheckParameters() -> No image dimensions provided !" << endl);
    return 1;
  }
  
  // Check verbosity
  if (m_verbose<0)
  {
    Cerr("***** vDeformation::CheckParameters() -> Wrong verbosity level provided !" << endl);
    return 1;
  }

  // Check number of basis functions
  if (m_nbTransformations <0)
  {
    Cerr("***** vDeformation::CheckParameters() -> Number of transformations in the deformation has not been initialized !" << endl);
    return 1;
  }
  
  // Check parameters of the child class (if this function is overloaded)
  if (CheckSpecificParameters())
  {
    Cerr("***** vDeformation::CheckParameters() -> An error occurred while checking parameters of the child dynamic class !" << endl);
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
/*
  \fn ApplyDeformationsToBackwardImage
  \param ap_Image : required to access the backward image and its deformation backup matrice
  \param a_fr : frame index
  \param a_defIdx : index of the deformation
  \brief Apply backward transformation of the backward image to the reference position
  \details Loop on frames
           Recover any potential data stored in the backup matrice m2p_defTmpBackwardImage
  \return 0 if success, positive value otherwise
*/
int vDeformation::ApplyDeformationsToBackwardImage(oImageSpace* ap_Image, int a_fr, int a_defIdx)
{
  if(m_verbose >= 2) Cout("vDeformation::ApplyDeformationsToBackwardImage ... " <<endl);

  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** vDeformation::ApplyDeformationsToBackwardImage() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  
  for(int bimg=0; bimg<ap_Image->GetNbBackwardImages(); bimg++)
    //for(int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
      for(int rimg=0; rimg<mp_ID->GetNbRespGates(); rimg++)
        for(int cimg=0; cimg<mp_ID->GetNbCardGates(); cimg++)
        {
          /*
                    int r;
          scanf("%d", &r);
          
        for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v=v+mp_ID->GetNbVoxX())
      cout << ap_Image->m6p_backwardImage[bimg][0][fr][rimg][cimg][v] << endl;
  
  
          scanf("%d", &r);*/
          
          // Perform backward deformation
          if(a_defIdx > 0 // a_defIdx == 0 -> reference position, no need for transformation
            && ApplyDeformations(ap_Image->m6p_backwardImage[bimg][0][a_fr][rimg][cimg], ap_Image->m6p_backwardImage[bimg][0][a_fr][rimg][cimg], BACKWARD_DEFORMATION, a_defIdx-1) )
          {
            Cerr("***** vDeformation::ApplyDeformationsToBackwardImage()-> An error occured while performing backward deformation of the backward image !" << endl);
            Cerr("*****                                                      frame index " << a_fr << " respiratory image index " << rimg<< " cardiac image index " << cimg<<  " !" << endl);
            return 1;
          }
          /*
                  for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v=v+mp_ID->GetNbVoxX())
      cout << ap_Image->m6p_backwardImage[bimg][0][fr][rimg][cimg][v] << endl;
  
  
          scanf("%d", &r);*/
          
          
          // Recover the content of the temporary backup deformation image to the backward image
          for(int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
            ap_Image->m6p_backwardImage[bimg][0][a_fr][rimg][cimg][v] += ap_Image->m5p_refDynBackwardImage[bimg][a_fr][rimg][cimg][v] ;
        }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ApplyDeformationsToHistoSensitivityImage
  \param ap_Image : required to access the backward image and its deformation backup matrice
  \param a_fr : frame index
  \param a_defIdx : index of the deformation
  \brief Apply backward transformations of the sensitivity image to the reference position (histogram mode)
  \details Loop on frames
          Recover any potential data stored in the backup matrice m4p_defTmpSensitivityImage
  \return 0 if success, positive value otherwise
*/
int vDeformation::ApplyDeformationsToHistoSensitivityImage(oImageSpace* ap_Image, int a_fr, int a_defIdx)
{
  if(m_verbose >= 2) Cout("vDeformation::ApplyDeformationsToHistoSensitivityImage ... " <<endl);
  
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** vDeformation::ApplyDeformationsToHistoSensitivityImage() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  
  //for(int fr=0; fr<mp_ID->GetNbTimeFrames(); fr++)
    for(int rimg=0; rimg<mp_ID->GetNbRespGates() ; rimg++)
      for(int cimg=0; cimg<mp_ID->GetNbCardGates() ; cimg++)
      {
        // Perform backward deformation
        if( a_defIdx > 0 // a_defIdx == 0 -> reference position, no need for transformation
           && ApplyDeformations(ap_Image->m5p_sensitivity[0][a_fr][rimg][cimg], ap_Image->m5p_sensitivity[0][a_fr][rimg][cimg], BACKWARD_DEFORMATION, a_defIdx-1) )
           {
             Cerr("***** vDeformation::ApplyDeformationsToHistoSensitivityImage()-> An error occured while performing backward deformation of the backward image !" << endl);
             Cerr("*****                                                            frame index " << a_fr << " respiratory image index " << rimg<< " cardiac image index " << cimg<<  " !" << endl);
             return 1;
           }
        
        // Recover the content of the temporary backup deformation image in the sensitivity image
        for(int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
          ap_Image->m5p_sensitivity[0][a_fr][rimg][cimg][v] += ap_Image->m4p_refDynSensitivityImage[a_fr][rimg][cimg][v];
      }

  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PerformDeformation
  \param ap_Image : required to access oImageSpace image matrices
  \param a_defIdx : index of the deformation
  \param a_fr : frame index
  \param a_rimg : respiratory image index
  \param a_cimg : cardiac image index
  \brief Apply deformations during reconstruction
  \details 1. Recover all the data of the multithreaded backward image matrices in the first one (thread index 0)
           2. Perform backward deformation of the backward image to the reference position with defIdx-1
           3. Add coefficients of the backward image matrice to the temporary backup image matrice & reset backward image
           4. Apply forward deformation of the forward image (backward deformation to the reference position with defIdx-1, and forward deformation with defIdx)
  \return 0 if success, positive value otherwise
*/
int vDeformation::PerformDeformation(oImageSpace* ap_Image, int a_defIdx, int a_fr, int a_rimg, int a_cimg)
{
  if(m_verbose >= 3) Cout("vDeformation::PerformDeformation, deformation #" << a_defIdx+1 << endl);
  
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** vDeformation::PerformDeformation() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  
  // REDUCE
  for (int bimg=0; bimg<ap_Image->GetNbBackwardImages(); bimg++)
  {
    for (int th=1 ; th<mp_ID->GetNbThreadsForProjection() ; th++) //TODO add loops
    {
      // Synchronization of the multi-threaded backwardImages inside the first image Reduce)
      for(int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
        ap_Image->m6p_backwardImage[bimg][0][a_fr][a_rimg][a_cimg][v] += ap_Image->m6p_backwardImage[bimg][th][a_fr][a_rimg][a_cimg][v];
    }

    // BACKWARD DEFORMATION of the backward image  //TODO add loops
    if (a_defIdx > 1 // a_defIdx == 1 -> already in reference position, no need for transformation
       && ApplyDeformations(ap_Image->m6p_backwardImage[bimg][0][a_fr][a_rimg][a_cimg], ap_Image->m6p_backwardImage[bimg][0][a_fr][a_rimg][a_cimg], BACKWARD_DEFORMATION, a_defIdx-2) )
    {
      Cerr("***** vDeformation::ApplyDeformations()-> An error occured while performing backward deformation of the backward image !" << endl);
      Cerr("*****                                     frame index " << a_fr << " respiratory image index " << a_rimg<< " cardiac image index " << a_cimg<<  " !" << endl);
      return 1;
    }

    // Store Backward deformation update coefficients in the temporary backup image  //TODO add loops
    for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
    {
      ap_Image->m5p_refDynBackwardImage[bimg][a_fr][a_rimg][a_cimg][v] += ap_Image->m6p_backwardImage[bimg][0][a_fr][a_rimg][a_cimg][v];
    }
  }
  
  // Reset the backward images (i.e set all voxels to the specific fr/rimg/cimg to 0)
  for (int bimg=0; bimg<ap_Image->GetNbBackwardImages(); bimg++)
    for(int th=0 ; th<mp_ID->GetNbThreadsForProjection() ; th++)
      for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
        ap_Image->m6p_backwardImage[bimg][th][a_fr][a_rimg][a_cimg][v] = 0.;

  /*
            int r;
          scanf("%d", &r);
                cout << "ApplyDef forward Image AVT" << endl;
        for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v=v+mp_ID->GetNbVoxX())
      if(ap_Image->m4p_forwardImage[fr][rimg][cimg][v] != 1)cout << ap_Image->m4p_forwardImage[fr][rimg][cimg][v] << endl;
  
  
          scanf("%d", &r);*/
  
  
  // FORWARD DEFORMATION of the Forward image
  /*
  if (ApplyDeformations(ap_Image->m4p_forwardImage[fr][rimg][cimg], ap_Image->m4p_forwardImage[fr][rimg][cimg], FORWARD_DEFORMATION, a_defIdx) )
  {
    Cerr("***** vDeformation::ApplyDeformations()-> An error occured while performing forward deformation of the forward image !" << endl);
    Cerr("*****                                     frame index " << fr << " respiratory image index " << rimg<< " cardiac image index " << cimg<<  " !" << endl);
    return 1;
  }
  */
  
  if (ApplyDeformations(ap_Image->m4p_refDynForwardImage[a_fr][a_rimg][a_cimg], ap_Image->m4p_forwardImage[a_fr][a_rimg][a_cimg], FORWARD_DEFORMATION, a_defIdx-1) )
  {
    Cerr("***** vDeformation::ApplyDeformations()-> An error occured while performing forward deformation of the forward image !" << endl);
    Cerr("*****                                     frame index " << a_fr << " respiratory image index " << a_rimg<< " cardiac image index " << a_cimg<<  " !" << endl);
    return 1;
  }
  /*
  cout << "ApplyDef forward Image APS" << endl;
            scanf("%d", &r);
          for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v=v+mp_ID->GetNbVoxX())
      if(ap_Image->m4p_forwardImage[fr][rimg][cimg][v] != 1)cout << ap_Image->m4p_forwardImage[fr][rimg][cimg][v] << endl;
  
  
          scanf("%d", &r);*/
          
  return 0;
}







int vDeformation::PerformDeformationBis(oImageSpace* ap_Image, int a_defIdx, int a_fr, int a_rimg, int a_cimg)
{
  if(m_verbose >= 3) Cout("vDeformation::PerformDeformationBis, deformation #" << a_defIdx+1 << endl);
  
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** vDeformation::PerformDeformation() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  

  
  if (ApplyDeformations(ap_Image->m4p_refDynForwardImage[a_fr][a_rimg][a_cimg], ap_Image->m4p_forwardImage[a_fr][a_rimg][a_cimg], FORWARD_DEFORMATION, a_defIdx-1) )
  {
    Cerr("***** vDeformation::ApplyDeformationsBis()-> An error occured while performing forward deformation of the forward image !" << endl);
    Cerr("*****                                     frame index " << a_fr << " respiratory image index " << a_rimg<< " cardiac image index " << a_cimg<<  " !" << endl);
    return 1;
  }

  return 0;
}











// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PerformHistoSensitivityDeformation
  \param ap_Image : required to access oImageSpace image matrices
  \param a_defIdx : index of the deformation
  \param fr : frame index
  \param rimg : respiratory image index
  \param cimg : cardiac image index
  \brief Apply deformations on the sensitivity image during reconstruction in histogram mode
  \details 1. Recover all the data of the multithreaded sensitivity image matrice in the first one (thread index 0)
           2. Perform backward deformation of the sensitivity image to the reference position with defIdx-1
           3. Add coefficients of the sensitivity image matrice to the temporary backup image matrice & reset sensitivity image
  \return 0 if success, positive value otherwise
*/
int vDeformation::PerformHistoSensitivityDeformation(oImageSpace* ap_Image, int a_defIdx, int fr, int rimg, int cimg)
{
  if(m_verbose >= 2) Cout("vDeformation::PerformHistoSensitivityDeformation ... " <<endl);
  
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** vDeformation::PerformHistoSensitivityDeformation() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  
  // REDUCE
  for (int th=1 ; th<mp_ID->GetNbThreadsForProjection() ; th++)
  {
    // Synchronisation of the multi-threaded sensitivity images inside the first image )
    for(int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
      ap_Image->m5p_sensitivity[0][fr][rimg][cimg][v] += ap_Image->m5p_sensitivity[th][fr][rimg][cimg][v];
  }
  
  // BACKWARD DEFORMATION of the sensitivity image
  if (a_defIdx > 1 // a_defIdx == 1 -> already in reference position, no need for transformation
     && ApplyDeformations(ap_Image->m5p_sensitivity[0][fr][rimg][cimg], ap_Image->m5p_sensitivity[0][fr][rimg][cimg], BACKWARD_DEFORMATION, a_defIdx-2) )
  {
    Cerr("***** vDeformation::PerformHistoSensitivityDeformation()-> An error occured while performing backward deformation of the sensitivity image !" << endl);
    Cerr("*****                                                      frame index " << fr << " respiratory image index " << rimg<< " cardiac image index " << cimg<<  " !" << endl);
    return 1;
  }
  
  // Store Backward deformation update coefficients in temporary image
  for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
  {
    ap_Image->m4p_refDynSensitivityImage[fr][rimg][cimg][v] += ap_Image->m5p_sensitivity[0][fr][rimg][cimg][v];
  }

  // Reset the sensitivity images (i.e set all voxels to the specific fr/rimg/cimg to 0) 
  for (int th=0; th<mp_ID->GetNbThreadsForProjection(); th++)
    for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
      ap_Image->m5p_sensitivity[th][fr][rimg][cimg][v] = 0.;

  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn PerformSensitivityDeformation
  \param ap_Image : required to access oImageSpace image matrices
  \param a_defDirection : a direction for the deformation to perform (forward or backward)
  \param a_defIdx : index of the deformation
  \param fr : frame index
  \param rg : respiratory gate index
  \param cg : cardiac gate index
  \brief Apply image deformations during sensitivity image generation for list-mode
  \details Depending on the deformation direction (forward or backward):
           Forward : Perform forward deformation of the forward image to the deformation index position
           Backward: Perform backward deformation of the backward image to the reference position
  \return 0 if success, positive value otherwise
*/
int vDeformation::PerformSensitivityDeformation(oImageSpace* ap_Image, int a_defDirection, int a_defIdx, int fr, int rg, int cg)
{
  if(m_verbose >= 2) Cout("vDeformation::PerformSensitivityDeformation ... " << a_defDirection << ";" <<  a_defIdx <<  endl);
  
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** vDeformation::PerformSensitivityDeformation() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  
  // Idx 0 == reference (no deformation)
  if(a_defIdx == 0)
    return 0;
    
  // tranformations are numbered from 0 to nbTransformations-1
  a_defIdx--;
  
  if (a_defDirection == FORWARD_DEFORMATION)
  {
    if (ApplyDeformations(ap_Image->m4p_forwardImage[fr][rg][cg], ap_Image->m4p_forwardImage[fr][rg][cg], a_defDirection, a_defIdx) )
    //if (ApplyDeformations(ap_Image->m4p_forwardImage[fr][0][cg], ap_Image->m4p_forwardImage[fr][a_defIdx][cg], a_defDirection, a_defIdx) )
    {
      Cerr("***** vDeformation::PerformSensitivityDeformation -> An error occured while performing forward deformation !" << endl);
      Cerr("*****                                                frame index " << fr << " respiratory gate index " << rg<< " cardiac gate index " << cg<<  " !" << endl);
      return 1;
    }
  }
  else if (a_defDirection == BACKWARD_DEFORMATION)
  {
    if (ApplyDeformations(ap_Image->m6p_backwardImage[0][0][fr][rg][cg], ap_Image->m6p_backwardImage[0][0][fr][rg][cg], a_defDirection, a_defIdx) )
    //if (ApplyDeformations(ap_Image->m6p_backwardImage[0][0][fr][a_defIdx][cg], ap_Image->m6p_backwardImage[0][0][fr][a_defIdx][cg], a_defDirection, a_defIdx) )
    {
      Cerr("***** vDeformation::PerformSensitivityDeformation -> An error occured while performing backward deformation !" << endl);
      Cerr("*****                                                frame index " << fr << " respiratory gate index " << rg<< " cardiac gate index " << cg<<  " !" << endl);
      return 1;
    }
  }
  else 
  {
    Cerr("***** vDeformation::PerformDeformation -> Unknown type of deformation !" << endl);
    return 1;
  }

  return 0;
}




/*
  \fn Tlerp
  \param ap_inputImage : input image matrix 
  \param ap_outputImage : output image matrix
  \param iov : index of the voxel to interpolate in the output image
  \param iiv : index of the input image central voxel for interpolation
  \param dx : x-axis difference between output voxel cartesian position after transformation and center of estimated vox position
  \param dy : y-axis difference between output voxel cartesian position after transformation and center of estimated vox position
  \param dz : z-axis difference between output voxel cartesian position after transformation and center of estimated vox position
  \brief This function performs a trilinear interpolation for a specific voxel
  \todo : perhaps use padded image in order to avoid if statements
  \return 0 if success, other value otherwise.
*/
int vDeformation::Tlerp(HPFLTNB *ap_inputImage, HPFLTNB *ap_outputImage, uint32_t iov, uint32_t iiv, FLTNB dX, FLTNB dY, FLTNB dZ)
{
  // Trilinear interpolation
  // Todo : Use padded image in order to avoid 'if' statements ?
  int incX = dX>0 ? 1 : -1;
  int incY = dY>0 ? mp_ID->GetNbVoxX() : -mp_ID->GetNbVoxX();
  int incZ = dZ>0 ? mp_ID->GetNbVoxXY() : -mp_ID->GetNbVoxXY();

  dX = (dX>=0) ? dX : -dX;
  dY = (dY>=0) ? dY : -dY;
  dZ = (dZ>=0) ? dZ : -dZ;
  
  uint32_t nb_tot_vox = (uint32_t)mp_ID->GetNbVoxXYZ();

  if ((iiv <  nb_tot_vox) 
   && (iiv >= 0)) 
    ap_outputImage[iov] += ap_inputImage[iiv] * (1.-dX)*(1.-dY)*(1.-dZ);

  if ((iiv+incX <  nb_tot_vox) 
   && (((int)iiv)+incX >= 0))
      ap_outputImage[iov] += ap_inputImage[iiv+incX] * dX*(1.-dY)*(1.-dZ);

  if ((iiv+incY <  nb_tot_vox )
  &&  (((int)iiv)+incY >= 0))
      ap_outputImage[iov] += ap_inputImage[iiv +incY] * (1.-dX)*dY*(1.-dZ);
                               
  if ((iiv+incX+incY <  nb_tot_vox) 
  && ( ((int)iiv)+incX+incY >= 0))
      ap_outputImage[iov] += ap_inputImage[iiv+incX+incY] * dX*dY*(1.-dZ);

  if ((iiv+incZ  <  nb_tot_vox) 
   && (((int)iiv)+incZ) >= 0)
      ap_outputImage[iov] += ap_inputImage[iiv+incZ] * (1.-dX)*(1.-dY)*dZ;

  if ((iiv+incX+incZ) <  nb_tot_vox  
  &&  (((int)iiv)+incX+incZ) >=0 )
      ap_outputImage[iov] += ap_inputImage[iiv+incX+incZ] * dX*(1.-dY)*dZ;

  if ((iiv+incY+incZ <  nb_tot_vox) 
   && (((int)iiv)+incY+incZ >= 0))
      ap_outputImage[iov] += ap_inputImage[iiv+incY+incZ] * (1.-dX)*dY*dZ;

  if ((iiv+incX+incY+incZ <  nb_tot_vox) 
  &&  (((int)iiv)+incX+incY+incZ >= 0))
      ap_outputImage[iov] += ap_inputImage[iiv+incX+incY+incZ] *dX*dY*dZ;

    
  return 0;
}
