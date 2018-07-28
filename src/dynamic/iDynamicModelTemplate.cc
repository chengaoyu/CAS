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
  \brief    Implementation of class iDynamicModelTemplate
*/

#include "iDynamicModelTemplate.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn iDynamicModelTemplate
  \brief Constructor of iDynamicModelTemplate. Simply set all data members to default values.
*/
iDynamicModelTemplate::iDynamicModelTemplate() : vDynamicModel() 
{
  m_nbTimeBF = 1; // Initialize the number of basis functions in the model
  m2p_parametricImages = NULL; // Should be allocated in Initialize() function
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ~iDynamicModelTemplate
  \brief Destructor of iDynamicModelTemplate
*/
iDynamicModelTemplate::~iDynamicModelTemplate() 
{
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ShowHelp
  \brief Print out specific help about the implementation of the model 
         and its initialization
*/
void iDynamicModelTemplate::ShowHelp()
{
  // ===================================================================
  // Here, display some help and guidance to how to use this dynamic model and what it does
  // ===================================================================
  cout << "This class is a template class dedicated to add your own dynamic model." << endl;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckConfigurationFile
  \param const string& a_configurationFile : ASCII file containing informations about a dynamic model
  \brief This function is used to read options from a configuration file.
  \return 0 if success, other value otherwise.
*/
int iDynamicModelTemplate::ReadAndCheckConfigurationFile(string a_fileOptions)
{
  if(m_verbose >=2) Cout("iDynamicModelTemplate::ReadAndCheckConfigurationFile ..."<< endl); 
  
  // ===================================================================
  // Implement here the reading of any options specific to this dynamic model 
  // (i.e : parameters or path to deformation files), through a configuration file
  // The ReadDataASCIIFile() functions could be helpful to recover data from a file
  // ===================================================================
    
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckOptionsList
  \param a_optionsList : a list of parameters separated by commas
  \brief This function is used to read parameters from a string.
  \return 0 if success, other value otherwise.
*/
int iDynamicModelTemplate::ReadAndCheckOptionsList(string a_listOptions)
{
  // ===================================================================
  // Implement here the reading of any options specific to this deformation model,
  // through a list of options separated by commas
  // The ReadStringOption() function could be helpful to parse the list of parameters in an array
  // ===================================================================
  
  // Normal end
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn CheckSpecificParameters
  \brief This function is used to check whether all member variables
         have been correctly initialized or not.
  \return 0 if success, positive value otherwise.
*/
int iDynamicModelTemplate::CheckSpecificParameters()
{
  // ===================================================================
  // Implement here checks over parameters which should be read using either
  // ReadAndCheckConfigurationFile() and ReadAndCheckOptionsList() functions
  // ===================================================================
    
  if(m_verbose >=2) Cout("iDynamicModelTemplate::CheckSpecificParameters ..."<< endl); 
  
  // Normal end
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn Initialize
  \brief This function is used to initialize the model parametric images and basis functions
  \return 0 if success, other value otherwise.
*/
int iDynamicModelTemplate::Initialize()
{
  if(m_verbose >=2) Cout("iDynamicModelTemplate::Initialize ..."<< endl); 


  // ===================================================================
  // Implement here the allocation/initialization of whatever member
  // variables specifically used by this deformation model
  // ===================================================================
  
  
  // Forbid initialization without check
  if (!m_checked)
  {
    Cerr("***** oDynamicModelManager::Initialize() -> Must call CheckParameters functions before Initialize() !" << endl);
    return 1;
  }
  
  // Normal end
  m_initialized = true;
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn EstimateModelParameters
  \param ap_ImageS : pointer to the ImageSpace
  \param a_ite : index of the actual iteration (not used)
  \param a_sset : index of the actual subset (not used)
  \brief Estimate parametric images
  \return 0 if success, other value otherwise.
*/
int iDynamicModelTemplate::EstimateModelParameters(oImageSpace* ap_ImageS, int a_ite, int a_sset) 
{
  if(m_verbose >=2) Cout("iDynamicModelTemplate::EstimateModelParameters ..." <<endl);

  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iDynamicModelTemplate::EstimateModelParameters() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  
  // ===================================================================
  // Implement here the parametric images estimate after the image update
  // which occurs at the end of the iteration (or subset)
  //
  //
  // The main image matrices are stored in the ImageSpace object, passed
  // in argument.
  // The main image contains 4 dimensions :
  // ap_ImageS->m4p_image[fr][rg][cg][v]
  // fr = time frames
  // rg = respiratory gates
  // cg = cardiac gates
  //  v = actual voxel of the 3D volume
  
  /* IMAGE DIMENSIONS :
  * For code efficiency and readability, the spatial index of a voxel is a cumulative 1D index. That is to say, given a voxel [indexX,indexY,indexZ],
  * its cumulative 1D index is computed by 'index = indexZ*nbVoxXY + indexY*nbVoxX + indexX'.
  *
  * The image dimensions can be recovered from the mp_ID class
  * Total number of voxels         : mp_ID->GetNbVoxXYZ()
  * Number of voxels in a slice    : mp_ID->GetNbVoxXY()
  * Number of voxels on the X-axis : mp_ID->GetNbVoxX()
  */
  
  // Any error should return a value >0.
    
  // ===================================================================
  
  return 0;
}
  
  
  
// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn EstimateImageWithModel
  \param ap_ImageS : pointer to the ImageSpace
  \param a_ite : index of the actual iteration (not used)
  \param a_sset : index of the actual subset (not used)
  \brief Estimate image using the model parametric images and basis functions
  \return 0 if success, other value otherwise.
*/
int iDynamicModelTemplate::EstimateImageWithModel(oImageSpace* ap_ImageS, int a_ite, int a_sset) 
{
  if(m_verbose >= 2) Cout("iDynamicModelTemplate::EstimateImageWithModel ... " <<endl);
  
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iDynamicModelTemplate::EstimateImageWithModel() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif


  // ===================================================================
  // This function can be used to generate the serie of dynamic images 
  // from the model, after estimation of the model parameters in EstimateModelParameters()
  // It is called right after EstimateModelParameters()
  //
  // The main image matrices are stored in the ImageSpace object, passed
  // in argument.
  // The main image contains 4 dimensions :
  // ap_ImageS->m4p_image[fr][rg][cg][v]
  // fr = time frames
  // rg = respiratory gates
  // cg = cardiac gates
  //  v = actual voxel of the 3D volume
  
  /* IMAGE DIMENSIONS :
  * For code efficiency and readability, the spatial index of a voxel is a cumulative 1D index. That is to say, given a voxel [indexX,indexY,indexZ],
  * its cumulative 1D index is computed by 'index = indexZ*nbVoxXY + indexY*nbVoxX + indexX'.
  *
  * The image dimensions can be recovered from the mp_ID class
  * Total number of voxels         : mp_ID->GetNbVoxXYZ()
  * Number of voxels in a slice    : mp_ID->GetNbVoxXY()
  * Number of voxels on the X-axis : mp_ID->GetNbVoxX()
  */
  
  // Any error should return a value >0.
    
  // ===================================================================
        
  return 0;
}





// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn SaveParametricImages
  \param ap_ImageS : pointer to the ImageSpace
  \param a_ite : index of the actual iteration
  \brief Write parametric images on disk if 'm_savePImgFlag' is enabled
  \return 0 if success, other value otherwise.
*/
int iDynamicModelTemplate::SaveParametricImages(int a_ite, int a_subset)
{
  if(m_verbose >=2) Cout("iDynamicModelTemplate::SaveParametricImages ..." <<endl);
  

  // ===================================================================
  // Implement here the parametric image output writting step, if needed
  // ===================================================================

  // Get the output manager which contains output directory
  sOutputManager* p_output_manager = sOutputManager::GetInstance();  

  // Interfile
  string path_to_image = p_output_manager->GetPathName() + p_output_manager->GetBaseName();
  
  // Add a suffix for iteration
  if (a_ite >= 0)
  {
    stringstream ss; ss << a_ite + 1;
    path_to_image.append("dynModel_it").append(ss.str());
  }
  

  
  // This function is dedicated to the output of parametric images
  if(IntfWriteImgDynCoeffFile(path_to_image, // path to the image (or metaheader) to write 
                              m2p_parametricImages, // parametric image to write.
                                                    // it must be a 2 dimensions array [nb of basis functions][nb of voxels]
                              mp_ID, // just the object containing image dimensions
                              m_nbTimeBF, // number of basis functions of the model
                              m_verbose) )
  {
    Cerr("***** iDynamicModelTemplate::SaveParametricImages()-> Error writing Interfile of output image !" << endl);  
    return 1;
  }
  
  return 0;
}


