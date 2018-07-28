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
  \brief    Implementation of class iDeformationTemplate
*/

#include "iDeformationTemplate.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn iDeformationTemplate
  \brief Constructor of iDeformationTemplate. Simply set all data members to default values.
*/
iDeformationTemplate::iDeformationTemplate() : vDeformation()
{
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ~iDeformationTemplate
  \brief Destructor of iDeformationTemplate. Free memory from all allocated tabs.
*/
iDeformationTemplate::~iDeformationTemplate() 
{
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn ShowHelp
  \brief This function is used to print out specific help about the deformation model and its options.
*/
void iDeformationTemplate::ShowHelp()
{
  // ===================================================================
  // Here, display some help and guidance to how to use this deformation and what it does
  // ===================================================================
  
  cout << "This class is a template class dedicated to add your own custom deformation model." << endl;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckConfigurationFile
  \param a_configurationFile
  \brief This function is an implementation of the pure virtual mother function. It is used to read options from a configuration file.
  \return 0 if success, other value otherwise.
*/
int iDeformationTemplate::ReadAndCheckConfigurationFile(const string& a_fileOptions)
{
  if(m_verbose >=2) Cout("iDeformationTemplate::ReadAndCheckConfigurationFile ..."<< endl); 
  
  // ===================================================================
  // Implement here the reading of any options specific to this deformation model 
  // (i.e : parameters or path to deformation files), through a configuration file
  // The ReadDataASCIIFile() functions could be helpful to recover data from a file
  // ===================================================================
  
  // Normal end
  return 0;  
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckOptionsList
  \param a_optionsList
  \brief This function is an implementation of the pure virtual mother function. It is used to read options from a list of options.
  \return 0 if success, other value otherwise.
*/
int iDeformationTemplate::ReadAndCheckOptionsList(const string& a_listOptions)
{
  if(m_verbose >=2) Cout("iDeformationTemplate::ReadAndCheckOptionsList ..."<< endl); 
  
  // ===================================================================
  // Implement here the reading of any options specific to this deformation model, through a list of options separated by commas
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
  \brief This function is an implementation of the pure virtual mother function. It is used to
         check parameters of the child deformation model before initialization.
  \return 0 if success, other value otherwise.
*/
int iDeformationTemplate::CheckSpecificParameters()
{
  if(m_verbose >=2) Cout("iDeformationTemplate::CheckSpecificParameters ..."<< endl); 

  // ===================================================================
  // Implement here checks over parameters which should be read using either
  // ReadAndCheckConfigurationFile() and ReadAndCheckOptionsList() functions
  // ===================================================================
  
  // Normal end
  m_checked = true;
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn Initialize
  \brief This function is an implementation of the pure virtual mother function. It is used to
         initialize specific stuff to the child deformation model.
  \return 0 if success, other value otherwise.
*/
int iDeformationTemplate::Initialize()
{
  if(m_verbose >=2) Cout("iDeformationTemplate::Initialize ..."<< endl); 

  // ===================================================================
  // Implement here the allocation/initialization of whatever member variables specifically used by this deformation model
  // ===================================================================
  
  if (!m_checked)
  {
    Cerr("***** iDeformationTemplate::Initialize() -> Parameters should be checked before Initialize() !" << endl);
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
  \fn ApplyDeformations
  \param ap_inputImage : input image to deform
  \param ap_outputImage : image in which the output of the deformation should be recovered
  \param a_direction : a direction for the deformation to perform (forward or backward)
  \param a_defIdx : index of the deformation
  \brief This function is an implementation of the pure virtual mother function. The actual deformation should be implemented here 
  \return 0 if success, other value otherwise.
*/
int iDeformationTemplate::ApplyDeformations(FLTNB* ap_inputImage, FLTNB* ap_outputImage, int a_direction, int a_defIdx)
{
  #ifdef CASTOR_VERBOSE
  if(m_verbose >=4) Cout("iDeformationTemplate::ApplyDeformations ..."<< endl); 
  #endif

  // ===================================================================
  // The deformation model should be implemented here, with the help of any private functions if required
  
  /* The 'a_defIdx' parameter defines the deformation index of the transformation
  * The 'a_direction' parameter is an integer which indicates the direction of the deformation to perform, i.e :
  * - FORWARD_DEFORMATION (from the reference position to the 'a_defIdx' position) 
  * - BACKWARD_DEFORMATION (from the 'a_defIdx' position to the reference position) 
  * The integers FORWARD_DEFORMATION & BACKWARD_DEFORMATION are macros defined in the beginning of vDeformation.hh
  */
  
  // The deformation should be applied to the ap_inputImage matrix, and the resulting image should be recovered in ap_outputImage matrix

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
  
  // Normal end
  return 0;
}
