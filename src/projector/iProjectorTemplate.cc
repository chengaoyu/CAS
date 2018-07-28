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
  \ingroup  projector
  \brief    Implementation of class iProjectorTemplate
*/

#include "iProjectorTemplate.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iProjectorTemplate::iProjectorTemplate() : vProjector()
{
  // Set all the data members to a default value

  // Also tell if the projector is compatible with SPECT attenuation correction. In order
  // to be so, all voxels contributing to a line must be strictly sorted with respect to
  // their distance to point 2 (the line must also go from point 1 to point 2 and not the
  // inverse)
  m_compatibleWithSPECTAttenuationCorrection = false;

  // Tell if the projector is compatible with compression. This means that for a given event,
  // multiple couples of indices are used to compute a centroid position for each line's end
  // point. Then the detection element indices in the projection line are set to -1 because
  // it is undefined. So inside the projection algorithm, if such indices must be used, then
  // the projector will not be compatible with compression.
  m_compatibleWithCompression = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iProjectorTemplate::~iProjectorTemplate()
{
  // Delete or free all structures allocated by this projector
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorTemplate::ReadConfigurationFile(const string& a_configurationFile)
{
  // Implement here the reading of any options specific to this projector, through a configuration file
  ;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorTemplate::ReadOptionsList(const string& a_optionsList)
{
  // Implement here the reading of any options specific to this projector, through a list of options separated by commas
  ;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iProjectorTemplate::ShowHelpSpecific()
{
  // Here, display some help and guidance to how to use this projector and what it does
  cout << "This projector is a template class dedicated to add your own custom projector." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorTemplate::CheckSpecificParameters()
{
  // Here, check that all parameters needed by this projector are allocated and have correct values
  ;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorTemplate::InitializeSpecific()
{
  // Implement here the initialization of whatever member variables specifically used by this projector
  ;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

INTNB iProjectorTemplate::EstimateMaxNumberOfVoxelsPerLine()
{
  // Implement here a way to precompute the estimated maximum number of voxels that will contribute to a line of response.
  // By default, it uses a buffer size corresponding to the total number of voxels of the image.
  // The idea is to optimize the RAM usage by providing a better estimate that suites the need of this projector.
  // If you do not have a better estimation, then you can remove this function from this class because it is already
  // implemented as is in the mother class.

  // Find the maximum number of voxels along a given dimension
  INTNB max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxXYZ();
  // Return the value
  return max_nb_voxels_in_dimension;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorTemplate::ProjectWithoutTOF(int a_direction, oProjectionLine* ap_ProjectionLine )
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iProjectorTemplate::ProjectWithoutTOF() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("iProjectorTemplate::Project without TOF -> Project line '" << ap_ProjectionLine << "' in " << direction << " direction" << endl);
  }
  #endif

  // --------------------------------------------------------------------------------------------------------------------------------------------
  // Please read the following information that will help implement your projector:

  // FLTNB is a macro defining the precision of the code (float, double, long double) that can be customized through some compilation options.
  // So please, DO NOT USE 'float' or 'double' keywords but USE INSTEAD 'FLTNB'.
  // Same for integers used to define image dimensions, DO NOT USE 'int' or 'long int' but USE INSTEAD 'INTNB'.

  // All 3D vectors of type FLTNB* or INTNB* carry the information in the following order: X then Y then Z.

  // The image dimensions can be accessed via some local copies of the parameters:
  //  - number of voxels: mp_nbVox[0] (along X), mp_nbVox[1] (along Y), mp_nbVox[2] (along Z),
  //  - size of voxels in mm: mp_voxSize[0] (along X), mp_voxSize[1] (along Y), mp_voxSize[2] (along Z),
  //  - half image dimensions in mm: mp_halfFOV[0] (along X), mp_halfFOV[1] (along Y), mp_halfFOV (along Z).

  // For code efficiency and readability, the spatial index of a voxel is a cumulative 1D index. That is to say, given a voxel [indexX,indexY,indexZ],
  // its cumulative 1D index is computed by 'index = indexZ*m_nbVoxXY + indexY*mp_nbVox[0] + indexX'.

  // All information that you may need about the line of response are embedded into the oProjectionLine object given as a parameter. So take a look
  // at this class to know how to get those information through some ap_ProjectionLine->GetXXX() functions.

  // The end points of the line are already computed by the vProjector with respect to the different options provided; e.g. mean depth of
  // interaction, actual position of interaction (POI), randomization of end points, etc. However, if you want to customize those end points,
  // take a look at what the vScanner and children can do through the use of some dedicated functions. If it cannot do what you want, consider adding
  // this function into the vScanner or children classes. The vScanner object can be accessed using the mp_Scanner member object of this class.

  // The projected line must go from point 1 to point 2 and voxel contributions by sorted in order to be compatible with SPECT attenuation correction.

  // Finally, to add the contribution of a given voxel to this projection line, simply use this instruction:
  // ap_ProjectionLine->AddVoxel(a_direction, my_index, my_weight), where 'my_index' is the spatial index of the voxel and 'my_weight' is its
  // associated weight (i.e. its contribution to the line).

  // Finally, remember that the mantra of CASToR is the genericity, so when you add some code, think about it twice in order to ensure that this
  // piece of code can be used by anyone in any context!

  // --------------------------------------------------------------------------------------------------------------------------------------------

  Cerr("***** iProjectorTemplate::ProjectWithoutTOF() -> Not yet implemented !" << endl);
  return 1;

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorTemplate::ProjectWithTOFPos(int a_Projector, oProjectionLine* ap_ProjectionLine)
{
  // Read the information in the ProjectWithoutTOF function to know the general guidelines.
  // This function implements a projection using a continuous TOF information = when using list-mode data.
  // The TOF resolution and measurement associated to the running event are accessible through the ap_ProjectionLine
  // parameter using some GetXXX() functions.

  Cerr("***** iProjectorTemplate::ProjectWithTOFPos() -> Not yet implemented !" << endl);
  return 1;

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorTemplate::ProjectWithTOFBin(int a_Projector, oProjectionLine* ap_ProjectionLine)
{
  // Read the information in the ProjectWithoutTOF function to know the general guidelines.
  // This function implements a projection using a binned TOF information = when using histogram data.
  // The number of TOF bins, TOF resolution, etc, are accessible through the ap_ProjectionLine using
  // some GetXXX() functions. This function is supposed to fill all TOF bins at once. To add voxel
  // contributions to a specific TOF bin, use the dedicated function ap_ProjectionLine->AddVoxelInTOFBin().
  // All forward and backward operations will be carried out later by the vOptimizer, automatically
  // managing all TOF bins.

  Cerr("***** iProjectorTemplate::ProjectWithTOFBin() -> Not yet implemented !" << endl);
  return 1;

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
