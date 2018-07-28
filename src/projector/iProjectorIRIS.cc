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
  \brief    Implementation of class iProjectorIRIS
*/

#include "iProjectorIRIS.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iProjectorIRIS::iProjectorIRIS() : vProjector()
{
  mp_pathToIDRFFiles = NULL;
  mp_IDRF = NULL;
  m2p_IDRF_CDFs = NULL;
  m_nbLinesPerLOR = -1;
  m_nVoxDepthIDRF = -1;
  m_nVoxTransaxialIDRF = -1;
  m_nVoxAxialIDRF = -1;
  m_nVoxXYZIDRF = -1;
  m_nBetaAnglesIDRF = -1;
  m_nAlphaAnglesIDRF = -1;
  m_sizeVoxDepthIDRF = -1.;
  m_sizeVoxTransaxialIDRF = -1.;
  m_sizeVoxAxialIDRF = -1.;
  m_stepBetaAnglesIDRF = -1.;
  m_stepAlphaAnglesIDRF = -1.;
  // This projector is not compatible with SPECT attenuation correction because
  // the voxels contributing to the line are not strictly ordered with respect to
  // their distance to point 2 (due to the use of multiple lines that are
  // stack one after the other)
  m_compatibleWithSPECTAttenuationCorrection = false;
  // This projector is not compatible with compression as it works only with the
  // detection element indices
  m_compatibleWithCompression = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iProjectorIRIS::~iProjectorIRIS()
{
  if (mp_pathToIDRFFiles) delete[] mp_pathToIDRFFiles;

  if (m_nAlphaAnglesIDRF>0 && m_nBetaAnglesIDRF>0 && m2p_IDRF_CDFs)
    for (int a=0 ; a<m_nAlphaAnglesIDRF*m_nBetaAnglesIDRF ; a++)
       if (m2p_IDRF_CDFs[a]) delete[] m2p_IDRF_CDFs[a];

  if (m2p_IDRF_CDFs) delete[] m2p_IDRF_CDFs;

  if (mp_IDRF) delete[] mp_IDRF;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIRIS::ReadConfigurationFile(const string& a_configurationFile)
{  
  if (ReadDataASCIIFile(a_configurationFile, "number lines per LOR", &m_nbLinesPerLOR, 1, 1) )
  {
    Cerr("***** iProjectorIRIS::ReadConfigurationFile() -> An error occurred while trying to read the number of lines per LOR from the IRIS configuration file at :" << a_configurationFile << " !" << endl);
    return 1;
  }

  if (ReadDataASCIIFile(a_configurationFile, "number voxels depth", &m_nVoxDepthIDRF, 1, 1) ||
      ReadDataASCIIFile(a_configurationFile, "number voxels transaxial", &m_nVoxTransaxialIDRF, 1, 1) ||
      ReadDataASCIIFile(a_configurationFile, "number voxels axial", &m_nVoxAxialIDRF, 1, 1) )
  {
    Cerr("***** iProjectorIRIS::ReadConfigurationFile() -> An error occurred while trying to read the IDRF volume dimensions (nb voxels) from the IRIS configuration file at :" << a_configurationFile << " !" << endl);
    return 1;
  }
 
  // Total number of voxels in the IDRF
  m_nVoxXYZIDRF = m_nVoxDepthIDRF*m_nVoxTransaxialIDRF*m_nVoxAxialIDRF;
  
  if (ReadDataASCIIFile(a_configurationFile, "size voxels depth", &m_sizeVoxDepthIDRF, 1, 1) ||
      ReadDataASCIIFile(a_configurationFile, "size voxels transaxial", &m_sizeVoxTransaxialIDRF, 1, 1) ||
      ReadDataASCIIFile(a_configurationFile, "size voxels axial", &m_sizeVoxAxialIDRF, 1, 1) )
  {
    Cerr("***** iProjectorIRIS::ReadConfigurationFile() -> An error occurred while trying to read the voxel sizes of the IDRF volumes from the IRIS configuration file at :" << a_configurationFile << " !" << endl);
    return 1;
  }

  if (ReadDataASCIIFile(a_configurationFile, "number beta angles", &m_nBetaAnglesIDRF, 1, 1) ||
      ReadDataASCIIFile(a_configurationFile, "number alpha angles", &m_nAlphaAnglesIDRF, 1, 1) )
  {
    Cerr("***** iProjectorIRIS::ReadConfigurationFile() -> An error occurred while trying to read the number of alpha/beta angles from the IRIS configuration file at :" << a_configurationFile << " !" << endl);
    return 1;
  }

  if (ReadDataASCIIFile(a_configurationFile, "step beta angles", &m_stepBetaAnglesIDRF, 1, 1) ||
      ReadDataASCIIFile(a_configurationFile, "step alpha angles", &m_stepAlphaAnglesIDRF, 1, 1) )
  {
    Cerr("***** iProjectorIRIS::ReadConfigurationFile() -> An error occurred while trying to read the alpha/beta angles steps from the IRIS configuration file at :" << a_configurationFile << " !" << endl);
    return 1;
  }
  
  // Convert angle steps in radians 
  m_stepAlphaAnglesIDRF *= M_PI/180.;
  m_stepBetaAnglesIDRF *= M_PI/180.;
  
  mp_pathToIDRFFiles = new string[m_nAlphaAnglesIDRF*m_nBetaAnglesIDRF];
  
  mp_IDRF = new float[m_nVoxXYZIDRF];
  m2p_IDRF_CDFs = new float*[m_nAlphaAnglesIDRF*m_nBetaAnglesIDRF];
  
  for(int a=0 ; a<m_nAlphaAnglesIDRF*m_nBetaAnglesIDRF ; a++)
    m2p_IDRF_CDFs[a] = new float[m_nVoxXYZIDRF];

  if (ReadDataASCIIFile(a_configurationFile, "IDRFs path", mp_pathToIDRFFiles, 1, m_nAlphaAnglesIDRF*m_nBetaAnglesIDRF, 1) )
  {
    Cerr("***** iProjectorIRIS::ReadConfigurationFile() -> An error occurred while trying to read the path to the IDRF files, from the IRIS configuration file at :" << a_configurationFile << " !" << endl);
    Cerr("*****                                            ("<< m_nAlphaAnglesIDRF*m_nBetaAnglesIDRF << " lines expected)" << endl);
    return 1;
  }
  
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIRIS::ReadOptionsList(const string& a_optionsList)
{
  // Configuration is possible only through a configuration file (no command-line options)
  Cerr("***** iProjectorIRIS::ReadOptionsList() -> This projector should be initialized with a configuration file, and not parameters (-proj IRIS:path_to_configuration_file) !" << endl);
  return 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iProjectorIRIS::ShowHelpSpecific()
{
  cout << "Multiray projector requiring Intrinsic Detector Response Functions (IDRFs) models" << endl;
  cout << "The projector will generate an user-defined number of lines using from pairs of random points generated using the IDRFs models" << endl;
  cout << "The lines will be rendered using the Incremental Siddon projector (incrementalSiddon)" << endl;
  cout << "This projector requires an initialization file and Monte-Carlo generated IDRFs. For more informations, please read [article]" << endl;
  cout << "INITIALIZATION FILE FIELDS :" << endl;
  cout << "number lines per LOR : number of lines generated to estimate the CDRF (int32)" << endl;
  cout << "number voxels depth : number of voxels in the depth direction in the IDRF volumes (int32)" << endl;
  cout << "number voxels transaxial : number of voxels in the transaxial direction in the IDRF volumes (int32)" << endl;
  cout << "number voxels axial : number of voxels in the axial direction in the IDRF volumes (int32)" << endl;
  cout << "size voxels depth : size of voxels in the depth direction of the IDRF volumes (float32)" << endl;
  cout << "size voxels transaxial : size of voxels in the transaxial direction of the IDRF volume (float32)" << endl;
  cout << "size voxels axial : size of voxels in the axial direction of the IDRF volume (float32)" << endl;
  cout << "number beta angles : number of axial angles in the IDRF volumes (int32)" << endl;
  cout << "number alpha angles : number of transaxial angles in the IDRF volumes (int32)" << endl;
  cout << "step beta angles : axial angular steps in the IDRF volumes (float32)" << endl;
  cout << "step alpha angles : transaxial angular steps in the IDRF volumes (float32)" << endl;
  cout << "IDRFs path : path to the IDRF .vol files (string)" << endl;
  cout << "           : each path should be entered on a separated line" << endl;
  cout << "                                                                                                                                     " << endl;
  cout << "Each .vol files should contain a header of 6 float32 fields :" << endl;
  cout << "number voxels axial ; number voxels transaxial ; number voxels depth ; size voxels axial ; size voxels transaxial ; size voxels depth" << endl;
  cout << "... followed by number voxels depth*number voxels transaxial*number voxels axial float32 elements corresponding to the IDRF coefficients" << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIRIS::CheckSpecificParameters()
{
  if(mp_pathToIDRFFiles==NULL)
  {
    Cerr("***** iProjectorIRIS::CheckSpecificParameters() -> path to IDRF files not initialized !" << endl);
    return 1;
  }
  if(mp_IDRF==NULL || m2p_IDRF_CDFs==NULL)
  {
    Cerr("***** iProjectorIRIS::CheckSpecificParameters() -> IDRF coefficients not initialized !" << endl);
    return 1;
  }
  if(m_nbLinesPerLOR<0)
  {
    Cerr("***** iProjectorIRIS::CheckSpecificParameters() -> number of lines per LOR not initialized !" << endl);
    return 1;
  }
  if(m_nVoxDepthIDRF<0 || m_nVoxTransaxialIDRF<0 || m_nVoxAxialIDRF<0 || m_nVoxXYZIDRF<0)
  {
    Cerr("***** iProjectorIRIS::CheckSpecificParameters() -> number of voxels for the IDRF volumes not initialized !" << endl);
    return 1;
  }
  if(m_sizeVoxDepthIDRF<0 || m_sizeVoxTransaxialIDRF<0 || m_sizeVoxAxialIDRF<0)
  {
    Cerr("***** iProjectorIRIS::CheckSpecificParameters() -> voxel sizes of the IDRF volumes not initialized !" << endl);
    return 1;
  }
  if(m_nBetaAnglesIDRF<0 || m_nAlphaAnglesIDRF<0)
  {
    Cerr("***** iProjectorIRIS::CheckSpecificParameters() -> number of angles for the IDRF volumes not initialized !" << endl);
    return 1;
  }
  if(m_stepBetaAnglesIDRF<0 || m_stepAlphaAnglesIDRF<0)
  {
    Cerr("***** iProjectorIRIS::CheckSpecificParameters() -> step of the IDRF angles not initialized !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIRIS::InitializeSpecific()
{
  // Verbose
  if (m_verbose>=2) Cout("iProjectorIRIS::InitializeSpecific() -> Use IRIS projector" << endl);

  // Read IDRFs coefficients from the raw data files, and write them on a single vector
  for(int a=0 ; a<m_nAlphaAnglesIDRF ; a++) // loop on alpha angles
    for(int b=0 ; b<m_nBetaAnglesIDRF ; b++) // loop on beta angles
    {
      // Open & check file existence
      ifstream IDRF_file(mp_pathToIDRFFiles[a*m_nBetaAnglesIDRF+b].c_str(), ios::binary| ios::in);
      if (!IDRF_file.is_open()) 
      {
        Cerr("***** iProjectorIRIS::InitializeSpecific() -> Error reading the IRDF .vol file at the path: '" << mp_pathToIDRFFiles[a*m_nBetaAnglesIDRF+b] << "'" << endl);
        return 1;
      }
      // Check file size consistency
      IDRF_file.seekg(0, ios::end);
      size_t size_in_bytes = IDRF_file.tellg();

      if(((size_t)((6+m_nVoxXYZIDRF)*sizeof(float))) != size_in_bytes)
      {
        Cerr("***** iProjectorIRIS::InitializeSpecific() -> Error : Size of the .vol " << mp_pathToIDRFFiles[a*m_nBetaAnglesIDRF+b] << " is not consistent with the information provided by the user configuration file!" << endl);
        Cerr("***** iProjectorIRIS::InitializeSpecific() -> Expected size : "<< (6+m_nVoxXYZIDRF)*sizeof(float) << endl);
        Cerr("***** iProjectorIRIS::InitializeSpecific() -> Actual size : "<< size_in_bytes << endl << endl);
        return 1;
      }
      // Check header of the volume 
      float vol_file_header[6];
      // Go back to first character
      IDRF_file.seekg(0);
      for(int i=0 ; i<6 ; i++)
        IDRF_file.read((char*)&vol_file_header[i], sizeof(float));
      // TODO : Take out header from .vol files. Just use the other checks and delete this condition
      if ((int)vol_file_header[0]!=m_nVoxAxialIDRF || 
          (int)vol_file_header[1]!=m_nVoxTransaxialIDRF || 
          (int)vol_file_header[2]!=m_nVoxDepthIDRF || 
          fabs(vol_file_header[3]-m_sizeVoxAxialIDRF) > 0.0001 || 
          fabs(vol_file_header[4]-m_sizeVoxTransaxialIDRF) > 0.0001 ||
          fabs(vol_file_header[5]-m_sizeVoxDepthIDRF) > 0.0001 ) 
          {
            Cerr("***** iProjectorIRIS::InitializeSpecific() -> Inconsistency between initialized and read header file in the .vol: '" << mp_pathToIDRFFiles[a*m_nBetaAnglesIDRF+b] << "'" << endl);
            return 1;
          }
      // Read raw data
      for(int i=0 ; i<m_nVoxXYZIDRF ; i++)
      {
        IDRF_file.read((char*)&vol_file_header[0], sizeof(float));
        mp_IDRF[i] = vol_file_header[0];
      }
      if(ComputeIDRF_CDF(a*m_nBetaAnglesIDRF + b) )
      {
        Cerr("***** iProjectorIRIS::InitializeSpecific() -> An error occurs while computing the coefficient of the following IDRF volume: alpha " << a*m_stepAlphaAnglesIDRF << " beta " << b*m_stepBetaAnglesIDRF << endl);
        return 1;
      }
      IDRF_file.close();
    }

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

INTNB iProjectorIRIS::EstimateMaxNumberOfVoxelsPerLine()
{
  // Find the maximum number of voxels along a given dimension
  INTNB max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxX();
  if (mp_ImageDimensionsAndQuantification->GetNbVoxY()>max_nb_voxels_in_dimension) max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxY();
  if (mp_ImageDimensionsAndQuantification->GetNbVoxZ()>max_nb_voxels_in_dimension) max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxZ();
  // We should have at most 4 voxels in a given plane, so multiply by 4
  // (note: this is not true however it ensures no overflow and is already quite optimized for RAM usage !)
  max_nb_voxels_in_dimension *= 4;
  // Finally multiply by the number of lines
  max_nb_voxels_in_dimension *= ((INTNB)m_nbLinesPerLOR);
  // Return the value
  return max_nb_voxels_in_dimension;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIRIS::ProjectWithoutTOF(int a_direction, oProjectionLine* ap_ProjectionLine )
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iProjectorIRIS::ProjectWithoutTOF() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("iProjectorIRIS::Project without TOF -> Project line '" << ap_ProjectionLine << "' in " << direction << " direction" << endl);
  }
  #endif

  // Get central positions of the crystals
  FLTNB* event1 = ap_ProjectionLine->GetPosition1();
  FLTNB* event2 = ap_ProjectionLine->GetPosition2();
  
  // Local variables to recover crystal central positions and orientations
  float x1, y1, z1, crystal_alpha1, x2, y2, z2, crystal_alpha2;
    
  x1 = event1[0];
  y1 = event1[1];
  z1 = event1[2];
  x2 = event2[0];
  y2 = event2[1];
  z2 = event2[2];

  // Get transaxial orientations of the two crystals // CHECK HERE
  FLTNB* orientation1 = ap_ProjectionLine->GetOrientation1();
  FLTNB* orientation2 = ap_ProjectionLine->GetOrientation2();

  // Compute transaxial angular orientation of the crystals
  crystal_alpha1 = atan2f(orientation1[0], orientation1[1]); //LUT_alpha[id1]; //HERE
  crystal_alpha2 = atan2f(orientation2[0], orientation2[1]); //LUT_alpha[id2]; //HERE

  // -------
  // Compute transaxial incident angles (alpha) between the LOR and the crystals
  // -------
  
  float alpha1, alpha2, alpha_lor;

  // Compute transaxial incident angle of the LOR
  alpha_lor = atan2f(x1-x2, y1-y2); // alpha_lor = atan2f( y1-y2 , x1-x2 ); //CHECK HERE

  // Compute incident angle on each crystal surfaces
  alpha1 = remainderf(alpha_lor-M_PI-crystal_alpha1, M_PI);
  if(alpha1>(M_PI/2.)) // Check if this operation and other surrounding should require HPFLTNB instead of float
    alpha1 -= M_PI;

  alpha2 = remainderf(alpha_lor-M_PI-crystal_alpha2, M_PI);
  if(alpha2>(M_PI/2.))
    alpha2 -= M_PI;

  //alpha1 = remainderf(alpha_lor-Pi-crystal_alpha1, M_PI_2);
  //alpha2 = remainderf(alpha_lor-crystal_alpha2, M_PI_2);
  
  
  
  // -------
  // Compute axial incident angles (beta) between the LOR and the crystals
  // -------
  
  float beta1, beta2, H;
  
  // Compute transaxial distance
  H = sqrtf((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));

  // Compute axial incident angle
  beta1 = atan2f( (z1-z2), H );

  // Normalize with Pi
  beta1 = remainderf(beta1, M_PI);
  if(beta1>(M_PI/2.))
      beta1 -= M_PI;

  beta2 = -beta1;
  beta2 = remainderf(beta2, M_PI);
  if(beta2>(M_PI/2.))
      beta2 -= M_PI;


  /*
  float r_to_deg = 180./M_PI;
  
  //cout << "pos1 X: "  << event1[0] << " Y: " << event1[1]  << " Z: " << event1[2] << " vs Cry: " << alpha1*r_to_deg << " ang " << crystal_alpha1*r_to_deg << " bet " << beta1*r_to_deg << endl;
  //cout << "pos2 X: "  << event2[0] << " Y: " << event2[1]  << " Z: " << event2[2] << " vs Cry: " << alpha2*r_to_deg << " ang " << crystal_alpha2*r_to_deg << " bet " << beta2*r_to_deg <<  endl;
  cout << "pos1 X: "  << event1[0] << " Y: " << event1[1]  << " Z: " << event1[2] << " vs Cry: " << alpha1*r_to_deg << " ang " << crystal_alpha1*r_to_deg << endl;
  cout << "pos2 X: "  << event2[0] << " Y: " << event2[1]  << " Z: " << event2[2] << " vs Cry: " << alpha2*r_to_deg << " ang " << crystal_alpha2*r_to_deg <<  endl;
  cout << ap_ProjectionLine->GetIndex1() << " : " << ap_ProjectionLine->GetIndex2() << " : " << alpha_lor*r_to_deg << endl;

  if(ap_ProjectionLine->GetIndex2()== 17556)
  {
  int r;
  scanf("%d" , &r);
  }
  
  if(ap_ProjectionLine->GetIndex1() == 0 
  && ap_ProjectionLine->GetIndex2() == 50)
  {
        cout << endl;
    cout << "pos1 X: "  << event1[0] << " Y: " << event1[1]  << " Z: " << event1[2] << " aTrs: " << alpha1*r_to_deg  << " aAxl: " << beta1*r_to_deg << endl;
    cout << "pos2 X: "  << event2[0] << " Y: " << event2[1]  << " Z: " << event2[2] << " aTrs: " << alpha2*r_to_deg  << " aAxl: " << beta2*r_to_deg << endl;
  }
  
  
  if(ap_ProjectionLine->GetIndex1() == 0 
  && ap_ProjectionLine->GetIndex2() == 308)
  {
        cout << endl;
    cout << "pos1 X: "  << event1[0] << " Y: " << event1[1]  << " Z: " << event1[2] << " aTrs: " << alpha1*r_to_deg  << " aAxl: " << beta1*r_to_deg << endl;
    cout << "pos2 X: "  << event2[0] << " Y: " << event2[1]  << " Z: " << event2[2] << " aTrs: " << alpha2*r_to_deg  << " aAxl: " << beta2*r_to_deg << endl;
  }

  if(ap_ProjectionLine->GetIndex1() == 0 
  && ap_ProjectionLine->GetIndex2() == 17298)
  {
        cout << endl;
    cout << "pos1 X: "  << event1[0] << " Y: " << event1[1]  << " Z: " << event1[2] << " aTrs: " << alpha1*r_to_deg  << " aAxl: " << beta1*r_to_deg << endl;
    cout << "pos2 X: "  << event2[0] << " Y: " << event2[1]  << " Z: " << event2[2] << " aTrs: " << alpha2*r_to_deg  << " aAxl: " << beta2*r_to_deg << endl;
  }

  if(ap_ProjectionLine->GetIndex1() == 0 
  && ap_ProjectionLine->GetIndex2() == 17556)
  {
    cout << endl;
    cout << "pos1 X: "  << event1[0] << " Y: " << event1[1]  << " Z: " << event1[2] << " aTrs: " << alpha1*r_to_deg  << " aAxl: " << beta1*r_to_deg << endl;
    cout << "pos2 X: "  << event2[0] << " Y: " << event2[1]  << " Z: " << event2[2] << " aTrs: " << alpha2*r_to_deg  << " aAxl: " << beta2*r_to_deg << endl;
  }
  */
  
  
  //////////////////////////////////////////// Draw LORs using Native Siddon   /////////////////////////////////////////////////
  for (int i_line=0; i_line<m_nbLinesPerLOR; ++i_line)
  {
    // Generate two random points using the IDRFs.
    float X1, Y1, Z1;
    float X2, Y2, Z2;
    float pos1[3], pos2[3];
    if (GenerateIRISRdmPos(pos1, alpha1, beta1) )
    {
      Cerr("***** iProjectorIRIS::ProjectWithoutTOF() -> An error occurred while trying to generate random position with the IRIS projector" << endl);
      return 1; 
    }
  
    X1 = pos1[0];
    Y1 = pos1[1];
    Z1 = pos1[2];
    
    if (GenerateIRISRdmPos(pos2, alpha2, beta2) )
    {
      Cerr("***** iProjectorIRIS::ProjectWithoutTOF() -> An error occurred while trying to generate random position with the IRIS projector" << endl);
      return 1; 
    }
    
    X2 = pos2[0];
    Y2 = pos2[1];
    Z2 = pos2[2];
    
    // Move the random points in the frame of reference of the scanner.
    // IRIS ref:
    // X = depth
    // Y = transaxial
    // Z = axial
    
    float tmp_X, tmp_Y;
    tmp_X = X1; tmp_Y = Y1;
    X1 = x1 + tmp_X*orientation1[0]-tmp_Y*orientation1[1]; //HERE
    Y1 = y1 + tmp_X*orientation1[1]+tmp_Y*orientation1[0]; //HERE
    Z1 += z1;
    /*
    cout << "pos1 X: "  << event1[0] << " Y: " << event1[1]  << " Z: " << event1[2] << endl ;
    cout << orientation1[0] << " ; " << orientation1[1] << endl;
    cout << "tmp_X1: " << tmp_X << " tmp_Y1 " << tmp_Y << " --- X1: " << X1 << " Y1 " << Y1 << endl;
    */
    
    tmp_X = X2; tmp_Y = Y2;
    X2 = x2 + tmp_X*orientation2[0]-tmp_Y*orientation2[1]; //HERE
    Y2 = y2 + tmp_X*orientation2[1]+tmp_Y*orientation2[0]; //HERE
    Z2 += z2;

/*
    cout << "pos2 X: "  << event2[0] << " Y: " << event2[1]  << " Z: " << event2[2] << endl;
    cout << orientation2[0] << " ; " << orientation2[1] << endl;
    cout << "tmp_X2: " << tmp_X << " tmp_Y2 " << tmp_Y << " --- X2: " << X2 << " Y2 " << Y2 << endl;
  int r;
  scanf("%d", &r);*/

    // SIDDON ALGORITHM
    
    // **************************************
    // STEP 1: LOR length calculation
    // **************************************
    HPFLTNB length_LOR = ap_ProjectionLine->GetLength();
  
    // **************************************
    // STEP 2: Compute entrance voxel indexes
    // **************************************
    HPFLTNB alphaFirst[] = { 0.0, 0.0, 0.0 };
    HPFLTNB alphaLast[] = { 0.0, 0.0, 0.0 };
  
    HPFLTNB alphaMin = 0.0f, alphaMax = 1.0f;
    HPFLTNB delta_pos[] = {
      X2 - X1, // event2[ 0 ] - event1[ 0 ]
      Y2 - Y1, // event2[ 1 ] - event1[ 1 ]
      Z2 - Z1  // event2[ 2 ] - event1[ 2 ]
    };
  
    // Computation of alphaMin et alphaMax values (entrance and exit point of the ray)
    for( int i = 0; i < 3; ++i )
    {
      if( delta_pos[ i ] != 0.0 )
      {
        alphaFirst[i] = (-mp_halfFOV[i] - event1[i]) / (delta_pos[ i ]);
        alphaLast[i]  = ( mp_halfFOV[i] - event1[i]) / (delta_pos[ i ]);
        alphaMin = (std::max)(alphaMin,(std::min)(alphaFirst[i],alphaLast[i]));
        alphaMax = (std::min)(alphaMax,(std::max)(alphaFirst[i],alphaLast[i]));
      }
    }
  
    // if alphaMax is less than or equal to alphaMin no intersection
    // and return an empty buffer
    if( alphaMax <= alphaMin ) return 0;
  
    // Now we have to find the indices of the particular plane
    // (iMin,iMax), (jMin,jMax), (kMin,kMax)
    int iMin = 0, iMax = 0;
    int jMin = 0, jMax = 0;
    int kMin = 0, kMax = 0;
  
    // For the X-axis
    if( delta_pos[ 0 ] > 0.0f )
    {
      iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMin * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
      iMax = ::floor( 1 + ( event1[ 0 ] + alphaMax * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
    }
    else if( delta_pos[ 0 ] < 0.0f )
    {
      iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMax * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
      iMax = ::floor( 1 + ( event1[ 0 ] + alphaMin * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
    }
    if( delta_pos[ 0 ] == 0 )
    {
      iMin = 1, iMax = 0;
    }
  
    // For the Y-axis
    if( delta_pos[ 1 ] > 0 )
    {
      jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMin * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
      jMax = ::floor( 1 + ( event1[ 1 ] + alphaMax * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
    }
    else if( delta_pos[ 1 ] < 0 )
    {
      jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMax * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
      jMax = ::floor( 1 + ( event1[ 1 ] + alphaMin * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
    }
    else if( delta_pos[ 1 ] == 0 )
    {
      jMin = 1, jMax = 0;
    }
  
    // For the Z-axis
    if( delta_pos[ 2 ] > 0 )
    {
      kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMin * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
      kMax = ::floor( 1 + ( event1[ 2 ] + alphaMax * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
    }
    else if( delta_pos[ 2 ] < 0 )
    {
      kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMax * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
      kMax = ::floor( 1 + ( event1[ 2 ] + alphaMin * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
    }
    else if( delta_pos[ 2 ] == 0 )
    {
      kMin = 1, kMax = 0;
    }
  
    // Computing the last term n number of intersection
    INTNB n = ( iMax - iMin + 1 ) + ( jMax - jMin + 1 )
      + ( kMax - kMin + 1 );
  
    // Array storing the first alpha in X, Y and Z
    HPFLTNB alpha_XYZ[ 3 ] = { 1.0f, 1.0f, 1.0f };
  
    // Computing the first alpha in X
    if( delta_pos[ 0 ] > 0 )
      alpha_XYZ[ 0 ] = ( ( (-mp_halfFOV[ 0 ]) + ( iMin - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];
    else if( delta_pos[ 0 ] < 0 )
      alpha_XYZ[ 0 ] = ( ( (-mp_halfFOV[ 0 ]) + ( iMax - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];
  
    // Computing the first alpha in Y
    if( delta_pos[ 1 ] > 0 )
      alpha_XYZ[ 1 ] = ( ( (-mp_halfFOV[ 1 ]) + ( jMin - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];
    else if( delta_pos[ 1 ] < 0 )
      alpha_XYZ[ 1 ] = ( ( (-mp_halfFOV[ 1 ]) + ( jMax - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];
  
    // Computing the first alpha in Z
    if( delta_pos[ 2 ] > 0 )
      alpha_XYZ[ 2 ] = ( ( (-mp_halfFOV[ 2 ]) + ( kMin - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];
    else if( delta_pos[ 2 ] < 0 )
      alpha_XYZ[ 2 ] = ( ( (-mp_halfFOV[ 2 ]) + ( kMax - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];
  
    // Computing the alpha updating
    HPFLTNB alpha_u[ 3 ] = {
      mp_sizeVox[0] / std::fabs( delta_pos[ 0 ] ),
      mp_sizeVox[1] / std::fabs( delta_pos[ 1 ] ),
      mp_sizeVox[2] / std::fabs( delta_pos[ 2 ] )
    };
  
    // Computing the index updating 
    INTNB index_u[ 3 ] = {
      (delta_pos[ 0 ] < 0) ? -1 : 1,
      (delta_pos[ 1 ] < 0) ? -1 : 1,
      (delta_pos[ 2 ] < 0) ? -1 : 1
    };
  
    // Check which alpha is the min/max and increment
    if( alpha_XYZ[ 0 ] == alphaMin ) alpha_XYZ[ 0 ] += alpha_u[ 0 ];
    if( alpha_XYZ[ 1 ] == alphaMin ) alpha_XYZ[ 1 ] += alpha_u[ 1 ];
    if( alpha_XYZ[ 2 ] == alphaMin ) alpha_XYZ[ 2 ] += alpha_u[ 2 ];
  
    // Computing the minimum value in the alpha_XYZ buffer
    HPFLTNB const min_alpha_XYZ = (std::min)( alpha_XYZ[ 0 ],
      (std::min)( alpha_XYZ[ 1 ], alpha_XYZ[ 2 ] ) );
  
    // Computing the first index of intersection
    HPFLTNB const alphaMid = ( min_alpha_XYZ + alphaMin ) / 2.0f;
    INTNB index_ijk[ 3 ] = {
      1 + (int)( ( event1[ 0 ] + alphaMid * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] ),
      1 + (int)( ( event1[ 1 ] + alphaMid * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] ),
      1 + (int)( ( event1[ 2 ] + alphaMid * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] )
    };
  
    INTNB const w = mp_nbVox[ 0 ];
    INTNB const wh = w * mp_nbVox[ 1 ];
  
    // Loop over the number of plane to cross
    HPFLTNB alpha_c = alphaMin;
    FLTNB coeff = 0.0f;
    INTNB numVox = 0;
    for( int nP = 0; nP < n - 1; ++nP )
    {
      if( ( alpha_XYZ[ 0 ] <= alpha_XYZ[ 1 ] )
       && ( alpha_XYZ[ 0 ] <= alpha_XYZ[ 2 ] ) )
      {
        // Storing values
        if( ( alpha_XYZ[ 0 ] >= alphaMin )
         && ( index_ijk[ 0 ] - 1 >= 0 )
         && ( index_ijk[ 0 ] - 1 <= mp_nbVox[ 0 ] - 1 )
         && ( index_ijk[ 1 ] - 1 >= 0 )
         && ( index_ijk[ 1 ] - 1 <= mp_nbVox[ 1 ] - 1 )
         && ( index_ijk[ 2 ] - 1 >= 0 )
         && ( index_ijk[ 2 ] - 1 <= mp_nbVox[ 2 ] - 1 ) )
        {
          coeff = ( alpha_XYZ[ 0 ] - alpha_c ) * length_LOR;
          numVox = ( index_ijk[ 0 ] - 1 ) + ( ( index_ijk[ 1 ] - 1 ) ) * w + ( ( index_ijk[ 2 ] - 1 ) ) * wh;
          ap_ProjectionLine->AddVoxel(a_direction, numVox, coeff/m_nbLinesPerLOR); // TODO Normalizations according to the number of lines are performed in these lines. Perhaps we should do this at a upper level (vProjector)
        }
  
        // Increment values
        alpha_c = alpha_XYZ[ 0 ];
        alpha_XYZ[ 0 ] += alpha_u[ 0 ];
        index_ijk[ 0 ] += index_u[ 0 ];
      }
      else if( ( alpha_XYZ[ 1 ] < alpha_XYZ[ 0 ] )
            && ( alpha_XYZ[ 1 ] <= alpha_XYZ[ 2 ] ) )
      {
        // Storing values
        if( ( alpha_XYZ[ 1 ] >= alphaMin )
         && ( index_ijk[ 0 ] - 1 >= 0 )
         && ( index_ijk[ 0 ] - 1 <= mp_nbVox[ 0 ] - 1 )
         && ( index_ijk[ 1 ] - 1 >= 0 )
         && ( index_ijk[ 1 ] - 1 <= mp_nbVox[ 1 ] - 1 )
         && ( index_ijk[ 2 ] - 1 >= 0 )
         && ( index_ijk[ 2 ] - 1 <= mp_nbVox[ 2 ] - 1 ) )
        {
          coeff = ( alpha_XYZ[ 1 ] - alpha_c ) * length_LOR;
          numVox = ( index_ijk[ 0 ] - 1 ) + ( ( index_ijk[ 1 ] - 1 ) ) * w + ( ( index_ijk[ 2 ] - 1 ) ) * wh;
          ap_ProjectionLine->AddVoxel(a_direction, numVox, coeff/m_nbLinesPerLOR);
        }
  
        // Increment values
        alpha_c = alpha_XYZ[ 1 ];
        alpha_XYZ[ 1 ] += alpha_u[ 1 ];
        index_ijk[ 1 ] += index_u[ 1 ];
      }
      else if( ( alpha_XYZ[ 2 ] < alpha_XYZ[ 0 ] )
            && ( alpha_XYZ[ 2 ] < alpha_XYZ[ 1 ] ) )
      {
        // Storing values
        if( ( alpha_XYZ[ 2 ] >= alphaMin )
         && ( index_ijk[ 0 ] - 1 >= 0 )
         && ( index_ijk[ 0 ] - 1 <= mp_nbVox[ 0 ] - 1 )
         && ( index_ijk[ 1 ] - 1 >= 0 )
         && ( index_ijk[ 1 ] - 1 <= mp_nbVox[ 1 ] - 1 )
         && ( index_ijk[ 2 ] - 1 >= 0 )
         && ( index_ijk[ 2 ] - 1 <= mp_nbVox[ 2 ] - 1 ) )
        {
          coeff = ( alpha_XYZ[ 2 ] - alpha_c ) * length_LOR;
          numVox = ( index_ijk[ 0 ] - 1 ) + ( ( index_ijk[ 1 ] - 1 ) ) * w + ( ( index_ijk[ 2 ] - 1 ) ) * wh;
          ap_ProjectionLine->AddVoxel(a_direction, numVox, coeff/m_nbLinesPerLOR);
        }
  
        // Increment values
        alpha_c = alpha_XYZ[ 2 ];
        alpha_XYZ[ 2 ] += alpha_u[ 2 ];
        index_ijk[ 2 ] += index_u[ 2 ];
      }
    }
    
  }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIRIS::ProjectWithTOFPos(int a_Projector, oProjectionLine* ap_ProjectionLine)
{
  Cerr("***** iProjectorIRIS::ProjectWithTOFPos() -> Not yet implemented !" << endl);
  return 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIRIS::ProjectWithTOFBin(int a_Projector, oProjectionLine* ap_ProjectionLine)
{
  Cerr("***** iProjectorIRIS::ProjectWithTOFBin() -> Not yet implemented !" << endl);
  return 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIRIS::ComputeIDRF_CDF(int a_angleId)
{  
  // Compute cumulated sum
  m2p_IDRF_CDFs[a_angleId][0] = mp_IDRF[0];
    
  for(int i = 1; i<m_nVoxXYZIDRF; ++i)
    m2p_IDRF_CDFs[a_angleId][i] = m2p_IDRF_CDFs[a_angleId][i-1] + mp_IDRF[i];
    
  // Normalization
  if(m2p_IDRF_CDFs[a_angleId][m_nVoxXYZIDRF-1] > 0) // check if not null
    for(int i = 0; i<m_nVoxXYZIDRF; ++i)
      m2p_IDRF_CDFs[a_angleId][i] /= m2p_IDRF_CDFs[a_angleId][m_nVoxXYZIDRF-1];
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIRIS::GenerateIRISRdmPos(float ap_generatedPos[3], float a_alpha, float a_beta)
{
  sRandomNumberGenerator* p_RNG = sRandomNumberGenerator::GetInstance(); 
  
  float rdm_pos[3];

  int alpha_id = (int) (fabsf(a_alpha)/m_stepAlphaAnglesIDRF);
  int beta_id  = (int) (fabsf(a_beta) /m_stepBetaAnglesIDRF);
  
  alpha_id = min(alpha_id, m_nAlphaAnglesIDRF-1);
  beta_id  = min(beta_id , m_nBetaAnglesIDRF -1);

  int id = FindGreaterValue(m2p_IDRF_CDFs[alpha_id*m_nBetaAnglesIDRF+beta_id], p_RNG->GenerateRdmNber(), m_nVoxXYZIDRF);
  
  // IDRFs indices sorting are axial/transaxial/DOI, DOI_id being the smallest index.
  int axial_id   =  id / (m_nVoxDepthIDRF*m_nVoxTransaxialIDRF);
  int trans_id   = (id-axial_id*m_nVoxTransaxialIDRF*m_nVoxDepthIDRF) / m_nVoxDepthIDRF;
  int DOI_id     =  id-axial_id*m_nVoxTransaxialIDRF*m_nVoxDepthIDRF - trans_id*m_nVoxDepthIDRF;

  // Use random number on each axis to generate a random position on the IDRF voxel
  rdm_pos[0] = (DOI_id   + p_RNG->GenerateRdmNber() - m_nVoxDepthIDRF/2.0f      )*m_sizeVoxDepthIDRF; // CHECK : should choose surface central position by default when using IRIS, or adress the difference here
  rdm_pos[1] = (trans_id + p_RNG->GenerateRdmNber() - m_nVoxTransaxialIDRF/2.0f )*m_sizeVoxTransaxialIDRF;
  rdm_pos[2] = (axial_id + p_RNG->GenerateRdmNber() - m_nVoxAxialIDRF/2.0f      )*m_sizeVoxAxialIDRF;

  if(a_alpha<0) rdm_pos[1] = -rdm_pos[1];
  if(a_beta<0)  rdm_pos[2] = -rdm_pos[2];

  ap_generatedPos[0] = rdm_pos[0];
  ap_generatedPos[1] = rdm_pos[1];
  ap_generatedPos[2] = rdm_pos[2];
  
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIRIS::FindGreaterValue(float *ap_val, float a_key, int a_maxValue, int a_minStart, int a_maxStart)
{
  int min = a_minStart, max, mid;
  
  max = a_maxStart ? a_maxStart : a_maxValue;
  
  int cpt=0;
  
  while (min < max)
  {
    mid = (min + max) >> 1; //(min+max)/2
    if (a_key > ap_val[mid])
    {
      min = mid + 1;
    } 
    else
    {
      max = mid;
    }
    
    cpt++;
  }
  
  return min;
}
