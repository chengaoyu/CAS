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
  \ingroup utils
  \brief This program generates a CASToR Look-Up-Table (LUT) containing geometry information related to each crystal of a PET system. \n
         The LUT contains x,y,z cartesian position of the center of each crystals, as well as their transaxial angular orientations. \n
         The LUT files is composed of two files : \n
         (1) A binary file containing for each crystal: Cartesian positions (X,Y,Z) of the center of the crystals, and Vector orientation (X,Y,Z), corresponding to the angular orientation of each crystal (see documentation) \n
         (2) A header file containing several informations about the system. \n
         CASToR algorithms will keep a similar crystal indexation as the LUT was written, and load it in RAM \n
         This script uses a GATE model of a General Electric PET system as example.
*/

#include "gVariables.hh"
#include "gOptions.hh"
#include "sOutputManager.hh"
#include "oMatrix.hh"


/*!
  \fn ShowHelp
  \param a_returnCode
  \brief Show usage
*/
void ShowHelp(int a_returnCode)
{
  cout << endl;
  cout << "Usage:  castor-PetScannerLutEx  -alias scanner_name  [settings]" << endl;
  cout << endl;
  cout << "[Input settings]:" << endl;
  cout << "  -alias scanner_name   : give the alias of the scanner for which the LUT will be generated (suggested template Modality-Constructor-Model (ex: PET_GE_DLS)" << endl;
  cout << "                          the resulting file will be written in the scanner repository (default : /config/scanner directory)" << endl;
  cout << endl;
  #ifdef CASTOR_VERSION
  cout << "  This program is part of the CASToR release version " << CASTOR_VERSION << "." << endl;
  cout << endl;
  #endif
  Exit(a_returnCode);
}


/*!
  \fn main
  \brief LUT Generator for PET scanner whose architecture consists in rotational sectors, blocs (modules/submodules) and crystals
*/
int main(int argc, char** argv)
{
  // No argument, then show help
  if (argc==1) ShowHelp(0);
 
  string scanner_name = "";
  string path_to_LUT = "";
  string path_to_headerLUT = "";
  ofstream LUT_file, header_LUT_file;
  
  // ============================================================================================================
  // Read command-line parameters
  // ============================================================================================================
  for (int i=1; i<argc; i++)
  {
    string option = (string)argv[i];

    if (option=="-h" || option=="--help" || option=="-help") ShowHelp(0);

    else if (option=="-alias")
    {
      if (i>=argc-1)
      {
        cerr << "***** castor-PetScannerLutEx :: Argument missing for option: " << option << endl;
        Exit(EXIT_FAILURE);
      }
      scanner_name = argv[i+1];
      string path_base = sOutputManager::GetInstance()->GetPathToConfigDir();

      path_base.append("scanner");
      path_to_LUT = path_base+OS_SEP;
      path_to_headerLUT = path_base+OS_SEP;
      path_to_LUT.append(scanner_name.c_str()).append(".lut");
      path_to_headerLUT.append(scanner_name.c_str()).append(".hscan");
      i++;
    }
    else
    {
      cerr << "***** castor-PetScannerLutEx :: Unknown option '" << option << "' !" << endl;
      Exit(EXIT_FAILURE);
    }
  }



  // ============================================================================================================
  // Checks options have been provided
  // ============================================================================================================

  // Checks
  // output directory
  if (scanner_name.empty() )
  {
    cerr << "***** castor-PetScannerLutEx :: Please provide an alias for this scanner !" << endl;
    ShowHelp(0);
    Exit(EXIT_FAILURE);
  }
  else
  {
    cout << endl << "Generating LUT with the following alias: " << scanner_name << "... " << endl << endl;
  }



  // ============================================================================================================
  // Input parameter declarations
  // ============================================================================================================
  int nb_rings;
  int nb_elts;
  FLTNBLUT min_trs_angle_diff;
  string scanner_modality;
  string description;
  FLTNBLUT angular_span;
            
  int *nb_rsectors_lyr,
      *nb_trans_mod_lyr,
      *nb_axial_mod_lyr,
      *nb_trans_submod_lyr,
      *nb_axial_submod_lyr,
      *nb_trans_crystal_lyr,
      *nb_axial_crystal_lyr,
      *nb_crystals_lyr;
              
  FLTNBLUT  *radius_lyr,
            *gap_trans_mod_lyr,
            *gap_axial_mod_lyr,
            *gap_trans_submod_lyr,
            *gap_axial_submod_lyr,
            *gap_trans_crystal_lyr,
            *gap_axial_crystal_lyr,
            *crystal_size_depth_lyr,  
            *crystal_size_trans_lyr,
            *crystal_size_axial_lyr,
            *mean_depth_of_interaction_lyr;

  
  // Layer-dependent variables
  int nbLayers = 1;
  nb_rsectors_lyr = new int[nbLayers];
  nb_trans_mod_lyr = new int[nbLayers];
  nb_axial_mod_lyr = new int[nbLayers];
  nb_trans_submod_lyr = new int[nbLayers];
  nb_axial_submod_lyr = new int[nbLayers];
  nb_trans_crystal_lyr = new int[nbLayers];
  nb_axial_crystal_lyr = new int[nbLayers];
  nb_crystals_lyr = new int[nbLayers];
  
  radius_lyr = new FLTNBLUT[nbLayers];
  
  gap_trans_mod_lyr = new FLTNBLUT[nbLayers];
  gap_axial_mod_lyr = new FLTNBLUT[nbLayers];
  gap_trans_submod_lyr = new FLTNBLUT[nbLayers];
  gap_axial_submod_lyr = new FLTNBLUT[nbLayers];
  gap_trans_crystal_lyr = new FLTNBLUT[nbLayers];
  gap_axial_crystal_lyr = new FLTNBLUT[nbLayers];
  crystal_size_depth_lyr = new FLTNBLUT[nbLayers];  
  crystal_size_trans_lyr = new FLTNBLUT[nbLayers];
  crystal_size_axial_lyr = new FLTNBLUT[nbLayers];
  mean_depth_of_interaction_lyr = new FLTNBLUT[nbLayers];


  // Initialize value of each scanner element according to the scanner system
  
  description = "User-made LUT of a GATE model of the GE DRX PET scanner system, generated by the castor-PetScannerLutEx script";
  scanner_modality = "PET";
  
  // System minimal transaxial angle difference between two scanner elements
  // to get a LOR 
  min_trs_angle_diff = 40.;
  
  // radius in mm (from isocenter to crystal surface) 
  radius_lyr[0]= 443; 
  
  // nb scanner elements
  nb_rsectors_lyr[0] = 70;
  nb_trans_mod_lyr[0] = 1;
  nb_axial_mod_lyr[0] = 4;
  nb_trans_submod_lyr[0] = 1;
  nb_axial_submod_lyr[0] = 1;
  nb_trans_crystal_lyr[0] = 9;
  nb_axial_crystal_lyr[0] = 6;


  // Gaps between scanner elements
  gap_trans_mod_lyr[0] = 0;
  gap_axial_mod_lyr[0] = 1.75;
  gap_trans_submod_lyr[0] = 0;
  gap_axial_submod_lyr[0] = 0;
  gap_trans_crystal_lyr[0] = 0.065;
  gap_axial_crystal_lyr[0] = 0.1;

  // crystal dimensions (mm)
  crystal_size_depth_lyr[0] = 30;  
  crystal_size_trans_lyr[0] = 4.230;
  crystal_size_axial_lyr[0] = 6.350;
  
  // mean depth of interaction in the crystal (mm) 
  // negative value means no depth of interaction
  mean_depth_of_interaction_lyr[0] = -1.;
  
  // angular span in degree = 360 by defaut
  // (rsectors will be uniformly positionned according to this value)
  angular_span = 360.;
  
  // Default reconstruction parameters for the scanner
  int default_dim_trans, default_dim_axial;
  default_dim_trans = 256;
  default_dim_axial = 47;
  
  // default field of view
  FLTNBLUT default_FOV_trans, default_FOV_axial; //mm
  default_FOV_trans = 700;
  default_FOV_axial = 153.69;


  // Z-shifts (rsector axial shift for each module)
      
  // Default initialization
  int nb_rsctr_axial_shift = 1;
  FLTNBLUT *rsctr_zshift;
  rsctr_zshift = new FLTNBLUT[nb_rsctr_axial_shift];
  
  // System contains z-shifts)
  if(nb_rsctr_axial_shift > 1)
  {
    for(int zs=0 ; zs<nb_rsctr_axial_shift ; zs++)
      rsctr_zshift[zs] = 0.; // Add specific z-shift values for each layer
  }
  // No z-shift, default initialization
  else
    rsctr_zshift[0] = 0.;


  // Compute the total number of elements
  nb_elts = 0;
  for(int lyr=0 ; lyr<nbLayers ; lyr++)
    nb_elts += nb_rsectors_lyr[lyr] * nb_trans_mod_lyr[lyr] * nb_axial_mod_lyr[lyr] 
                                    * nb_trans_submod_lyr[lyr] * nb_axial_submod_lyr[lyr] 
                                    * nb_trans_crystal_lyr[lyr] * nb_axial_crystal_lyr[lyr];

  // Cumulative number of crystal (in order to keep track of crystal index if nb_layer>1)
  uint32_t nb_cry_cur = 0;
  
  // Number of crystal in this layer
  uint32_t nb_cry_in_layer = 0;
  
  // Loop on layers. 
  // Write the crytal LUT elements successively for each layer
  for(int lyr=0 ; lyr<nbLayers ; lyr++)
  {
    // Default position for the first rsector in CASToR is directly above isocenter
    // The gate model of this scanner locates the first rsector on the right side of the scanner
    // This variable allows to provides the position of the first rsector in comparison with
    // CASToR convention
    // This is required to recover angular orientations of the crystals
    FLTNBLUT rsector_first_angle = 90;
    
    int nb_rsectors = nb_rsectors_lyr[lyr],
        nb_trans_mod = nb_trans_mod_lyr[lyr],
        nb_axial_mod = nb_axial_mod_lyr[lyr],
        nb_trans_submod = nb_trans_submod_lyr[lyr],
        nb_axial_submod = nb_axial_submod_lyr[lyr],
        nb_trans_crystal = nb_trans_crystal_lyr[lyr],
        nb_axial_crystal = nb_axial_crystal_lyr[lyr];
  
    FLTNBLUT  radius = radius_lyr[lyr],
              gap_trans_mod = gap_trans_mod_lyr[lyr],
              gap_axial_mod = gap_axial_mod_lyr[lyr],
              gap_trans_submod = gap_trans_submod_lyr[lyr],
              gap_axial_submod = gap_axial_submod_lyr[lyr],
              gap_trans_crystal = gap_trans_crystal_lyr[lyr],
              gap_axial_crystal = gap_axial_crystal_lyr[lyr],
              crystal_size_trans = crystal_size_trans_lyr[lyr],
              crystal_size_axial = crystal_size_axial_lyr[lyr],
              crystal_size_depth = crystal_size_depth_lyr[lyr];
  
    // Compute system element sizes
    FLTNBLUT size_trans_submod = nb_trans_crystal*crystal_size_trans + (nb_trans_crystal-1)*gap_trans_crystal;
    FLTNBLUT size_axial_submod = nb_axial_crystal*crystal_size_axial + (nb_axial_crystal-1)*gap_axial_crystal;
    FLTNBLUT size_trans_mod = nb_trans_submod*size_trans_submod + (nb_trans_submod-1)*gap_trans_submod;
    FLTNBLUT size_axial_mod = nb_axial_submod*size_axial_submod + (nb_axial_submod-1)*gap_axial_submod;
  
    int nb_mod = nb_axial_mod*nb_trans_mod;
    int nb_submod = nb_axial_submod*nb_trans_submod;
    int nb_crystal = nb_trans_crystal*nb_axial_crystal;
  
    nb_cry_in_layer = nb_rsectors
                    * nb_mod
                    * nb_submod
                    * nb_crystal;
    
    nb_crystals_lyr[lyr] = nb_cry_in_layer;
    
    nb_rings = nb_rsectors
             * nb_trans_mod
             * nb_trans_submod
             * nb_trans_crystal;
                    
    int number_crystals_in_ring = nb_crystals_lyr[lyr]/nb_rings;
  
    // Variables gathering the LUT elements
    FLTNBLUT* crystal_positionX = new FLTNBLUT[ nb_crystals_lyr[lyr] ];
    FLTNBLUT* crystal_positionY = new FLTNBLUT[ nb_crystals_lyr[lyr] ];
    FLTNBLUT* crystal_positionZ = new FLTNBLUT[ nb_crystals_lyr[lyr] ];
    FLTNBLUT* crystal_orientationX = new FLTNBLUT[ nb_crystals_lyr[lyr] ];
    FLTNBLUT* crystal_orientationY = new FLTNBLUT[ nb_crystals_lyr[lyr] ];
    FLTNBLUT* crystal_orientationZ = new FLTNBLUT[ nb_crystals_lyr[lyr] ];
  
  
    // ============================================================================================================
    // Main part of the program: Generate the LUT
    // ============================================================================================================

    // Loop to nb_rsectors+1. crystal_center[0] will be used to gather position of the reference rsector (directly above isocenter)
    oMatrix *****crystal_center = new oMatrix ****[nb_rsectors];
  
    for(int i = 0; i < nb_rsectors+1 ; i++)
    {
      crystal_center[i] = new oMatrix ***[nb_axial_mod*nb_trans_mod];
  
      for (int j = 0; j<nb_axial_mod*nb_trans_mod; j++)
      {
        crystal_center[i][j] = new oMatrix **[nb_axial_submod*nb_trans_submod];
  
        for (int k = 0; k<nb_axial_submod*nb_trans_submod; k++)
        {
          crystal_center[i][j][k] = new oMatrix*[nb_axial_crystal*nb_trans_crystal];
  
          for (int l = 0; l<nb_axial_crystal*nb_trans_crystal; l++)
            crystal_center[i][j][k][l]  = new oMatrix(3,1);
        }
      }
    }
    
  
    // ============================================================================================================
    // Generation of the rotation matrix allowing to compute the position of all the rsectors. 
    // ============================================================================================================
    oMatrix** rotation_mtx = new oMatrix*[nb_rsectors];
  
    for(int i=0; i<nb_rsectors; i++)
      rotation_mtx[i] = new oMatrix(3,3);

    FLTNBLUT angular_span_rad = angular_span*M_PI/180.;
    for (int i = 0; i<nb_rsectors; i++)
    {
      FLTNBLUT angle = remainderf((FLTNB)i*angular_span_rad/((FLTNB)nb_rsectors), 2.*M_PI);
      
      rotation_mtx[i]->SetMatriceElt(0,0,cos(angle) );
      rotation_mtx[i]->SetMatriceElt(1,0,-sin(angle) );
      rotation_mtx[i]->SetMatriceElt(2,0,0);
      rotation_mtx[i]->SetMatriceElt(0,1,sin(angle) );
      rotation_mtx[i]->SetMatriceElt(1,1,cos(angle) );
      rotation_mtx[i]->SetMatriceElt(2,1,0);
      rotation_mtx[i]->SetMatriceElt(0,2,0);
      rotation_mtx[i]->SetMatriceElt(1,2,0);
      rotation_mtx[i]->SetMatriceElt(2,2,1);
    }



    // ============================================================================================================
    // Compute scanner elements positions for the first rsector 
    // (For the example scanner, it is located and centered on the right side of isocenter in the GATE model
    // ============================================================================================================
    
    for (int i=0; i < nb_mod ; i++)
    {
      // Define the transaxial and axial edge start positions for the rsector
      FLTNBLUT y_start_m = (nb_trans_mod*size_trans_mod + (nb_trans_mod-1)*gap_trans_mod) / 2;
      FLTNBLUT z_start_m = -(nb_axial_mod*size_axial_mod + (nb_axial_mod-1)*gap_axial_mod) / 2 ;
  
      // Define the transaxial and axial edge start positions for the i-Module in the rsector. 
      // Enumeration starting with the transaxial modules.
      y_start_m -= (i%nb_trans_mod) *  (size_trans_mod + gap_trans_mod);
      z_start_m += int(i/nb_trans_mod) * (size_axial_mod + gap_axial_mod);
  
      for (int j=0 ; j < nb_submod ; j++)
      {
        FLTNBLUT y_start_sm = y_start_m;
        FLTNBLUT z_start_sm = z_start_m;
        
        y_start_sm -= (j%nb_trans_submod) *  (size_trans_submod + gap_trans_submod);
        z_start_sm += int(j/nb_trans_submod) * (size_axial_submod + gap_axial_submod);
         
        for (int k=0 ; k < nb_crystal ; k++) 
        {
          // Define the transaxial and axial center positions for the j-SubModule (crystal) i-Module of the rsector.
          // Enumeration starting with the transaxial submodules.
          FLTNBLUT Xcrist = radius + crystal_size_depth/2;
          FLTNBLUT Ycrist = y_start_sm - (k%nb_trans_crystal) * (crystal_size_trans + gap_trans_crystal) - crystal_size_trans/2; 
          FLTNBLUT Zcrist = z_start_sm + int(k/nb_trans_crystal) * (crystal_size_axial + gap_axial_crystal) + crystal_size_axial/2;
          
          crystal_center[0][i][j][k]->SetMatriceElt(0,0,Xcrist);
          crystal_center[0][i][j][k]->SetMatriceElt(1,0,Ycrist);
          crystal_center[0][i][j][k]->SetMatriceElt(2,0,Zcrist);
        }
      }
    }
  

    // ============================================================================================================
    // Loop over all the other rsectors.
    // Mandatory informations about crystals are recovered in :
    // crystal_position(X,Y,Z) : Cartesian positions of the center of the crystals
    // crystal_orientation(X,Y,Z) : Vector orientation recovering crystal angles 
    //                             (usually identical for all crystals inside the same rsector)
    // ============================================================================================================

    for (int rs=0 ; rs<nb_rsectors ; rs++)
    {

      // Compute angle for orientation vector for this rsector
      FLTNBLUT rsector_first_angle_rad = rsector_first_angle*M_PI/180.;  
      FLTNBLUT orientation_angle = remainderf(rsector_first_angle_rad + (FLTNB)rs*angular_span_rad/((FLTNB)nb_rsectors), 2.*M_PI);

      for (int j=0 ; j<nb_mod ; j++)
        for (int k=0 ; k<nb_submod ; k++)
          for (int l=0 ; l<nb_crystal ; l++)
          {
            // crystal indexation
            int cryID = int(j/nb_trans_mod)*nb_axial_submod*nb_axial_crystal*number_crystals_in_ring // = nb indexed crystals in the rings covered by the previous (axial) modules
                      + int(k/nb_trans_submod)*nb_axial_crystal*number_crystals_in_ring // = nb indexed crystals in the rings covered by the previous (axial) submodules
                      + int(l/nb_trans_crystal)*number_crystals_in_ring // = nb indexed crystals in the rings covered by the previous (axial) crystals
                      + rs*nb_trans_mod*nb_trans_submod*nb_trans_crystal // = nb indexed crystals in the previous rsectors
                      + j/nb_axial_mod*nb_trans_submod*nb_trans_crystal // = nb indexed crystals in the previous modules
                      + k/nb_axial_submod*nb_trans_crystal // = nb indexed crystals in the previous submodules
                      + l%nb_trans_crystal // previous crystals in the submodule
                      + nb_cry_cur; // = number of crystals already indexed if lyr>0

            rotation_mtx[rs]->Multiplication(crystal_center[0][j][k][l], crystal_center[rs+1][j][k][l]);
            crystal_positionX[cryID]  = crystal_center[rs+1][j][k][l]->GetMatriceElt(0,0);
            crystal_positionY[cryID]  = crystal_center[rs+1][j][k][l]->GetMatriceElt(1,0);
            crystal_positionZ[cryID]  = crystal_center[rs+1][j][k][l]->GetMatriceElt(2,0);
            crystal_positionZ[cryID] += rsctr_zshift[rs%nb_rsctr_axial_shift];

            crystal_orientationX[cryID] = cos(orientation_angle);
            crystal_orientationY[cryID] = sin(orientation_angle);
            crystal_orientationZ[cryID] = 0;
          }
    }

    // Update nb of crystal for systems with (required if nb_layers>1)
    nb_cry_cur += nb_crystals_lyr[lyr];
    
    // ============================================================================================================
    // Write the binary LUT file (append if nb_layers > 1)
    // ============================================================================================================

    cout << ">>> Start writing binary LUT file for layer #" << lyr << "..." << endl;
    // Write binary file, overwrite mode
    if(lyr == 0)
      LUT_file.open(path_to_LUT.c_str(), ios::binary | ios::out);
    // Append otherwise (data for crystals from the next layer)
    else
      LUT_file.open(path_to_LUT.c_str(), ios::binary | ios::out | ios::app);
    
    
    // Write the corresponding crystal parameters in the LUT
    for(int i=0 ; i<nb_crystals_lyr[lyr] ; i++)
    {
      LUT_file.write(reinterpret_cast<char*>(&crystal_positionX[i]), sizeof(FLTNBLUT));
      LUT_file.write(reinterpret_cast<char*>(&crystal_positionY[i]), sizeof(FLTNBLUT));
      LUT_file.write(reinterpret_cast<char*>(&crystal_positionZ[i]), sizeof(FLTNBLUT));
      LUT_file.write(reinterpret_cast<char*>(&crystal_orientationX[i]), sizeof(FLTNBLUT));
      LUT_file.write(reinterpret_cast<char*>(&crystal_orientationY[i]), sizeof(FLTNBLUT));
      LUT_file.write(reinterpret_cast<char*>(&crystal_orientationZ[i]), sizeof(FLTNBLUT));
    }

    LUT_file.close(); 
    cout << ">>> Binary LUT writing OK" << endl;
    
    // ============================================================================================================
    // Free memory
    // ============================================================================================================
  
    // Delete objects
  
    for (int i = 0; i < nb_rsectors ; i++)
     for (int j = 0; j<nb_axial_mod*nb_trans_mod; j++)
      for (int k = 0; k<nb_axial_submod*nb_trans_submod; k++)
       for (int l = 0; l<nb_axial_crystal*nb_trans_crystal; l++)
         delete crystal_center[i][j][k][l];
  
    for(int i = 0; i < nb_rsectors ; i++)
     for (int j = 0; j<nb_axial_mod*nb_trans_mod; j++)
      for (int k = 0; k<nb_axial_submod*nb_trans_submod; k++)
        delete[] crystal_center[i][j][k];

    for(int i = 0; i < nb_rsectors ; i++)
     for (int j = 0; j<nb_axial_mod*nb_trans_mod; j++)
       delete[] crystal_center[i][j];

    for(int i = 0; i < nb_rsectors ; i++)
    {
      delete[] crystal_center[i];
      delete rotation_mtx[i];
    }

    delete[] crystal_center;
    delete[] rotation_mtx;
    delete[] crystal_positionX;
    delete[] crystal_positionY;
    delete[] crystal_positionZ;
    delete[] crystal_orientationX;
    delete[] crystal_orientationY;
    delete[] crystal_orientationZ;
    
  } // end of loop on layers
  
  
  // ============================================================================================================
  // Write header file
  // ============================================================================================================

  cout << ">>> Start writing header LUT file..." << endl;
  header_LUT_file.open(path_to_headerLUT.c_str(), ios::out); 

  header_LUT_file << "scanner name:" << "    " << scanner_name << endl; 
  header_LUT_file << "modality:" << "    " << scanner_modality << endl; 
  
  header_LUT_file << "scanner radius:" << "    " << radius_lyr[0];
  for (int lyr=1 ; lyr<nbLayers ; lyr++) 
    header_LUT_file << "," << radius_lyr[lyr] ; header_LUT_file << endl;
    
  header_LUT_file << "number of rings in scanner:" << "    " << nb_rings << endl;
  header_LUT_file << "number of elements:" << "    " << nb_elts << endl; 
  header_LUT_file << "number of layers:" << "    " << nbLayers << endl;
  header_LUT_file << "number of crystals in layer(s):" << "    " << nb_crystals_lyr[0];
  for (int lyr=1 ; lyr<nbLayers ; lyr++) 
    header_LUT_file << ","<< nb_crystals_lyr[lyr] ; header_LUT_file << endl;
    
  header_LUT_file << "crystals size depth:" << "    " << crystal_size_depth_lyr[0];
  for (int lyr=1 ; lyr<nbLayers ; lyr++) 
    header_LUT_file << ","<< crystal_size_depth_lyr[lyr] ; header_LUT_file << endl;
    
  header_LUT_file << "crystals size transaxial:" << "    " << crystal_size_trans_lyr[0];  
  for (int lyr=1 ; lyr<nbLayers ; lyr++) 
    header_LUT_file << ","<< crystal_size_trans_lyr[lyr] ; header_LUT_file << endl;
    
  header_LUT_file << "crystals size axial:" << "    " << crystal_size_axial_lyr[0]; 
  for (int lyr=1 ; lyr<nbLayers ; lyr++) 
    header_LUT_file << ","<< crystal_size_axial_lyr[lyr] ; header_LUT_file << endl;
    
 
  //default reconstruction parameters
  header_LUT_file << "voxels number transaxial:" << "    " << default_dim_trans << endl;  
  header_LUT_file << "voxels number axial:" << "    " << default_dim_axial << endl; 

  header_LUT_file << "field of view transaxial:" << "    " << default_FOV_trans << endl;  
  header_LUT_file << "field of view axial:" << "    " << default_FOV_axial << endl; 
  
  header_LUT_file << "min angle difference:" << "    " << min_trs_angle_diff << " #deg" << endl;
  
  header_LUT_file << "mean depth of interaction:" << "    " << mean_depth_of_interaction_lyr[0];
  for (int lyr=1 ; lyr<nbLayers ; lyr++) 
    header_LUT_file << ","<< mean_depth_of_interaction_lyr[lyr] ;
  header_LUT_file << " #optional (default value : center of crystal ). Input value must correspond to the distance from the crystal surface, or negative value if default" << endl;
  
  header_LUT_file << "description:" << "    " << description << endl;
  
  cout << ">>> Header LUT file writing OK" << endl;
    
  // Free memory
  delete rsctr_zshift;
  delete nb_rsectors_lyr;
  delete nb_trans_mod_lyr;
  delete nb_axial_mod_lyr;
  delete nb_trans_submod_lyr;
  delete nb_axial_submod_lyr;
  delete nb_trans_crystal_lyr;
  delete nb_axial_crystal_lyr;

  delete radius_lyr;
  delete gap_trans_mod_lyr;
  delete gap_axial_mod_lyr;
  delete gap_trans_submod_lyr;
  delete gap_axial_submod_lyr;
  delete gap_trans_crystal_lyr;
  delete gap_axial_crystal_lyr;
  delete crystal_size_depth_lyr;  
  delete crystal_size_trans_lyr;
  delete crystal_size_axial_lyr;

  cout << "Binary file has been created in: " << path_to_LUT << endl;
  cout << "Header file has been created in: " << path_to_headerLUT << endl << endl;
  cout << "End of LUT generation" << endl << endl;
  
  return EXIT_SUCCESS;
}

