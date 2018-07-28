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
  \ingroup management
  \brief   This file gathers various function dedicated to data conversion in order to convert various type of GATE files (macro file, root datafile) to the CASToR file formats.
*/

#ifndef UTILS_HH
#define UTILS_HH 1

#include "gVariables.hh"
#include "gOptions.hh"
#include "iDataFilePET.hh"
#include "sOutputManager.hh"
#include "sScannerManager.hh"
#include "oInterfileIO.hh"

#ifdef CASTOR_ROOT
  #ifdef _WIN32
    #include "Windows4Root.h"
  #endif
  #include "TROOT.h"
  #include "TApplication.h"
  #include "TGClient.h"
  #include "TCanvas.h"
  #include "TSystem.h"
  #include "TTree.h"
  #include "TBranch.h"
  #include "TFile.h"
#endif

/**
 * @defgroup GATE_SYS_TYPE GATE system type
 *   \brief Variables related to the system type used in a GATE macro file \n
 *          Only cylindricalPET, ecat ans SPECThead are supported in this script \n
 *          Defined in gDataConversionUtilities.hh
 */
/*@{*/
/** Variable corresponding to an unknown system (=1, default) */
#define GATE_SYS_UNKNOWN -1
/** Variable corresponding to a cylindrical PET system (=0) */
#define GATE_SYS_CYLINDRICAL 0
/** Variable corresponding to an ecat system (=1) */
#define GATE_SYS_ECAT 1
/** Variable corresponding to a SPECThead system (=2) */
#define GATE_SYS_SPECT 2
/** @} */


#define GATE_NB_MAX_LAYERS 4 // Max number of layer in GATE = 4


/*!
  \fn CheckGATECommand(const string& a_key, const string& a_line)
  \param a_key: string containing a GATE command to recover
  \param a_line : string containing the line to check
  \brief Check if the line contains the provided GATE command. In this case,
         parse the line and returns the values in a vector of strings
  \details Values are converted in mm if 'cm' is found in the line 
  \return the vector of strings containing the elements of the line.
*/
vector<string> CheckGATECommand(const string& a_key, const string& a_line);

/*!
  \fn Split(string a_line)
  \param a_line : string to split
  \brief Split the line provided in parameter into a vector of strings (separator is blankspace)
  \return vector of string containing the splitted elements of the line
*/
vector<string> Split(string a_line);


/*!
  \fn ConvertValuesTomm(vector<string>& values)
  \param ap_v : vector of strings
  \brief Check if the vector of strings passed in parameter contains the 'cm' unit
         In this case, convert all its values in mm
*/
void ConvertValuesTomm(vector<string>& values);


/*!
  \fn toString(T a_val)
  \param a_val : value to convert
  \brief Convert a value of any type into string
  \return the value in parameter converted in string
*/
template <class T>
string toString(T a_val)
{
  stringstream ss;
  ss << a_val;
  return ss.str();
}


/*!
  \fn WriteVector(ofstream& file, const string& a_key, vector <T> a_vals)
  \param file : the output file
  \param a_key : key to write
  \param a_vals : vector containing the key values
  \brief Write the key and its values in the file provided in parameter
  \return 0 if success, positive value otherwise
*/
template <typename T>
int WriteVector(ofstream& file, const string& a_key, vector <T> a_vals);

/*!
  \fn WriteVector(ofstream& file, const string& a_key, vector <string> a_vals)
  \param file : the output file
  \param a_key : key to write
  \param a_vals : vector containing the key values
  \brief Write the key and its values in the file provided in parameter
  \return 0 if success, positive value otherwise
*/
int WriteVector(ofstream& file, const string& a_key, vector <string> a_vals);

/*!
  \fn WriteVector(ofstream& file, const string& a_key, vector <vector<string> > a_vals)
  \param file : the output file
  \param a_key : key to write
  \param a_vals : vector containing the key values
  \brief Write the key and its values contained in a 2 level vector of strings in the file provided in parameter
  \return 0 if success, positive value otherwise
*/
int WriteVector(ofstream& file, const string& a_key, vector <vector<string> > a_vals);



/*!
  \fn GetGATESystemType(const string& a_pathMac)
  \param a_pathMac : path to a GATE macro file
  \brief Read a GATE macro file and identify the system type from the 'gate/systems/' command lines
  \return system type, as described by the macro GATE_SYS_TYPE, -1 if error
*/
int GetGATESystemType(const string& a_pathMac);




/*!
  \fn GetGATEMacFiles(const string& a_pathMac, vector<string> &ap_pathToMacFiles)
  \param a_pathMac : path to a GATE main macro file
  \param ap_pathToMacFiles : array containing the paths of macro files
  \brief Extract the paths to each macro file contained in the main macro file.     
  \return 0 if success, positive value otherwise
*/
int GetGATEMacFiles(const string& a_pathMac, vector<string> &ap_pathToMacFiles);




/*!
  \fn GetGATEAliasesCylindrical(vector<string>  path_mac_files,
                                string&         rsector_name,
                                string&         module_name,
                                string&         submodule_name,
                                string&         crystal_name,
                                vector<string>& layers_name ,
                                int             vb );
  \param path_mac_files
  \param rsector_name
  \param module_name
  \param submodule_name
  \param crystal_name
  \param layers_name
  \brief Loop over a list of path to GATE macro files passed in parameter 
         to recover aliases of the different parts of a cylindricalPET
  \return 0 if success, positive value otherwise
*/
int GetGATEAliasesCylindrical(vector<string>  path_mac_files,
                              string&         rsector_name,
                              string&         module_name,
                              string&         submodule_name,
                              string&         crystal_name,
                              vector<string>& layers_name ,
                              int             vb );
                              
                              
                              

/*!
  \fn GetGATEAliasesEcat(vector<string> path_mac_files,
                         string&        block_name,
                         string&        crystal_name,
                         int            vb )
  \param path_mac_files
  \param block_name
  \param crystal_name
  \param vb
  \brief Loop over a list of path to GATE macro files passed in parameter 
         to recover aliases of the different parts of an ecat system
  \return 0 if success, positive value otherwise
*/
int GetGATEAliasesEcat(vector<string> path_mac_files,
                       string&        block_name,
                       string&        crystal_name,
                       int            vb );
                       
                       
                       

/*!
  \fn GetGATEAliasesSPECT(vector<string> path_mac_files,
                          string&        base_name,
                          string&        crystal_name,
                          string&        pixel_name,
                           int            vb )
  \param path_mac_files
  \param base_name
  \param crystal_name
  \param pixel_name
  \param vb
  \brief Loop over a list of path to GATE macro files passed in parameter 
         to recover aliases of the different parts of an ecat system
  \return 0 if success, positive value otherwise
*/
int GetGATEAliasesSPECT(vector<string> path_mac_files,
                        string&        base_name,
                        string&        crystal_name,
                        string&        pixel_name,
                        int            vb );
                        
                        
                        
                        
/*!
  \fn ConvertIDecat(int32_t nBlocksPerRing, 
                    int32_t nBlocksLine, 
                    int32_t nCrystalsTransaxial, 
                    int32_t nCrystalsAxial, 
                    int32_t crystalID, 
                    int32_t blockID)
  \param nBlocksPerRing
  \param nBlocksLine
  \param nCrystalsTransaxial
  \param nCrystalsAxial
  \param crystalID
  \param blockID
  \brief Compute a CASToR crystal index of a GATE ecat system from its indexes (block/crystal) and the system layout
  \return CASToR crystal index
*/
uint32_t ConvertIDecat(int32_t nBlocksPerRing, 
                       int32_t nBlocksLine, 
                       int32_t nCrystalsTransaxial, 
                       int32_t nCrystalsAxial, 
                       int32_t crystalID, 
                       int32_t blockID);





/*!
  \fn ConvertIDSPECTRoot1( int32_t  a_headID,
                           float_t  a_rotAngle,
                           float_t  a_angStep,
                           uint32_t a_nProjectionsByHead)
  \param a_headID : head index as provided in the root file
  \param a_rotAngle : rotation angle (deg) as provided in the root file
  \param a_angStep : angular step (deg) computed from the macro file
  \param a_nProjectionsByHead : total number of projections for each head
  \brief Compute a CASToR projection index of a GATE SPECThead system
  \return CASToR projection index for a SPECThead system
*/
uint32_t ConvertIDSPECTRoot1( int32_t  a_headID,
                              float_t  a_rotAngle,
                              float_t  a_angStep,
                              uint32_t a_nProjectionsByHead);
                              
                              
                              
/*!
  \fn ConvertIDSPECTRoot2( uint32_t a_nbSimulatedPixels,
                           uint32_t a_nPixTrs, 
                           uint32_t a_nPixAxl, 
                            int32_t a_headID,
                            int32_t a_crystalID,
                            int32_t a_pixelID,
                              float_t a_rotAngle,
                              float_t a_headAngPitch,
                              float_t a_crystalSizeAxl,
                              float_t a_crystalSizeTrs,
                              float_t a_gPosX,
                              float_t a_gPosY,
                              float_t a_gPosZ)
  \param a_nbSimulatedPixels
  \param a_nPixTrs
  \param a_nPixAxl
  \param a_headID
  \param a_crystalID
  \param a_pixelID
  \param a_rotAngle
  \param a_headAngPitch
  \param a_crystalSizeAxl
  \param a_crystalSizeTrs
  \param a_gPosX
  \param a_gPosY
  \param a_gPosZ
  \brief Compute a CASToR crystal index of a GATE SPECThead system
  \return CASToR crystal index. Return a higher number than the max index if error
*/
uint32_t ConvertIDSPECTRoot2( uint32_t a_nbSimulatedPixels,
                              uint32_t a_nPixTrs, 
                              uint32_t a_nPixAxl, 
                              int32_t a_headID,
                              int32_t a_crystalID,
                              int32_t a_pixelID,
                              float_t a_rotAngle,
                              float_t a_headAngPitch,
                              float_t a_crystalSizeAxl,
                              float_t a_crystalSizeTrs,
                              float_t a_gPosX,
                              float_t a_gPosY,
                              float_t a_gPosZ);
                             
                             
                             
                             
/*!
  \fn ConvertIDcylindrical(uint32_t  nRsectorsAngPos,
                           uint32_t  nRsectorsAxial,
                           int       a_rsectorIdOrder,
                           uint32_t  nModulesTransaxial,
                           uint32_t  nModulesAxial,
                           uint32_t  nSubmodulesTransaxial,
                           uint32_t  nSubmodulesAxial,
                           uint32_t  nCrystalsTransaxial,
                           uint32_t  nCrystalsAxial, 
                            uint8_t  nLayers,
                           uint32_t* nCrystalPerLayer, 
                    vector<uint32_t> nLayersRptTransaxial,
                    vector<uint32_t> nLayersRptAxial,
                            int32_t  layerID,
                            int32_t  crystalID,
                            int32_t  submoduleID,
                            int32_t  moduleID,
                            int32_t  rsectorID)
  \param nRsectorsAngPos
  \param nRsectorsAxial
  \param a_rsectorIdOrder
  \param nModulesTransaxial
  \param nModulesAxial
  \param nSubmodulesTransaxial
  \param nSubmodulesAxial
  \param nCrystalsTransaxial
  \param nCrystalsAxial
  \param nLayers
  \param nb_crystal_per_layer
  \param nLayersRptTransaxial
  \param nLayersRptAxial
  \param layerID
  \param crystalID
  \param submoduleID
  \param moduleID
  \param rsectorID
  \brief Compute a CASToR crystal index of a GATE cylindricalPET system from its indexes (rsector/module/submodule/crystal) and the system layout
  \return CASToR crystal index
*/
uint32_t ConvertIDcylindrical(uint32_t  nRsectorsAngPos,
                              uint32_t  nRsectorsAxial, 
                              int       a_rsectorIdOrder,
                              uint32_t  nModulesTransaxial,
                              uint32_t  nModulesAxial,
                              uint32_t  nSubmodulesTransaxial,
                              uint32_t  nSubmodulesAxial,
                              uint32_t  nCrystalsTransaxial, 
                              uint32_t  nCrystalsAxial, 
                              uint8_t   nLayers,
                              uint32_t* nCrystalPerLayer, 
                       vector<uint32_t> nLayersRptTransaxial,
                       vector<uint32_t> nLayersRptAxial,
                               int32_t  layerID,
                               int32_t  crystalID, 
                               int32_t  submoduleID, 
                               int32_t  moduleID, 
                               int32_t  rsectorID);

/*!
  \fn ComputeKindGATEEvent(uint32_t eventID1, uint32_t eventID2,
                           int comptonPhantom1, int comptonPhantom2, 
                           int rayleighPhantom1, int rayleighPhantom2)
  \param eventID1
  \param eventID2
  \param comptonPhantom1
  \param comptonPhantom2
  \param rayleighPhantom1
  \param rayleighPhantom2
  \brief Determine kind of a given coincidence event, from its attributes 
  \return kind of the coincidence : unknown (=0, default), true(=1), single scat(=2), multiple scat(=3), random(=4)) (default value =0)
*/
int ComputeKindGATEEvent(uint32_t eventID1, uint32_t eventID2,
                              int comptonPhantom1, int comptonPhantom2, 
                              int rayleighPhantom1, int rayleighPhantom2);



/*!
  \fn ReadMacECAT(string a_pathMac,
              uint32_t &nCrystalsTot,
              uint32_t &nCrystalsAxial,
              uint32_t &nCrystalsTransaxial,
              uint32_t &nBlocksLine,
              uint32_t &nBlocksPerRing,
              uint32_t &start_time_ms,
              uint32_t &duration_ms,
                   int vb);
  \param a_pathMac : path to a main macro file
  \param nCrystalsAxial : nb of axial crystals
  \param nCrystalsTransaxial : nb of transaxial crystals
  \param nBlocksLine : nb of axial blocks
  \param nBlocksPerRing : nb of blocks per ring
  \param start_time_ms : acquisition start time converted in ms
  \param duration_ms : acquisition duration converted in ms
  \param vb : verbosity
  \brief Recover informations about the scanner element of an ECAT system and acquisition duration, from a GATE macro file
  \return 0 if success, positive value otherwise
*/
int ReadMacECAT(string a_pathMac,
              uint32_t &nCrystalsTot,
              uint32_t &nCrystalsAxial,
              uint32_t &nCrystalsTransaxial,
              uint32_t &nBlocksLine,
              uint32_t &nBlocksPerRing,
              uint32_t &start_time_ms,
              uint32_t &duration_ms,
                   int vb);


/*!
  \fn ReadMacSPECT( string   a_pathMac,
                    float_t  &distToDetector,
                    uint32_t &nHeads,
                    uint32_t &nPixAxl,
                    uint32_t &nPixTrs,
                    float_t  &crystalSizeAxl,
                    float_t  &crystalSizeTrs,
                    uint32_t &nProjectionsTot,
                    uint32_t &nProjectionsByHead,
                    float_t  &head1stAngle,
                    float_t  &headAngPitch,
                    float_t  &headAngStepDeg,
                    int      &headRotDirection,
                    uint32_t &start_time_ms,
                    uint32_t &duration_ms,
                         int vb)
  \param a_pathMac : path to a macro file
  \param distToDetector : distance between center of rotation and detector surface
  \param nHeads : nb of cameras
  \param nPixAxl : nb of transaxial pixels
  \param nPixTrs : nb of axial pixels
  \param crystalSizeAxl : crystal axial dimension
  \param crystalSizeTrs : crystal transaxial dimension
  \param nProjectionsTot : total number of projections (cumulated over each head)
  \param nProjectionsByHead : total number of projections for each head
  \param head1stAngle : head first angle
  \param headAngPitch : angular pitch between heads
  \param headAngStepDeg : angular step between each projection
  \param headRotDirection : rotation direction of the head (CW or CCW)
  \param start_time_ms : acquisition start time converted in ms
  \param duration_ms : acquisition duration converted in ms
  \param vb : verbosity
  \brief Recover informations about the scanner element of an ECAT system, and acquisition duration, from a GATE macro file
  \return 0 if success, positive value otherwise
*/
int ReadMacSPECT( string   a_pathMac,
                  float_t  &distToDetector,
                  uint32_t &nHeads,
                  uint32_t &nPixAxl,
                  uint32_t &nPixTrs,
                  float_t  &crystalSizeAxl,
                  float_t  &crystalSizeTrs,
                  uint32_t &nProjectionsTot,
                  uint32_t &nProjectionsByHead,
                  float_t  &head1stAngle,
                  float_t  &headAngPitch,
                  float_t  &headAngStepDeg,
                  int      &headRotDirection,
                  uint32_t &start_time_ms,
                  uint32_t &duration_ms,
                       int vb);
                       
                       
                       


/*!
  \fn ReadIntfSPECT(string      a_pathIntf,
                    float_t     &ap_distToDetector,
                    uint32_t    &ap_nHeads,
                    uint32_t    &ap_nPixAxl,
                    uint32_t    &ap_nPixTrs,
                    float_t     &ap_crystalSizeAxl,
                    float_t     &ap_crystalSizeTrs,
                    uint32_t    &ap_nProjectionsTot,
                    uint32_t    &ap_nProjectionsByHead,
                    float_t     &ap_head1stAngle,
                    float_t     &ap_headAngPitch,
                    float_t     &headAngStepDeg,
                    int         &ap_headRotDirection,
                    uint32_t    &ap_start_time_ms,
                    uint32_t    &ap_duration_ms,
                    int         vb)
  \param a_pathIntf : path to the interfile header
  \param distToDetector : distance between center of rotation and detector surface
  \param nHeads : nb of cameras
  \param nPixAxl : nb of transaxial pixels
  \param nPixTrs : nb of axial pixels
  \param crystalSizeAxl : crystal axial dimension
  \param crystalSizeTrs : crystal transaxial dimension
  \param head1stAngle : head first angle
  \param headAngPitch : angular pitch between heads
  \param headRotSpeed : angle between projection
  \param headRotDirection : rotation direction of the head (CW or CCW)
  \param start_time_ms : acquisition start time converted in ms
  \param duration_ms : acquisition duration converted in ms
  \param Intf_fiels : Interfile key structure to initialize with projection image parameters
  \param vb : verbosity
  \brief Recover informations about the scanner element of an ECAT system, and acquisition duration, from a GATE macro file
  \return 0 if success, positive value otherwise
*/
int ReadIntfSPECT(string      a_pathIntf,
                  float_t     &ap_distToDetector,
                  uint32_t    &ap_nHeads,
                  uint32_t    &ap_nPixAxl,
                  uint32_t    &ap_nPixTrs,
                  float_t     &ap_crystalSizeAxl,
                  float_t     &ap_crystalSizeTrs,
                  uint32_t    &ap_nProjectionsTot,
                  uint32_t    &ap_nProjectionsByHead,
                  float_t     &ap_head1stAngle,
                  float_t     &ap_headAngPitch,
                  float_t     &headAngStepDeg,
                  int         &ap_headRotDirection,
                  uint32_t    &ap_start_time_ms,
                  uint32_t    &ap_duration_ms,
                  int         vb);
                  
                  
                  
                   
/*!
  \fn ReadMacCylindrical( string a_pathMac,
                      uint8_t  &nLayers,
                      uint32_t *nb_crystal_per_layer,
                      uint32_t &nCrystalsTot,
                      uint32_t &nCrystalsAxial,
                      uint32_t &nCrystalsTransaxial, 
              vector<uint32_t> &nLayersRptAxial,
              vector<uint32_t> &nLayersRptTransaxial,
                      uint32_t &nSubmodulesAxial, 
                      uint32_t &nSubmodulesTransaxial,
                      uint32_t &nModulesAxial, 
                      uint32_t &nModulesTransaxial,
                      uint32_t &nRsectorsAxial,
                      uint32_t &nRsectorsAngPos, 
                      int      &rsector_id_order,
                      uint32_t &start_time_ms, 
                      uint32_t &duration_ms,
                           int vb);
  \param a_pathMac : path to a main macro file
  \param nLayers : nb of crystal layers
  \param nCrystalsAxial : nb of axial crystals (in a submodule)
  \param nCrystalsTransaxial : nb of transaxial crystals (in a submodule)
  \param nCrystalLayersAxial : Array containing the number of axial crystals in each layer
  \param nCrystalLayersTransaxial : Array containing the number of transaxial crystals in each layer
  \param nSubmodulesAxial : nb of axial submodules (in a module)
  \param nSubmodulesTransaxial : nb of transaxial submodules (in a module)
  \param nModulesAxial : nb of axial modules (in a rsector)
  \param nModulesTransaxial : nb of transaxial modules (in a rsector)
  \param nRsectorsAxial : nb of axial rsectors
  \param nRsectorsAngPos : nb of rsectors per ring
  \param rsector_id_order : ordering of rsector id
  \param start_time_ms : acquisition start time converted in ms
  \param duration_ms : acquisition duration converted in ms
  \param vb : verbosity
  \brief Recover informations about the scanner element of a cylindricalPET system and acquisition duration, from a GATE macro file
  \return 0 if success, positive value otherwise
*/
int ReadMacCylindrical( string a_pathMac,
                      uint8_t  &nLayers,
                      uint32_t *nb_crystal_per_layer,
                      uint32_t &nCrystalsTot,
                      uint32_t &nCrystalsAxial,
                      uint32_t &nCrystalsTransaxial, 
              vector<uint32_t> &nLayersRptAxial,
              vector<uint32_t> &nLayersRptTransaxial,
                      uint32_t &nSubmodulesAxial, 
                      uint32_t &nSubmodulesTransaxial,
                      uint32_t &nModulesAxial, 
                      uint32_t &nModulesTransaxial,
                      uint32_t &nRsectorsAxial,
                      uint32_t &nRsectorsAngPos, 
                      int      &rsector_id_order,
                      uint32_t &start_time_ms, 
                      uint32_t &duration_ms,
                           int vb);



/*!
  \fn      CreateGeomWithECAT(string a_pathMac, string a_pathGeom)
  \param   a_pathMac : string containing the path to a GATE macro file
  \param   a_pathGeom : string containing the path to a CASToR output geom file
  \brief   Read a GATE macro file containing the description of an ecat system, and convert it to a geom file
  \return  0 if success, positive value otherwise
*/
int CreateGeomWithECAT(string a_pathMac, string a_pathGeom);




/*!
  \fn      CreateGeomWithCylindrical(string a_pathMac, string a_pathGeom)
  \param   a_pathMac : string containing the path to a GATE macro file
  \param   a_pathGeom : string containing the path to a CASToR output geom file
  \brief   Read a GATE macro file containing the description of a cylindricalPET system, and convert it to a geom file
  \return  0 if success, positive value otherwise
*/
int CreateGeomWithCylindrical(string a_pathMac, string a_pathGeom);




/*!
  \fn      CreateGeomWithSPECT(string a_pathMac, string a_pathGeom)
  \param   a_pathMac : string containing the path to a GATE macro file
  \param   a_pathGeom : string containing the path to a CASToR output geom file
  \brief   Read a GATE macro file containing the description of a SPECThead system, and convert it to a geom file
  \return  0 if success, positive value otherwise
*/
int CreateGeomWithSPECT(string a_pathMac, string a_pathGeom);

#endif
