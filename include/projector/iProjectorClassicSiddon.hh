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
  \brief    Declaration of class iProjectorClassicSiddon
*/

#ifndef IPROJECTORCLASSICSIDDON_HH
#define IPROJECTORCLASSICSIDDON_HH 1

#include "gVariables.hh"
#include "sAddonManager.hh"
#include "vProjector.hh"

/*!
  \class   iProjectorClassicSiddon
  \brief   This class is a child of the vProjector class implementing the classic Siddon ray tracer
  \details This class implements the classic Siddon algorithm which is a ray-tracer algorithm calculating
           the exact length of a line through each voxel it crosses. Reference: R. L. Siddon, "Fast calculation
           of the exact radiological path for a three-dimensional CT array", Med. Phys., vol. 12, pp. 252-5, 1985.
*/
class iProjectorClassicSiddon : public vProjector
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public iProjectorClassicSiddon::iProjectorClassicSiddon()
      \brief   The constructor of iProjectorClassicSiddon
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    iProjectorClassicSiddon();
    /*!
      \fn      public iProjectorClassicSiddon::~iProjectorClassicSiddon()
      \brief   The destructor of iProjectorClassicSiddon
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built by this class.
    */
    ~iProjectorClassicSiddon();


  // -------------------------------------------------------------------
  // Public member functions
  public:
    // Function for automatic insertion (put the class name as the parameter and do not add semi-column at the end of the line)
    FUNCTION_PROJECTOR( iProjectorClassicSiddon )
    /*!
      \fn      public int iProjectorClassicSiddon::ReadConfigurationFile()
      \param   const string& a_configurationFile
      \brief   A function used to read options from a configuration file
      \details This function implements the reading of all options associated to the child projector, from
               a configuration file. It is the implementation of the pure virtual function inherited
               from the abstract class vProjector. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadConfigurationFile(const string& a_configurationFile);
    /*!
      \fn      public int iProjectorClassicSiddon::ReadOptionsList()
      \param   const string& a_configurationFile
      \brief   A function used to read options from a list of options
      \details This function implements the reading of all options associated to the child projector, from
               a list of options. It is the implementation of the pure virtual function inherited
               from the abstract class vProjector. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadOptionsList(const string& a_optionsList);
    /*!
      \fn      public INTNB iProjectorNativeSiddon::EstimateMaxNumberOfVoxelsPerLine()
      \brief   This function is used to compute and provide an estimate of the maximum number of voxels that could
               contribute to a projected line.
      \details This function is an overloaded implementation of the virtual mother function. It is used to compute
               and provide an estimate of the maximum number of voxels that could contribute to a projected line.
      \return  The estimate of the maximum number of voxels contributing to a line.
    */
    INTNB EstimateMaxNumberOfVoxelsPerLine();


  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      private void iProjectorClassicSiddon::ShowHelpSpecific()
      \brief   A function used to show help about the child projector
      \details This function must describe what the module does and how to use it. It describes in
               details the different parameters of the projector, and how to set them through the use
               of a configuration file or a list of options. It is the implementation of the pure
               virtual function inherited from the abstract class vProjector.
    */
    void ShowHelpSpecific();
    /*!
      \fn      private int iProjectorClassicSiddon::CheckSpecificParameters()
      \brief   A private function used to check the parameters settings specific to the child projector
      \details This function is used to check that all parameters specific to the projector are correctly set
               within allowed values. It is called by the CheckParameters() function of the mother class.
               It is the implementation of the pure virtual function inherited from the abstract mother
               class vProjector.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      private int iProjectorClassicSiddon::InitializeSpecific()
      \brief   This function is used to initialize specific stuff to the child projector.
      \details It is called by the public Initialize() function from the mother.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */
    int InitializeSpecific();
    /*!
      \fn      private int iProjectorClassicSiddon::ProjectWithoutTOF()
      \param   int a_direction
      \param   oProjectionLine* ap_ProjectionLine
      \brief   A function to project without TOF.
      \details Projects the provided line following the provided direction, without TOF. It fills the provided
               oProjectionLine. It is an implementation of the pure virtual function from the mother class.
      \return  An integer reflecting the projection status; 0 if no problem, another value otherwise.
    */
    int ProjectWithoutTOF( int a_direction, oProjectionLine* ap_ProjectionLine );
    /*!
      \fn      private int iProjectorClassicSiddon::ProjectWithTOFPos()
      \param   int a_direction
      \param   oProjectionLine* ap_ProjectionLine
      \brief   A function to project with TOF continuous information.
      \details Projects the provided line following the provided direction, with TOF described as a continuous
               measurement. It fills the provided oProjectionLine. It is an implementation of the pure virtual
               function from the mother class.
      \return  An integer reflecting the projection status; 0 if no problem, another value otherwise.
    */
    int ProjectWithTOFPos( int a_direction, oProjectionLine* ap_ProjectionLine );
    /*!
      \fn      private int iProjectorClassicSiddon::ProjectWithTOFBin()
      \param   int a_direction
      \param   oProjectionLine* ap_ProjectionLine
      \brief   A function to project with TOF binned information.
      \details Projects the provided line following the provided direction, with TOF information describe as a
               histogram bin. It fills the provided oProjectionLine. It is an implementation of the pure virtual
               function from the mother class.
      \return  An integer reflecting the projection status; 0 if no problem, another value otherwise.
    */
    int ProjectWithTOFBin( int a_direction, oProjectionLine* ap_ProjectionLine );
};


// Class for automatic insertion (set here the visible projector's name as the first parameter,
// put the class name as the second parameter and do NOT add semi-colon at the end of the line)
CLASS_PROJECTOR( classicSiddon, iProjectorClassicSiddon )

#endif

