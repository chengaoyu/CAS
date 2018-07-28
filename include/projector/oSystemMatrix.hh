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
  \brief    Declaration of class oSystemMatrix
*/

#ifndef OSYSTEMMATRIX_HH
#define OSYSTEMMATRIX_HH 1

#include "gVariables.hh"
#include "oProjectionLine.hh"
#include "vScanner.hh"

/**
 * \ingroup projector
 * @defgroup SYSTEM_MATRIX_KEYWORD System matrix keyword
 *
 *    \brief Keyword for a pre-computed and loaded system matrix \n
 *           Defined in oSystemMatrix.hh
 */
 /*@{*/
/** String constant corresponding to a key word used as a projector
    name for pre-computed and loaded system matrix instead of using
    an on-the-fly projector */
#define SYSTEM_MATRIX_KEYWORD "matrix" 
/** @} */


/*!
  \class   oSystemMatrix
  \brief   This class is designed to manage pre-computed system matrices
  \details This class is basically a container for pre-computed system matrices
           able to read/write oProjectionLines. It can be used during the
           reconstruction process by the oProjectorManager.
           Everything needs to be implemented !
*/
class oSystemMatrix
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public oSystemMatrix::oSystemMatrix()
      \brief   The constructor of oSystemMatrix
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    oSystemMatrix();
    /*!
      \fn      virtual public oSystemMatrix::~oSystemMatrix()
      \brief   The destructor of oSystemMatrix
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built by this class.
               It is virtual, so that it is automatically called when a child object is deleted.
    */
    ~oSystemMatrix();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      public int vProjector::Project()
      \param   int a_direction
      \param   oProjectionLine* ap_ProjectionLine
      \param   int* ap_index1
      \param   int* ap_index2
      \param   int a_nbIndices
      \brief   A function use to computed the projection elements with respect to the provided parameters
      \details This function is used to fill the provided oProjectionLine with the system matrix elements associated
               to the provided indices. Probably do this by copying the pointers of voxel indices and weights tabs
               stored in maps into the oProjectionLine.
      \return  An integer reflecting the projection status; 0 if no problem, another value otherwise.
    */
    int Project(int a_direction, oProjectionLine* ap_ProjectionLine, uint32_t* ap_index1, uint32_t* ap_index2, int a_nbIndices);


  // -------------------------------------------------------------------
  // Private member functions
  private:


  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      inline public bool oSystemMatrix::GetCompatibilityWithSPECTAttenuationCorrection()
      \return  m_compatibleWithSPECTAttenuationCorrection
    */
    inline bool GetCompatibilityWithSPECTAttenuationCorrection()
           {return m_compatibleWithSPECTAttenuationCorrection;}


  // -------------------------------------------------------------------
  // Data members
  private:
    map<pair<int,int>,INTNB>  mp_nbVoxels;            /*!< The map associating pairs of indices to the number of contributing voxels */
    map<pair<int,int>,INTNB*> m2p_voxelIndicesMap;    /*!< The map associating pairs of indices to the list of voxels indices */
    map<pair<int,int>,FLTNB*> m2p_voxelWeightsMap;    /*!< The map associating pairs of indices to the list of voxels weights*/
//    vScanner* mp_Scanner;                             /*!< The pointer to the vScanner in use */
    bool m_compatibleWithSPECTAttenuationCorrection;  /*!< Boolean that says if the projector is compatible with SPECT attenuation correction */
};

#endif
