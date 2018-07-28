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
  \ingroup  scanner
  \brief    Declaration of class oMatrix
*/

#ifndef OMATRIX_HH
#define OMATRIX_HH 1

#include "gVariables.hh"
#include "sOutputManager.hh"

/*!
  \class   oMatrix
  \brief   Structure designed for basic matrices operations.
  \details This structure is mostly used in scanner classes, for geometrical rotation purposes
*/
class oMatrix
{
  // Constructor & Destructor
  public:
    /*!
      \brief   oMatrix constructor. 
               Initialize the member variables to their default values.
    */
    oMatrix();
    /*!
      \brief   oMatrix destructor. 
               Free memory of the oMatrix object.
    */
    ~oMatrix();
    /*!
      \param   nl : a number of lines
      \param   nc : a number of colons
      \brief   oMatrix constructor. 
               Instanciate a matrix structure with the number of lines and colons provided in parameters
    */
    oMatrix(uint16_t nl, uint16_t nc);


  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      void oMatrix::Allocate()
      \param   nl: a number of lines
      \param   nc: a number of colons
      \brief   Instanciate a Matrix structure with the number of lines and colons provided in parameters
    */
    void Allocate(uint16_t nl, uint16_t nc);
    /*!
      \fn      int oMatrix::SetMatriceElt()
      \param   l: a line index
      \param   c: a colon index
      \param   a_val : a value to initialize the matrix element with
      \brief   Set the matrix element corresponding to the argument indices with the provided value.
      \return  0 if success, positive value otherwise
    */
    int SetMatriceElt(uint16_t l, uint16_t c, FLTNB a_val);
    /*!
      \fn      FLTNB oMatrix::GetMatriceElt()
      \param   l: a line index
      \param   c: a colon index
      \return  the matrix element value corresponding to the argument indices provided in arguments.
      \todo    return error if unmatching number of lines/colons with the matrix initialization ?
    */
    FLTNB GetMatriceElt(uint16_t l,uint16_t c);
    /*!
      \fn      FLTNB** oMatrix::GetMtx()
      \return  pointer to the matrix array.
    */
    FLTNB** GetMtx() {return m2p_Mat; }
    /*!
      \fn      int oMatrix::Multiplication()
      \param   ap_Mtx : input matrix to multiply to the matrix from which the function is called
      \param   ap_MtxResult : output matrix containing the result of the multiplication
      \brief   Multiply the member matrix with the matrix provided in 1st parameter
               Return the result in the matric provided in 2nd parameter
      \return  0 if success, positive value otherwise
    */
    int Multiplication(oMatrix *ap_Mtx, oMatrix *ap_MtxResult);
    /*!
      \fn      int oMatrix::Transpose()
      \param   ap_MtxResult : output matrix containing the result of the transposition
      \brief   Transpose the elements of the matrix
      \return  0 if success, positive value otherwise
    */
    int Transpose(oMatrix *a_MtxResult);
    /*!
      \fn      int oMatrix::Inverse()
      \param   ap_MtxResult : output matrix containing the inverse matrix
      \brief   Inverse the matrix on which the function is called
               An error is returned if the matrix is not square,
               or if there is a nil component on the diagonal 
               The resulting matrix object must be different than the input matrix
      \return  0 if success, positive value otherwise
    */
    int Inverse(oMatrix *ap_MtxResult);
    /*!
      \fn      int oMatrix::SetXRotMtx()
      \param   ang : angle in degree
      \brief   Set a (3,3) X-axis rotation matrix using the provided angle
      \return 0 if sucess, positive value otherwise
    */
    int SetXRotMtx(FLTNB ang);
    
    /*!
      \fn      int oMatrix::SetYRotMtx()
      \param   ang : angle in degree
      \brief   Set a (3,3) Y-axis rotation matrix using the provided angle
      \return 0 if sucess, positive value otherwise
    */
    int SetYRotMtx(FLTNB ang);
    
    /*!
      \fn      int oMatrix::SetZRotMtx()
      \param   ang : angle in degree
      \brief   Set a (3,3) Z-axis rotation matrix using the provided angle
      \return 0 if sucess, positive value otherwise
    */
    int SetZRotMtx(FLTNB ang);


    /*!
      \fn      int oMatrix::Describe()
      \brief   Display the element of the matrix
    */
    void Describe();

  // -------------------------------------------------------------------
  // Data members. (Nb lines/colons set public in order to optimize Multiplication function)
  public:
    uint16_t m_lin;     /*!< Number of lines in the matrix. Default =0 */
    uint16_t m_col;     /*!< Number of colons in the matrix. Default =0 */
  private:
    FLTNB **m2p_Mat; /*!< 2D pointer containing the matrix elements. Default =NULL */
};

#endif
