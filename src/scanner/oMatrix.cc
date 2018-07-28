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
  \brief    Implementation of class oMatrix
*/

#include "oMatrix.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \brief oMatrix constructor. 
         Initialize the member variables to their default values.
*/
oMatrix::oMatrix() 
{
  m_lin = 0;
  m_col = 0;

  m2p_Mat = NULL;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \param nl : a number of lines
  \param nc : a number of colons
  \brief oMatrix constructor. 
         Instanciate a Matrix structure with the number of lines and colons provided in parameters
*/
oMatrix::oMatrix(uint16_t nl, uint16_t nc) 
{
  m_lin = nl;
  m_col = nc;
  m2p_Mat = new FLTNB *[nl];

  for(uint16_t l=0 ; l<nl ; l++)
    m2p_Mat[l] = new FLTNB[nc];
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \brief oMatrix destructor. 
         Free memory of the oMatrix object.
*/
oMatrix::~oMatrix() 
{
  for(uint16_t l=0 ; l<m_lin ; l++)
    if(m2p_Mat[l]) delete[] m2p_Mat[l];

  if(m2p_Mat) delete[] m2p_Mat;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn Allocate
  \param nl : a number of lines
  \param nc : a number of colons
  \brief Instanciate a Matrix structure with the number of lines and colons provided in parameters
*/
void oMatrix::Allocate(uint16_t nl, uint16_t nc) 
{
  // No verbosity in oMatrix structure. Show this only if CASTOR_DEBUG is enabled (for bug tracking)
  #ifdef CASTOR_DEBUG
  Cout("oMatrix::Allocate() ...");
  #endif
  
  // Free memory in case the matrix had already been allocated
  if(m2p_Mat != NULL)
  {
    for(uint16_t l=0 ; l<m_lin ; l++)
      if(m2p_Mat[l]) delete[] m2p_Mat[l];

    delete[] m2p_Mat;
  }
  
  m_lin = nl;
  m_col = nc;
  m2p_Mat = new FLTNB *[nl];

  for(uint16_t l=0 ; l<nl ; l++)
    m2p_Mat[l] = new FLTNB[nc];
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn SetMatriceElt
  \param l : a line index
  \param c : a colon index
  \param a_val : a value to initialize the matrix element with
  \brief set the matrix element corresponding to the argument indices with the provided value.
  \return 0 if success, positive value otherwise
*/
int oMatrix::SetMatriceElt(uint16_t l, uint16_t c, FLTNB a_val) 
{
  if(l>=m_lin || c>=m_col)
  {
    Cerr("***** oMatrix::SetMatriceElt()-> Nb of (lin,col) ("<<l+1<<","<<c+1<<") in parameters ");
    Cerr("> to the number of (lin,col) of this matrix ("<<m_lin<<","<<m_col<<") !" << endl);
    return 1;
  }
  
  m2p_Mat[l][c] = a_val;
  
  return 0;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn GetMatriceElt
  \param l : a line index
  \param c : a colon index
  \return the matrix element value corresponding to the argument indices provided in arguments.
  \todo return error if unmatching number of lines/colons with the matrix initialization ?
*/
FLTNB oMatrix::GetMatriceElt(uint16_t l,uint16_t c) 
{
  return m2p_Mat[l][c];
}
    
    
    
// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn Multiplication
  \param ap_Mtx : input matrix to multiply to the matrix from which the function is called
  \param ap_MtxResult : output matrix containing the result of the multiplication
  \brief Multiply the member matrix with the matrix provided in 1st parameter
         Return the result in the matric provided in 2nd parameter
  \return 0 if success, positive value otherwise
*/
int oMatrix::Multiplication(oMatrix *ap_Mtx, oMatrix *ap_MtxResult)
{  
  HPFLTNB sum=0;

  if ( m_col != ap_Mtx->m_lin ) 
  {
    Cerr("***** oMatrix::Multiplication()-> Not matching number of colons and lines of the two matrices !");
    return 1;
  }
  else if (ap_MtxResult == NULL)
  {
    Cerr("***** oMatrix::Multiplication()-> The resulting matrix has not been allocated !");
    return 1;
  }
  else
  {
    for ( uint16_t tl = 0; tl < m_lin; tl++ ) 
      for ( uint16_t c = 0; c < ap_Mtx->m_col; c++ ) 
      {
        for ( uint16_t l = 0; l < ap_Mtx->m_lin; l++ ) 
        {
          sum += GetMatriceElt(tl,l) * ap_Mtx->GetMatriceElt(l,c);
        }
        ap_MtxResult->SetMatriceElt(tl,c,sum);
        sum = 0;
      }
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      int oMatrix::Transpose()
  \param   ap_MtxResult : output matrix containing the result of the transposition
  \brief   Transpose the elements of the matrix
  \return  0 if success, positive value otherwise
*/
int oMatrix::Transpose(oMatrix *ap_Mtx)
{  
  if ( this->m_col != ap_Mtx->m_lin
    || this->m_lin != ap_Mtx->m_col ) 
  {
    Cerr("***** oMatrix::Transpose()-> Not matching number of colons <-> lines between the input and transpose matrices !");
    return 1;
  }
  if (ap_Mtx == NULL)
  {
    Cerr("***** oMatrix::Transpose()-> The resulting matrix has not been allocated !");
    return 1;
  }

  uint16_t lsize = m_lin>m_col ? m_lin : m_col;
     
  for ( uint16_t l = 0; l < lsize; l++ ) 
    for ( uint16_t c = l; c < lsize; c++ ) 
    {
      HPFLTNB tmp = 0.;
      
      if( l < m_lin
       && c < m_col )
        tmp = m2p_Mat[l][c];
        
      if(l < ap_Mtx->m_lin 
      && c < ap_Mtx->m_col
      && l < m_col
      && c < m_lin )
        ap_Mtx->SetMatriceElt( l , c , m2p_Mat[c][l] );
        
      if(l < ap_Mtx->m_col 
      && c < ap_Mtx->m_lin
      && l < m_lin
      && c < m_col )
        ap_Mtx->SetMatriceElt( c , l , tmp );
    }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      int oMatrix::Inverse()
  \param   ap_MtxResult : output matrix containing the inverse matrix
  \brief   Inverse the matrix on which the function is called
           An error is returned if the matrix is not square,
           or if there is a nil component on the diagonal 
           The resulting matrix object must be different than the input matrix
  \return  0 if success, positive value otherwise
*/
int oMatrix::Inverse(oMatrix *ap_Mtx)
{
  // Check matching nb colons/lines
  if ( (this->m_col   != this->m_lin)
    || (ap_Mtx->m_col != ap_Mtx->m_lin) ) 
  {
    Cerr("***** oMatrix::Inverse()-> Not matching number of colons <-> lines between the input and transpose matrices !");
    return 1;
  }
  // Check if resulting matrix has been allocated
  if (ap_Mtx == NULL)
  {
    Cerr("***** oMatrix::Inverse()-> The resulting matrix has not been allocated !");
    return 1;
  }
  // Check if resulting matrix is the same as this matrix
  // Throw error as the resulting matrix will overwrite the components of the input matrix in the process
  if (this == ap_Mtx)
  {
    Cerr("***** oMatrix::Inverse()-> The resulting matrix has not been allocated !");
    return 1;
  }
  
  for ( uint16_t l = 0; l < m_lin; l++ ) 
    for ( uint16_t c = 0; c < m_col; c++ ) 
    {
      if (m2p_Mat[l][c] == 0)
      {
        Cerr("***** oMatrix::Inverse()-> One element on the diagonal is nil !");
        return 1;
      }
    }

  oMatrix* p_tMtx = new oMatrix(m_lin, m_col);
  
  
  // Init
  for( uint16_t l=0 ; l<m_lin ; l++ )
    for ( uint16_t c = 0; c < m_col; c++ )
    {
      p_tMtx->SetMatriceElt( l, c, 0. );
      ap_Mtx->SetMatriceElt( l, c, m2p_Mat[l][c] );
    }
      
  FLTNB** pMtx = ap_Mtx->GetMtx();
      
  for( uint16_t l=0 ; l<m_lin ; l++ )
  {
    for ( uint16_t c = 0; c < m_col; c++ ) 
    {
      p_tMtx->SetMatriceElt(l, l, 1./pMtx[l][l]);
      
      if(c != l)
        p_tMtx->SetMatriceElt(l, c, -pMtx[l][c] / pMtx[l][l]);
      
      for( uint16_t r=0 ; r<m_lin ; r++ )
      {
        if(r!=l)
          p_tMtx->SetMatriceElt(r, l, pMtx[r][l]/pMtx[l][l]);
        if(c!=l && r!=l)
          p_tMtx->SetMatriceElt(r, c, pMtx[r][c] - pMtx[l][c]*pMtx[r][l]/pMtx[l][l]);
      }
      
    }
  
    for(int ll=0 ; ll<m_lin ; ll++)
      for(int cc=0 ; cc<m_col ; cc++)
        pMtx[ll][cc] = p_tMtx->GetMatriceElt(ll,cc) ;

  }
  
  delete p_tMtx;  
    
  return 0;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      int oMatrix::SetXRotMtx()
  \param   ang : angle in degree
  \brief   Set a (3,3) X-axis rotation matrix using the provided angle
  \return 0 if sucess, positive value otherwise
*/
int oMatrix::SetXRotMtx(FLTNB ang)
{
  if(m_lin !=3 || m_col != 3)
  {
    Cerr("***** oMatrix::SetXRotMtx()-> The matrix must be (3,3) in order to use this function !");
    return 1;
  }
  
  SetMatriceElt(0,0, 1 );
  SetMatriceElt(0,1, 0 );
  SetMatriceElt(0,2, 0 );
  SetMatriceElt(1,0, 0 );
  SetMatriceElt(1,1, cos(ang));
  SetMatriceElt(1,2, -sin(ang));
  SetMatriceElt(2,0, 0 );
  SetMatriceElt(2,1, sin(ang));
  SetMatriceElt(2,2, cos(ang)) ;
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      int oMatrix::SetYRotMtx()
  \param   ang : angle in degree
  \brief   Set a (3,3) Y-axis rotation matrix using the provided angle
  \return 0 if sucess, positive value otherwise
*/
int oMatrix::SetYRotMtx(FLTNB ang)
{
  if(m_lin !=3 || m_col != 3)
  {
    Cerr("***** oMatrix::SetYRotMtx()-> The matrix must be (3,3) in order to use this function !");
    return 1;
  }
  
  SetMatriceElt(0,0, cos(ang));
  SetMatriceElt(0,1, 0 );
  SetMatriceElt(0,2, sin(ang));
  SetMatriceElt(1,0, 0 );
  SetMatriceElt(1,1, 1 );
  SetMatriceElt(1,2, 0 );
  SetMatriceElt(2,0, -sin(ang));
  SetMatriceElt(2,1, 0 );
  SetMatriceElt(2,2, cos(ang));
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      int oMatrix::SetZRotMtx()
  \param   ang : angle in degree
  \brief   Set a (3,3) Z-axis rotation matrix using the provided angle
  \return 0 if sucess, positive value otherwise
*/
int oMatrix::SetZRotMtx(FLTNB ang)
{
  if(m_lin !=3 || m_col != 3)
  {
    Cerr("***** oMatrix::SetZRotMtx()-> The matrix must be (3,3) in order to use this function !");
    return 1;
  }
  
  SetMatriceElt(0,0, cos(ang) );
  SetMatriceElt(0,1, -sin(ang) );
  SetMatriceElt(0,2, 0);
  SetMatriceElt(1,0, sin(ang));
  SetMatriceElt(1,1, cos(ang) );
  SetMatriceElt(1,2, 0);
  SetMatriceElt(2,0, 0);
  SetMatriceElt(2,1, 0);
  SetMatriceElt(2,2, 1);
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      int oMatrix::Describe()
  \brief   Display the element of the matrix
*/
void oMatrix::Describe()
{  
  for ( uint16_t l = 0; l < m_lin; l++ ) 
  {
    for ( uint16_t c = 0; c < m_col; c++ ) 
      cout << m2p_Mat[l][c] << " ; ";

    cout << endl;
  }
}
