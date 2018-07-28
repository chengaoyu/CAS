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
  \brief    Implementation of class oSystemMatrix
*/

#include "oSystemMatrix.hh"
#include "vScanner.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oSystemMatrix::oSystemMatrix()
{
  // Not yet implemented
  Cerr("***** oSystemMatrix::oSystemMatrix() -> Not yet implemented !" << endl);
  Exit(EXIT_FAILURE);
  // Default all members
  m_compatibleWithSPECTAttenuationCorrection = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oSystemMatrix::~oSystemMatrix() {}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oSystemMatrix::Project(int a_projector, oProjectionLine* ap_ProjectionLine, uint32_t* ap_index1, uint32_t* ap_index2, int a_nbIndices)
{
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
