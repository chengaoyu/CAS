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
  \ingroup datafile

  \brief Implementation of class iEventListPET
*/

#include "iEventListPET.hh"
#include "vDataFile.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventListPET::iEventListPET() : iEventPET()
{
  m_dataType = TYPE_PET;
  m_dataSpec = SPEC_EMISSION;
  m_dataMode = MODE_LIST;
  mp_POI1[0] = 0.;
  mp_POI1[1] = 0.;
  mp_POI1[2] = -1.;
  mp_POI2[0] = 0.;
  mp_POI2[1] = 0.;
  mp_POI2[2] = -1.;
  m_kind = KIND_UNKNOWN;
  m_eventScatRate = 0.;
  m_TOFMeasurement = 0.;
  m_eventValue = 1.;
  m_nbLines = 1;
  m_spatialTOFMeasurementRange = 0.;
  m_hasTOFinfo = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iEventListPET::~iEventListPET() {}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iEventListPET::SetEventValue(int a_bin, FLTNBDATA a_value) 
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  Cout("***** iEventListPET::SetEventValue() -> Trying to set the value of a list mode event !");
  Exit(EXIT_FAILURE);
} 

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iEventListPET::Describe()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  Cout("iEventListPET::Describe() -> Display contents" << endl);
  Cout("Time: " << m_timeInMs << " ms" << endl);
  Cout("Number of lines: " << m_nbLines << endl);
  for (uint16_t l=0; l<m_nbLines; l++) Cout("  --> ID1: " << mp_ID1[l] << " | ID2: " << mp_ID2[l] << endl);
  Cout("Random rate: " << m_eventRdmRate << endl);
  Cout("Normalization factor: " << m_eventNormFactor << endl);
  Cout("ACF: " << m_atnCorrFactor << endl);
//  Cout("kind: " << m_kind << endl);
  Cout("Scatter rate: " << m_eventScatRate << endl);
  Cout("TOF measurement: " << m_TOFMeasurement << " ps" << endl);
  Cout("Spatial TOF measurement range: " << m_spatialTOFMeasurementRange << " mm" << endl);
  Cout("POI1 (x ; y ; z ) = " << mp_POI1[0] <<" ; " << mp_POI1[1] <<" ; " << mp_POI1[2] <<" ; " << endl);
  Cout("POI2 (x ; y ; z ) = " << mp_POI2[0] <<" ; " << mp_POI2[1] <<" ; " << mp_POI2[2] <<" ; " << endl);
  Cout(flush);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iEventListPET::GetAdditiveCorrections(int a_bin)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  // If TOF, convert the LOR random rate into the random rate matching a single event with continuous TOF measurement
  // (division by the total spatial range of TOF measurements)
  return m_eventScatRate + (m_hasTOFinfo ? (m_eventRdmRate/m_spatialTOFMeasurementRange) : m_eventRdmRate); 
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
