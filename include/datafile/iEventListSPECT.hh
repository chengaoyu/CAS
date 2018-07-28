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
  \ingroup  datafile
  \brief    Declaration of class iEventListSPECT
*/

#ifndef IEVENTLISTSPECT_HH
#define IEVENTLISTSPECT_HH 1

#include "iEventSPECT.hh"

/*!
  \class   iEventListSPECT
  \brief   Inherit from iEventSPECT. Class for SPECT list-mode events 
  \details It manages data and functions specific to list mode SPECT.
*/
class iEventListSPECT : public iEventSPECT
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   iEventListSPECT constructor. 
               Initialize the member variables to their default values.
    */
    iEventListSPECT();
    /*!
      \brief   iEventListSPECT destructor. 
    */
    ~iEventListSPECT();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      inline int iEventListSPECT::AllocateSpecificData()
      \brief   Function allowing the allocation of specific data. Return 0 by default for iEventListSPECT
      \return  0 is success, positive value otherwise
    */
    inline int AllocateSpecificData()
           {return 0;}
    /*!
      \fn      void iEventListSPECT::Describe()
      \brief   This function can be used to get a description of the event printed out
    */
    void Describe();

  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      inline uint8_t iEventListSPECT::GetKind() 
      \return  the kind of event
    */
    inline uint8_t GetKind()
           {return m_kind;}
    /*!
      \fn      inline FLTNB* iEventListSPECT::GetPOI()
      \return  the pointer containing the 3 indices for the Point Of Interaction of the crystal
    */
    inline FLTNB* GetPOI()
           {return mp_POI;}
    /*!
      \fn      inline FLTNB* iEventListSPECT::GetPOI()
      \param   a_axis
      \return  the POI value for the specific axis
    */
    inline FLTNB GetPOI(int a_axis)
           {return mp_POI[a_axis];}
    /*!
      \fn      inline FLTNB iEventListSPECT::GetEventValue()
      \return  1 as default for a list-mode Event (a_bin is dedicated to histogram mode, and ignored for list-mode)
    */
    inline FLTNB GetEventValue(int a_bin)
           {return 1.;}
    /*!
      \fn      inline void iEventListSPECT::SetKind()
      \param   a_value
      \brief   Set the kind of coincidence
    */
    inline void SetKind(uint8_t a_value)
           {m_kind = a_value;}
    /*!
      \fn      inline void iEventListSPECT::SetPOI()
      \param   a_axis
      \param   a_value
      \brief   Initialize the POI of the crystal with a value for the specific axis.
    */
    inline void SetPOI(int a_axis, FLTNBDATA a_value)
           {mp_POI[a_axis] = (FLTNB)a_value;}
    /*!
      \fn      void iEventListSPECT::SetEventValue()
      \param   a_bin
      \param   a_value
      \brief   Throw a warning (depending of verbosity) as the event value of a list-mode Event should be equal to 1 and not modified
    */
    void SetEventValue(int a_bin, FLTNBDATA a_value);
    /*!
      \fn      virtual INTNB vEvent::GetNbValueBins()
      \brief   Get the number of event value bins
      \return  Single value so 1 bin
    */
    inline INTNB GetNbValueBins()
           {return 1;}

  // -------------------------------------------------------------------
  // Private member functions
  private:

  // -------------------------------------------------------------------
  // Data members
  private:
    uint8_t m_kind;  /*!< Event type : unknown (=0), true(=1), single scat(=2), multiple scat(=3)) Default value =0  */
    FLTNB mp_POI[3]; /*!< Position of interaction in the crystal along each axis (mm). Default value =0.0;0.0;-1.0 */
};

#endif
