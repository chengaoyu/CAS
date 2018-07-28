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
  \brief    Declaration of class iEventHistoPET
*/

#ifndef IEVENTHISTOPET_HH
#define IEVENTHISTOPET_HH 1

#include "iEventPET.hh"

/*!
  \class   iEventHistoPET
  \brief   Inherit from iEventPET. Class for PET histogram mode events 
  \details It manages data and functions specific to histo mode PET.
*/
class iEventHistoPET : public iEventPET
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   iEventHistoPET constructor. 
               Initialize the member variables to their default values.
    */
    iEventHistoPET();
    /*!
      \brief   iEventHistoPET destructor. 
    */
    ~iEventHistoPET();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      int iEventHistoPET::AllocateSpecificData()
      \brief   Function allowing the allocation of specific data. \n
               Instantiate and initialize the mp_eventValue and mp_eventScatIntensity arrays depending of the number of TOF bins
      \return  0 is success, positive value otherwise
    */
    int AllocateSpecificData();
    /*!
      \fn      void iEventHistoPET::Describe()
      \brief   This function can be used to get a description of the event printed out
    */
    void Describe();

  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      inline FLTNB iEventHistoPET::GetEventScatRate()
      \param   a_bin (0 if noTOF)
      \return  the scatter correction rate in 1/s for this list-mode event for the specific TOF bin
    */
    inline FLTNB GetEventScatRate(int a_bin)
           {return mp_eventScatRate[a_bin];}
    /*!
      \fn      inline uint16_t iEventHistoPET::GetEventNbTOFBins()
      \return  the number of TOF bins in the Event
    */
    inline uint16_t GetEventNbTOFBins()
           {return m_eventNbTOFBins;}
    /*!
      \fn      inline FLTNB iEventHistoPET::GetEventValue()
      \param   a_bin (0 if noTOF)
      \return  the event value corresponding to the specific TOF bin passed as parameter
    */
    inline FLTNB GetEventValue(int a_bin)
           {return mp_eventValue[a_bin];}
    /*!
      \fn      inline FLTNB iEventHistoPET::GetEventScatRate()
      \param   a_bin (0 if noTOF)
      \param   a_value
      \return  the scatter correction rate in 1/s for this histo-mode event for the specific TOF bin
    */
    inline void SetScatterRate(int a_bin, FLTNBDATA a_value)
           {mp_eventScatRate[a_bin] = (FLTNB)a_value;}
    /*!
      \fn      inline FLTNB iEventHistoPET::GetEventScatRate()
      \param   a_value
      \return  the scatter correction rate in 1/s for this histo-mode event for the specific TOF bin
    */
    inline void SetEventNbTOFBins(uint16_t a_value)
           {m_eventNbTOFBins = a_value;}
    /*!
      \fn      inline void iEventHistoPET::SetEventValue()
      \param   a_bin (0 if noTOF)
      \param   a_value
      \brief   Cast the FLTNBDATA value passed in parameters in FLTNB, and use it to set the event value of the specific TOF bin
    */
    inline void SetEventValue(int a_bin, FLTNBDATA a_value)
           {mp_eventValue[a_bin] = (FLTNB)a_value;}
    /*!
      \fn      virtual INTNB vEvent::GetNbValueBins()
      \brief   Get the number of event value bins
      \return  Number of TOF bins here
    */
    inline INTNB GetNbValueBins()
           {return m_eventNbTOFBins;}
    /*!
      \fn      FLTNB iEventHistoPET::GetAdditiveCorrections()
      \return  the sum of additive correction terms, summed for all the TOF bins
    */
    FLTNB GetAdditiveCorrections(int a_bin);
    /*!
      \fn      inline void iEventHistoPET::MultiplyAdditiveCorrections()
      \param   FLTNB a_factor
      \brief   Divide additive corrections by the provided factor (scatters and randoms)
    */
    inline void MultiplyAdditiveCorrections(FLTNB a_factor)
           {m_eventRdmRate *= a_factor; for (uint16_t tof=0; tof<m_eventNbTOFBins; tof++) mp_eventScatRate[tof] *= a_factor;}

  // -------------------------------------------------------------------
  // Private member functions
  private:

  // -------------------------------------------------------------------
  // Data members
  private:
    FLTNB* mp_eventValue;      /*!< Pointer containing the amount of data in each potential TOF bin. Default value =1.0 */
    FLTNB* mp_eventScatRate;   /*!< Pointer containing the scatter correction term (as a rate in s-1) for each potential TOF bin. Default value =0.0 */
    uint16_t m_eventNbTOFBins; /*!< Number of TOF bins in the Event. Default value =1 */
};

#endif
