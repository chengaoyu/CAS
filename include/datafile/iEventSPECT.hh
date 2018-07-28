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
  \brief    Declaration of class iEventSPECT
*/

#ifndef IEVENTSPECT_HH
#define IEVENTSPECT_HH 1

#include "vEvent.hh"

/*!
  \class   iEventSPECT
  \brief   Inherit from vEvent. Main SPECT class for the Event objects
  \details This class is designed to be an abstract class that should not be used on its own; only its children are used. \n
           It manages data and functions common to both class dedicated to histogram and list mode SPECT.
*/
class iEventSPECT : public vEvent
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   iEventSPECT constructor. 
               Initialize the member variables to their default values.
    */
    iEventSPECT();
    /*!
      \brief   iEventSPECT destructor. 
    */
    virtual ~iEventSPECT();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      virtual int iEventSPECT::AllocateSpecificData() = 0
      \brief   Pure virtual function implemented in the child classes, dedicated to the allocation of specific data in the child classes
      \return  0 is success, positive value otherwise
    */
    virtual int AllocateSpecificData() = 0;
    /*!
      \fn      virtual void iEventSPECT::Describe()
      \brief   This function can be used to get a description of the event printed out
    */
    virtual void Describe();

  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      inline FLTNB iEventSPECT::GetNormFactor()
      \return  the normalization term
    */
    inline FLTNB GetNormFactor()
           {return m_eventNormFactor;}
    /*!
      \fn      inline FLTNB iEventSPECT::GetEventScatRate()
      \return  the correction term for scatters as a rate in 1/s
    */
    inline FLTNB GetEventScatRate()
           {return m_eventScatRate;}
    /*!
      \fn      inline void iEventSPECT::SetNormalizationFactor()
      \param   a normalization term
      \brief   Cast the FLTNBDATA value passed as parameter in FLTNB, and set it to the normalization term
    */
    inline void SetNormalizationFactor(FLTNBDATA a_value)
           {m_eventNormFactor = (FLTNB)a_value;}
    /*!
      \fn      inline void iEventSPECT::SetScatterRate()
      \param   a_value
      \brief   Cast the FLTNBDATA value passed as parameter in FLTNB, and set it to the scatter correction rate
    */
    inline void SetScatterRate(FLTNBDATA a_value)
           {m_eventScatRate = (FLTNB)a_value;}
    /*!
      \fn      inline FLTNB iEventSPECT::GetAdditiveCorrections()
      \param   a bin (0 if only one line)
      \return  Just return the scatter rate (the TOF bin parameter is ignored for SPECT)
    */
    inline FLTNB GetAdditiveCorrections(int a_bin)
           {return m_eventScatRate;}
    /*!
      \fn      inline FLTNB iEventSPECT::GetMultiplicativeCorrections()
      \return  the product of the multiplicative corrections terms for this event
    */
    inline FLTNB GetMultiplicativeCorrections()
           {return m_eventNormFactor;}
    /*!
      \fn      virtual FLTNB iEventSPECT::GetEventValue() = 0
      \param   a bin (0 if only one line)
      \brief   Pure virtual function implemented in the child classes
      \return  the value of the event
    */
    virtual FLTNB GetEventValue(int a_bin) = 0;
    /*!
      \fn      virtual void iEventSPECT::SetEventValue() = 0
      \param   a bin (0 if only one line)
      \param   a_value
      \brief   Set the event value, this is a pure virtual function implemented in the child classes
    */
    virtual void SetEventValue(int a_bin, FLTNBDATA a_value) = 0;
    /*!
      \fn      virtual INTNB vEvent::GetNbValueBins() = 0
      \brief   Get the number of event value bins
      \return  the number of value bins
    */
    virtual INTNB GetNbValueBins() = 0;
    /*!
      \fn      inline void iEventSPECT::MultiplyAdditiveCorrections()
      \param   FLTNB a_factor
      \brief   Divide additive corrections by the provided factor (scatters)
    */
    inline void MultiplyAdditiveCorrections(FLTNB a_factor)
           {m_eventScatRate *= a_factor;}

  // -------------------------------------------------------------------
  // Private member functions
  private:

  // -------------------------------------------------------------------
  // Data members
  protected:
    FLTNB m_eventScatRate;   /*!< Correction term for scatter rate (unit: s-1). Default value =0.0 */
    FLTNB m_eventNormFactor; /*!< Normalization term. Default value =1.0 */
};

#endif
