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
  \brief    Declaration of class iEventNorm
*/

#ifndef IEVENTNORM_HH
#define IEVENTNORM_HH 1

#include "vEvent.hh"

/*!
  \class   iEventNorm
  \brief   Inherit from vEvent. Used for normalization events for sensitivity computation
  \details This class is designed to represent a normalization event. It is used for sensitivity computation. \n
           A normalization datafile will be interpreted as a collection of iEventNorm events where each event \n
           contains a list of considered LORs with associated normalization factors.
*/
class iEventNorm : public vEvent
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \brief   iEventNorm constructor. 
               Initialize the member variables to their default values.
    */
    iEventNorm();
    /*!
      \brief   iEventNorm destructor
    */
    ~iEventNorm();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      int iEventNorm::AllocateSpecificData()
      \brief   Inherited and pure virtual from vEvent
      \details Just do nothing here
      \return  0 is success, positive value otherwise
    */
    inline int AllocateSpecificData() {return 0;}
    /*!
      \fn      void iEventNorm::Describe()
      \brief   This function can be used to get a description of the event printed out
    */
    void Describe();

  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      inline void iEventNorm::SetNormalizationFactor(FLTNBDATA a_value)
      \param   a normalization term
      \brief   Cast the FLTNBDATA value passed as parameter in FLTNB, and set it to the normalization term
    */
    inline void SetNormalizationFactor(FLTNBDATA a_value)
           {m_normalizationFactor = (FLTNB)a_value;}
    /*!
      \fn      inline void iEventNorm::SetAttenuationCorrectionFactor(FLTNBDATA a_value)
      \param   an attenuation correction factor
      \brief   Cast the FLTNBDATA value passed as parameter in FLTNB, and set it to the attenuation correction factor
    */
    inline void SetAttenuationCorrectionFactor(FLTNBDATA a_value)
           {m_attenuationCorrectionFactor = (FLTNB)a_value;}
    /*!
      \fn      FLTNB iEventNorm::GetMultiplicativeCorrections()
      \return  The normalization factor times the attenuation correction factor
    */
    inline FLTNB GetMultiplicativeCorrections()
           {return m_normalizationFactor * m_attenuationCorrectionFactor;}
    /*!
      \fn      FLTNB iEventNorm::GetNormalizationFactor()
      \return  The normalization factor
    */
    inline FLTNB GetNormalizationFactor()
           {return m_normalizationFactor;}
    /*!
      \fn      FLTNB iEventNorm::GetAttenuationCorrectionFactor()
      \return  The attenuation correction factor
    */
    inline FLTNB GetAttenuationCorrectionFactor()
           {return m_attenuationCorrectionFactor;}

  // -------------------------------------------------------------------
  // Public functions that should not be used, so they throw an error and crash
  public:
    /*!
      \fn      FLTNB iEventNorm::GetEventValue()
      \param   a bin (0 if noTOF)
      \brief   Not used, so throw an error and Exit
      \return  -1., but crash before
    */
    FLTNB GetEventValue(int a_bin);
    /*!
      \fn      FLTNB iEventNorm::GetAdditiveCorrections()
      \param   a bin (0 if noTOF)
      \brief   Not used, so throw an error and Exit
      \return  -1., but crash before
    */
    FLTNB GetAdditiveCorrections(int a_bin);
    /*!
      \fn      void iEventNorm::SetEventValue()
      \param   a bin (ignored)
      \param   a_value
      \brief   Not used, so throw an error and Exit
    */
    void SetEventValue(int a_bin, FLTNBDATA a_value);
    /*!
      \fn      virtual INTNB iEventNorm::GetNbValueBins()
      \brief   Not used, so throw an error and Exit
      \return  -1., but crash before
    */
    INTNB GetNbValueBins();
    /*!
      \fn      void iEventNorm::MultiplyAdditiveCorrections()
      \param   FLTNB a_factor
      \brief   Not used, so throw an error and Exit
    */
    void MultiplyAdditiveCorrections(FLTNB a_factor);

  // -------------------------------------------------------------------
  // Private member functions
  private:

  // -------------------------------------------------------------------
  // Data members
  protected:
    FLTNB m_normalizationFactor;         /*!< Normalization correction term. Default value = 1.0 */
    FLTNB m_attenuationCorrectionFactor; /*!< Attenuation correction factor. Default value = 1.0 */
};

#endif
