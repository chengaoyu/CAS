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
  \ingroup  optimizer
  \brief    Declaration of class vPenalty
*/

#ifndef VPENALTY_HH
#define VPENALTY_HH 1

#include "gVariables.hh"
#include "oImageDimensionsAndQuantification.hh"

/*!
  \class   vPenalty
  \brief   This class is designed to generically described any penalty applied to MAP algorithms
  \details This class is an abstract one, in the sense that it cannot be used on its own
           because several pure virtual functions belong to it. Its children are
           implementations of actual penalties. Everywhere in the code, this parent class
           should be used instead of any of its children.
           Nothing is yet implemented. To be designed.
*/
class vPenalty
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public vPenalty::vPenalty()
      \brief   The constructor of vPenalty
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    vPenalty();
    /*!
      \fn      virtual public vPenalty::~vPenalty()
      \brief   The destructor of vPenalty
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built by this class.
               It is virtual, so that it is automatically called when a child object is deleted.
    */
    virtual ~vPenalty();


  // -------------------------------------------------------------------
  // Public member functions
  public:


  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      public inline void vPenalty::SetVerbose()
      \param   int a_verboseLevel
      \brief   Set the verbose level.
    */
    inline void SetVerbose(int a_verbose)
           {m_verbose = a_verbose;}
    /*!
      \fn      public inline void vPenalty::SetImageDimensionsAndQuantification()
      \param   oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification
      \brief   Set the pointer to the image dimensions in use.
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
           {mp_ImageDimensionsAndQuantification = ap_ImageDimensionsAndQuantification;}
    /*
    inline int GetPenaltyEnergyFunctionDerivativesOrder()
           {return m_penaltyEnergyFunctionDerivativesOrder;}
    */


  // -------------------------------------------------------------------
  // Private member functions
  private:
  

  // -------------------------------------------------------------------
  // Data members
  protected:
    // Image dimensions
    oImageDimensionsAndQuantification* 
      mp_ImageDimensionsAndQuantification;         /*!< The pointer to the image dimensions and quantification object */
    // Verbosity
    int m_verbose;                                 /*!< The verbose level */
    // Order of penalty energy function derivatives
    //int m_penaltyEnergyFunctionDerivativesOrder; /*!< The derivative order of the penalty function */
};


// ---------------------------------------------------------------------
// Part of code that manages the auto declaration of children classes
// ---------------------------------------------------------------------

// Macro for the function that creates the object
#define FUNCTION_PENALTY(CLASS) \
  static vPenalty *make_penalty() { return new CLASS(); };

// Macro for the class that links the appropriate function to the map of objects
#define CLASS_PENALTY(NAME,CLASS)                                                          \
  class NAME##PenaltyCreator                                                               \
  {                                                                                        \
    public:                                                                                \
      NAME##PenaltyCreator()                                                               \
        { sAddonManager::GetInstance()->mp_listOfPenalties[#NAME] = CLASS::make_penalty; } \
  };                                                                                       \
  static NAME##PenaltyCreator PenaltyCreator##NAME;

#endif
