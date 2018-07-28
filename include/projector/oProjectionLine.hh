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
  \brief    Declaration of class oProjectionLine
*/

#ifndef OPROJECTIONLINE_HH
#define OPROJECTIONLINE_HH 1

#include "gVariables.hh"
#include "oImageDimensionsAndQuantification.hh"

/**
 * @defgroup COMPUTATION_STRATEGY Computation strategy
 *
 *    \brief Organization strategy for system matrix elements \n
 *           Defined in oProjectionLine.hh
 * @{
 */
/** Constant corresponding to an image-based system matrix elements storage:
    The voxels weights are added in a matrix representing the whole image, so
    the addition of a new line to the previous ones is straightforward only
    by adding the weights to the corresponding voxels. As it is insanely long,
    it can possibly be used for example with extremely complex projectors that
    makes use of huge number of ray tracings for a single event, where the list
    of contributions can become longer than the number of voxels in the image.
    This computation strategy is obviously not compatible with SPECT reconstruction
    with attenuation correction. */
#define IMAGE_COMPUTATION_STRATEGY 1
/** Constant corresponding to a fixed-size list storage of system matrix elements:
    The voxels are added one by one in two separated lists, one containing voxel
    indices and the other voxel weights. When a voxel is added to the oProjectionLine,
    it is simply pilled-up to the list. The list has a fixed size which is provided
    by the EstimateMaxNumberOfVoxelsPerLine() function from the vProjector class.
    There are no ckecks at all for possible buffer overflows. */
#define FIXED_LIST_COMPUTATION_STRATEGY 2
/** Constant corresponding to an adaptative-size list storage of system matrix elements:
    This is the same as the fixed-size strategy except that the size can be upgraded if
    the current number of contributing voxels exceed the list's size. The first allocated
    size corresponds to the diagonal of the image. */
#define ADAPTATIVE_LIST_COMPUTATION_STRATEGY 3
/** @} */


/**
 * @defgroup PROJECTION_DIRECTION Projection direction
 *
 *    \brief Direction of projection (Forward, backward) \n
 *           Defined in oProjectionLine.hh
 * @{
 */
/** Store in the FORWARD direction */
#define FORWARD 0
/** Store in the BACKWARD direction */
#define BACKWARD 1
/** @} */

class vProjector;


/*!
  \class   oProjectionLine
  \brief   This class is designed to manage and store system matrix elements associated to a vEvent.
  \details This class is basically a container for system matrix elements associated to a vEvent.
           It can use different storage and computation strategies for this role. It contains
           the voxel contributions to the LOR associated to a vEvent. It manages TOF bins.
*/
class oProjectionLine
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public oProjectionLine::oProjectionLine()
      \brief   The constructor of oProjectionLine
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    oProjectionLine();
    /*!
      \fn      public oProjectionLine::~oProjectionLine()
      \brief   The destructor of oProjectionLine
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built by this class.
    */
    ~oProjectionLine();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      public int oProjectionLine::CheckParameters()
      \brief   A function used to check the parameters settings
      \details This function does not take any parameter and is used to check that all mandatory
               members were correctly parameterized.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckParameters();
    /*!
      \fn      public int oProjectionLine::Initialize()
      \brief   A function used to initialize a bunch of stuff after parameters have been checked.
      \details It allocates all tables depending on the computation strategy.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */
    int Initialize();
    /*!
      \fn      public void oProjectionLine::ComputeLineLength()
      \brief   Simply compute and update the m_length using the associated mp_position1 and mp_position2.
    */
    void ComputeLineLength();
    /*!
      \fn      public bool oProjectionLine::NotEmptyLine()
      \brief   This function is used to know if the line contains any voxel contribution.
      \return  True if any voxel contribution, false if empty.
    */
    bool NotEmptyLine();
    /*!
      \fn      public void oProjectionLine::Reset()
      \brief   Reset length and all the voxel indices and weights tabs.
    */
    void Reset();
    /*!
      \fn      public void oProjectionLine::ApplyOffset()
      \brief   Apply the offset of oImageDimensionsAndQuantification to the mp_position1 and mp_position2.
    */
    void ApplyOffset();
    /*!
      \fn      public void oProjectionLine::ApplyBedOffset()
      \brief   Apply the bed offset of m_bedOffset to the mp_position1 and mp_position2.
    */
    void ApplyBedOffset();
    /*!
      \fn      public INTNB oProjectionLine::GetVoxelIndex()
      \param   int a_direction
      \param   int a_TOFBin
      \param   INTNB a_voxelInLine
      \brief   This function is used to get the contributing voxel index of the provided direction, TOF bin and voxel rank.
      \return  The voxel index.
    */
    INTNB GetVoxelIndex(int a_direction, int a_TOFBin, INTNB a_voxelInLine);
    /*!
      \fn      public void oProjectionLine::AddVoxelInTOFBin()
      \param   int a_direction
      \param   int a_TOFBin
      \param   INTNB a_voxelIndex
      \param   FLTNB a_voxelWeight
      \brief   This function is used to add a voxel contribution to the line and provided TOF bin.
    */
    void AddVoxelInTOFBin(int a_direction, int a_TOFBin, INTNB a_voxelIndice, FLTNB a_voxelWeight);
    /*!
      \fn      public void oProjectionLine::AddVoxelAllTOFBins()
      \param   int a_direction
      \param   INTNB a_voxelIndex
      \param   FLTNB a_voxelWeight
      \param   HPFLTNB* a_tofWeights : TOF weightings for all TOF bins for this voxel
      \param   INTNB a_tofBinFirst : first relevant TOF bin for this voxel
      \param   INTNB a_tofBinLast : last relevant TOF bin for this voxel
      \brief   Add a voxel contribution to the line for all relevant TOF bins
    */
    void AddVoxelAllTOFBins(int a_direction, INTNB a_voxelIndex, FLTNB a_voxelWeight, HPFLTNB* a_tofWeights, INTNB a_tofBinFirst, INTNB a_tofBinLast);
    /*!
      \fn      public void oProjectionLine::AddVoxel()
      \param   int a_direction
      \param   INTNB a_voxelIndice
      \param   FLTNB a_voxelWeight
      \brief   This function is used to add a voxel contribution to the line, assuming TOF bin 0 (i.e. no TOF).
    */
    void AddVoxel(int a_direction, INTNB a_voxelIndice, FLTNB a_voxelWeight);
    /*!
      \fn      public FLTNB oProjectionLine::ForwardProject()
      \param   FLTNB* ap_image = NULL
      \brief   Simply forward projects the provided image if not null, or else 1, for the current TOF bin.
      \details It forward projects for the current TOF bin the provided image if not NULL, or else it assumes a uniform image of 1.
               It also applies the inverse of the multiplicative correction term, and assumes it is not zero.
      \return  The value of the forward projection
    */
    FLTNB ForwardProject(FLTNB* ap_image = NULL);
    /*!
      \fn      public FLTNB oProjectionLine::ForwardProjectWithSPECTAttenuation()
      \param   FLTNB* ap_attenuation
      \param   FLTNB* ap_image = NULL
      \brief   Forward projects the provided image for the current TOF bin with an inner loop on the attenuation (for SPECT).
      \details It forward projects the provided image if not null, or else it assumes a uniform image of 1. It does it for the current
               TOF bin (for genericity purpose) even if this function should only be called for SPECT. The order of the voxels contributions
               is assumed to be from the outside to the detector, so it only works with the list computation strategy. It also applies the
               inverse of the multiplicative correction term, and assumes it is not zero.
      \return  The value of the forward projection
    */
    FLTNB ForwardProjectWithSPECTAttenuation(FLTNB* ap_attenuation, FLTNB* ap_image = NULL);
    /*!
      \fn      public void oProjectionLine::BackwardProject()
      \param   FLTNB* ap_image
      \param   FLTNB a_value
      \brief   Simply backward projects the provided value inside the provided image, for the current TOF bin.
      \details It backward projects for the current TOF bin the provided value inside the provided image.
               It also applies the inverse of the multiplicative correction term before the backward projection, and assumes it is not zero.
    */
    void BackwardProject(FLTNB* ap_image, FLTNB a_value);
    /*!
      \fn      public void oProjectionLine::BackwardProjectWithSPECTAttenuation()
      \param   FLTNB* ap_attenuation
      \param   FLTNB* ap_image
      \param   FLTNB a_value
      \brief   Backward project the provided value inside the provided image with an inner loop on the attenuation (for SPECT).
      \details It backward projects the provided into the provided image. It does it for the current TOF bin (for genericity purpose)
               even if this function should only be called for SPECT. The order of the voxels contributions is assumed to be from
               the outside to the detector, so it only works with the list computation strategy. It also applies the inverse of the
               multiplicative correction term before the backward projection, and assumes it is not zero.
    */
    void BackwardProjectWithSPECTAttenuation(FLTNB* ap_attenuation, FLTNB* ap_image, FLTNB a_value);
    /*!
      \fn      public FLTNB oProjectionLine::ComputeLineIntegral()
      \param   int a_direction
      \brief   It simply computes the sum of all voxels contributions following the provided direction
      \return  The sum of all voxels contributions
    */
    FLTNB ComputeLineIntegral(int a_direction);


  // -------------------------------------------------------------------
  // Get functions
  public:
    /*!
      \fn      public inline FLTNB oProjectionLine::GetVoxelWeights()
      \param   int a_direction
      \param   int a_TOFBin
      \param   INTNB a_voxelInLine
      \brief   This function is used to get the contributing voxel weight of the provided direction, TOF bin and voxel rank.
      \return  The voxel weight.
    */
    inline FLTNB GetVoxelWeights(int a_direction, int a_TOFBin, INTNB a_voxelInLine)
           {return m3p_voxelWeights[a_direction][a_TOFBin][a_voxelInLine];}
    /*!
      \fn      public inline INTNB oProjectionLine::GetCurrentNbVoxels()
      \param   int a_direction
      \param   int a_TOFBin
      \brief   This function is used to get the current number of contributing voxels to the line.
      \return  The current number of contributing voxels.
    */
    inline INTNB GetCurrentNbVoxels (int a_direction, int a_TOFBin)
           {return m2p_currentNbVoxels[a_direction][a_TOFBin];}
    /*!
      \fn      public inline int oProjectionLine::GetNbTOFBins()
      \brief   This function is used to get the number of TOF bins in use.
      \return  The number of TOF bins.
    */
    inline int GetNbTOFBins()
           {return m_nbTOFBins;}
    /*!
      \fn      public inline FLTNB oProjectionLine::GetLength()
      \brief   This function is used to get the length of the line.
      \return  The length.
    */
    inline FLTNB GetLength()
           {return m_length;}
    /*!
      \fn      public inline int oProjectionLine::GetComputationStrategy()
      \brief   This function is used to get the computation strategy.
      \return  The computation strategy.
    */
    inline int GetComputationStrategy()
           {return m_computationStrategy;}
    /*!
      \fn      public inline FLTNB* oProjectionLine::GetPosition1()
      \brief   This function is used to get the pointer to the mp_position1 (3-values tab).
      \return  The position of point 1.
    */
    inline FLTNB* GetPosition1()
           {return mp_position1;}
    /*!
      \fn      public inline FLTNB* oProjectionLine::GetPosition2()
      \brief   This function is used to get the pointer to the mp_position2 (3-values tab).
      \return  The position of point 2.
    */
    inline FLTNB* GetPosition2()
           {return mp_position2;}
    /*!
      \fn      public inline FLTNB* oProjectionLine::GetBufferPosition1()
      \brief   This function is used to get the pointer to the mp_bufferPosition1 (3-values tab).
      \return  The buffer position of point 1.
    */
    inline FLTNB* GetBufferPosition1()
           {return mp_bufferPosition1;}
    /*!
      \fn      public inline FLTNB* oProjectionLine::GetBufferPosition2()
      \brief   This function is used to get the pointer to the mp_bufferPosition2 (3-values tab).
      \return  The buffer position of point 2.
    */
    inline FLTNB* GetBufferPosition2()
           {return mp_bufferPosition2;}
    /*!
      \fn      public inline FLTNB* oProjectionLine::GetOrientation1()
      \brief   This function is used to get the pointer to the mp_orientation1 (3-values tab).
      \return  The orientation of point 1.
    */
    inline FLTNB* GetOrientation1()
           {return mp_orientation1;}
    /*!
      \fn      public inline FLTNB* oProjectionLine::GetOrientation2()
      \brief   This function is used to get the pointer to the mp_orientation2 (3-values tab).
      \return  The orientation of point 2.
    */
    inline FLTNB* GetOrientation2()
           {return mp_orientation2;}
    /*!
      \fn      public inline FLTNB* oProjectionLine::GetBufferOrientation1()
      \brief   This function is used to get the pointer to the mp_bufferOrientation1 (3-values tab).
      \return  The buffer orientation of point 1.
    */
    inline FLTNB* GetBufferOrientation1()
           {return mp_bufferOrientation1;}
    /*!
      \fn      public inline FLTNB* oProjectionLine::GetBufferOrientation2()
      \brief   This function is used to get the pointer to the mp_bufferOrientation2 (3-values tab).
      \return  The buffer orientation of point 2.
    */
    inline FLTNB* GetBufferOrientation2()
           {return mp_bufferOrientation2;}
    /*!
      \fn      public inline FLTNB* oProjectionLine::GetPOI1()
      \brief   This function is used to get the pointer to POI of point 1 (3-values tab).
      \return  POI of point 1.
    */
    inline FLTNB* GetPOI1()
           {return mp_POI1;}
    /*!
      \fn      public inline FLTNB* oProjectionLine::GetPOI2()
      \brief   This function is used to get the pointer to POI of point 2 (3-values tab).
      \return  POI of point 2 mp_POI2
    */
    inline FLTNB* GetPOI2()
           {return mp_POI2;}
    /*!
      \fn      public inline int oProjectionLine::GetTOFResolution()
      \brief   This function is used to get the TOF resolution.
      \return  TOF resolution.
    */
    inline FLTNB GetTOFResolution()
           {return m_TOFResolution;}
    /*!
      \fn      public inline int oProjectionLine::GetTOFMeasurement()
      \brief   This function is used to get the TOF measurement.
      \return  TOF measurement.
    */
    inline FLTNB GetTOFMeasurement()
           {return m_TOFMeasurement;}
    /*!
      \fn      public inline int oProjectionLine::GetTOFBinSize()
      \brief   This function is used to get the size in ps of a TOF bin.
      \return  TOF bin size in ps.
    */
    inline FLTNB GetTOFBinSize()
           {return m_TOFBinSize;}
    /*!
      \fn      public inline int oProjectionLine::GetIndex1()
      \brief   This function is used to get the index associated to point 1.
      \return  Index of point 1.
    */
    inline int GetIndex1()
           {return m_index1;}
    /*!
      \fn      public inline int oProjectionLine::GetIndex2()
      \brief   This function is used to get the index associated to point 2.
      \return  Index of point 2.
    */
    inline int GetIndex2()
           {return m_index2;}
    /*!
      \fn      public inline int oProjectionLine::GetThreadNumber()
      \brief   This function is used to get the thread number associated to this line.
      \return  The associated thread number.
    */
    inline int GetThreadNumber()
           {return m_threadNumber;}
    /*!
      \fn      public inline FLTNB oProjectionLine::GetBedOffset()
      \brief   This function is used to get the axial bed offset associated to this line.
      \return  The associated bed offset
    */
    inline FLTNB GetBedOffset()
           {return m_bedOffset;}


  // -------------------------------------------------------------------
  // Set functions
  public:
    /*!
      \fn      public inline void oProjectionLine::SetLength()
      \param   FLTNB a_length
      \brief   This function is used to set the length of the line.
    */
    inline void SetLength(FLTNB a_length)
           {m_length = a_length;}
    /*!
      \fn      public inline void oProjectionLine::SetPOI1()
      \param   FLTNB* ap_POI1
      \brief   This function is used to set the POI of point 1.
    */
    inline void SetPOI1(FLTNB* ap_POI1)
           {mp_POI1 = ap_POI1;}
    /*!
      \fn      public inline void oProjectionLine::SetPOI2()
      \param   FLTNB* ap_POI2
      \brief   This function is used to set the POI of point 2.
    */
    inline void SetPOI2(FLTNB* ap_POI2)
           {mp_POI2 = ap_POI2;}
    /*!
      \fn      public inline void oProjectionLine::SetTOFMeasurement()
      \param   FLTNB a_TOFMeasurement
      \brief   This function is used to set the TOF measurement associated to the line.
    */
    inline void SetTOFMeasurement(FLTNB a_TOFMeasurement)
           {m_TOFMeasurement = a_TOFMeasurement;}
    /*!
      \fn      public inline void oProjectionLine::SetTOFBinSize()
      \param   FLTNB a_TOFBinSize
      \brief   This function is used to set the size of a TOF bin in ps.
    */
    inline void SetTOFBinSize(FLTNB a_TOFBinSize)
           {m_TOFBinSize = a_TOFBinSize;}
    /*!
      \fn      public inline void oProjectionLine::SetIndex1()
      \param   int a_index1
      \brief   This function is used to set the index m_index1 of point 1.
    */
    inline void SetIndex1(int a_index1)
           {m_index1 = a_index1;}
    /*!
      \fn      public inline void oProjectionLine::SetIndex2()
      \param   int a_index2
      \brief   This function is used to set the index m_index1 of point 2.
    */
    inline void SetIndex2(int a_index2)
           {m_index2 = a_index2;}
    /*!
      \fn      public inline void oProjectionLine::SetNbTOFBins()
      \param   int a_nbTOFBins
      \brief   This function is used to set the number of TOF bins in use.
    */
    inline void SetNbTOFBins(int a_nbTOFBins)
           {m_nbTOFBins = a_nbTOFBins;}
    /*!
      \fn      public inline void oProjectionLine::SetCurrentTOFBin()
      \param   int a_TOFBin
      \brief   This function is used to set the current TOF bin that is used.
    */
    inline void SetCurrentTOFBin(int a_TOFBin)
           {m_currentTOFBin = a_TOFBin;}
    /*!
      \fn      public inline void oProjectionLine::SetMatchedProjectors()
      \param   bool a_UseMatchedProjectors
      \brief   This function is used to set the boolean that says if we use matched projectors.
    */
    inline void SetMatchedProjectors(bool a_UseMatchedProjectors)
           {m_useMatchedProjectors = a_UseMatchedProjectors;}
    /*!
      \fn      public inline void oProjectionLine::SetTOFResolution()
      \param   FLTNB a_TOFResolution
      \brief   This function is used to set the TOF resolution in use.
    */
    inline void SetTOFResolution(FLTNB a_TOFResolution)
           {m_TOFResolution = a_TOFResolution;}
    /*!
      \fn      public inline void oProjectionLine::SetPOIResolution()
      \param   FLTNB* ap_POIResolution
      \brief   This function is used to set the POI resolution along the 3 axes.
    */
    inline void SetPOIResolution(FLTNB* ap_POIResolution)
           {mp_POIResolution = ap_POIResolution;}
    /*!
      \fn      public inline void oProjectionLine::SetImageDimensionsAndQuantification()
      \param   oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification
      \brief   This function is used to set the pointer to the oImageDimensionsAndQuantification in use.
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
           {mp_ImageDimensionsAndQuantification = ap_ImageDimensionsAndQuantification;}
    /*!
      \fn      public inline void oProjectionLine::SetComputationStrategy()
      \param   int a_computationStrategy
      \brief   This function is used to set the computation strategy in use.
    */
    inline void SetComputationStrategy(int a_computationStrategy)
           {m_computationStrategy = a_computationStrategy;}
    /*!
      \fn      public inline void oProjectionLine::SetForwardProjector()
      \param   vProjector* ap_Projector
      \brief   This function is used to set the pointer to the forward projector.
    */
    inline void SetForwardProjector(vProjector* ap_Projector)
           {mp_ForwardProjector = ap_Projector;}
    /*!
      \fn      public inline void oProjectionLine::SetBackwardProjector()
      \param   vProjector* ap_Projector
      \brief   This function is used to set the pointer to the backward projector.
    */
    inline void SetBackwardProjector(vProjector* ap_Projector)
           {mp_BackwardProjector = ap_Projector;}
    /*!
      \fn      public inline void oProjectionLine::SetThreadNumber()
      \param   int a_threadNumber
      \brief   This function is used to set the thread number of this particular line.
    */
    inline void SetThreadNumber(int a_threadNumber)
           {m_threadNumber = a_threadNumber;}
    /*!
      \fn      public inline void oProjectionLine::SetMultiplicativeCorrection()
      \param   FLTNB a_multiplicativeCorrection
      \brief   This function is used to set the multiplicative correction to be applied during forward and backward projections.
    */
    inline void SetMultiplicativeCorrection(FLTNB a_multiplicativeCorrection)
           {m_multiplicativeCorrection = a_multiplicativeCorrection;}
    /*!
      \fn      public inline void oProjectionLine::SetVerbose()
      \param   int a_verbose
      \brief   This function is used to set the verbose level.
    */
    inline void SetVerbose(int a_verbose)
           {m_verbose = a_verbose;}
    /*!
      \fn      public inline void oProjectionLine::SetBedOffset()
      \param   FLTNB a_bedOffset
      \brief   This function is used to set the bed offset
    */
    inline void SetBedOffset(FLTNB a_bedOffset)
           {m_bedOffset = a_bedOffset;}


  // -------------------------------------------------------------------
  // Data members
  private:

    // Verbose level
    int m_verbose;                         /*!< The verbose level */
    // Has been checked ?
    bool m_checked;                        /*!< Boolean that says if the parameters were checked or not */
    // Has been initialized ?
    bool m_initialized;                    /*!< Boolean that says if the manager was initialized or not */

    // ---------------------------------------------------------------------------------------------------
    // Common stuff
    // ---------------------------------------------------------------------------------------------------

    // The thread number associated to this projection line
    int m_threadNumber;                    /*!< Thread number associated to this projection line */
    // The current multiplicative correction factor (that will be applied during forward and backward projections)
    FLTNB m_multiplicativeCorrection;      /*!< Multiplicative correction factor that will be applied during forward and backward projections */
    // Image dimensions
    oImageDimensionsAndQuantification* 
      mp_ImageDimensionsAndQuantification; /*!< Pointer to the oImageDimensionsAndQuantification object in use */
    // The computation strategy for multiple lines (see comments above)
    int m_computationStrategy;             /*!< Integer defining the computation/storage strategy of the contributing voxels */
    // The bed offset (when multiple bed positions are reconstructed simultaneously
    FLTNB m_bedOffset;                     /*!< Bed axial offset when reconstructing multiple bed positions at once, in mm */
    // The number of TOF bins, TOF resolution and measurement
    int m_nbTOFBins;                       /*!< The number of TOF bins in use */
    FLTNB m_TOFResolution;                 /*!< The TOF resolution in use */
    FLTNB m_TOFMeasurement;                /*!< The current TOF measurement of the event */
    int m_currentTOFBin;                   /*!< The current TOF bin in use, can be used for simple to calls to forward or backward projections */
    FLTNB m_TOFBinSize;                    /*!< The size in ps of a single TOF bin */
    // The POI and its resolution along the 3 axis
    FLTNB* mp_POI1;                        /*!< The current POI of point 1 of the event (along the 3 axes) */
    FLTNB* mp_POI2;                        /*!< The current POI of point 2 of the event (along the 3 axes) */
    FLTNB* mp_POIResolution;               /*!< The POI resolution in use (along the 3 axes) */
    // This is the length of the line
    FLTNB m_length;                        /*!< The current length of the line */
    // These are the positions and orientations of the two end points (the buffer ones are used only as buffer when compression)
    FLTNB* mp_position1;                   /*!< The current position of point 1 (along the 3 axes) */
    FLTNB* mp_position2;                   /*!< The current position of point 2 (along the 3 axes) */
    FLTNB* mp_bufferPosition1;             /*!< A buffer for position of point 1 (used when compression) */
    FLTNB* mp_bufferPosition2;             /*!< A buffer for position of point 2 (used when compression) */
    FLTNB* mp_orientation1;                /*!< The current orientation of point 1 (along the 3 axes) */
    FLTNB* mp_orientation2;                /*!< The current orientation of point 2 (along the 3 axes) */
    FLTNB* mp_bufferOrientation1;          /*!< A buffer for orientation of point 1 (used when compression) */
    FLTNB* mp_bufferOrientation2;          /*!< A buffer for orientation of point 2 (used when compression) */
    // These are the indices associated to the two end points (in case of compression, set to -1)
    int m_index1;                          /*!< The current index of point 1 (associated to the vScanner in use) */
    int m_index2;                          /*!< The current index of point 2 (associated to the vScanner in use) */
    // The rest of the data members can be different for forward and backward operations.
    // However only one operation can be used for simple forward or backward needs.
    // This is managed by this boolean flag that said if we use matched projectors or not.
    bool m_useMatchedProjectors;           /*!< Boolean that says if we use matched projectors (forward = backward) */
    vProjector* mp_ForwardProjector;       /*!< Pointer to the forward projector in use */
    vProjector* mp_BackwardProjector;      /*!< Pointer to the backward projector in use */

    // ------------------------------------------------------------------------------------------------------
    // For the rest below, the first pointer is used to discriminate between forward and backward projectors,
    // while the second pointer discriminates between TOF bins
    // ------------------------------------------------------------------------------------------------------

    // The allocated number of voxels corresponds to the maximum number of contributing
    // voxels a line of the ProjectionLine can handle. Whereas the current number of voxels
    // is the number of voxels that a line is currently using, when using the LIST_COMPUTATION
    // strategy.
    INTNB** m2p_allocatedNbVoxels;         /*!< Number of allocated voxels for each direction and each TOF bin */
    INTNB** m2p_currentNbVoxels;           /*!< Current number of voxels for each direction and each TOF bin */
    // This contains the voxels indices for each lines, when using the both LIST_COMPUTATION strategies
    INTNB*** m3p_voxelIndices;             /*!< List of contributing voxel indices for each direction and each TOF bin */
    // This contains the voxels weights for both computation strategies
    FLTNB*** m3p_voxelWeights;             /*!< List of contributing voxel weights for each direction and each TOF bin */
};

#endif
