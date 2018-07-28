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
  \brief   Implementation of file to memory mapping
  \details This class was taken from Stephan Brumme website and
           adapted to the CASToR project with almost no modifications.
           Many thanks to him!
           Copyright (c) 2013 Stephan Brumme. All rights reserved.
           see http://create.stephan-brumme.com/disclaimer.html
*/

#ifndef OMEMORYMAPPED_HH
#define OMEMORYMAPPED_HH 1

// define fixed size integer types
#if defined(_WIN32)
#include <inttypes.h>
typedef unsigned __int64 uint64_t;
#else
#include <stdint.h>
#endif

#include "gVariables.hh"
#include "sOutputManager.hh"

/// Portable read-only memory mapping (Windows and Linux)
/** Filesize limited by size_t, usually 2^32 or 2^64 */
class oMemoryMapped
{
  public:
    /// tweak performance
    enum CacheHint
    {
      Normal,         ///< good overall performance
      SequentialScan, ///< read file only once with few seeks
      RandomAccess    ///< jump around
    };

    /// how much should be mappend
    enum MapRange
    {
      WholeFile = 0   ///< everything ... be careful when file is larger than memory
    };

    /// do nothing, must use open()
    oMemoryMapped();
    /// open file, mappedBytes = 0 maps the whole file
    oMemoryMapped(const std::string& filename, size_t mappedBytes = WholeFile, CacheHint hint = Normal);
    /// close file (see close() )
    ~oMemoryMapped();

    /// open file, mappedBytes = 0 maps the whole file
    int Open(const std::string& filename, size_t mappedBytes = WholeFile, CacheHint hint = Normal);
    /// close file
    void Close();

    /// access position, no range checking (faster)
    unsigned char operator[](size_t offset) const;
    /// access position, including range checking
    unsigned char at        (size_t offset) const;

    /// raw access
    const unsigned char* GetData() const;

    /// true, if file successfully opened
    bool IsValid() const;

    /// get file size
    uint64_t size() const;
    /// get number of actually mapped bytes
    size_t   mappedSize() const;

    /// replace mapping by a new one of the same file, offset MUST be a multiple of the page size
    int Remap(uint64_t offset, size_t mappedBytes);

  private:
    /// don't copy object
    oMemoryMapped(const oMemoryMapped&);
    /// don't copy object
    oMemoryMapped& operator=(const oMemoryMapped&);

    /// get OS page size (for remap)
    static int GetPageSize();

    /// file name
    std::string _filename;
    /// file size
    uint64_t    _filesize;
    /// caching strategy
    CacheHint   _hint;
    /// mapped size
    size_t      _mappedBytes;

    /// define handle
#if defined(_WIN32) || defined(CASTOR_USE_MINGW)
    typedef void* FileHandle;
    /// Windows handle to memory mapping of _file
    void*       _mappedFile;
#else
    typedef int   FileHandle;
#endif

    /// file handle
    FileHandle  _file;
    /// pointer to the file contents mapped into memory
    void*       _mappedView;
};

#endif
