#ifndef GAUGE_FIELD_HDF5_H
#define GAUGE_FIELD_HDF5_H

#include "latt_field.h"
#include <string>
#include <vector>

namespace LaphEnv {

// *******************************************************************************
// *                                                                             *
// *   Output of a gauge configuration (e.g. the smeared gauge field) to an      *
// *   HDF5 file, with each timeslice in its own group, so that a reader can     *
// *   load single timeslices (pylaph reads these files, and chroma_laph's       *
// *   single-group CONVERT_GAUGE_FIELD_TO_HDF5 files, with                      *
// *   extract_spatial_links.py):                                                *
// *                                                                             *
// *      /GaugeField            group                                           *
// *         attribute "Header"  variable-length string: XML header, root tag    *
// *                             <GaugeFieldMetaData>, which must contain        *
// *                             <LatticeDimensions>Nx Ny Nz Nt</...>            *
// *      /GaugeField/t<t>       one group per timeslice t = 0,1,...,Nt-1        *
// *         re, im              1-d float64 datasets of length 4*Nx*Ny*Nz*9,    *
// *                             real/imaginary parts of the links on that       *
// *                             timeslice in (mu,x,y,z,c1,c2) order, c2 varying *
// *                             fastest; mu=0,1,2,3 is x,y,z,t, and each color  *
// *                             matrix is row major (c1 = row, c2 = column).    *
// *                                                                             *
// *   The data are always written in double precision, whatever the host       *
// *   precision of the LattField.  In a parallel build (ARCH_PARALLEL), HDF5    *
// *   must be built with parallel support: all MPI ranks must call              *
// *   writeGaugeFieldHDF5, and each writes its own sublattice into the one      *
// *   file collectively.  Any existing file of the same name is overwritten.    *
// *                                                                             *
// *******************************************************************************

   // true if "filename" ends in ".hdf5" or ".h5" (case insensitive)
bool isHDF5FileName(const std::string& filename);

   // "gauge" must hold LayoutInfo::Ndim ColorMatrix fields (one per direction)
void writeGaugeFieldHDF5(const std::string& filename,
                         const std::vector<LattField>& gauge,
                         const std::string& header_xml);

}
#endif
