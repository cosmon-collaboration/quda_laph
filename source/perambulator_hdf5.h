#ifndef PERAMBULATOR_HDF5_H
#define PERAMBULATOR_HDF5_H

#include <complex>
#include <string>
#include <vector>

namespace LaphEnv {

// *******************************************************************************
// *                                                                             *
// *   Output of quark perambulators to one HDF5 file, in the layout written by  *
// *   doPerambulatorConvert (last_laph, inline_convert_perambulators.cc):       *
// *   (sparse-grid perambulators: doPerambulatorConvertSG, with the group       *
// *   /PerambulatorSparseGridData, header <PerambulatorHandlerSparseGrid-       *
// *   DataFile>, and nrows = 3*NGrid rows n = color + 3*grid point instead of   *
// *   the nEv sink Laph eigvecs: length 3*NGrid*nEv, element n*nEv+m)           *
// *                                                                             *
// *      /PerambulatorData            group                                     *
// *         attribute "Header"        variable-length string: the XML header    *
// *                                   <PerambulatorHandlerDataFile> of          *
// *                                   PerambulatorHandler::getHeader            *
// *      /PerambulatorData/srcTime<t0>_snkTime<t>/srcSpin<a>/snkSpin<b>/re, im  *
// *                                   1-d float64 datasets of length nEv*nEv:   *
// *                                   element n*nEv+m is the real/imaginary     *
// *                                   part of the perambulator for sink Laph    *
// *                                   eigvec n and source Laph eigvec m; the    *
// *                                   spins a, b run from 1 to Nspin.           *
// *                                                                             *
// *   The datasets are created filled with NaN, and each source eigvec (one    *
// *   column m) is written when computed, so entries not (yet) computed are     *
// *   NaN.  Existing files are extended, so that the columns can be computed    *
// *   in several runs; "queryColumn" tells if a column is complete.             *
// *                                                                             *
// *   All HDF5 I/O is done by the primary MPI rank.  All ranks must make the    *
// *   same calls (the results of open and queryColumn are broadcast); the data  *
// *   passed to putColumns are only used on the primary rank.                   *
// *                                                                             *
// *******************************************************************************

class PerambulatorHDF5Writer
{
   long long m_file;        // hid_t (primary rank only)
   std::string m_filename;
   std::string m_group;     // "/PerambulatorData" or "/PerambulatorSparseGridData"
   int m_nrows;             // sink index range (nEv, or 3*NGrid)
   int m_nev;               // source Laph eigvecs (columns)

 public:

   PerambulatorHDF5Writer() : m_file(-1), m_nrows(0), m_nev(0) {}
   ~PerambulatorHDF5Writer() {close();}

      // Opens "filename" for writing perambulators under group "group" as
      // nrows x nev matrices ("nev" source Laph eigvecs).  If "overwrite" or
      // the file does not exist, a new file is created with "header_xml" as
      // the header of "group", and an empty string returned.  Otherwise the
      // existing file is opened, and its header returned (on all ranks) for
      // the caller to check.
   std::string open(const std::string& filename, const std::string& group,
                    const std::string& header_xml, int nrows, int nev, bool overwrite);

   bool isOpen() const {return !m_filename.empty();}

   const std::string& getFileName() const {return m_filename;}

      // True (on all ranks) if column "srcev" of the matrices for source
      // time/spin "src_time", "src_spin", every sink time in "snk_times" and
      // every sink spin 1..nspin exists and contains no NaN.
   bool queryColumn(int src_time, int src_spin, const std::vector<int>& snk_times,
                    int nspin, int srcev);

      // Writes columns "srcevs" of the matrix for (src_time, src_spin,
      // snk_time, snk_spin), creating it if needed: "data" holds nev*srcevs.size()
      // values, data[n*srcevs.size()+k] for sink eigvec n and source eigvec
      // srcevs[k].
   void putColumns(int src_time, int src_spin, int snk_time, int snk_spin,
                   const std::vector<int>& srcevs,
                   const std::vector<std::complex<double>>& data);

   void flush();

   void close();

 private:

   PerambulatorHDF5Writer(const PerambulatorHDF5Writer&);
   PerambulatorHDF5Writer& operator=(const PerambulatorHDF5Writer&);
};

}
#endif
