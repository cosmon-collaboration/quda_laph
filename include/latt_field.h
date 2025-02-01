#ifndef FIELD_CONTAINER_H
#define FIELD_CONTAINER_H

#include <cstddef>
#include <vector>

namespace LaphEnv {

// **************************************************************************
// *                                                                        *
// *  Class "LattField" is a simple container for a lattice field on the    *
// *  hosts (CPUs).  The type of lattice field is specified by the possible *
// *  values in the FieldSiteType enum.  Each type differs by how many      *
// *  complex values the field has at each site of the lattice.  The        *
// *  possibilities are                                                     *
// *                                                                        *
// *      Complex, ColorMatrix, ColorVector, ColorSpinVector, Unknown.      *
// *                                                                        *
// *  The global constants "FieldNspin" and "FieldNcolor" are also          *
// *  defined and initialized here. The precision of the floating point     *
// *  numbers in the field are taken from the QudaInfo::cpu_prec value.     *
// *                                                                        *
// *  The field itself is stored as a vector<char>.  The constructor sizes  *
// *  this array based on the FieldSiteType given and the precision of      *
// *  floating point numbers.  The "reset" member resizes the vector if     *
// *  necessary.  This class is mainly useful for file I/O of lattice       *
// *  fields and simple storage on the hosts.                               *
// *                                                                        *
// *  QDPXX format is assumed: lattice sites are assumed to be even-odd     *
// *  checkerboard (even/odd sites have x+y+z+t even/odd), with ordering    *
// *  (x,y,z,t) with t varying most slowly for each parity; SU3 color       *
// *  matrix is all 9 complex components row major; ColorSpinVector has     *
// *  the color index varying fastest.                                      *
// *                                                                        *
// **************************************************************************

typedef enum FieldSiteType_s {
  Complex,
  ColorMatrix,
  ColorVector,
  ColorSpinVector,
  Unknown
} FieldSiteType;

// global variables for field

const int FieldNspin = 4;
const int FieldNcolor = 3;

// ***********************************************************************

class LattField {

  std::vector<char> m_data;

  FieldSiteType m_sitetype;

  size_t m_site_elems; // how many complex numbers at each site

  size_t m_site_bytes; // number of bytes at each site

  static int cpu_prec_bytes;

public:
  LattField(FieldSiteType sitetype = Unknown);

  // caution with copy: memory usage
  LattField(const LattField &in)
      : m_data(in.m_data), m_sitetype(in.m_sitetype),
        m_site_elems(in.m_site_elems), m_site_bytes(in.m_site_bytes) {}

  LattField &operator=(const LattField &in) {
    m_data = in.m_data;
    m_sitetype = in.m_sitetype;
    m_site_elems = in.m_site_elems;
    m_site_bytes = in.m_site_bytes;
    return *this;
  }

  size_t bytesPerSite() const { return m_site_bytes; }

  size_t elemsPerSite() const { return m_site_elems; }

  static size_t bytesPerWord() // bytes of complex
  {
    return cpu_prec_bytes;
  }

  FieldSiteType getFieldSiteType() const { return m_sitetype; }

  std::vector<char> &getDataRef() { return m_data; }

  const std::vector<char> &getDataConstRef() const { return m_data; }

  const char *getDataConstPtr() const { return m_data.data(); }

  char *getDataPtr() { return m_data.data(); }

  LattField &reset(FieldSiteType sitetype);

  LattField &clear();

  // Next two routines are useful for testing/checks (slow though); gets
  // and puts data into the right location according to the lattice site
  // "coord"; even/odd checkboard with x,y,z,t column major

  std::vector<char> getSiteData(const std::vector<int> &latt_coord) const;

  void putSiteData(const std::vector<int> &latt_coord,
                   const std::vector<char> &siteData);

  static int get_cpu_prec_bytes() { return cpu_prec_bytes; }

private:
  void do_resize();

  void calc_site_elems();

  void applyFermionTemporalAntiPeriodic();

  LattField &reset_by_precision(FieldSiteType sitetype, char prec);

  LattField &reset_by_bytes_per_site(int bytes_per_site);

  void to_single_precision(std::vector<char> &buffer);

  void to_double_precision(std::vector<char> &buffer);

  LattField &to_quda_precision();

  friend class IOFMHandler;

  friend class GaugeCERNConfigReader;

  friend class QudaInfo;

  friend class GaugeConfigurationHandler;
};
} // namespace LaphEnv
#endif
