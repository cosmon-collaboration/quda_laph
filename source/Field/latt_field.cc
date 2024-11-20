#include "latt_field.h"
#include "laph_stdio.h"
#include "layout_info.h"
#include "utils.h"
#include <complex>
#include <cstring>

using namespace std;

namespace LaphEnv {

int LattField::cpu_prec_bytes = 2 * sizeof(double); // default value

LattField::LattField(FieldSiteType sitetype) : m_sitetype(sitetype) {
  do_resize();
}

LattField &LattField::reset(FieldSiteType sitetype) {
  m_sitetype = sitetype;
  do_resize();
  return *this;
}

LattField &LattField::clear() {
  m_sitetype = FieldSiteType::Unknown;
  do_resize();
  return *this;
}

void LattField::do_resize() {
  calc_site_elems();
  m_site_bytes = cpu_prec_bytes * m_site_elems;
  if (m_site_bytes > 0)
    m_data.resize(LayoutInfo::getRankLatticeNumSites() * m_site_bytes);
  else
    m_data.clear();
}

void LattField::calc_site_elems() {
  switch (m_sitetype) {
  case Complex:
    m_site_elems = 1;
    break;
  case ColorMatrix:
    m_site_elems = FieldNcolor * FieldNcolor;
    break;
  case ColorVector:
    m_site_elems = FieldNcolor;
    break;
  case ColorSpinVector:
    m_site_elems = FieldNspin * FieldNcolor;
    break;
  default:
    m_site_elems = 0;
    ;
    break;
  }
}

//  This is private.  Useful for constructing a lattice field that has
//  a precision that DIFFERS from that in QudaInfo.  This is sometimes
//  needed when reading from file, since the data in the file may not
//  have the same precision as that requested in QudaInfo.
//  "prec" should be 'S' or 'D' for single or double precision.

LattField &LattField::reset_by_precision(const FieldSiteType sitetype,
                                         const char prec) {
  m_sitetype = sitetype;
  calc_site_elems();
  if (prec == 'S') {
    m_site_bytes = 2 * sizeof(float) * m_site_elems;
  } else if (prec == 'D') {
    m_site_bytes = 2 * sizeof(double) * m_site_elems;
  } else {
    throw(std::invalid_argument("Invalid precision in LattField"));
  }
  if (m_site_bytes > 0) {
    m_data.resize(LayoutInfo::getRankLatticeNumSites() * m_site_bytes);
  } else {
    m_data.clear();
  }
  return *this;
}

//  This is used by the IOHandlerFM.  It resets the storage
//  space to accommodate the field as stored in file.  The precision
//  may or may not match the current quda precision.  The IOHandlerFM
//  read must then convert the precision as needed.

LattField &LattField::reset_by_bytes_per_site(const int bytes_per_site) {
  m_site_bytes = bytes_per_site;
  switch (m_site_bytes) {
  case 8:
  case 16:
    m_site_elems = 1;
    m_sitetype = Complex;
    break;
  case 24:
  case 48:
    m_site_elems = 3;
    m_sitetype = ColorVector;
    break;
  case 72:
  case 144:
    m_site_elems = 9;
    m_sitetype = ColorMatrix;
    break;
  case 96:
  case 192:
    m_site_elems = 12;
    m_sitetype = ColorSpinVector;
    break;
  default:
    errorLaph(make_str("Invalid resizing of LattField with bytes_per_site = ",
                       bytes_per_site));
    break;
  }
  if (m_site_bytes > 0) {
    m_data.resize(LayoutInfo::getRankLatticeNumSites() * m_site_bytes);
  } else {
    m_data.clear();
  }
  return *this;
}

std::vector<char>
LattField::getSiteData(const std::vector<int> &latt_coords) const {
  std::vector<char> buffer(m_site_bytes);
  int send_rank, rank_site_linear_index;
  LayoutInfo::getCommInfoFromLatticeCoords(latt_coords, send_rank,
                                           rank_site_linear_index);
  if (send_rank == LayoutInfo::getMyRank()) {
    int offset = m_site_bytes * rank_site_linear_index;
    std::memcpy(buffer.data(), &m_data[offset], m_site_bytes);
  }
#if ARCH_PARALLEL
  comm_broadcast(buffer.data(), m_site_bytes, send_rank);
#endif
  return buffer;
}

void LattField::putSiteData(const std::vector<int> &latt_coords,
                            const std::vector<char> &siteData) {
  if (siteData.size() != m_site_bytes) {
    throw(std::runtime_error("Invalid site data to put into LattField"));
  }
  int recv_rank, rank_site_linear_index;
  LayoutInfo::getCommInfoFromLatticeCoords(latt_coords, recv_rank,
                                           rank_site_linear_index);
  if (recv_rank == LayoutInfo::getMyRank()) {
    int offset = m_site_bytes * rank_site_linear_index;
    std::memcpy(&m_data[offset], siteData.data(), m_site_bytes);
  }
}

// This converts to the precision requested by Quda (useful
// when reading field that was stored in a different
// precision).

LattField &LattField::to_quda_precision() {
  size_t cmplx_bytes = m_site_bytes / m_site_elems;
  if (cmplx_bytes == size_t(cpu_prec_bytes)) {
    return *this;
  } // nothing to do
  else if (cmplx_bytes ==
           size_t(2 * cpu_prec_bytes)) { // convert double to single
    vector<char> buffer;
    to_single_precision(buffer);
    m_data = std::move(buffer);
    m_site_bytes /= 2;
  } else if ((2 * cmplx_bytes) ==
             size_t(cpu_prec_bytes)) { // convert single to double (warn!)
    vector<char> buffer;
    to_double_precision(buffer);
    m_data = std::move(buffer);
    m_site_bytes *= 2;
  }
  return *this;
}

void LattField::to_single_precision(vector<char> &buffer) {
  buffer.resize(m_data.size() / 2);
  int n = 2 * LayoutInfo::getRankLatticeNumSites() *
          m_site_elems; // number of reals
  float *dest = reinterpret_cast<float *>(buffer.data());
  double *src = reinterpret_cast<double *>(m_data.data());
  for (int k = 0; k < n; ++k, ++dest, ++src) {
    *dest = *src;
  }
}

void LattField::to_double_precision(vector<char> &buffer) {
  printLaph(
      "WARNING: converting lattice field from single to double precision");
  buffer.resize(2 * m_data.size());
  int n = 2 * LayoutInfo::getRankLatticeNumSites() *
          m_site_elems; // number of reals
  double *dest = reinterpret_cast<double *>(buffer.data());
  float *src = reinterpret_cast<float *>(m_data.data());
  for (int k = 0; k < n; ++k, ++dest, ++src) {
    *dest = *src;
  }
}

void LattField::applyFermionTemporalAntiPeriodic() {
  const int Tdir = LayoutInfo::Ndim - 1;
  if (LayoutInfo::getMyCommCoords()[Tdir] !=
      (LayoutInfo::getCommNumPartitions()[Tdir] - 1)) {
    return;
  }
  const int start_parity = LayoutInfo::getMyStartParity();
  const int npsites = LayoutInfo::getRankLatticeNumSites() / 2;
  const vector<int> &N = LayoutInfo::getRankLattExtents();
  const bool dp = (bytesPerWord() == sizeof(complex<double>));
  const int nreal = 2 * elemsPerSite();
  const int t = N[Tdir] - 1;
  for (int z = 0; z < N[2]; ++z)
    for (int y = 0; y < N[1]; ++y)
      for (int x = 0; x < N[0]; ++x) {
        const int sindex =
            nreal * ((x + N[0] * (y + N[1] * (z + N[2] * t))) / 2 +
                     npsites * ((start_parity + x + y + z + t) % 2));
        if (dp) {
          double *ptr = reinterpret_cast<double *>(m_data.data()) + sindex;
          for (int s = 0; s < nreal; ++s, ++ptr)
            *ptr = -(*ptr);
        } else {
          float *ptr = reinterpret_cast<float *>(m_data.data()) + sindex;
          for (int s = 0; s < nreal; ++s, ++ptr)
            *ptr = -(*ptr);
        }
      }
}
} // namespace LaphEnv
