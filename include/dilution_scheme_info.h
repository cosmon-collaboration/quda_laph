#ifndef LAPH_DILUTION_SCHEME_H
#define LAPH_DILUTION_SCHEME_H

#include "xml_handler.h"

namespace LaphEnv {

// *******************************************************************
// *                                                                 *
// *  ZN noise "rho" is introduced depending on time, spin, and      *
// *  Laplacian eigenvector index:                                   *
// *                                                                 *
// *           rho[ time_index, spin_index, laph_index ]             *
// *                                                                 *
// *  Dilution projectors are matrices in this space, but we choose  *
// *  matrices that are outer products of time dilutions, spin       *
// *  dilutions, and Laph eigenvector index dilutions:               *
// *                                                                 *
// *   P[ t',s',l' ; t,s,l ] = P_t[ t'; t] P_s[ s'; s ] P_l[ l'; l ] *
// *                                                                 *
// *  Each projector can have one of the following forms (N is the   *
// *  number of values that the indices i,j can take, and A denotes  *
// *  the different projectors in each scheme):                      *
// *                                                                 *
// *            No dilution:      A = 0                              *
// *                                                                 *
// *                P(A)[i,j] = delta(i,j)                           *
// *                                                                 *
// *          Full dilution:      A = 0,1,...,N-1                    *
// *                                                                 *
// *                P(A)[i,j] = delta(i,j) delta(A,i)                *
// *                                                                 *
// *       Block-K dilution:      A = 0,1,...,K-1                    *
// *                                                                 *
// *                P(A)[i,j] = delta(i,j) delta(A, floor(K*i/N))    *
// *                                                                 *
// *   Interlace-K dilution:      A = 0,1,...,K-1                    *
// *                                                                 *
// *                P(A)[i,j] = delta(i,j) delta(A, i mod K)         *
// *                                                                 *
// *   Note that all off-diagonal elements are always zero in the    *
// *   above projectors, and in each case, some of the diagonal      *
// *   elements are unity while the rest are zero.  Given this       *
// *   special form, a projector of index A can be fully described   *
// *   by a list of the diagonal entries that are "on" or set to     *
// *   unity.  This is the representation that is used below.        *
// *                                                                 *
// *                                                                 *
// *   Objects of class "DilutionSchemeInfo" store identifying info  *
// *   about a Laph dilution scheme, which involves time dilution,   *
// *   spin dilution, and Laph eigenvector dilution.  XML input      *
// *   format for a Laph dilution scheme is:                         *
// *                                                                 *
// *     <LaphDilutionScheme>                                        *
// *        <TimeDilution>                                           *
// *           <DilutionType> full </DilutionType>                   *
// *        </TimeDilution>                                          *
// *        <SpinDilution>                                           *
// *           <DilutionType> none </DilutionType>                   *
// *        </SpinDilution>                                          *
// *        <EigvecDilution>                                         *
// *           <DilutionType> block </DilutionType>                  *
// *           <NumberProjectors> 4 </NumberProjectors>              *
// *        </EigvecDilution>                                        *
// *     </LaphDilutionScheme>                                       *
// *                                                                 *
// *   Spin, Laph eigenvector, and time dilution can be of type      *
// *        "full", "none", "block", or "interlace".                 *
// *   If "block" or "interlace" is specified, a tag                 *
// *   specifying the number of projectors must be included          *
// *                                                                 *
// *                                                                 *
// *   The members most likely to be used are                        *
// *                                                                 *
// *     XMLHandler xlm_in(...);                                     *
// *     DilutionSchemeInfo dil(xml_in);  // construct from xml      *
// *                                                                 *
// *     DilutionSchemeInfo dil2(...);                               *
// *     dil.checkEqual(dil2);  // if dil != dil2 throw exception    *
// *                                                                 *
// *     dil.match(dil2);  // source-sink match up                   *
// *                                                                 *
// *     string xml_out = dil.output();     // xml output            *
// *                                                                 *
// *                                                                 *
// *******************************************************************

class DilutionSchemeInfo {

  int spinDilutionType;   //  0 = none, 1 = full
  int eigvecDilutionType; //  x  (x>=2) = block with no. projectors x
  int timeDilutionType;   // -y  (y>=2) = interlace with no. projectors y

public:
  DilutionSchemeInfo(const XMLHandler &xml_in);

  DilutionSchemeInfo(); // constructor for no dilution

  DilutionSchemeInfo(const DilutionSchemeInfo &in);

  DilutionSchemeInfo &operator=(const DilutionSchemeInfo &in);

  ~DilutionSchemeInfo() {}

  bool isFullTimeDilution() const { return (timeDilutionType == 1); }

  void checkEqual(const DilutionSchemeInfo &in) const;

  bool operator==(const DilutionSchemeInfo &in) const;

  bool operator!=(const DilutionSchemeInfo &in) const;

  bool operator<(const DilutionSchemeInfo &in) const;

  // output functions

  std::string output(int indent = 0) const; // XML output

  void output(XMLHandler &xmlout) const; // write header

  //  return true if current dilution scheme can be changed
  //  (undiluted) to scheme "newscheme"

  bool canUndilute(const DilutionSchemeInfo &newscheme) const;

private:
  void assign(int spin_dil_type, int eigvec_dil_type, int time_dil_type = 1);
  void assign_from_reader(XMLHandler &xml_in);
  void dil_in(XMLHandler &xml_in, const std::string &path, int &DilType);
  std::string dil_out(int indent, int DilType, bool out_nproj = false) const;
  void dil_out(XMLHandler &xmlout, int DilType, bool out_nproj = false) const;
  bool can_undilute(int dilorig, int dilnew) const;

  friend class DilutionHandler;
};

// **************************************************
} // namespace LaphEnv
#endif
