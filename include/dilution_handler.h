#ifndef DILUTION_HANDLER_H
#define DILUTION_HANDLER_H

#include "dilution_scheme_info.h"
#include "field_smearing_info.h"

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
// *   NOTE: Dirac spin indices are zero-based here:  spin=0,1,2,3   *
// *   Elsewhere in the chroma_laph code, spin is unit-based.        *
// *                                                                 *
// *                                                                 *
// *******************************************************************

class DilutionHandler {

  DilutionSchemeInfo *dilPtr;
  int Textent; //, minTime, maxTime;
  int nEigvecs;
  int nSpinProjectors;
  int nEigvecProjectors;
  int nTimeProjectors;
  uint Nspin;

  std::vector<std::list<int>> spinProjs; // holds the spin dilution projectors
  std::vector<std::list<int>>
      eigvecProjs;                       // holds the eigvec dilution projectors
  std::vector<std::list<int>> timeProjs; // holds the time dilution projectors

  std::vector<int> spin_proj_indices;
  std::vector<int> eigvec_proj_indices;

  std::vector<int> which_time_proj;
  std::vector<int> which_spin_proj;
  std::vector<int> which_eigvec_proj;

  // prevent copying

  DilutionHandler(const DilutionHandler &in);
  DilutionHandler &operator=(const DilutionHandler &in);

public:
  DilutionHandler();

  DilutionHandler(const DilutionSchemeInfo &dilScheme,
                  const QuarkSmearingInfo &qSmear,
                  bool UpperSpinComponentsOnly = false);

  void setInfo(const DilutionSchemeInfo &dilScheme,
               const QuarkSmearingInfo &qSmear,
               bool UpperSpinComponentsOnly = false);

  ~DilutionHandler();

  void clear();

  // access to the info

  bool isInfoSet() const;

  const DilutionSchemeInfo &getDilutionSchemeInfo() const;

  int getNumberOfSpinEigvecProjectors() const;

  int getNumberOfTimeProjectors() const;

  int getNumberOfSpinProjectors() const;

  int getNumberOfEigvecProjectors() const;

  int getSpinProjectorIndex(int spineigvec_index) const;

  int getEigvecProjectorIndex(int spineigvec_index) const;

  int getTimeProjectorIndex(int time_val) const;

  const std::list<int> &getOnSpinIndices(int spineigvec_index) const;

  const std::list<int> &getOnEigvecIndices(int spineigvec_index) const;

  const std::list<int> &getSourceOnEigvecIndices(int eigvec_index) const;

  const std::list<int> &getOnTimes(int time_proj_index) const;

  bool isOnSpin(int spineigvec_index, int spin_val) const;

  bool isOnEigvec(int spineigvec_index, int eigvec_index) const;

  bool isOnTime(int time_proj_index, int time_val) const;

  bool isValidTimeProjectorIndex(int time_proj_index) const;

  bool isValidSpinEigvecProjectorIndex(int spineigvec_index) const;

  bool isFullTimeDilution() const;

private:
  void set_info(const DilutionSchemeInfo &dilScheme,
                const QuarkSmearingInfo &qSmear, bool UpperSpinComponentsOnly);

  void setProjectorMasks(std::vector<std::list<int>> &projs, int dil_type,
                         int nBasis, std::vector<int> &projind, int &nproj);

  void check_info_set(const std::string &name) const;

  void check_valid_spineig(int spineigvec_index) const;

  void check_valid_timeproj(int timeprojindex) const;
};

// **************************************************
} // namespace LaphEnv
#endif
