#include "laph_noise_info.h"
#include "laph_stdio.h"

using namespace std;

namespace LaphEnv {

LaphNoiseInfo::LaphNoiseInfo() : store(4) {}

// XMLHandler constructor

LaphNoiseInfo::LaphNoiseInfo(const XMLHandler &xml_in) {
  XMLHandler xmlr(xml_in);
  xml_tag_assert(xmlr, "LaphNoiseInfo");
  extract_info_from_reader(xmlr);
}

LaphNoiseInfo::LaphNoiseInfo(int znGroup, int seed) {
  try {
    check_assignment(znGroup, seed);
    encode(znGroup, seed);
  } catch (const std::exception &msg) {
    throw(std::runtime_error("Invalid LaphNoiseInfo construction"));
  }
}

void LaphNoiseInfo::extract_info_from_reader(XMLHandler &xml_in) {
  int znGroup, seed;
  xmlread(xml_in, "ZNGroup", znGroup, "LaphNoiseInfo");
  xmlread(xml_in, "Seed", seed, "LaphNoiseInfo");
  try {
    check_assignment(znGroup, seed);
    encode(znGroup, seed);
  } catch (const std::exception &msg) {
    xmlreadfail(xml_in, "LaphNoiseInfo");
  }
}

void LaphNoiseInfo::check_assignment(int znGroup, int seed) {
  if ((znGroup != 4) && (znGroup != 8) && (znGroup != 32) && (znGroup != 1)) {
    printLaph("improper initialization of LaphNoiseInfo");
    printLaph("ZNGroup must have integer value 4, 8, or 32");
    throw(std::runtime_error("error"));
  }
  if (znGroup == 1) {
    printLaph("Warning: znGroup = 1 so debugging ONLY");
  }
  if ((seed < 0) || (seed > 65535)) {
    printLaph("improper initialization of LaphNoiseInfo");
    printLaph("seed must have value 0 to 65535");
    throw(std::runtime_error("error"));
  }
}

void LaphNoiseInfo::encode(int znGroup, int seed) {
  store = (unsigned int)seed;
  store <<= 6;
  store |= (unsigned int)znGroup;
}

// *************************************************************

// copy constructor

LaphNoiseInfo::LaphNoiseInfo(const LaphNoiseInfo &in) : store(in.store) {}

LaphNoiseInfo &LaphNoiseInfo::operator=(const LaphNoiseInfo &in) {
  store = in.store;
  return *this;
}

void LaphNoiseInfo::checkEqual(const LaphNoiseInfo &in) const {
  if (store != in.store)
    throw(std::invalid_argument("LaphNoiseInfo does not checkEqual"));
}

bool LaphNoiseInfo::operator==(const LaphNoiseInfo &in) const {
  return (store == in.store);
}

bool LaphNoiseInfo::operator!=(const LaphNoiseInfo &in) const {
  return (store != in.store);
}

bool LaphNoiseInfo::operator<(const LaphNoiseInfo &in) const {
  return (store < in.store);
}

// **********************************************************

unsigned long LaphNoiseInfo::getSeed(const GaugeConfigurationInfo &G) const {
  int traj_num = G.getConfigNumber();
  if ((traj_num < 0) || (traj_num > 65535)) {
    std::cerr << " Error in LaphNoiseInfo::getSeed; unsupported HMC"
              << " trajectory number...limited range 0..65535" << std::endl;
    throw(std::invalid_argument("error"));
  }
  unsigned long m_seed = getSeed();
  unsigned long currTraj = traj_num;
  unsigned long rngseed = m_seed & 0xFFUL; // 8 least sig bits of m_seed
  rngseed =
      (rngseed << 16) | (currTraj & 0xFF00UL);   // 8 most sig bits of currTraj
  rngseed = (rngseed << 8) | (m_seed & 0xFF00UL) // 8 most sig bits of m_seed
            | (currTraj & 0xFFUL);               // 8 least sig bits of currTraj
  return rngseed;
}

string LaphNoiseInfo::output(int indent) const {
  XMLHandler xmlh;
  output(xmlh);
  return xmlh.output(indent);
}

string LaphNoiseInfo::str() const {
  XMLHandler xmlh;
  output(xmlh);
  return xmlh.str();
}

void LaphNoiseInfo::output(XMLHandler &xmlout) const {
  xmlout.set_root("LaphNoiseInfo");
  xmlout.put_child("ZNGroup", make_string(getZNGroup()));
  xmlout.put_child("Seed", make_string(getSeed()));
}

// *************************************************************
} // namespace LaphEnv
