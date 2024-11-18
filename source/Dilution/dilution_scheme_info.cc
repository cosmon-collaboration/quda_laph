#include "dilution_scheme_info.h"
#include "laph_stdio.h"
#include "multi_compare.h"

using namespace std;

namespace LaphEnv {

DilutionSchemeInfo::DilutionSchemeInfo(const XMLHandler &xmlin) {
  XMLHandler xml_in(xmlin);
  xml_tag_assert(xml_in, "LaphDilutionScheme", "DilutionSchemeInfo");
  XMLHandler xmlr(xml_in, "LaphDilutionScheme");
  assign_from_reader(xmlr);
}

void DilutionSchemeInfo::assign_from_reader(XMLHandler &xmlr) {
  int time_dil_type = 0, spin_dil_type = 0, eigvec_dil_type = 0;
  try {
    dil_in(xmlr, "TimeDilution", time_dil_type);
    dil_in(xmlr, "SpinDilution", spin_dil_type);
    dil_in(xmlr, "EigvecDilution", eigvec_dil_type);
  } catch (const std::exception &err) {
    xmlreadfail(xmlr, "DilutionSchemeInfo",
                "could not initialize DilutionSchemeInfo from XML input");
  }
  try {
    assign(spin_dil_type, eigvec_dil_type, time_dil_type);
  } catch (const std::exception &xp) {
    errorLaph("abort");
  }
}

DilutionSchemeInfo::DilutionSchemeInfo() {
  spinDilutionType = 0;
  eigvecDilutionType = 0;
  timeDilutionType = 0;
}

// copy constructor

DilutionSchemeInfo::DilutionSchemeInfo(const DilutionSchemeInfo &in)
    : spinDilutionType(in.spinDilutionType),
      eigvecDilutionType(in.eigvecDilutionType),
      timeDilutionType(in.timeDilutionType) {}

DilutionSchemeInfo &
DilutionSchemeInfo::operator=(const DilutionSchemeInfo &in) {
  spinDilutionType = in.spinDilutionType;
  eigvecDilutionType = in.eigvecDilutionType;
  timeDilutionType = in.timeDilutionType;
  return *this;
}

void DilutionSchemeInfo::assign(const int spin_dil_type,
				const int eigvec_dil_type,
                                const int time_dil_type) {
  try {
    if ((time_dil_type == -1) || (spin_dil_type == -1) ||
        (eigvec_dil_type == -1))
      throw(std::invalid_argument("dilution types cannot have value -1"));
    if ((spin_dil_type > 2) || (spin_dil_type < -2))
      throw(std::invalid_argument(
          "spin dilution type must have value -2, 0, 1, 2"));
  } catch (const std::exception &errmsg) {
    std::cerr << "Invalid DilutionSchemeInfo assigment:" << std::endl;
    std::cerr << "   ..." << errmsg.what() << std::endl;
    throw(std::invalid_argument("error"));
  }

  timeDilutionType = time_dil_type;
  spinDilutionType = spin_dil_type;
  eigvecDilutionType = eigvec_dil_type;
}

void DilutionSchemeInfo::checkEqual(const DilutionSchemeInfo &in) const {
  if ((spinDilutionType != in.spinDilutionType) ||
      (eigvecDilutionType != in.eigvecDilutionType) ||
      (timeDilutionType != in.timeDilutionType)) {
    std::cerr << "DilutionSchemeInfo checkEqual failed" << std::endl;
    std::cerr << "LHS:" << std::endl
              << output() << std::endl
              << "RHS:" << std::endl
              << in.output() << std::endl;
    throw(std::invalid_argument(
        "DilutionSchemeInfo does not checkEqual...abort"));
  }
}

bool DilutionSchemeInfo::operator==(const DilutionSchemeInfo &in) const {
  return multiEqual(spinDilutionType, in.spinDilutionType, eigvecDilutionType,
                    in.eigvecDilutionType, timeDilutionType,
                    in.timeDilutionType);
}

bool DilutionSchemeInfo::operator!=(const DilutionSchemeInfo &in) const {
  return multiNotEqual(spinDilutionType, in.spinDilutionType,
                       eigvecDilutionType, in.eigvecDilutionType,
                       timeDilutionType, in.timeDilutionType);
}

bool DilutionSchemeInfo::operator<(const DilutionSchemeInfo &in) const {
  return multiLessThan(spinDilutionType, in.spinDilutionType,
                       eigvecDilutionType, in.eigvecDilutionType,
                       timeDilutionType, in.timeDilutionType);
}

string DilutionSchemeInfo::output(int indent) const {
  XMLHandler xmlh;
  output(xmlh);
  return xmlh.output(indent);
}

void DilutionSchemeInfo::output(XMLHandler &xmlout) const {
  xmlout.set_root("LaphDilutionScheme");
  xmlout.put_child("TimeDilution");
  xmlout.put_child("SpinDilution");
  xmlout.put_child("EigvecDilution");
  xmlout.seek_root();
  xmlout.seek_first_child();
  dil_out(xmlout, timeDilutionType, true);
  xmlout.seek_next_sibling();
  dil_out(xmlout, spinDilutionType, false);
  xmlout.seek_next_sibling();
  dil_out(xmlout, eigvecDilutionType, true);
}

//  return true if current dilution scheme can be changed
//  (undiluted) to scheme "newscheme"

bool DilutionSchemeInfo::canUndilute(
    const DilutionSchemeInfo &newscheme) const {
  return ((can_undilute(spinDilutionType, newscheme.spinDilutionType)) &&
          (can_undilute(eigvecDilutionType, newscheme.eigvecDilutionType)) &&
          (can_undilute(timeDilutionType, newscheme.timeDilutionType)));
}

void DilutionSchemeInfo::dil_in(XMLHandler &xml_in, const string &tagname,
                                int &DilType) {
  string tmp;
  DilType = 0;
  try {
    XMLHandler xmlr(xml_in, tagname);
    xmlread(xmlr, "DilutionType", tmp, "DilutionSchemeInfo");
    tmp = tidyString(tmp);
    if (tmp == "none") {
      DilType = 0;
    } else if (tmp == "full") {
      DilType = 1;
    } else if (tmp == "block") {
      DilType = 2;
      int nproj = 0;
      if (xmlr.count("NumberProjectors") == 1) {
        xmlread(xmlr, "NumberProjectors", nproj, "DilutionSchemeInfo");
      }
      if (nproj > 1)
        DilType = nproj;
      else
        throw(std::runtime_error("invalid number of block projectors"));
    } else if (tmp == "interlace") {
      DilType = -2;
      int nproj = 0;
      if (xmlr.count("NumberProjectors") == 1) {
        xmlread(xmlr, "NumberProjectors", nproj, "DilutionSchemeInfo");
      }
      if (nproj > 1)
        DilType = -nproj;
      else
        throw(std::runtime_error("invalid number of interlace projectors"));
    } else {
      throw(std::runtime_error("invalid read"));
    }
  } catch (const std::exception &errstr) {
    std::cerr << "invalid DilutionSchemeInfo read from XML" << errstr.what()
              << endl;
    throw(std::runtime_error("error"));
  }
}

void DilutionSchemeInfo::dil_out(XMLHandler &xmlout, int DilType,
                                 bool out_nproj) const {
  string dtype;
  if (DilType == 0)
    dtype = "none";
  else if (DilType == 1)
    dtype = "full";
  else if (DilType > 1)
    dtype = "block";
  else if (DilType < -1)
    dtype = "interlace";
  XMLHandler xmld;
  xmld.set_root("DilutionType", dtype);
  if (out_nproj) {
    if (DilType > 1)
      xmlout.put_child("NumberProjectors", make_string(DilType));
    else if (DilType < -1)
      xmlout.put_child("NumberProjectors", make_string(-DilType));
  }
  xmlout.put_child(xmld);
}

bool DilutionSchemeInfo::can_undilute(int dilorig, int dilnew) const {
  if (dilorig == 1)
    return true;
  if (dilnew == 1)
    return false;
  if (dilorig == 0)
    return (dilnew == 0);
  if ((dilorig >= 2) && (dilnew < 0))
    return false;
  if ((dilorig <= -2) && (dilnew > 0))
    return false;
  const int dorig = (dilorig >= 2) ? dilorig : -dilorig;
  const int dnew = (dilnew >= 2) ? dilnew : -dilnew;
  if (dorig >= 2) {
    if (dnew == 0)
      return true;
    else if (dnew == 1)
      return false;
    if ((dorig % dnew) == 0)
      return true;
    return false;
  }
  return false;
}
} // namespace LaphEnv
