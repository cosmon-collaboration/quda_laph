#include "field_smearing_info.h"

using namespace std;

namespace LaphEnv {

// *************************************************************

// XMLHandler constructor

GluonSmearingInfo::GluonSmearingInfo(const XMLHandler &xmlin) {
  XMLHandler xml_in(xmlin);
  xml_tag_assert(xml_in, "GluonStoutSmearingInfo", "GluonSmearingInfo");
  XMLHandler xmlr(xml_in, "GluonStoutSmearingInfo");
  extract_info_from_reader(xmlr);
}

void GluonSmearingInfo::extract_info_from_reader(XMLHandler &xml_in) {
  xmlread(xml_in, "LinkIterations", linkIterations, "GluonSmearingInfo");
  xmlread(xml_in, "LinkStapleWeight", linkStapleWeight, "GluonSmearingInfo");
  if ((linkIterations < 0) || (linkStapleWeight < 0.0)) {
    xmlreadfail(xml_in, "GluonSmearingInfo",
                "invalid smearing scheme parameters in GluonSmearingInfo");
  }
}

// *************************************************************

// copy constructor

GluonSmearingInfo::GluonSmearingInfo(const GluonSmearingInfo &in)
    : linkIterations(in.linkIterations), linkStapleWeight(in.linkStapleWeight) {
}

GluonSmearingInfo &GluonSmearingInfo::operator=(const GluonSmearingInfo &in) {
  linkIterations = in.linkIterations;
  linkStapleWeight = in.linkStapleWeight;
  return *this;
}

void GluonSmearingInfo::checkEqual(const GluonSmearingInfo &in) const {
  if ((linkIterations != in.linkIterations) ||
      (abs(linkStapleWeight - in.linkStapleWeight) > 1e-12)) {
    std::cerr << "GluonSmearingInfo checkEqual failed" << std::endl;
    std::cerr << "LHS:" << std::endl
              << output() << std::endl
              << "RHS:" << std::endl
              << in.output() << std::endl;
    throw(std::invalid_argument("GluonSmearingInfo checkEqual failed..."));
  }
}

bool GluonSmearingInfo::operator==(const GluonSmearingInfo &in) const {
  return ((linkIterations == in.linkIterations) &&
          (abs(linkStapleWeight - in.linkStapleWeight) < 1e-12));
}

string GluonSmearingInfo::output(int indent) const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

void GluonSmearingInfo::output(XMLHandler &xmlout) const {
  xmlout.set_root("GluonStoutSmearingInfo");
  xmlout.put_child("LinkIterations", make_string(linkIterations));
  xmlout.put_child("LinkStapleWeight", make_string(linkStapleWeight));
}

void GluonSmearingInfo::setQudaGaugeSmearParam(
    QudaGaugeSmearParam &gauge_smear_param) const {
  gauge_smear_param.smear_type = QUDA_GAUGE_SMEAR_STOUT;
  gauge_smear_param.n_steps = linkIterations;
  gauge_smear_param.rho = linkStapleWeight;
  gauge_smear_param.meas_interval = std::max(linkIterations, 1);
  gauge_smear_param.struct_size = sizeof(gauge_smear_param);
}

// *************************************************************

// XMLHandler constructor

QuarkSmearingInfo::QuarkSmearingInfo(const XMLHandler &xmlin) {
  XMLHandler xml_in(xmlin);
  xml_tag_assert(xml_in, "QuarkLaphSmearingInfo", "QuarkSmearingInfo");
  XMLHandler xmlr(xml_in, "QuarkLaphSmearingInfo");
  extract_info_from_reader(xmlr);
}

void QuarkSmearingInfo::extract_info_from_reader(XMLHandler &xml_in) {
  xmlread(xml_in, "LaphSigmaCutoff", laphSigma, "QuarkSmearingInfo");
  xmlread(xml_in, "NumberLaphEigvecs", laphNumEigvecs, "QuarkSmearingInfo");
  if ((laphNumEigvecs < 1) || (laphSigma <= 0)) {
    throw(std::invalid_argument(
        "error: invalid smearing scheme parameters in QuarkSmearingInfo"));
  }
}

// ************************************************************

// copy constructor

QuarkSmearingInfo::QuarkSmearingInfo(const QuarkSmearingInfo &in)
    : laphNumEigvecs(in.laphNumEigvecs), laphSigma(in.laphSigma) {}

QuarkSmearingInfo &QuarkSmearingInfo::operator=(const QuarkSmearingInfo &in) {
  laphNumEigvecs = in.laphNumEigvecs;
  laphSigma = in.laphSigma;
  return *this;
}

void QuarkSmearingInfo::increaseUpdate(const QuarkSmearingInfo &in) {
  if (in.laphNumEigvecs > laphNumEigvecs) {
    laphNumEigvecs = in.laphNumEigvecs;
    laphSigma = in.laphSigma;
  }
}

void QuarkSmearingInfo::checkEqual(const QuarkSmearingInfo &in) const {
  if ((abs(laphSigma - in.laphSigma) > 1e-12) ||
      (laphNumEigvecs != in.laphNumEigvecs)) {
    std::cerr << "QuarkSmearingInfo checkEqual failed" << endl;
    std::cerr << "LHS:" << endl
              << output() << endl
              << "RHS:" << endl
              << in.output() << endl;
    throw(std::invalid_argument("QuarkSmearingInfo checkEqual failed..."));
  }
}

bool QuarkSmearingInfo::operator==(const QuarkSmearingInfo &in) const {
  return ((abs(laphSigma - in.laphSigma) < 1e-12) &&
          (laphNumEigvecs == in.laphNumEigvecs));
}

void QuarkSmearingInfo::checkOK(const QuarkSmearingInfo &in) const {
  if (laphNumEigvecs > in.laphNumEigvecs) {
    std::cerr << "QuarkSmearingInfo checkOK failed" << endl;
    std::cerr << "LHS:" << endl
              << output() << endl
              << "RHS:" << endl
              << in.output() << endl;
    throw(std::invalid_argument("QuarkSmearingInfo checkOK failed..."));
  }
}

string QuarkSmearingInfo::output(int indent) const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

void QuarkSmearingInfo::output(XMLHandler &xmlout) const {
  xmlout.set_root("QuarkLaphSmearingInfo");
  xmlout.put_child("LaphSigmaCutoff", make_string(laphSigma));
  xmlout.put_child("NumberLaphEigvecs", make_string(laphNumEigvecs));
}

// *************************************************************
} // namespace LaphEnv
