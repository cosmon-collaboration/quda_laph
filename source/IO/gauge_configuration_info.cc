#include "gauge_configuration_info.h"
#include "laph_stdio.h"
#include "layout_info.h"
#include "quda_info.h"
#include <sstream>

namespace LaphEnv {

GaugeConfigurationInfo::GaugeConfigurationInfo(const XMLHandler &xmlin) {
  XMLHandler xml_in(xmlin);
  set_info(xml_in);
}

void GaugeConfigurationInfo::set_info(XMLHandler &xml_in) {
  xml_tag_assert(xml_in, "GaugeConfigurationInfo");
  XMLHandler xmlg(xml_in, "GaugeConfigurationInfo");
  xmlread(xmlg, "EnsembleName", ensemble_name, "GaugeConfigurationInfo");
  xmlread(xmlg, "FileFormat", file_format, "GaugeConfigurationInfo");
  check_valid(file_format, "FileFormat",
              {"CERN", "CLS", "SZIN_SP", "SZIN_DP", "USQCD_SP", "USQCD_DP"});
  config_type = "WilsonImproved";
  xmlreadif(xmlg, "ConfigType", config_type, "GaugeConfigurationInfo");
  check_valid(config_type, "ConfigType", {"Wilson", "WilsonImproved"});
  xmlread(xmlg, "FileName", file_name, "GaugeConfigurationInfo");
  xmlread(xmlg, "ConfigNumber", config_num, "GaugeConfigurationInfo");
  mc_chain_num = 0;
  xmlreadif(xmlg, "MarkovChainNumber", mc_chain_num, "GaugeConfigurationInfo");
  nom_id = "default_gauge_id";
  xmlreadif(xmlg, "NOMId", nom_id, "GaugeConfigurationInfo");
  gluon_anisotropy = 1.0;
  xmlreadif(xmlg, "GluonAnisotropy", gluon_anisotropy,
            "GaugeConfigurationInfo");
  std::string response = "antiperiodic";
  xmlreadif(xmlg, "FermionTimeBC", response, "GaugeConfigurationInfo");
  check_valid(response, "FermionTimeBC", {"antiperiodic", "periodic"});
  fermion_time_bc = (response == "antiperiodic") ? 'A' : 'P';
}

void GaugeConfigurationInfo::check_valid(const std::string &sdata,
					 const std::string &tag,
                                         const std::vector<std::string> &allowed) {
  bool ok = false;
  for (std::vector<std::string>::const_iterator it = allowed.begin(); it != allowed.end();
       ++it) {
    if (sdata == (*it)) {
      ok = true;
      break;
    }
  }
  if (!ok) {
    std::string mesg("Invalid <");
    mesg += tag;
    mesg += "> ";
    mesg += sdata;
    mesg += " in GaugeConfigurationInfo";
    throw(std::invalid_argument(mesg));
  }
}

// copy constructor

GaugeConfigurationInfo::GaugeConfigurationInfo(
    const GaugeConfigurationInfo &rhs)
    : ensemble_name(rhs.ensemble_name), file_format(rhs.file_format),
      config_type(rhs.config_type), file_name(rhs.file_name),
      config_num(rhs.config_num), mc_chain_num(rhs.mc_chain_num),
      nom_id(rhs.nom_id), gluon_anisotropy(rhs.gluon_anisotropy),
      fermion_time_bc(rhs.fermion_time_bc) {}

GaugeConfigurationInfo &
GaugeConfigurationInfo::operator=(const GaugeConfigurationInfo &rhs) {
  ensemble_name = rhs.ensemble_name;
  file_format = rhs.file_format;
  config_type = rhs.config_type;
  file_name = rhs.file_name;
  config_num = rhs.config_num;
  mc_chain_num = rhs.mc_chain_num;
  nom_id = rhs.nom_id;
  gluon_anisotropy = rhs.gluon_anisotropy;
  fermion_time_bc = rhs.fermion_time_bc;
  return *this;
}

void GaugeConfigurationInfo::checkEqual(
    const GaugeConfigurationInfo &rhs) const {
  if (!((*this) == rhs)) {
    printLaph("GaugeConfigurationInfo checkEqual failed");
    printLaph(make_strf("LHS:\n%s\nRHS:\n%s\n", output(), rhs.output()));
    throw(std::invalid_argument("GaugeConfigurationInfo checkEqual failed"));
  }
}

bool GaugeConfigurationInfo::operator==(
    const GaugeConfigurationInfo &rhs) const {
  return ((ensemble_name == rhs.ensemble_name) &&
          (config_type == rhs.config_type) && (config_num == rhs.config_num) &&
          (mc_chain_num == rhs.mc_chain_num) &&
          (std::abs(gluon_anisotropy - rhs.gluon_anisotropy) < 1e-12) &&
          (fermion_time_bc == rhs.fermion_time_bc));
}

std::string GaugeConfigurationInfo::output(const int indent) const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

void GaugeConfigurationInfo::output(XMLHandler &xmlout) const {
  xmlout.set_root("GaugeConfigurationInfo");
  xmlout.put_child("EnsembleName", ensemble_name);
  xmlout.put_child("FileFormat", file_format);
  xmlout.put_child("ConfigType", config_type);
  xmlout.put_child("FileName", file_name);
  xmlout.put_child("ConfigNumber", make_string(config_num));
  xmlout.put_child("MarkovChainNumber", make_string(mc_chain_num));
  xmlout.put_child("NOMId", nom_id);
  xmlout.put_child("GluonAnisotropy", make_string(gluon_anisotropy));
  std::string fermtimebc((fermion_time_bc == 'A') ? "antiperiodic" : "periodic");
  xmlout.put_child("FermionTimeBC", fermtimebc);
}

void GaugeConfigurationInfo::setQudaGaugeParam(
    QudaGaugeParam &gauge_param) const {
  gauge_param.type = QUDA_WILSON_LINKS;

  gauge_param.X[0] = LayoutInfo::getRankLattExtents()[0];
  gauge_param.X[1] = LayoutInfo::getRankLattExtents()[1];
  gauge_param.X[2] = LayoutInfo::getRankLattExtents()[2];
  gauge_param.X[3] = LayoutInfo::getRankLattExtents()[3];

  gauge_param.cpu_prec = QudaInfo::cpu_prec;
  gauge_param.cuda_prec = QudaInfo::cuda_prec;
  gauge_param.cuda_prec_sloppy = QudaInfo::cuda_prec_sloppy;
  gauge_param.cuda_prec_precondition = QudaInfo::cuda_prec;
  gauge_param.cuda_prec_eigensolver = QudaInfo::cuda_prec;

  gauge_param.reconstruct = QudaInfo::link_recon;
  gauge_param.reconstruct_sloppy = QudaInfo::link_recon;
  gauge_param.reconstruct_precondition = QudaInfo::link_recon;
  gauge_param.reconstruct_eigensolver = QudaInfo::link_recon;
  gauge_param.reconstruct_refinement_sloppy = QudaInfo::link_recon_sloppy;

  gauge_param.anisotropy = gluon_anisotropy;
  gauge_param.tadpole_coeff = 1.0; // used only by staggered

  gauge_param.ga_pad = 0;
  gauge_param.mom_ga_pad = 0;
  gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;

  gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER;
  gauge_param.t_boundary =
      (fermion_time_bc == 'A') ? QUDA_ANTI_PERIODIC_T : QUDA_PERIODIC_T;

  gauge_param.cuda_prec_sloppy = QudaInfo::cuda_prec_sloppy;
  gauge_param.cuda_prec_precondition = QudaInfo::cuda_prec_sloppy;
  gauge_param.cuda_prec_eigensolver = QudaInfo::cuda_prec_sloppy;
  gauge_param.cuda_prec_refinement_sloppy = QudaInfo::cuda_prec_sloppy;

  gauge_param.reconstruct_sloppy = QudaInfo::link_recon_sloppy;
  gauge_param.reconstruct_precondition = QudaInfo::link_recon_sloppy;
  gauge_param.reconstruct_eigensolver = QudaInfo::link_recon_sloppy;
  gauge_param.reconstruct_refinement_sloppy = QudaInfo::link_recon_sloppy;

  int pad_size = 0;
  // For multi-GPU, ga_pad must be large enough to store a time-slice
#ifdef ARCH_PARALLEL
  const int x_face_size = gauge_param.X[1] * gauge_param.X[2] * gauge_param.X[3] / 2;
  const int y_face_size = gauge_param.X[0] * gauge_param.X[2] * gauge_param.X[3] / 2;
  const int z_face_size = gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[3] / 2;
  const int t_face_size = gauge_param.X[0] * gauge_param.X[1] * gauge_param.X[2] / 2;
  pad_size = std::max(std::max(x_face_size, y_face_size),
                      std::max(z_face_size, t_face_size));
#endif
  gauge_param.ga_pad = pad_size;

  gauge_param.struct_size = sizeof(gauge_param);
}
} // namespace LaphEnv
