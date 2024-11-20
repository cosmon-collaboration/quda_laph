#include "quda_info.h"
#include "laph_stdio.h"
#include "latt_field.h"

using namespace std;

namespace LaphEnv {

//   Expects the following XML:
//    <QudaInfo>
//       <CPUPrecision>double</CPUPrecision> (or single)
//       <CUDAPrecision>double</CUDAPrecision> (or single)
//       <CUDASloppyPrecision>double</CUDASloppyPrecision> (or single)
//    </QudaInfo>

void QudaInfo::init(const XMLHandler &xml_in, const bool echo) {
  XMLHandler xmlr(xml_in, "QudaInfo");
  string response;
  if (xmlreadif(xmlr, "CPUPrecision", response, "QudaInfo")) {
    if (response == "single") {
      cpu_prec = QUDA_SINGLE_PRECISION;
      cpu_prec_bytes = 2 * sizeof(float);
    } else if (response == "double") {
      cpu_prec = QUDA_DOUBLE_PRECISION;
      cpu_prec_bytes = 2 * sizeof(double);
    } else {
      throw(std::invalid_argument("Unsupported CPUPrecision"));
    }
  }
  LattField::cpu_prec_bytes = cpu_prec_bytes;
  if (xmlreadif(xmlr, "CUDAPrecision", response, "QudaInfo")) {
    if (response == "single") {
      cuda_prec = QUDA_SINGLE_PRECISION;
    } else if (response == "double") {
      cuda_prec = QUDA_DOUBLE_PRECISION;
    } else {
      throw(std::invalid_argument("Unsupported CUDAPrecision"));
    }
  }
  if (xmlreadif(xmlr, "CUDASloppyPrecision", response, "QudaInfo")) {
    if (response == "single") {
      cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
    } else if (response == "double") {
      cuda_prec_sloppy = QUDA_DOUBLE_PRECISION;
    } else {
      throw(std::invalid_argument("Unsupported CUDAPrecision"));
    }
  }
  if (echo) {
    XMLHandler xmlout("QudaInfo");
    xmlout.put_child("CPUPrecision",
                     (cpu_prec == QUDA_DOUBLE_PRECISION) ? "double" : "single");
    xmlout.put_child("CUDAPrecision", (cuda_prec == QUDA_DOUBLE_PRECISION)
                                          ? "double"
                                          : "single");
    xmlout.put_child("CUDASloppyPrecision",
                     (cuda_prec_sloppy == QUDA_DOUBLE_PRECISION) ? "double"
                                                                 : "single");
    printLaph("Quda Precisions:");
    printLaph(xmlout.output());
  }
}

#ifdef ARCH_PARALLEL
int QudaInfo::device_ordinal = -1;
#else
int QudaInfo::device_ordinal = 0;
#endif

// default values

QudaPrecision QudaInfo::cpu_prec = QUDA_DOUBLE_PRECISION;

int QudaInfo::cpu_prec_bytes = 2 * sizeof(double);

QudaPrecision QudaInfo::cuda_prec = QUDA_DOUBLE_PRECISION;

QudaPrecision QudaInfo::cuda_prec_sloppy = QUDA_DOUBLE_PRECISION;

QudaReconstructType QudaInfo::link_recon = QUDA_RECONSTRUCT_NO;

QudaReconstructType QudaInfo::link_recon_sloppy = QUDA_RECONSTRUCT_NO;

bool QudaInfo::gauge_config_on_device = false;

bool QudaInfo::smeared_gauge_on_device = false;

bool QudaInfo::clover_on_device = false;

// remove all quantities from memory on the device

void QudaInfo::clearDevice() {
  if (gauge_config_on_device) {
    freeGaugeQuda();
    gauge_config_on_device = false;
  }
  if (smeared_gauge_on_device) {
    freeGaugeSmearedQuda();
    smeared_gauge_on_device = false;
  }
  if (clover_on_device) {
    freeCloverQuda();
    clover_on_device = false;
  }
}

void QudaInfo::clearDeviceGaugeConfiguration() {
  if (gauge_config_on_device) {
    freeGaugeQuda();
    gauge_config_on_device = false;
  }
}

void QudaInfo::clearDeviceSmearedGaugeConfiguration() {
  if (smeared_gauge_on_device) {
    freeGaugeSmearedQuda();
    smeared_gauge_on_device = false;
  }
}

void QudaInfo::clearDeviceCloverField() {
  if (clover_on_device) {
    freeCloverQuda();
    clover_on_device = false;
  }
}
} // namespace LaphEnv
