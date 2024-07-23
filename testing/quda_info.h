#ifndef QUDA_PARAMS_H
#define QUDA_PARAMS_H

#include "xml_handler.h"
#include "quda.h"

namespace LaphEnv {

// ***********************************************************************
// *                                                                     *
// *   Class "QudaInfo" is a singleton which stores and makes            *
// *   available global information needed by Quda. Various              *
// *   quantities, such as the gauge field, quark sinks, will need       *
// *   complicated parameter set ups.  The data members in this          *
// *   class will be necessary to initialize such structures.            *
// *   Some of the information is set by default and cannot be changed.  *
// *   The class also tries to keep track of whether various quantities  *
// *   are currently on the device or not.  At the end of program        *
// *   execution, this class clears the device.                          *
// *                                                                     *
// *   The input XML should have the following form (all are             *
// *   optional since double precision is default):                      *
// *                                                                     *
// *   <QudaInfo>                                                        *
// *      <CPUPrecision>double</CPUPrecision> (or single)                *
// *      <CUDAPrecision>double</CUDAPrecision> (or single)              *
// *      <CUDASloppyPrecision>double</CUDASloppyPrecision> (or single)  *
// *   </QudaInfo>                                                       *
// *                                                                     *
// ***********************************************************************


class QudaInfo
{

   QudaInfo() {}

   ~QudaInfo()
    {clearDevice();}

   static bool gauge_config_on_device;

   static bool smeared_gauge_on_device;

   static bool clover_on_device;
   
   static bool mg_preconditioner;

   static int device_ordinal;

     // data members

   static QudaPrecision cpu_prec;
   
   static int cpu_prec_bytes;  // bytes in a complex number with precision cpu_prec

   static QudaPrecision cuda_prec;

   static QudaPrecision cuda_prec_sloppy;

   static QudaReconstructType link_recon;

   static QudaReconstructType link_recon_sloppy;


 public:

   QudaInfo(const QudaInfo &in) = delete;           // no copy constructor

   QudaInfo(QudaInfo &in) = delete;                 // no copy constructor

   QudaInfo& operator=(const QudaInfo &in) = delete; // not assignable

   QudaInfo& operator=(QudaInfo &in) = delete;       // not assignable
    
   static void init(const XMLHandler& xml_in, bool echo=true);
   
   static int getDeviceOrdinal() {return device_ordinal;}
 
   static void clearDevice();

   static QudaPrecision get_cpu_prec() {return cpu_prec;}
   
   static int get_cpu_prec_bytes() {return cpu_prec_bytes;} 

   static QudaPrecision get_cuda_prec() {return cuda_prec;}

   static QudaPrecision get_cuda_prec_sloppy() {return cuda_prec_sloppy;}

   static QudaReconstructType get_link_recon() {return link_recon;}

   static QudaReconstructType get_link_recon_sloppy() {return link_recon_sloppy;}


   friend class GaugeConfigurationInfo;
   friend class GaugeConfigurationHandler;
   friend class GluonSmearingHandler;
   friend class QuarkHandler;

};


//extern QudaReconstructType link_recon_sloppy;
//extern QudaReconstructType link_recon_precondition;
//extern QudaReconstructType link_recon_eigensolver;
//extern QudaPrecision prec_sloppy;
//extern QudaPrecision prec_refinement_sloppy;
//extern QudaPrecision prec_precondition;
//extern QudaPrecision prec_eigensolver;
//extern QudaPrecision prec_null;
//extern QudaPrecision prec_ritz;
//extern QudaVerbosity verbosity;

// *************************************************************
}
#endif
