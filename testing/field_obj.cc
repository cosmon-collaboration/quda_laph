#include "field_obj.h"
#include "layout_info.h"
#include <complex>

using namespace std;


namespace LaphEnv {

// *********************************************************************************


LatticeFieldObj::LatticeFieldObj(SiteType sitetype, void *data, bool double_precision)
                : m_sitetype(sitetype), m_data_ptr(data), m_dbleprec(double_precision)
{
 if (m_sitetype==Complex) m_site_elems=1;
 else if (m_sitetype==ColorMatrix) m_site_elems=Ncolor*Ncolor;
 else if (m_sitetype==ColorVector) m_site_elems=Ncolor;
 else if (m_sitetype==ColorSpinVector) m_site_elems=Nspin*Ncolor;
}


const int LatticeFieldObj::Nspin=4;

const int LatticeFieldObj::Ncolor=3;


size_t LatticeFieldObj::bytesPerSite() const
{
 return (m_dbleprec ? sizeof(complex<double>) : sizeof(complex<float>))*m_site_elems;
}


size_t LatticeFieldObj::bytesPerWord() const
{
 return (m_dbleprec ? sizeof(complex<double>) : sizeof(complex<float>));
}




// *********************************************************************************
}
