#ifndef MULTI_LOOPER_H
#define MULTI_LOOPER_H
#include <list>
#include <set>
#include <vector>
#include <complex>
#include <string>
#include <utility>
#include <stdexcept>

// *******************************************************************************
// *                                                                             *
// *    The classes "MultiIntLooper" and "MultiIntLooperWR" (with restrictions)  *
// *    are defined in this file.  These are used in the CorrelatorHandler       *
// *    compute routines to loop over dilution indices and noise choices.        *
// *                                                                             *
// *    The following utility routines are also defined here:                    *
// *          "make_list",  "make_vector",   "cmplx_same"                        *
// *                                                                             *
// *******************************************************************************

  //   Simple routine to help make a list of objects

template <typename T>
std::list<T> make_list()
{
 std::list<T> result;
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1)
{
 std::list<T> result;
 result.push_back(v1);
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2)
{
 std::list<T> result;
 result.push_back(v1);
 result.push_back(v2);
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3)
{
 std::list<T> result;
 result.push_back(v1);
 result.push_back(v2);
 result.push_back(v3);
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3, const T& v4)
{
 std::list<T> result;
 result.push_back(v1);
 result.push_back(v2);
 result.push_back(v3);
 result.push_back(v4);
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3, const T& v4,
                       const T& v5)
{
 std::list<T> result;
 result.push_back(v1);
 result.push_back(v2);
 result.push_back(v3);
 result.push_back(v4);
 result.push_back(v5);
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3, const T& v4,
                       const T& v5, const T& v6)
{
 std::list<T> result;
 result.push_back(v1);
 result.push_back(v2);
 result.push_back(v3);
 result.push_back(v4);
 result.push_back(v5);
 result.push_back(v6);
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3, const T& v4,
                       const T& v5, const T& v6, const T& v7)
{
 std::list<T> result;
 result.push_back(v1);
 result.push_back(v2);
 result.push_back(v3);
 result.push_back(v4);
 result.push_back(v5);
 result.push_back(v6);
 result.push_back(v7);
 return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3, const T& v4,
                       const T& v5, const T& v6, const T& v7, const T& v8)
{
 std::list<T> result;
 result.push_back(v1);
 result.push_back(v2);
 result.push_back(v3);
 result.push_back(v4);
 result.push_back(v5);
 result.push_back(v6);
 result.push_back(v7);
 result.push_back(v8);
 return result;
}

  //   Simple routine to help make a vector of objects

template <typename T>
std::vector<T> make_vector()
{
 std::vector<T> result;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1)
{
 std::vector<T> result(1);
 result[0]=v1;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2)
{
 std::vector<T> result(2);
 result[0]=v1;
 result[1]=v2;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3)
{
 std::vector<T> result(3);
 result[0]=v1;
 result[1]=v2;
 result[2]=v3;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3, const T& v4)
{
 std::vector<T> result(4);
 result[0]=v1;
 result[1]=v2;
 result[2]=v3;
 result[3]=v4;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3, const T& v4,
                           const T& v5)
{
 std::vector<T> result(5);
 result[0]=v1;
 result[1]=v2;
 result[2]=v3;
 result[3]=v4;
 result[4]=v5;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3, const T& v4,
                           const T& v5, const T& v6)
{
 std::vector<T> result(6);
 result[0]=v1;
 result[1]=v2;
 result[2]=v3;
 result[3]=v4;
 result[4]=v5;
 result[5]=v6;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3, const T& v4,
                           const T& v5, const T& v6, const T& v7)
{
 std::vector<T> result(7);
 result[0]=v1;
 result[1]=v2;
 result[2]=v3;
 result[3]=v4;
 result[4]=v5;
 result[5]=v6;
 result[6]=v7;
 return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3, const T& v4,
                           const T& v5, const T& v6, const T& v7, const T& v8)
{
 std::vector<T> result(8);
 result[0]=v1;
 result[1]=v2;
 result[2]=v3;
 result[3]=v4;
 result[4]=v5;
 result[5]=v6;
 result[6]=v7;
 result[7]=v8;
 return result;
}


template<typename T>
std::vector<T> concatenate(const std::vector<T>& A, const std::vector<T>& B)
{
 std::vector<T> AB(A);
 AB.insert(AB.end(),B.begin(),B.end());
 return AB;
}


  // ********************************************************************************
  // *                                                                              *
  // *    Multi-dimensional looper over unsigned integers.  For dimension "k",      *
  // *    its index varies from 0..upper[k]-1.  The number of dimensions            *
  // *    and the upper limits must be specified in the constructor.  Each          *
  // *    dimension must have value 1 or larger (zero values are not allowed).      *
  // *                                                                              *
  // *    Usage:  MultiIntLooper v(3,4);   // row value 0,1,2   column = 0,1,2,3    *
  // *            for (v.start();v.notdone();++v){....}                             *
  // *     loop order = v(0,0), v(1,0), v(2,0), v(0,1), v(1,1), v(2,1), ...         *
  // *                                                                              *
  // ********************************************************************************

class MultiIntLooper
{

     unsigned int m_dim;
     std::vector<unsigned int> m_upper;
     std::vector<unsigned int> m_current;
     bool m_notdone;
     
  public:

     MultiIntLooper(const std::vector<unsigned int> upper_limits) 
          : m_dim(upper_limits.size()), m_upper(upper_limits), 
            m_current(upper_limits.size(),0), m_notdone(true)
            {checkUpper();}

     MultiIntLooper(unsigned int upper1) 
          : m_dim(1), m_upper(make_vector(upper1)), m_current(1,0), m_notdone(true)
            {checkUpper();}
 
     MultiIntLooper(int upper1) 
          : m_dim(1), m_upper(make_vector(unsigner(upper1))), m_current(1,0), m_notdone(true)
            {checkUpper();}

     MultiIntLooper(unsigned int upper1, unsigned int upper2) 
          : m_dim(2), m_upper(make_vector(upper1,upper2)), m_current(2,0), m_notdone(true)
            {checkUpper();}

     MultiIntLooper(int upper1, int upper2) 
          : m_dim(2), m_upper(make_vector(unsigner(upper1),unsigner(upper2))), 
            m_current(2,0), m_notdone(true)  {checkUpper();}

     MultiIntLooper(unsigned int upper1, unsigned int upper2, 
                    unsigned int upper3) 
          : m_dim(3), m_upper(make_vector(upper1,upper2,upper3)), m_current(3,0), m_notdone(true)
            {checkUpper();}

     MultiIntLooper(int upper1, int upper2, int upper3) 
          : m_dim(3), m_upper(make_vector(unsigner(upper1),unsigner(upper2),unsigner(upper3))), 
            m_current(3,0), m_notdone(true)  {checkUpper();}

     MultiIntLooper(unsigned int upper1, unsigned int upper2, 
                    unsigned int upper3, unsigned int upper4) 
          : m_dim(4), m_upper(make_vector(upper1,upper2,upper3,upper4)), m_current(4,0), 
            m_notdone(true)  {checkUpper();}

     MultiIntLooper(int upper1, int upper2, int upper3, int upper4) 
          : m_dim(4), m_upper(make_vector(unsigner(upper1),unsigner(upper2),unsigner(upper3),
            unsigner(upper4))), m_current(4,0), m_notdone(true)
            {checkUpper();}

     ~MultiIntLooper(){}

     const MultiIntLooper& start()
      {
       for (unsigned int k=0;k<m_dim;++k) m_current[k]=0;
       m_notdone=true;
       return *this;}

     bool notdone() const {return m_notdone;}
       
     const MultiIntLooper& operator++()        // prefix
      {return increment();}

     unsigned int numDimensions() const {return m_dim;}
       
     unsigned int operator[](unsigned int k) const
      {return m_current[k];}

     const std::vector<unsigned int>& get_current() const {return m_current;}


  private:

     const MultiIntLooper& increment()
      {
       if (!m_notdone) return *this;
       ++m_current[0];
       unsigned int j=0;
       while (m_current[j]==m_upper[j]){
          m_current[j]=0;
          if (j==(m_dim-1)){
             m_notdone=false; break;}
          ++m_current[++j];}
       return *this;}

     void checkUpper()
      {
       if (m_upper.size()==0) throw(std::invalid_argument("Bad MultiIntLooper: no dimensions"));
       for (unsigned int k=0;k<m_dim;++k) 
          if (m_upper[k]<1) throw(std::invalid_argument("Bad MultiIntLooper upper limits"));
      }
     
     unsigned int unsigner(int k)
      {
       return (k>=0) ? (unsigned int)(k) : 0;
      }
      
};




  // ********************************************************************************
  // *                                                                              *
  // *    Multi-dimensional looper over unsigned integers.  For dimension "k",      *
  // *    its index varies from 0..upper[k]-1.  The number of dimensions            *
  // *    and the upper limits must be specified in the constructor.  Each          *
  // *    dimension must have value 1 or larger (zero values are not allowed).      *
  // *    The "WR" here means "with restrictions".  This looper will avoid any      *
  // *    instances in which any two indices have the same values, except for       *
  // *    those pairs explicitly allowed by an "allowEqual" call.                   *
  // *                                                                              *
  // *    Usage:  MultiIntLooperWR v(3,4);   // row value 0,1,2   column = 0,1,2,3  *
  // *            v.allowEqual(0,1);          // index 0 and 1 can be equal         *
  // *            for (v.start();v.notdone();++v){....}                             *
  // *     loop order = v(0,0), v(1,0), v(2,0), v(0,1), v(1,1), v(2,1), ...         *
  // *                                                                              *
  // ********************************************************************************

class MultiIntLooperWR
{

     unsigned int m_dim;
     std::vector<unsigned int> m_upper;
     std::vector<unsigned int> m_current;
     bool m_notdone;
     std::set<std::pair<unsigned int,unsigned int> > m_allowEqual;

  public:

     MultiIntLooperWR(const std::vector<unsigned int> upper_limits) 
          : m_dim(upper_limits.size()), m_upper(upper_limits), 
            m_current(upper_limits.size(),0), m_notdone(true)  {checkUpper();}

     MultiIntLooperWR(unsigned int upper1) 
          : m_dim(1), m_upper(make_vector(upper1)), m_current(1,0), m_notdone(true)
            {checkUpper();}
 
     MultiIntLooperWR(int upper1) 
          : m_dim(1), m_upper(make_vector(unsigner(upper1))), m_current(1,0), m_notdone(true)
            {checkUpper();}

     MultiIntLooperWR(unsigned int upper1, unsigned int upper2) 
          : m_dim(2), m_upper(make_vector(upper1,upper2)), m_current(2,0), m_notdone(true)
            {checkUpper();}

     MultiIntLooperWR(int upper1, int upper2) 
          : m_dim(2), m_upper(make_vector(unsigner(upper1),unsigner(upper2))), 
            m_current(2,0), m_notdone(true) {checkUpper();}

     MultiIntLooperWR(unsigned int upper1, unsigned int upper2, unsigned int upper3) 
          : m_dim(3), m_upper(make_vector(upper1,upper2,upper3)), m_current(3,0), m_notdone(true)
            {checkUpper();}

     MultiIntLooperWR(int upper1, int upper2, int upper3) 
          : m_dim(3), m_upper(make_vector(unsigner(upper1),unsigner(upper2),unsigner(upper3))), 
            m_current(3,0), m_notdone(true) {checkUpper();}

     MultiIntLooperWR(unsigned int upper1, unsigned int upper2, 
                    unsigned int upper3, unsigned int upper4) 
          : m_dim(4), m_upper(make_vector(upper1,upper2,upper3,upper4)), m_current(4,0), 
            m_notdone(true)  {checkUpper();}

     MultiIntLooperWR(int upper1, int upper2, int upper3, int upper4) 
          : m_dim(4), m_upper(make_vector(unsigner(upper1),unsigner(upper2),unsigner(upper3),
            unsigner(upper4))), m_current(4,0), m_notdone(true)
            {checkUpper();}

     MultiIntLooperWR& allowEqual(int index1, int index2)
      {return allowEqual((unsigned int)(index1),(unsigned int)(index2));}
        

     MultiIntLooperWR& allowEqual(unsigned int index1, int unsigned index2)
      {
       if ((index1<index2)&&(index2<m_dim)) m_allowEqual.insert(std::make_pair(index1,index2));
       else if ((index2>index1)&&(index1<m_dim)) m_allowEqual.insert(std::make_pair(index2,index1));
       return *this;
      }

     ~MultiIntLooperWR(){}

     const MultiIntLooperWR& start()
      {
       for (unsigned int k=0;k<m_dim;++k) m_current[k]=0;
       m_notdone=true;
       if (m_dim<2) return *this;
       if (!allowed()) operator++();
       return *this;}

     bool notdone() const {return m_notdone;}
       
     const MultiIntLooperWR& operator++()        // prefix
      {
       increment();
       while ((!allowed())&&(m_notdone)) increment();
       return *this;
      }

     unsigned int numDimensions() const {return m_dim;}
       
     unsigned int operator[](unsigned int k) const
      {return m_current[k];}

     const std::vector<unsigned int>& get_current() const {return m_current;}

  private:

     const MultiIntLooperWR& increment()
      {
       if (!m_notdone) return *this;
       ++m_current[0];
       unsigned int j=0;
       while (m_current[j]==m_upper[j]){
          m_current[j]=0;
          if (j==(m_dim-1)){
             m_notdone=false; break;}
          ++m_current[++j];}
       return *this;}

     bool allowed()
      {
       for (unsigned int j=0;j<(m_dim-1);++j)
       for (unsigned int k=j+1;k<m_dim;++k)
          if ((m_current[j]==m_current[k])&&(m_allowEqual.count(std::make_pair(j,k))==0)) return false;
       return true;}
     
     void checkUpper()
      {
       if (m_upper.size()==0) throw(std::invalid_argument("Bad MultiIntLooperWR: no dimensions"));
       for (unsigned int k=0;k<m_dim;++k) 
          if (m_upper[k]<1) throw(std::invalid_argument("Bad MultiIntLooperWR upper limits"));
      }
     
     unsigned int unsigner(int k)
      {
       return (k>=0) ? (unsigned int)(k) : 0;
      }
      
};


   // ****************************************************************************

inline bool float_same(const float& x1, const float& x2)
{
 float abstol=1e-12, reltol=1e-4;
 return (std::abs(x1-x2)<=(abstol+reltol*std::abs(x2)));
}

inline bool cmplx_same(const std::complex<float>& z1, const std::complex<float>& z2)
{
 return (float_same(real(z1),real(z2)) && float_same(imag(z1),imag(z2)));
}

inline bool dble_same(const double& x1, const double& x2)
{
 double abstol=1e-15, reltol=1e-12;
 return (std::abs(x1-x2)<=(abstol+reltol*std::abs(x2)));
}

inline bool cmplx_same(const std::complex<double>& z1, const std::complex<double>& z2)
{
 return (dble_same(real(z1),real(z2)) && dble_same(imag(z1),imag(z2)));
}


   // ****************************************************************************

template <typename T>
std::ostream& operator<<(std::ostream& output, std::vector<T> const& values)
{
 output<< "[";
 for (auto const& value : values){
    output << value << ", ";}
 output<< "]";
 return output;
}

template <typename T>
std::ostream& operator<<(std::ostream& output, std::list<T> const& values)
{
 output<< "[";
 for (auto const& value : values){
    output << value << ", ";}
 output<< "]";
 return output;
}

   // ****************************************************************************

#endif
