#include <string>
#include <sstream>
#include <iostream>
#include <complex>
#include <iomanip>
#include "laph_stdio.h"

     //  make a string from various components; useful for printing out
     //  informational messages or error messages
/*
std::string make_str(const std::string& mesg);

template <typename T1>
std::string make_str(const std::string& startmesg, const T1& arg1)
{
 std::stringstream sout;
 sout << startmesg << arg1;
 return sout.str();
}

template <typename T1, typename T2>
std::string make_str(const std::string& startmesg, const T1& arg1, const T2& arg2)
{
 std::stringstream sout;
 sout << startmesg << arg1 << arg2;
 return sout.str();
}

template <typename T1, typename T2, typename T3>
std::string make_str(const std::string& startmesg, const T1& arg1, const T2& arg2, const T3& arg3)
{
 std::stringstream sout;
 sout << startmesg << arg1 << arg2 << arg3;
 return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4>
std::string make_str(const std::string& startmesg, const T1& arg1, const T2& arg2, const T3& arg3,
                     const T4& arg4)
{
 std::stringstream sout;
 sout << startmesg << arg1 << arg2 << arg3 << arg4;
 return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5>
std::string make_str(const std::string& startmesg, const T1& arg1, const T2& arg2, const T3& arg3,
                     const T4& arg4, const T5& arg5)
{
 std::stringstream sout;
 sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5;
 return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
std::string make_str(const std::string& startmesg, const T1& arg1, const T2& arg2, const T3& arg3,
                     const T4& arg4, const T5& arg5, const T6& arg6)
{
 std::stringstream sout;
 sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5 << arg6;
 return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
std::string make_str(const std::string& startmesg, const T1& arg1, const T2& arg2, const T3& arg3,
                     const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7)
{
 std::stringstream sout;
 sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5 << arg6 << arg7;
 return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
std::string make_str(const std::string& startmesg, const T1& arg1, const T2& arg2, const T3& arg3,
                     const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8)
{
 std::stringstream sout;
 sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5 << arg6 << arg7 << arg8;
 return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8,
          typename T9>
std::string make_str(const std::string& startmesg, const T1& arg1, const T2& arg2, const T3& arg3,
                     const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8,
                     const T9& arg9)
{
 std::stringstream sout;
 sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5 << arg6 << arg7 << arg8 << arg9;
 return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8,
          typename T9, typename T10>
std::string make_str(const std::string& startmesg, const T1& arg1, const T2& arg2, const T3& arg3,
                     const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8,
                     const T9& arg9, const T10& arg10)
{
 std::stringstream sout;
 sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5 << arg6 << arg7 << arg8 << arg9 << arg10;
 return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8,
          typename T9, typename T10, typename T11>
std::string make_str(const std::string& startmesg, const T1& arg1, const T2& arg2, const T3& arg3,
                     const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8,
                     const T9& arg9, const T10& arg10, const T11& arg11)
{
 std::stringstream sout;
 sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5 << arg6 << arg7 << arg8 << arg9 << arg10 << arg11;
 return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8,
          typename T9, typename T10, typename T11, typename T12>
std::string make_str(const std::string& startmesg, const T1& arg1, const T2& arg2, const T3& arg3,
                     const T4& arg4, const T5& arg5, const T6& arg6, const T7& arg7, const T8& arg8,
                     const T9& arg9, const T10& arg10, const T11& arg11, const T12& arg12)
{
 std::stringstream sout;
 sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5 << arg6 << arg7 << arg8 << arg9 
      << arg10 << arg11 << arg12;
 return sout.str();
}


int to_c_str(int arg) { return arg; }

uint to_c_str(uint arg) { return arg; }

long to_c_str(long arg) { return arg; }

unsigned long to_c_str(unsigned long arg) { return arg; }

short to_c_str(short arg) { return arg; }

unsigned short to_c_str(unsigned short arg) { return arg; }

char to_c_str(char arg) { return arg; }

bool to_c_str(bool arg) { return arg; }

float to_c_str(float arg) { return arg; }

double to_c_str(double arg) { return arg; }

const char* to_c_str(const char* arg) { return arg; }

const char* to_c_str(const std::string& arg)
{
 return arg.c_str();
}

template <typename T1>
std::string make_strf(const std::string& fmt, const T1& arg1)
{
 char buffer[1024];
 sprintf(buffer,fmt.c_str(),to_c_str(arg1));
 std::string result(buffer);
 if (result.length()>=1024){
    throw(std::invalid_argument("Buffer overflow in make_strf"));}
 return result;
}

template <typename T1, typename T2>
std::string make_strf(const std::string& fmt, const T1& arg1, const T2& arg2)
{
 char buffer[1024];
 sprintf(buffer,fmt.c_str(),to_c_str(arg1),to_c_str(arg2));
 std::string result(buffer);
 if (result.length()>=1024){
    throw(std::invalid_argument("Buffer overflow in make_strf"));}
 return result;
}
*/
// *************************************************************

using namespace std;

using namespace LaphEnv;
/*
std::string make_str(const std::string& mesg)
{
 return mesg;
}

*/
int main(int argc, char** argv)
{

 std::string result;
 char token1='G';
 int token2=-78;
 uint token3=55;
 bool token4=true;
 double token5=-9.4325/3.0;
 float token6=4.321;
 std::complex<double> token7(1.543,-2.216);
 std::complex<float> token8(-5.4,9.2);
 std::string token9("astring");
 long token10=6;
 short token11=4;
 char token12[]="Cstring";

 cout << make_str("here I am")<<endl;

 cout << make_str("here I am with 1 arg: ",token1)<<endl;
 cout << make_str("here I am with 1 arg: ",token2)<<endl;
 cout << make_str("here I am with 1 arg: ",token3)<<endl;
 cout << make_str("here I am with 1 arg: ",token4)<<endl;
 cout << make_str("here I am with 1 arg: ",token5)<<endl;
 cout << make_str("here I am with 1 arg: ",token6)<<endl;
 cout << make_str("here I am with 1 arg: ",token7)<<endl;
 cout << make_str("here I am with 1 arg: ",token8)<<endl;
 cout << make_str("here I am with 1 arg: ",token9)<<endl;
 cout << make_str("here I am with 1 arg: ",token10)<<endl;
 cout << make_str("here I am with 1 arg: ",token11)<<endl;
 cout << make_str("here I am with 1 arg: ",token12)<<endl;

 cout << make_str("here I am with 2 arg: ",token1 ,token1 )<<endl;
 cout << make_str("here I am with 2 arg: ",token2 ,token1 )<<endl;
 cout << make_str("here I am with 2 arg: ",token3 ,token1 )<<endl;
 cout << make_str("here I am with 2 arg: ",token4 ,token1 )<<endl;
 cout << make_str("here I am with 2 arg: ",token5 ,token1 )<<endl;
 cout << make_str("here I am with 2 arg: ",token6 ,token1 )<<endl;
 cout << make_str("here I am with 2 arg: ",token7 ,token1 )<<endl;
 cout << make_str("here I am with 2 arg: ",token8 ,token1 )<<endl;
 cout << make_str("here I am with 2 arg: ",token9 ,token1 )<<endl;
 cout << make_str("here I am with 2 arg: ",token10,token1 )<<endl;
 cout << make_str("here I am with 2 arg: ",token11,token1 )<<endl;
 cout << make_str("here I am with 2 arg: ",token12,token1 )<<endl;
 cout << make_str("here I am with 2 arg: ",token1 ,token2 )<<endl;
 cout << make_str("here I am with 2 arg: ",token2 ,token2 )<<endl;
 cout << make_str("here I am with 2 arg: ",token3 ,token2 )<<endl;
 cout << make_str("here I am with 2 arg: ",token4 ,token2 )<<endl;
 cout << make_str("here I am with 2 arg: ",token5 ,token2 )<<endl;
 cout << make_str("here I am with 2 arg: ",token6 ,token2 )<<endl;
 cout << make_str("here I am with 2 arg: ",token7 ,token2 )<<endl;
 cout << make_str("here I am with 2 arg: ",token8 ,token2 )<<endl;
 cout << make_str("here I am with 2 arg: ",token9 ,token2 )<<endl;
 cout << make_str("here I am with 2 arg: ",token10,token2 )<<endl;
 cout << make_str("here I am with 2 arg: ",token11,token2 )<<endl;
 cout << make_str("here I am with 2 arg: ",token12,token2 )<<endl;
 cout << make_str("here I am with 2 arg: ",token1 ,token4 )<<endl;
 cout << make_str("here I am with 2 arg: ",token2 ,token4 )<<endl;
 cout << make_str("here I am with 2 arg: ",token3 ,token4 )<<endl;
 cout << make_str("here I am with 2 arg: ",token4 ,token4 )<<endl;
 cout << make_str("here I am with 2 arg: ",token5 ,token4 )<<endl;
 cout << make_str("here I am with 2 arg: ",token6 ,token4 )<<endl;
 cout << make_str("here I am with 2 arg: ",token7 ,token4 )<<endl;
 cout << make_str("here I am with 2 arg: ",token8 ,token4 )<<endl;
 cout << make_str("here I am with 2 arg: ",token9 ,token4 )<<endl;
 cout << make_str("here I am with 2 arg: ",token10,token4 )<<endl;
 cout << make_str("here I am with 2 arg: ",token11,token4 )<<endl;
 cout << make_str("here I am with 2 arg: ",token12,token4 )<<endl;
 cout << make_str("here I am with 2 arg: ",token1 ,token6 )<<endl;
 cout << make_str("here I am with 2 arg: ",token2 ,token6 )<<endl;
 cout << make_str("here I am with 2 arg: ",token3 ,token6 )<<endl;
 cout << make_str("here I am with 2 arg: ",token4 ,token6 )<<endl;
 cout << make_str("here I am with 2 arg: ",token5 ,token6 )<<endl;
 cout << make_str("here I am with 2 arg: ",token6 ,token6 )<<endl;
 cout << make_str("here I am with 2 arg: ",token7 ,token6 )<<endl;
 cout << make_str("here I am with 2 arg: ",token8 ,token6 )<<endl;
 cout << make_str("here I am with 2 arg: ",token9 ,token6 )<<endl;
 cout << make_str("here I am with 2 arg: ",token10,token6 )<<endl;
 cout << make_str("here I am with 2 arg: ",token11,token6 )<<endl;
 cout << make_str("here I am with 2 arg: ",token12,token6 )<<endl;
 cout << make_str("here I am with 2 arg: ",token1 ,token7 )<<endl;
 cout << make_str("here I am with 2 arg: ",token2 ,token7 )<<endl;
 cout << make_str("here I am with 2 arg: ",token3 ,token7 )<<endl;
 cout << make_str("here I am with 2 arg: ",token4 ,token7 )<<endl;
 cout << make_str("here I am with 2 arg: ",token5 ,token7 )<<endl;
 cout << make_str("here I am with 2 arg: ",token6 ,token7 )<<endl;
 cout << make_str("here I am with 2 arg: ",token7 ,token7 )<<endl;
 cout << make_str("here I am with 2 arg: ",token8 ,token7 )<<endl;
 cout << make_str("here I am with 2 arg: ",token9 ,token7 )<<endl;
 cout << make_str("here I am with 2 arg: ",token10,token7 )<<endl;
 cout << make_str("here I am with 2 arg: ",token11,token7 )<<endl;
 cout << make_str("here I am with 2 arg: ",token12,token7 )<<endl;

 cout << make_str("here I am with 3 arg: ",token3," ",token7 )<<endl;
 cout << make_str("here I am with 3 arg: ",token3,token5,token7 )<<endl;
 cout << make_str("here I am with 3 arg: ",token3,token4,token7 )<<endl;
 cout << make_str("here I am with 3 arg: ",token3," ",token7 )<<endl;

 cout << make_str("here I am with 4 arg: ",token1,token3,token5,token7 )<<endl;

 cout << make_str("here I am with 5 arg: ",token1," and ",token3," or ",token5 )<<endl;

 cout << make_str("here I am with 6 arg: ",token1," and ",token3," or ",std::setprecision(12),token5 )<<endl;

 cout << make_str("here I am with 8 arg: ",token1," and ",token3," or ",token5,
              " with ",token9, " bye" )<<endl;

 cout << make_str("No strings:")<<endl;
 cout << make_str("1 string: "," arg1")<<endl;
 cout << make_str("2 strings: "," arg1"," arg2")<<endl;
 cout << make_str("3 strings: "," arg1"," arg2"," arg3")<<endl;
 cout << make_str("4 strings: "," arg1"," arg2"," arg3"," arg4")<<endl;
 cout << make_str("5 strings: "," arg1"," arg2"," arg3"," arg4"," arg5")<<endl;
 cout << make_str("6 strings: "," arg1"," arg2"," arg3"," arg4"," arg5"," arg6")<<endl;
 cout << make_str("7 strings: "," arg1"," arg2"," arg3"," arg4"," arg5"," arg6",
                  " arg7")<<endl;
 cout << make_str("8 strings: "," arg1"," arg2"," arg3"," arg4"," arg5"," arg6",
                  " arg7"," arg8")<<endl;
 cout << make_str("9 strings: "," arg1"," arg2"," arg3"," arg4"," arg5"," arg6",
                  " arg7"," arg8"," arg9")<<endl;
 cout << make_str("10 strings: "," arg1"," arg2"," arg3"," arg4"," arg5"," arg6",
                  " arg7"," arg8"," arg9"," arg10")<<endl;
 cout << make_str("11 strings: "," arg1"," arg2"," arg3"," arg4"," arg5"," arg6",
                  " arg7"," arg8"," arg9"," arg10"," arg11")<<endl;
 cout << make_str("12 strings: "," arg1"," arg2"," arg3"," arg4"," arg5"," arg6",
                  " arg7"," arg8"," arg9"," arg10"," arg11"," arg12")<<endl<<endl;

 cout << make_strf("Formatted string %c",token1)<<endl;
 cout << make_strf("Formatted string %d",token2)<<endl;
 cout << make_strf("Formatted string %d",token3)<<endl;
 cout << make_strf("Formatted string %d",token4)<<endl;
 cout << make_strf("Formatted string %12.9f",token5)<<endl;
 cout << make_strf("Formatted string %16.12f",token6)<<endl;
 cout << make_strf("Formatted string (%16.12f, %16.12f)",real(token7),imag(token7))<<endl;
 cout << make_strf("Formatted string (%16.12f, %16.12f)",real(token8),imag(token8))<<endl;
 cout << make_strf("Formatted string %s",token9)<<endl;
 cout << make_strf("Formatted string %d",token10)<<endl;
 cout << make_strf("Formatted string %d",token11)<<endl;
 cout << make_strf("Formatted string %s",token12)<<endl;
 
 string argA(" <argA> ");
 string argB(" <argB> ");
 string argC(" <argC> ");
 cout << make_strf("Formatted string %s %s %s",argA,argB,argC)<<endl;
 cout << make_strf("Formatted string %d %16.9f %s",token2,token5,argA)<<endl;
 cout << make_strf("Formatted string %s %s %s %s",argA,argB,argC,argA)<<endl;
 cout << make_strf("Formatted string %d %16.9f %s %d",token2,token5,argA,token2+6)<<endl;
 cout << make_strf("Formatted string %s %s %s %s %s",argA,argB,argC,argA,argB)<<endl;
 cout << make_strf("Formatted string %d %16.9f %s %d %s",token2,token5,argA,token2+6,argB)<<endl;
 cout << make_strf("Formatted string %s %s %s %s %s %s",argA,argB,argC,argA,argB,argC)<<endl;
 cout << make_strf("Formatted string %d %16.9f %s %d %s %e",token2,token5,argA,token2+6,argB,token5*3.4)<<endl;
 cout << make_strf("Formatted string %s %s %s %s %s %s %s",argA,argB,argC,argA,argB,argC,argA)<<endl;
 cout << make_strf("Formatted string %d %16.9f %s %d %s %e %d",token2,token5,argA,token2+6,argB,token5*3.4,token2-77)<<endl;
 cout << make_strf("Formatted string %s %s %s %s %s %s %s %s",argA,argB,argC,argA,argB,argC,argA,argB)<<endl;
 cout << make_strf("Formatted string %d %16.9f %s %d %s %e %d %s",token2,token5,argA,token2+6,argB,token5*3.4,token2-77,argC)<<endl;

 return 0;
}
