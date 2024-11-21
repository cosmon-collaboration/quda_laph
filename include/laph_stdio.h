#ifndef LAPH_STDIO_H
#define LAPH_STDIO_H

#include <iomanip>
#include <sstream>
#include <string>

// *********************************************************************
// *                                                                   *
// *  This file provides the functions                                 *
// *                                                                   *
// *      printLaph(const std::string& mesg)                           *
// *      errorLaph(const std::string& mesg, bool abort=true)          *
// *      make_str(const std::string& startmesg, ...)                  *
// *      make_strf(const std::string& format, ...)                    *
// *                                                                   *
// *********************************************************************

// *********************************************************************
// *                                                                   *
// *      printLaph(const std::string& mesg)                           *
// *                                                                   *
// *  This function prints a message to standard output only from      *
// *  the primary MPI rank.  An end-of-line character is inserted      *
// *  at the end of the message.  The user should use the              *
// *  make_str(...) or make_strf(...) to make the message.             *
// *  Typical usage:                                                   *
// *                                                                   *
// *      double x=5.6;                                                *
// *      printLaph(make_str("Value of x is ",x));                     *
// *                                                                   *
// *********************************************************************

// *********************************************************************
// *                                                                   *
// *      errorLaph(const std::string& mesg, bool abort=true)          *
// *                                                                   *
// *  This function prints a message to standard output from every     *
// *  MPI rank.  An end-of-line character is inserted at the end of    *
// *  the message.  The user should use the make_str(...) or           *
// *  make_strf(...) to make the message.  If "abort" is true, then    *
// *  the program aborts due to the error condition.                   *
// *  Typical usage:                                                   *
// *                                                                   *
// *      double x=-5.6;                                               *
// *      if (x<0.0){                                                  *
// *         errorLaph(make_str("Bad value of x is ",x),true);}        *
// *                                                                   *
// *********************************************************************

// *********************************************************************
// *                                                                   *
// *      make_str(mesg,t1,t2,....)                                    *
// *                                                                   *
// *  This routine returns a string formed by concatenating the        *
// *  strings associated with mesg, t1, t2, ... in order.  A           *
// *  stringstream is used.  The objects in t1, t2, ... have to have   *
// *  an operator<< defined.  For formatting the strings, io           *
// *  manipulators can be used.  Currently, there is a maximum of 8    *
// *  arguments after the original "mesg" string. A usage example      *
// *  would be                                                         *
// *                                                                   *
// *      string mesg("This is start of message k = ");                *
// *      int k=6;                                                     *
// *      float x=4.232;                                               *
// *      make_str(mesg,k," and x = ",std::setprecision(12),x);        *
// *                                                                   *
// *********************************************************************

// *********************************************************************
// *                                                                   *
// *      make_strf(fmt,t1,t2,...)                                     *
// *                                                                   *
// *  This routine returns a string using sprintf, except that t1,t2,..*
// *  can be C++ strings: they don't have to be null-terminanted       *
// *  C-style strings. This eliminates the annoying need to use        *
// *  .c_str() for strings. Currently, there is a maximum of 12        *
// *  arguments after the "fmt" string. In case you forget how to      *
// *  form the format string, reminders are below (end of line is "\n")*
// *  For complex number, you need to use the real and imaginary parts *
// *  separately as floats.                                            *
// *                                                                   *
// *  Each format specifier must have the form                         *
// *         %[flags][width][.precision]specifier                      *
// *  Specifier characters:                                            *
// *     %d or %i -> integer (signed or unsigned)                      *
// *     %u -> unsigned integer                                        *
// *     %f -> floating point                                          *
// *     %e -> scientific notation float                               *
// *     %g -> shortest of %f or %e                                    *
// *     %c -> character                                               *
// *     %s -> string                                                  *
// *  Flags:                                                           *
// *     -  left-justify (right is default)                            *
// *     +  print sign (prints + for positive)                         *
// *  Width:                                                           *
// *     minimum width of field in characters                          *
// *  Precision:                                                       *
// *     number of decimal places                                      *
// *                                                                   *
// *********************************************************************

namespace LaphEnv {

//  make a string from various components; useful for printing out
//  informational messages or error messages

std::string make_str(const std::string &mesg);

template <typename T1>
std::string make_str(const std::string &startmesg, const T1 &arg1) {
  std::stringstream sout;
  sout << startmesg << arg1;
  return sout.str();
}

template <typename T1, typename T2>
std::string make_str(const std::string &startmesg, const T1 &arg1,
                     const T2 &arg2) {
  std::stringstream sout;
  sout << startmesg << arg1 << arg2;
  return sout.str();
}

template <typename T1, typename T2, typename T3>
std::string make_str(const std::string &startmesg, const T1 &arg1,
                     const T2 &arg2, const T3 &arg3) {
  std::stringstream sout;
  sout << startmesg << arg1 << arg2 << arg3;
  return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4>
std::string make_str(const std::string &startmesg, const T1 &arg1,
                     const T2 &arg2, const T3 &arg3, const T4 &arg4) {
  std::stringstream sout;
  sout << startmesg << arg1 << arg2 << arg3 << arg4;
  return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5>
std::string make_str(const std::string &startmesg, const T1 &arg1,
                     const T2 &arg2, const T3 &arg3, const T4 &arg4,
                     const T5 &arg5) {
  std::stringstream sout;
  sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5;
  return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6>
std::string make_str(const std::string &startmesg, const T1 &arg1,
                     const T2 &arg2, const T3 &arg3, const T4 &arg4,
                     const T5 &arg5, const T6 &arg6) {
  std::stringstream sout;
  sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5 << arg6;
  return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7>
std::string make_str(const std::string &startmesg, const T1 &arg1,
                     const T2 &arg2, const T3 &arg3, const T4 &arg4,
                     const T5 &arg5, const T6 &arg6, const T7 &arg7) {
  std::stringstream sout;
  sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5 << arg6 << arg7;
  return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7, typename T8>
std::string make_str(const std::string &startmesg, const T1 &arg1,
                     const T2 &arg2, const T3 &arg3, const T4 &arg4,
                     const T5 &arg5, const T6 &arg6, const T7 &arg7,
                     const T8 &arg8) {
  std::stringstream sout;
  sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5 << arg6 << arg7
       << arg8;
  return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7, typename T8, typename T9>
std::string make_str(const std::string &startmesg, const T1 &arg1,
                     const T2 &arg2, const T3 &arg3, const T4 &arg4,
                     const T5 &arg5, const T6 &arg6, const T7 &arg7,
                     const T8 &arg8, const T9 &arg9) {
  std::stringstream sout;
  sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5 << arg6 << arg7
       << arg8 << arg9;
  return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7, typename T8, typename T9, typename T10>
std::string make_str(const std::string &startmesg, const T1 &arg1,
                     const T2 &arg2, const T3 &arg3, const T4 &arg4,
                     const T5 &arg5, const T6 &arg6, const T7 &arg7,
                     const T8 &arg8, const T9 &arg9, const T10 &arg10) {
  std::stringstream sout;
  sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5 << arg6 << arg7
       << arg8 << arg9 << arg10;
  return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7, typename T8, typename T9, typename T10,
          typename T11>
std::string make_str(const std::string &startmesg, const T1 &arg1,
                     const T2 &arg2, const T3 &arg3, const T4 &arg4,
                     const T5 &arg5, const T6 &arg6, const T7 &arg7,
                     const T8 &arg8, const T9 &arg9, const T10 &arg10,
                     const T11 &arg11) {
  std::stringstream sout;
  sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5 << arg6 << arg7
       << arg8 << arg9 << arg10 << arg11;
  return sout.str();
}

template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7, typename T8, typename T9, typename T10,
          typename T11, typename T12>
std::string make_str(const std::string &startmesg, const T1 &arg1,
                     const T2 &arg2, const T3 &arg3, const T4 &arg4,
                     const T5 &arg5, const T6 &arg6, const T7 &arg7,
                     const T8 &arg8, const T9 &arg9, const T10 &arg10,
                     const T11 &arg11, const T12 &arg12) {
  std::stringstream sout;
  sout << startmesg << arg1 << arg2 << arg3 << arg4 << arg5 << arg6 << arg7
       << arg8 << arg9 << arg10 << arg11 << arg12;
  return sout.str();
}

inline int to_c_str(int arg) { return arg; }

inline uint to_c_str(uint arg) { return arg; }

inline long to_c_str(long arg) { return arg; }

inline unsigned long to_c_str(unsigned long arg) { return arg; }

inline short to_c_str(short arg) { return arg; }

inline unsigned short to_c_str(unsigned short arg) { return arg; }

inline char to_c_str(char arg) { return arg; }

inline bool to_c_str(bool arg) { return arg; }

inline float to_c_str(float arg) { return arg; }

inline double to_c_str(double arg) { return arg; }

inline const char *to_c_str(const char *arg) { return arg; }

inline const char *to_c_str(const std::string &arg) { return arg.c_str(); }

std::string charstar_to_string(const char *buffer, int checklength);

template <typename T1>
std::string make_strf(const std::string &fmt, const T1 &arg1) {
  const int bufsize = 1024;
  char buffer[bufsize];
  sprintf(buffer, fmt.c_str(), to_c_str(arg1));
  return charstar_to_string(buffer, bufsize);
}

template <typename T1, typename T2>
std::string make_strf(const std::string &fmt, const T1 &arg1, const T2 &arg2) {
  const int bufsize = 1024;
  char buffer[bufsize];
  sprintf(buffer, fmt.c_str(), to_c_str(arg1), to_c_str(arg2));
  return charstar_to_string(buffer, bufsize);
}

template <typename T1, typename T2, typename T3>
std::string make_strf(const std::string &fmt, const T1 &arg1, const T2 &arg2,
                      const T3 &arg3) {
  const int bufsize = 1024;
  char buffer[bufsize];
  sprintf(buffer, fmt.c_str(), to_c_str(arg1), to_c_str(arg2), to_c_str(arg3));
  return charstar_to_string(buffer, bufsize);
}

template <typename T1, typename T2, typename T3, typename T4>
std::string make_strf(const std::string &fmt, const T1 &arg1, const T2 &arg2,
                      const T3 &arg3, const T4 &arg4) {
  const int bufsize = 1024;
  char buffer[bufsize];
  sprintf(buffer, fmt.c_str(), to_c_str(arg1), to_c_str(arg2), to_c_str(arg3),
          to_c_str(arg4));
  return charstar_to_string(buffer, bufsize);
}

template <typename T1, typename T2, typename T3, typename T4, typename T5>
std::string make_strf(const std::string &fmt, const T1 &arg1, const T2 &arg2,
                      const T3 &arg3, const T4 &arg4, const T5 &arg5) {
  const int bufsize = 1024;
  char buffer[bufsize];
  sprintf(buffer, fmt.c_str(), to_c_str(arg1), to_c_str(arg2), to_c_str(arg3),
          to_c_str(arg4), to_c_str(arg5));
  return charstar_to_string(buffer, bufsize);
}

template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6>
std::string make_strf(const std::string &fmt, const T1 &arg1, const T2 &arg2,
                      const T3 &arg3, const T4 &arg4, const T5 &arg5,
                      const T6 &arg6) {
  const int bufsize = 1024;
  char buffer[bufsize];
  sprintf(buffer, fmt.c_str(), to_c_str(arg1), to_c_str(arg2), to_c_str(arg3),
          to_c_str(arg4), to_c_str(arg5), to_c_str(arg6));
  return charstar_to_string(buffer, bufsize);
}

template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7>
std::string make_strf(const std::string &fmt, const T1 &arg1, const T2 &arg2,
                      const T3 &arg3, const T4 &arg4, const T5 &arg5,
                      const T6 &arg6, const T7 &arg7) {
  const int bufsize = 1024;
  char buffer[bufsize];
  sprintf(buffer, fmt.c_str(), to_c_str(arg1), to_c_str(arg2), to_c_str(arg3),
          to_c_str(arg4), to_c_str(arg5), to_c_str(arg6), to_c_str(arg7));
  return charstar_to_string(buffer, bufsize);
}

template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7, typename T8>
std::string make_strf(const std::string &fmt, const T1 &arg1, const T2 &arg2,
                      const T3 &arg3, const T4 &arg4, const T5 &arg5,
                      const T6 &arg6, const T7 &arg7, const T8 &arg8) {
  const int bufsize = 1024;
  char buffer[bufsize];
  sprintf(buffer, fmt.c_str(), to_c_str(arg1), to_c_str(arg2), to_c_str(arg3),
          to_c_str(arg4), to_c_str(arg5), to_c_str(arg6), to_c_str(arg7),
          to_c_str(arg8));
  return charstar_to_string(buffer, bufsize);
}

// *************************************************************

void printLaph(const std::string &mesg);

void errorLaph(const std::string &mesg, bool abort = true);

// *************************************************************
} // namespace LaphEnv
#endif
