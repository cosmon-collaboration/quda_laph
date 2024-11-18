#include "laph_stdio.h"
#include "utils.h"
#include <iostream>

namespace LaphEnv {

std::string make_str(const std::string &mesg) { return mesg; }

std::string charstar_to_string(const char *buffer, int checklength) {
  std::string result(buffer);
  if (int(result.length()) >= checklength) {
    throw(std::invalid_argument("Buffer overflow in make_strf"));
  }
  return result;
}

void printLaph(const std::string &mesg) {
  if (isPrimaryRank()) {
    std::cout << mesg << std::endl;
    std::cout.flush();
  }
}

void errorLaph(const std::string &mesg, bool abort) {
#ifdef ARCH_PARALLEL
  std::string errmsg("Error from rank ");
  errmsg += make_string(comm_rank());
  errmsg += ": " + mesg;
#else
  std::string errmsg("Error: ");
  errmsg += mesg;
#endif
  std::cout << errmsg << std::endl;
  if (abort) {
    laph_abort();
  }
}
} // namespace LaphEnv
