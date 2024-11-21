#include "init_quda_laph.h"

using namespace std;

int main(int argc, char **argv) {

  XMLHandler xml_info;
  const int ntasks = init_quda_laph(argc, argv, xml_info);

  run_tasks(xml_info, ntasks);

  finalize();

  return 0;
}
