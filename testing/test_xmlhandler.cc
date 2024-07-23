#ifdef TESTING
#include <string>
#include "task_tests.h"
#include "xml_handler.h"
#include "util_quda.h"
#include "laph_stdio.h"
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

namespace QLTestEnv {

// ************************************************

void testXMLHandler(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestXMLHandler")==0)
 return;
 
 printLaph("Running TestXMLHandler\n\n");

 string xmlstr("<TopTag>");
 xmlstr+=string(" <Pet>Dog</Pet>");
 xmlstr+=string(" <Pet>Rabbit</Pet>");
 xmlstr+=string(" <Pet>Frog</Pet>");
 xmlstr+=string(" <Pet>Cat</Pet>");
 xmlstr+=string(" <Pet>Slug</Pet>");
 xmlstr+=string(" <Tag1>");
 xmlstr+=string("   <Fruit>Apple</Fruit>");
 xmlstr+=string("   <Fruit>Strawberry</Fruit>");
 xmlstr+=string(" </Tag1>");
 xmlstr+=string(" <Tag2/>");
 xmlstr+=string("</TopTag>");
 
 XMLHandler xmls; xmls.set_from_string(xmlstr);
 
#ifdef MPI_COMMS
 int nrank=quda::comm_size(); 
 int thisrank=quda::comm_rank();
#else
 int nrank=1;
 int thisrank=0;
#endif
 for (int k=0;k<nrank;++k){
    if (thisrank==k){
       ostringstream fname;
       fname << "xmlcontent" << k<<".xml";
       ofstream fout(fname.str());
       fout << xmls.output()<<endl;
       fout.close();}}

 int fvals=1;
 XMLHandler xmltest;
 xmltest.set_from_file("xmlcontent0.xml");  // This tests if gets communicated to all ranks
 xmltest.seek_first_child();
 xmltest.seek_next_sibling();
 xmltest.seek_next_sibling();
 string nodename=xmltest.get_node_name();
 printLaph(make_strf("node name is %s\n",nodename));
 xmltest.seek_first_child();
 string nodeval=xmltest.get_node_value();
 printLaph(make_strf("node val is %s\n",nodeval));
 if (!((nodename=="Pet")&&(nodeval=="Frog"))){
    fvals=0;}

#ifdef MPI_COMMS
 quda::comm_allreduce_int(fvals);
#endif

 bool success=(fvals==nrank);
 if (success) printLaph("ALL TESTS PASSED!!\n");
 else         printLaph("Some tests FAILED\n");

}

// ***********************************************
}
#endif
