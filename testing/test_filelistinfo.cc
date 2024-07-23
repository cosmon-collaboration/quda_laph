#include "task_tests.h"
#include "xml_handler.h"
#include "laph_stdio.h"
#include "filelist_info.h"
#include "layout_info.h"
#include <fstream>

using namespace std;
using namespace LaphEnv;


namespace QLTestEnv {

// ************************************************

void testFileListInfo(XMLHandler& xml_in)
{
 if (xml_tag_count(xml_in,"TestFileListInfo")==0)
 return;

 XMLHandler xmlrr(xml_in,"TestFileListInfo");
 printLaph( "In testFileListInfo...\n");

 XMLHandler xmlr(xmlrr,"First");

 FileListInfo files1(xmlr);
 {XMLHandler xmlout;
 files1.output(xmlout);
 printLaph(make_strf("files1:\n %s",xmlout.str()));}

 {XMLHandler xmlout;
 files1.setOverWrite();
 files1.output(xmlout);
 printLaph(make_strf("files1:\n %s",xmlout.str()));}

 {XMLHandler xmlout;
 files1.setNoOverWrite();
 files1.output(xmlout);
 printLaph(make_strf("files1:\n %s",xmlout.str()));}

 FileListInfo files2("stubber",0,5,true);
 {XMLHandler xmlout;
 files2.output(xmlout);
 printLaph(make_strf("files2:\n %s",xmlout.str()));}
 
 stringstream oss;
 oss << "fouter"<<LayoutInfo::getMyRank()<<".dat";
 string fname(oss.str());
 ofstream fout(fname.c_str());
 if (!fout) std::cout << "could not open file "<<fname<<endl;
 fout << files1.getFileStub()<<endl;
 fout << files1.getMaxFileNumber()<<endl;
 fout << files1.getMinFileNumber()<<endl;
 fout.close();
}


// ***********************************************
}
