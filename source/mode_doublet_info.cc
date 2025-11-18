#include "mode_doublet_info.h"
#include "multi_compare.h"

using namespace std;

namespace LaphEnv {

  // "current" element of xml_in should be <ModeDoubletInfo>

ModeDoubletInfo::ModeDoubletInfo(XMLHandler& xml_in)
{
	xml_tag_assert(xml_in,"ModeDoublet","ModeDoubletInfo");
	XMLHandler xmlr(xml_in, "ModeDoublet");

	int momx, momy, momz, dispLength, ddir, dq; 
	dispLength=0; ddir=0, dq=0;

	vector<int> p;
	xmlread(xmlr, "Momentum", p, "ModeDoubletInfo");
	if (p.size()!=3){
		throw(invalid_argument("Bad mode doublet Momentum"));}
	momx=p[0];
	momy=p[1];
	momz=p[2];

	if (xml_tag_count(xmlr,"DispLength")==1) 
		xmlread(xmlr, "DispLength", dispLength, "ModeDoubletInfo");
	if ((dispLength<0)||(dispLength>63)){
		throw(invalid_argument("Displength must be postive. only < 64 is currently supported."));}
	if (xml_tag_count(xmlr,"DispDir")==1) 
		xmlread(xmlr, "DispDir", ddir, "ModeDoubletInfo");
	if ((ddir!=0)&&(abs(ddir)!=1)&&(abs(ddir)!=2)&&(abs(ddir)!=3)) {
		throw(invalid_argument("DispDir must be -1,-2,-3,0,1,2,3"));}
	if (((ddir!=0)&&(dispLength==0))||((ddir==0)&&(dispLength!=0))) {
		throw(invalid_argument("DispDir must be non-zero iff DisplacementLength is non-zero"));}

	//Do the encoding
	icode=dq;
	icode<<=1; icode|=(momx<0)?1:0; 
	icode<<=6; icode|=abs(momx);
	icode<<=1; icode|=(momy<0)?1:0;
	icode<<=6; icode|=abs(momy);
	icode<<=1; icode|=(momz<0)?1:0;
	icode<<=6; icode|=abs(momz);
	icode<<=6; icode|=dispLength;
	icode<<=1; icode|=(ddir<0)?1:0;
	icode<<=2; icode|=abs(ddir);

}

int ModeDoubletInfo::getDisplacementLength() const
{
 return ((icode>>3)&0x3F);
}

int ModeDoubletInfo::getDisplacementDir() const
{
	int tmp = icode & 0x3u;
  if (((icode>>2)&0x1u)==1) tmp=-tmp;
  return tmp;
}


Momentum ModeDoubletInfo::getMomentum() const
{
 unsigned int tmp=(icode>>9);
 int pz=tmp & 0x3Fu;
 if (((tmp>>6)&0x1u)==1) pz=-pz;
 tmp>>=7;
 int py=tmp & 0x3Fu;
 if (((tmp>>6)&0x1u)==1) py=-py;
 tmp>>=7;
 int px=tmp & 0x3Fu;
 if (((tmp>>6)&0x1u)==1) px=-px;
 return Momentum(px,py,pz);
}

int ModeDoubletInfo::getXMomentum() const
{
 unsigned int tmp=(icode>>23);
 int res=tmp & 0x3Fu;
 if (((tmp>>6)&0x1u)==1) return -res;
 return res;
}
  
int ModeDoubletInfo::getYMomentum() const
{
 unsigned int tmp=(icode>>16);
 int res=tmp & 0x3Fu;
 if (((tmp>>6)&0x1u)==1) return -res;
 return res;
}

int ModeDoubletInfo::getZMomentum() const
{
 unsigned int tmp=(icode>>9);
 int res=tmp & 0x3Fu;
 if (((tmp>>6)&0x1u)==1) return -res;
 return res;
}

bool ModeDoubletInfo::operator==(const ModeDoubletInfo& rhs) const
{
 return (icode==rhs.icode);
}

bool ModeDoubletInfo::operator!=(const ModeDoubletInfo& rhs) const
{
 return (icode!=rhs.icode);
}

bool ModeDoubletInfo::operator<(const ModeDoubletInfo& rhs) const
{
 return (icode<rhs.icode);
}

string ModeDoubletInfo::output(int indent) const {
  XMLHandler xmlout;
	output(xmlout); 
	return xmlout.output(indent);
}

void ModeDoubletInfo::output(XMLHandler& xmlout) const {
 xmlout.set_root("ModeDoublet");
 vector<int> p(3);  p[0]=getXMomentum(); 
 p[1]=getYMomentum(); p[2]=getZMomentum();
 xmlout.put_child("Momentum",make_string(p));
 xmlout.put_child("DispLength",make_string(getDisplacementLength()));
 xmlout.put_child("DispDir",make_string(getDisplacementDir()));
}


 // ********************************************************

   //  useful routine for reading multiple ModeDoubletInfo;
   //  also reads name of coefficient file

void readModeDoublets(XMLHandler& xml_in, 
                         set<ModeDoubletInfo>& mdoubs)
{
 if (xml_tag_count(xml_in,"ModeDoublets")!=1){
    throw(invalid_argument("XML input to CreateModeDoublets must have a single <ModeDoublets> ... </ModeDoublets>"));}
 mdoubs.clear();

 XMLHandler xmlrd(xml_in,"ModeDoublets");
 list<XMLHandler> xmlcs(xmlrd.find("ModeDoublet"));
 for (list<XMLHandler>::iterator it=xmlcs.begin();it!=xmlcs.end();++it){
    mdoubs.insert(ModeDoubletInfo(*it));}
}

// **************************************************************
  }
