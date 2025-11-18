#include "mode_triplet_info.h"
#include "multi_compare.h"

using namespace std;

namespace LaphEnv {

  // "current" element of xml_in should be <ModeTripletInfo>

ModeTripletInfo::ModeTripletInfo(XMLHandler& xml_in)
{
	xml_tag_assert(xml_in,"ModeTriplet","ModeTripletInfo");
	XMLHandler xmlr(xml_in, "ModeTriplet");


	int momx, momy, momz, dq, dispLength, ddir; 
	dispLength=0; dq=0; ddir=0;

	vector<int> p;
	xmlread(xmlr, "Momentum", p, "ModeTripletInfo");
	if (p.size()!=3){
		throw(invalid_argument("Bad mode triplet momentum"));}
	momx=p[0];
	momy=p[1];
	momz=p[2];

	if (xml_tag_count(xmlr,"DisplacedQuark")==1) 
		xmlread(xmlr, "DisplacedQuark", dq, "ModeTripletInfo");

	if ((dq!=0)&&(dq!=1)&&(dq!=2)&&(dq!=3)){
		throw(invalid_argument("Bad displaced quark, must be 0,1,2,3"));}

	if (xml_tag_count(xmlr,"DispLength")==1) 
		xmlread(xmlr, "DispLength", dispLength, "ModeTripletInfo");
	if ((dispLength<0)||(dispLength>63)){
		throw(invalid_argument("Displength must be postive. only < 64 is currently supported."));}
	if (xml_tag_count(xmlr,"DispDir")==1) 
		xmlread(xmlr, "DispDir", ddir, "ModeTripletInfo");
	if ((ddir!=0)&&(abs(ddir)!=1)&&(abs(ddir)!=2)&&(abs(ddir)!=3)) {
		throw(invalid_argument("DispDir must be -1,-2,-3,0,1,2,3"));}
	if (((ddir!=0)&&(dq==0))||((ddir==0)&&(dq!=0))) {
		throw(invalid_argument("DispDir must be non-zero iff DisplacedQuark is non-zero"));}
	if (((ddir!=0)&&(dispLength==0))||((ddir==0)&&(dispLength!=0))) {
		throw(invalid_argument("DispDir must be non-zero iff DisplacementLength is non-zero"));}
	if (((dq!=0)&&(dispLength==0))||((dq==0)&&(dispLength!=0))) {
		throw(invalid_argument("DisplacedQuark must be non-zero iff DisplacementLength is non-zero"));}

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

int ModeTripletInfo::getDisplacementLength() const
{
 return ((icode>>3)&0x3F);
}

int ModeTripletInfo::getDisplacementDir() const
{
	int tmp = icode & 0x3u;
  if (((icode>>2)&0x1u)==1) tmp=-tmp;
  return tmp;
}


int ModeTripletInfo::getDisplacedQuark() const
{
 return (icode>>30);
}


Momentum ModeTripletInfo::getMomentum() const
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

int ModeTripletInfo::getXMomentum() const
{
 unsigned int tmp=(icode>>23);
 int res=tmp & 0x3Fu;
 if (((tmp>>6)&0x1u)==1) return -res;
 return res;
}
  
int ModeTripletInfo::getYMomentum() const
{
 unsigned int tmp=(icode>>16);
 int res=tmp & 0x3Fu;
 if (((tmp>>6)&0x1u)==1) return -res;
 return res;
}

int ModeTripletInfo::getZMomentum() const
{
 unsigned int tmp=(icode>>9);
 int res=tmp & 0x3Fu;
 if (((tmp>>6)&0x1u)==1) return -res;
 return res;
}

bool ModeTripletInfo::operator==(const ModeTripletInfo& rhs) const
{
 return (icode==rhs.icode);
}

bool ModeTripletInfo::operator!=(const ModeTripletInfo& rhs) const
{
 return (icode!=rhs.icode);
}

bool ModeTripletInfo::operator<(const ModeTripletInfo& rhs) const
{
 return (icode<rhs.icode);
}

string ModeTripletInfo::output(int indent) const {
	XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

void ModeTripletInfo::output(XMLHandler& xmlout) const {
 xmlout.set_root("ModeTriplet");
 vector<int> p(3);  p[0]=getXMomentum(); 
 p[1]=getYMomentum(); p[2]=getZMomentum();
 xmlout.put_child("Momentum",make_string(p));
 xmlout.put_child("DisplacedQuark",make_string(getDisplacedQuark()));
 xmlout.put_child("DispLength",make_string(getDisplacementLength()));
 xmlout.put_child("DispDir",make_string(getDisplacementDir()));
}


 // ********************************************************

   //  useful routine for reading multiple ModeTripletInfo;
   //  also reads name of coefficient file

void readModeTriplets(XMLHandler& xml_in, 
                         set<ModeTripletInfo>& mtrips)
{
 if (xml_tag_count(xml_in,"ModeTriplets")!=1){
    throw(invalid_argument("XML input to CreateModeTriplets must have a single <ModeTriplets> ... </ModeTriplets>"));}
 mtrips.clear();

 XMLHandler xmlrd(xml_in,"ModeTriplets");
 list<XMLHandler> xmlcs(xmlrd.find("ModeTriplet"));
 for (list<XMLHandler>::iterator it=xmlcs.begin();it!=xmlcs.end();++it){
    mtrips.insert(ModeTripletInfo(*it));}
}

// **************************************************************
}
