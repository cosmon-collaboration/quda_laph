#include "verbosity_info.h"
#include "laph_stdio.h"


using namespace std;

namespace LaphEnv {

 // *************************************************************


void Verbosity::set_info(XMLHandler& xmlin)
{
 string tagvalue;
 xmlread(xmlin,"Verbosity",tagvalue,"QudaLaph");
 set_info(tagvalue);
}

void Verbosity::set_info(const std::string& tagvalue)
{
 if (tagvalue=="none"){
    level=QUDA_SILENT;}
 else if (tagvalue=="low"){
    level=QUDA_SUMMARIZE;}
 else if ((tagvalue=="medium")||(tagvalue=="med")){
    level=QUDA_VERBOSE;}
 else if ((tagvalue=="high")||(tagvalue=="full")){
    level=QUDA_DEBUG_VERBOSE;}
 else{
    errorLaph("Invalid input string for verbosity");}
}

void Verbosity::set_info(int ivalue)
{
 if (ivalue<=0){
    level=QUDA_SILENT;}
 else if (ivalue==1){
    level=QUDA_SUMMARIZE;}
 else if (ivalue==2){
    level=QUDA_VERBOSE;}
 else{
    level=QUDA_DEBUG_VERBOSE;}
}

void Verbosity::output(XMLHandler& xmlout) const
{
 string value;
 if (level==QUDA_SILENT) value="none";
 else if (level==QUDA_SUMMARIZE) value="low";
 else if (level==QUDA_VERBOSE) value="medium";
 else value="high";
 xmlout.set_root("Verbosity",value);
}

std::string Verbosity::str() const
{
 XMLHandler xmlout;
 output(xmlout);
 return xmlout.str();
}

std::string Verbosity::message() const
{
 if (level==QUDA_SILENT){
    return string("Verbosity set to none: SILENT");}
 else if (level==QUDA_SUMMARIZE){
    return string("Verbosity set to low: SUMMARIZE");}
 else if (level==QUDA_VERBOSE){
    return string("Verbosity set to medium: VERBOSE");}
 else{
    return string("Verbosity set to high: DEBUG_VERBOSE");}
}


void xml_read_if(XMLHandler& xmlin, Verbosity& v)
{
 if (xml_tag_count(xmlin,"Verbosity")==1){
    Verbosity vv(xmlin);
    v=vv;}
}

// ***************************************************************
}

