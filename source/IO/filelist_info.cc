#include "filelist_info.h"
#include "laph_stdio.h"
#include "named_obj_map.h"
#include "utils.h"

using namespace std;

namespace LaphEnv {

FileListInfo::FileListInfo(const XMLHandler &xml_in) {
  XMLHandler xmlr(xml_in, "FileListInfo");
  set_info(xmlr);
}

FileListInfo::FileListInfo(const XMLHandler &xml_in, const string &outertag) {
  XMLHandler xmlq(xml_in, outertag);
  XMLHandler xmlr(xmlq, "FileListInfo");
  set_info(xmlr);
}

FileListInfo::FileListInfo(const std::string &stub, int min_suffix,
                           int max_suffix, bool over_write) {
  try {
    set_info(stub, min_suffix, max_suffix, over_write);
  } catch (const std::exception &msg) {
    errorLaph(make_strf("invalid FileListInfo construction\n %s", msg.what()));
  }
}

void FileListInfo::set_info(const XMLHandler &xmlin) {
  XMLHandler xmlr(xmlin);
  string stub;
  xmlread(xmlr, "FileNameStub", stub, "FileListInfo");
  int min_suffix = 0;
  if (xml_tag_count(xmlr, "MinFileNumber") == 1)
    xmlread(xmlr, "MinFileNumber", min_suffix, "FileListInfo");
  int max_suffix = 0;
  if (xml_tag_count(xmlr, "MaxFileNumber") == 1)
    xmlread(xmlr, "MaxFileNumber", max_suffix, "FileListInfo");
  bool overwrite = false; // protect mode
  if (xml_tag_count(xmlr, "FileMode") == 1) {
    string fmode;
    xmlread(xmlr, "FileMode", fmode, "FileListInfo");
    fmode = tidyString(fmode);
    if (fmode == "overwrite")
      overwrite = true;
  }
  try {
    set_info(stub, min_suffix, max_suffix, overwrite);
  } catch (const std::exception &msg) {
    xmlreadfail(xmlr, "FileListInfo", msg.what());
  }
}

void FileListInfo::set_info(const std::string &stub, int min_suffix,
                            int max_suffix, bool over_write) {
  m_file_stub = tidyString(stub);
  if (m_file_stub.empty()) {
    throw(std::invalid_argument("Blank file name in FileListInfo"));
  }
  if ((m_file_stub.find("NOM_") == 0) && (m_file_stub.length() == 4)) {
    throw(std::invalid_argument("Blank NOM name in FileListInfo"));
  }
  m_max_file_number = max_suffix;
  m_overwrite_mode = over_write; // protect mode
  m_min_file_number = min_suffix;
  if ((m_min_file_number > m_max_file_number) || (m_min_file_number < 0)) {
    throw(std::invalid_argument(
        "minimum file number > maximum file number in FileListInfo"));
  }
}

FileListInfo::FileListInfo(const FileListInfo &fin)
    : m_file_stub(fin.m_file_stub), m_max_file_number(fin.m_max_file_number),
      m_min_file_number(fin.m_min_file_number),
      m_overwrite_mode(fin.m_overwrite_mode) {}

FileListInfo &FileListInfo::operator=(const FileListInfo &fin) {
  m_file_stub = fin.m_file_stub;
  m_max_file_number = fin.m_max_file_number;
  m_min_file_number = fin.m_min_file_number;
  m_overwrite_mode = fin.m_overwrite_mode;
  return *this;
}

std::string FileListInfo::getFileName(int suffix) const {
  stringstream fs;
  fs << m_file_stub << "." << suffix;
  return fs.str();
}

int FileListInfo::getFirstAvailableSuffix(bool global_mode) const {
  for (int suffix = m_min_file_number; suffix <= m_max_file_number; suffix++) {
    string filename = getFileName(suffix);
    if (!fileExists(filename))
      return suffix;
  }
  printLaph("no suffix numbers are available for writing\n ... increase "
            "maxFilenumber");
  throw(std::invalid_argument("error"));
}

int FileListInfo::getFirstAvailableSuffixNOM() const {
  for (int suffix = m_min_file_number; suffix <= m_max_file_number; suffix++) {
    string filename = getFileName(suffix);
    if (!NamedObjMap::query(filename))
      return suffix;
  }
  printLaph(
      "no suffix numbers are available for inserting into TheNamedObjMap");
  printLaph(" ... increase maxFilenumber");
  throw(std::invalid_argument("error"));
}

bool FileListInfo::operator==(const FileListInfo &in) const {
  return ((m_file_stub == in.m_file_stub) &&
          (m_max_file_number == in.m_max_file_number) &&
          (m_min_file_number == in.m_min_file_number));
}

void FileListInfo::output(XMLHandler &xmlout) const {
  xmlout.set_root("FileListInfo");
  xmlout.put_child("FileNameStub", m_file_stub);
  xmlout.put_child("MinFileNumber", make_string(m_min_file_number));
  xmlout.put_child("MaxFileNumber", make_string(m_max_file_number));
  if (m_overwrite_mode)
    xmlout.put_child("FileMode", "overwrite");
  else
    xmlout.put_child("FileMode", "protect");
}

string FileListInfo::str() const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.str();
}
} // namespace LaphEnv
