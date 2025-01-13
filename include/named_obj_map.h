#ifndef NAMED_OBJ_H
#define NAMED_OBJ_H

#include <map>

namespace LaphEnv {

//                                NamedObjMap *
//   This is a named object map that allows the user to have persistent data in
// 
//   memory between the different tasks.  The ID keys for the map must be of *
//   type string, and the values can be of any type.  Information about the *
//   data can be stored in an XMLHandler which accompanies each data object. *
//   "NamedObjMap" is a singleton. *
// 
//   Insertion into the NamedObjMap does a COPY on the heap, so any original *
//   data structures or variables used for insertion can be destroyed without *
//   affected the data in the NamedObjMap.  All "get" routines return a *
//   reference to the data. *
//
//   In order to accommodate all different kinds of values to store in the *
//   the named object map, the map stores base pointers.  If the user attempts
// 
//   to get from the map into an object of the wrong type, the dynamic_cast *
//   of the base class pointer/reference to the derived class will throw a *
//   "bad_cast" exception. *
// 
//   To create or insert a new entry of class T value into NamedObjMap with *
//   string key "KeyName", use one of the following methods: *
// 
//       T& var = NamedObjMap::insert<T>("KeyName"); *
//           --> inserts data with default constructor *
//           --> fatal error if key is already in the map *
//           --> can then use "var" to access the data *
// 
//       T data (value to copy into map) *
//       T& var = NamedObjMap::insert<T>("KeyName",data); *
//           --> copies "data: into the map *
//           --> fatal error if key is already in the map *
//           --> can then use "var" to access the data *
// 
//       T data (value to copy into map) *
//       XMLHandler xmlinfo  (info about the data to copy into map) *
//       T& var = NamedObjMap::insert<T>("KeyName",data,xmlinfo); *
//           --> copies "data: into the map *
//           --> copies xmlinfo into the map *
//           --> fatal error if key is already in the map *
//           --> can then use "var" to access the data *
// 
//   To check if a particular object of name "KeyName" is already in *
//   TheNamedObjMap, use *
// 
//       NamedObjMap::query("KeyName") -> returns bool *
// 
//   To get the set of all key names in the map, use *
// 
//       set<string> ids_in_map=NamedObjMap::getIds(); *
// 
//   Accessing entries in the map are done using *
// 
//       T& var = NamedObjMap::getData<T>("KeyName"); *
// --> fatal error if not in the map so check first *
// 
//       T* varptr; *
//       NamedObjMap::getData("KeyName",varptr); *
// 
//       const T* constvarptr; *
//       NamedObjMap::getData("KeyName",constvarptr); *
// 
//       XMLHandler xmlinfo; *
//       NamedObjMap::getDataAndXMLInfo("KeyName",varptr,xmlinfo); *
// 
//       NamedObjMap::getDataAndXMLInfo("KeyName",constvarptr,xmlinfo); *
//  *
//       NamedObjMap::getXMLInfo("KeyName",xmlinfo); *
// 
//   Once you have a reference or a pointer to the data in the map, *
//   you can access or change the data.  To change the XML info, *
//   use *
//       NamedObjMap::setXMLInfo("KeyName",xmlinfo); *
// 
//   To remove an entry in the map, use *
// 
//       NamedObjMap::erase("KeyName"); *
//           --> deletes the data and metadata too *
// 
//   To remove all entries and clear the map, use *
// 
//       NamedObjMap::clear(); *
// 

class NamedObjBase {
  NamedObjBase() {}
  virtual void setXMLInfo(XMLHandler &xml) = 0;
  virtual void getXMLInfo(XMLHandler &xml) const = 0;
  virtual ~NamedObjBase() {}

  friend class NamedObjMap;
  template <typename T> friend class NamedObj;
};

template <typename T> class NamedObj : public NamedObjBase {

  std::unique_ptr<T> data; // automatically deletes upon destructor called
  std::string xml_info;

  NamedObj() : data(new T) {}

  NamedObj(const T &val) : data(new T(val)) {}

  NamedObj(const T &val, XMLHandler &xml_in)
      : data(new T(val)), xml_info(xml_in.output_current()) {}
  ~NamedObj() {}

  virtual void setXMLInfo(XMLHandler &xml) { xml_info = xml.output_current(); }

  virtual void getXMLInfo(XMLHandler &xml) const {
    if (xml_info.empty()) {
      xml.clear();
      return;
    }
    xml.set_from_string(xml_info);
  }

  virtual T &getData() { return *data; }

  virtual const T &getData() const { return *data; }

  friend class NamedObjMap;
};

//  Defined as a singleton (there can only be one instance of a NamedObjMap)

class NamedObjMap {

  static std::map<std::string, NamedObjBase *> the_map;

  typedef std::map<std::string, NamedObjBase *>::iterator nom_iterator;

  typedef std::map<std::string, NamedObjBase *>::const_iterator
      nom_const_iterator;

  NamedObjMap() {}

  ~NamedObjMap() { clear(); }

public:
  NamedObjMap(const NamedObjMap &obj) = delete; // disallow copies

  NamedObjMap &operator=(const NamedObjMap &) = delete; // disallow assignment

  static void clear() {
    for (std::map<std::string, NamedObjBase *>::iterator it = the_map.begin();
         it != the_map.end(); ++it) {
      delete it->second;
    }
    the_map.clear();
  }

  template <typename T> static T &insert(const std::string &id) {
    if (the_map.find(id) != the_map.end()) {
      throw(std::runtime_error(
          "Attempt to NamedObjMap::insert an id already in map"));
    }
    NamedObj<T> *newTdata = new NamedObj<T>();
    NamedObjBase *newdata = dynamic_cast<NamedObjBase *>(newTdata);
    if (newdata == 0) {
      throw(std::runtime_error("Attempt to NamedObjMap::insert failed"));
    }
    std::pair<nom_iterator, bool> result =
        the_map.insert(make_pair(id, newdata));
    if (!result.second) {
      throw(std::runtime_error("Attempt to NamedObjMap::insert failed"));
    }
    return newTdata->getData();
  }

  template <typename T> static T &insert(const std::string &id, const T &data) {
    if (the_map.find(id) != the_map.end()) {
      throw(std::runtime_error(
          "Attempt to NamedObjMap::insert an id already in map"));
    }
    NamedObj<T> *newTdata = new NamedObj<T>(data);
    NamedObjBase *newdata = dynamic_cast<NamedObjBase *>(newTdata);
    if (newdata == 0) {
      throw(std::runtime_error("Attempt to NamedObjMap::insert failed"));
    }
    std::pair<nom_iterator, bool> result =
        the_map.insert(make_pair(id, newdata));
    if (!result.second) {
      throw(std::runtime_error("Attempt to NamedObjMap::insert failed"));
    }
    return newTdata->getData();
  }

  template <typename T>
  static T &insert(const std::string &id, const T &data, XMLHandler &xmlinfo) {
    if (the_map.find(id) != the_map.end()) {
      throw(std::runtime_error(
          "Attempt to NamedObjMap::insert an id already in map"));
    }
    NamedObj<T> *newTdata = new NamedObj<T>(data, xmlinfo);
    NamedObjBase *newdata = dynamic_cast<NamedObjBase *>(newTdata);
    if (newdata == 0) {
      throw(std::runtime_error("Attempt to NamedObjMap::insert failed"));
    }
    std::pair<nom_iterator, bool> result =
        the_map.insert(make_pair(id, newdata));
    if (!result.second) {
      throw(std::runtime_error("Attempt to NamedObjMap::insert failed"));
    }
    return newTdata->getData();
  }

  static bool query(const std::string &id) {
    return (the_map.find(id) != the_map.end()) ? true : false;
  }

  static void erase(const std::string &id) {
    nom_iterator it = the_map.find(id);
    if (it != the_map.end()) {
      delete it->second;
      the_map.erase(it);
    }
  }

  static std::set<std::string> getIds() {
    std::set<std::string> ids;
    for (nom_const_iterator it = the_map.begin(); it != the_map.end(); it++) {
      ids.insert(it->first);
    }
    return ids;
  }

  template <typename T> static T &getData(const std::string &id) {
    return dynamic_cast<NamedObj<T> &>(get(id)).getData();
  }

  template <typename T> static void getData(const std::string &id, T *&data) {
    data = &(dynamic_cast<NamedObj<T> &>(get(id)).getData());
  }

  template <typename T>
  static void getData(const std::string &id, const T *&data) {
    data = &(dynamic_cast<const NamedObj<T> &>(get(id)).getData());
  }

  template <typename T>
  static void getDataAndXMLInfo(const std::string &id, T *&data,
                                XMLHandler &xmlinfo) {
    NamedObjBase &ref = get(id);
    ref.getXMLInfo(xmlinfo);
    data = &(dynamic_cast<NamedObj<T> &>(ref).getData());
  }

  template <typename T>
  static void getDataAndXMLInfo(const std::string &id, const T *&data,
                                XMLHandler &xmlinfo) {
    const NamedObjBase &ref = get(id);
    ref.getXMLInfo(xmlinfo);
    data = &(dynamic_cast<const NamedObj<T> &>(ref).getData());
  }

  static void getXMLInfo(const std::string &id, XMLHandler &xmlinfo) {
    const NamedObjBase &ref = get(id);
    ref.getXMLInfo(xmlinfo);
  }

  static void setXMLInfo(const std::string &id, XMLHandler &xmlinfo) {
    NamedObjBase &ref = get(id);
    ref.setXMLInfo(xmlinfo);
  }

private:
  static NamedObjBase &get(const std::string &id) {
    nom_iterator it = the_map.find(id);
    if (it == the_map.end()) {
      throw(std::runtime_error(std::string("NamedObjMap::get : unknown id = ") +
                               id));
    }
    return *(it->second);
  }
};
} // namespace LaphEnv
#endif
