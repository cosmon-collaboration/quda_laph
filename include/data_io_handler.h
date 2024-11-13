#ifndef LAPH_DATA_IO_HANDLER_H
#define LAPH_DATA_IO_HANDLER_H

#include "io_map.h"
#include "filelist_info.h"
#include "xml_handler.h"
#include "utils.h"
#include "laph_stdio.h"
#include "named_obj_map.h"
#include <set>
#include <sstream>

namespace LaphEnv {

 // *********************************************************************************
 // *                                                                               *
 // *   The classes defined in this file that are important are                     *
 // *                                                                               *
 // *         DataPutHandlerMFO   (multi-file or NamedObjMap)                       *
 // *         DataGetHandlerMFO   (multi-file or NamedObjMap)                       *
 // *         DataPutHandlerSFO   (single file or NamedObjMap)                      *
 // *         DataGetHandlerSFO   (single file or NamedObjMap)                      *
 // *                                                                               *
 // *         DataPutHandlerMF   (multi-file)                                       *
 // *         DataGetHandlerMF   (multi-file)                                       *
 // *         DataPutHandlerSF   (single file)                                      *
 // *         DataGetHandlerSF   (single file)                                      *
 // *                                                                               *
 // *         DataPutHandlerMFNOM   (to NamedObjMap)                                *
 // *         DataGetHandlerMFNOM   (from NamedObjMap)                              *
 // *         DataPutHandlerSFNOM   (to NamedObjMap)                                *
 // *         DataGetHandlerSFNOM   (from NamedObjMap)                              *
 // *                                                                               *
 // *   In the "MFO" classes above, if the stub of the FileListInfo begins with     *
 // *   "NOM_", then NamedObjMap is used; otherwise, files are used.                *
 // *                                                                               *
 // *********************************************************************************

 
 // *********************************************************************************
 // *                                                                               *
 // *   The classes defined in this file that are important are                     *
 // *                                                                               *
 // *         DataPutHandlerMF   (multi-file)                                       *
 // *         DataGetHandlerMF   (multi-file)                                       *
 // *         DataPutHandlerSF   (single file)                                      *
 // *         DataGetHandlerSF   (single file)                                      *
 // *                                                                               *
 // *                                                                               *
 // *   The "DataPutHandlerMF" and "DataGetHandlerMF" (multi-file) classes handle   *
 // *   inserting data into IOMap records in files specified in a FileListInfo and  *
 // *   subsequently reading the data.  The file key determines which file to       *
 // *   access, and the record key determines which IOMap record in a given file    *
 // *   to access.  The file key and record key are combined into a storage key     *
 // *   which is used in an internal map to store the data in memory.  Do not use   *
 // *   the "put" and "get" classes simultaneously;  use "put" to build up the data *
 // *   in the files, then destroy the insert object; after this, a "get" object    *
 // *   can be created to access the  data in the files.  Objects of these classes  *
 // *   can be created globally (all nodes) or locally (one node).                  *
 // *                                                                               *
 // *   Objects of these "put" handlers always assume an "updating" mode. Existing  *
 // *   files are never erased, and new files are created as needed.  New records   *
 // *   are added to the files.  If the key of a record to be put already exists    *
 // *   in a file, the put will only occur if "overwrite" is specified AND the      *
 // *   size of the data to be put does not exceed the size of the data already in  *
 // *   the file for that key.  For the multi-file handler, the "overwrite" bool    *
 // *   must be given in the FileListInfo; for the single-file handler, the         *
 // *   constructor takes an explicit bool parameter.                               *
 // *                                                                               *
 // *   To use these classes, one needs the following ingredients:                  *
 // *                                                                               *
 // *    (a) a FileListInfo object and an IOHandler file ID string                  *
 // *                                                                               *
 // *    (b) a handler class "H" that has members                                   *
 // *                                                                               *
 // *            bool checkHeader(XMLHandler& xmlr, int suffix) **MF**              *
 // *            bool checkHeader(XMLHandler& xmlr)             **SF**              *
 // *                                                                               *
 // *          that checks that the header information in the file is good,         *
 // *          returning a boolean value, and                                       *
 // *                                                                               *
 // *            void writeHeader(XMLHandler&,const F&,int)   **MF**                *
 // *            void writeHeader(XMLHandler&)                **SF**                *
 // *                                                                               *
 // *        where "F" is the file key type, that writes out the header string      *
 // *        for each file                                                          *
 // *                                                                               *
 // *    (c) both classes need to extract the file key from the header string, so   *
 // *          the file key type must have a constructor that takes only an         *
 // *          XMLHandler and an output(XMLHandler&) member; you need to specify the*
 // *          header tag that the file key will be enclosed in; the file key must  *
 // *          also have a "<" operator and an "==" operator defined                *
 // *                                                                               *
 // *    (d) the record key class must have all of the features of an IOMap key:    *
 // *                                                                               *
 // *       -- since used in a C++ map, a less than operator must be define         *
 // *              const K& K::operator<(const K& rhs);                             *
 // *       -- a numbytes(ioh,K) function, where ioh is an IOHandler, must be       *
 // *           defined to give the number of bytes each key occupies in an         *
 // *           IOHandler file                                                      *
 // *       -- a copy constructor K(const K& in) must be defined                    *
 // *           (a default constructor is not needed)                               *
 // *       -- a multi_read(ioh, vector<K>&,n) must be defined to read n keys       *
 // *       -- a multi_write(ioh, const vector<K>&) must be defined                 *
 // *                                                                               *
 // *    (e) the value type must have the following features:                       *
 // *                                                                               *
 // *       -- a write(ioh, const V&) must be defined (ioh is an IOHandler object)  *
 // *       -- a read(ioh, V&) must be defined                                      *
 // *       -- a numbytes(ioh,V) must be defined giving number of bytes occupied    *
 // *            by V in an IOHandler file                                          *
 // *                                                                               *
 // *                                                                               *
 // *   Usage:  (inserting)                                                         *
 // *                                                                               *
 // *     FileListInfo files;     // specify overwrite mode in here                 *
 // *     string fid("id string");                                                  *
 // *     string HeaderTag("HeaderTag");                                            *
 // *     SomeHandler H;  //  R=record type, F=file type, D=data type               *
 // *     bool globalmode=true;  // false is local mode                             *
 // *                                                                               *
 // *     DataPutHandlerMF<SomeHandler,F,R,D>                                       *
 // *               DIH(H,files,fid,HeaderTag,globalmode);                          *
 // *                                                                               *
 // *        // best to insert data having same file key                            *
 // *        // in sequence to minimize file open/close                             *
 // *                                                                               *
 // *     F fkey; R rkey; D data;                                                   *
 // *     DIH.openFile(fkey);                                                       *
 // *     DIH.putData(rkey,data);                                                   *
 // *     R rkey2; D data2;                                                         *
 // *     DIH.putData(rkey2,data2);                                                 *
 // *     R rkey3; D data3;                                                         *
 // *     DIH.queryData(fkey,rkey3); // check if exists                             *
 // *                                                                               *
 // *        // other members                                                       *
 // *                                                                               *
 // *     DIH.getFileListInfo();                                                    *
 // *     DIH.setOverWrite();                                                       *
 // *     DIH.setNoOverWrite();                                                     *
 // *     DIH.close();    // manual close current file                              *
 // *                                                                               *
 // *                                                                               *
 // *   Usage:  (reading)                                                           *
 // *                                                                               *
 // *     FileListInfo files;                                                       *
 // *     string fid("id string");                                                  *
 // *     string HeaderTag("HeaderTag");                                            *
 // *     SomeHandler H;  //  R=record type, F=file type, D=data type               *
 // *     bool globalmode=true;  // false is local mode                             *
 // *                                                                               *
 // *     DataGetHandlerMF<SomeHandler,F,R,D>                                       *
 // *               DRH(H,files,fid,HeaderTag,globalmode);                          *
 // *                                                                               *
 // *     F fkey; R rkey; D data;                                                   *
 // *     const D& data=DRH.getData(fkey,rkey);                                     *
 // *     bool flag1=DRH.queryData(fkey,rkey);                                      *
 // *     bool flag2=DRH.queryData(fkey);                                           *
 // *                                                                               *
 // *     DRH.remove(fkey,rkey);                                                    *
 // *     DRH.clearData();                                                          *
 // *                                                                               *
 // *        // other members                                                       *
 // *                                                                               *
 // *     DRH.getFileListInfo();                                                    *
 // *     XMLHandler xmlout;                                                        *
 // *     DRH.getFileMap(xmlout);                                                   *
 // *     list<pair<F,list<R> > > keys=DRH.getKeys();                               *
 // *     DRH.outputKeys(xmlout);                                                   *
 // *                                                                               *
 // *                                                                               *
 // *   The "DataPutHandlerSF" and "DataGetHandlerSF" classes are similar to their  *
 // *   "MF" (multi-file) counterpart but only one single file (SF) is handled.     *
 // *                                                                               *
 // *                                                                               *
 // *********************************************************************************
 
     //  Purely virtual base classes

template <typename H, typename F, typename R, typename D>
class DataGetHandlerBaseMF
{
 public:
    virtual ~DataGetHandlerBaseMF() {};
    virtual const FileListInfo& getFileListInfo() const = 0;
    virtual bool isGlobal() const = 0;
    virtual bool isLocal() const = 0;
    virtual bool queryData(const F& fkey, const R& rkey) = 0;
    virtual bool queryFile(const F& fkey) = 0;
    virtual const D& getData(const F& fkey, const R& rkey) = 0;
    virtual void removeData(const F& fkey, const R& rkey) = 0;
    virtual void removeData(const F& fkey) = 0;
    virtual void clearData() = 0;
    virtual void getFileMap(XMLHandler& xmlout) const = 0;
    virtual std::map<int,F> getSuffixMap() const = 0;
    virtual std::set<F> getFileKeys() const = 0;
    virtual std::set<R> getKeys(const F& fkey) = 0;
    virtual void outputKeys(XMLHandler& xmlout) = 0;
    virtual void removeDataMEM(const F& fkey, const R& rkey) = 0;
    virtual void removeDataMEM(const F& fkey) = 0;
    virtual void clearDataMEM() = 0;
    virtual void closeAll() = 0;  // files closed, but data retained
};

template <typename H, typename F, typename R, typename D>
class DataPutHandlerBaseMF
{
 public:
    virtual ~DataPutHandlerBaseMF() {}
    virtual void setOverWrite() = 0;
    virtual void setNoOverWrite() = 0;
    virtual const FileListInfo& getFileListInfo() const = 0;
    virtual bool isGlobal() const = 0;
    virtual bool isLocal() const = 0;
    virtual void open(const F& fkey) = 0;
    virtual void putData(const R& rkey, const D& data) = 0;  // insert into current file
    virtual void flush() = 0;   // flush current file
    virtual void close() = 0;   // close current file
    virtual void putData(const F& fkey, const R& rkey, const D& data) = 0;
    virtual void flush(const F& fkey) = 0;
    virtual void close(const F& fkey) = 0;
    virtual void flushAll() = 0;
    virtual void closeAll() = 0;
    virtual bool queryData(const F& fkey, const R& rkey) = 0;
    virtual bool queryData(const R& rkey) = 0;  // query in current open file
    virtual bool queryFile(const F& fkey) = 0;
    virtual void getFileMap(XMLHandler& xmlout) const = 0;
    virtual std::map<int,F> getSuffixMap() const = 0;
    virtual std::set<F> getFileKeys() const = 0;
};

template <typename H, typename R, typename D>
class DataGetHandlerBaseSF
{
 public:
    virtual ~DataGetHandlerBaseSF() {};
    virtual std::string getFileName() const = 0;
    virtual bool isGlobal() const = 0;
    virtual bool isLocal() const = 0;
    virtual bool queryData(const R& rkey) = 0;
    virtual const D& getData(const R& rkey) = 0;
    virtual void removeData(const R& rkey) = 0;
    virtual void clearData() = 0;
    virtual std::set<R> getKeys() = 0;
    virtual void outputKeys(XMLHandler& xmlout) = 0;
    virtual void removeDataMEM(const R& rkey) = 0;
    virtual void clearDataMEM() = 0;
};

template <typename H, typename R, typename D>
class DataPutHandlerBaseSF
{
 public:
    virtual ~DataPutHandlerBaseSF() {}
    virtual std::string getFileName() const = 0;
    virtual bool isGlobal() const = 0;
    virtual bool isLocal() const = 0;
    virtual void putData(const R& rkey, const D& data) = 0;  
    virtual void flush() = 0;  
    virtual bool queryData(const R& rkey) = 0;  
};

 
 
   // **************************************************************
   // *                                                            *
   // *                      DataGetHandlerMF                      *
   // *                                                            *
   // **************************************************************


   // "get" class handles the internal storage in "m_storage"
   // which is a map of pointers.  The use of pointers is efficient
   // whenever each data structure is fairly large (100 elements or 
   // more) since only pointers get copied.  
   
   // "fileMap" is a map associating a file key to a FileMapValue, 
   // which contains a suffix and an IOMap pointer.  Upon construction, 
   // "fileMap" is assigned by opening each file one by one, reading the 
   // header string, and extracting the file key.  None of the files is
   // left open.  As data is accessed, the files are opened and left
   // open.  While a file is open, the IOMap keeps all of the record
   // keys in memory.
   
   // "queryData" and "getData" open files as needed, and leave them open.
   // "removeData" removes the data from memory, but the files are left open.
   // "clearData" wipes the data from memory and closes all files.
   // "removeDataMEM" is same as "removeData", except that the data
   // from NamedObjMap is removed too.


template <typename H, typename F, typename R, typename D>
class DataGetHandlerMF : public DataGetHandlerBaseMF<H,F,R,D>
{

    struct StorageKey
    {
      F fkey;
      R rkey;
      StorageKey(const F& in_fkey, const R& in_rkey) : fkey(in_fkey), rkey(in_rkey) {}
      StorageKey(const StorageKey& in)  : fkey(in.fkey), rkey(in.rkey) {}
      StorageKey& operator=(const StorageKey& in)
        {fkey=in.fkey; rkey=in.rkey; return *this;}
      bool operator<(const StorageKey& rhs) const
        {return ((fkey<rhs.fkey) || ((fkey==rhs.fkey)&&(rkey<rhs.rkey)));}
    };

    struct FileMapValue
    {
      int suffix;
      IOMap<R,D> *fptr;
      FileMapValue(int in_suff) : suffix(in_suff), fptr(0) {}
      ~FileMapValue() {delete fptr;}
    };

    typedef std::map<StorageKey,D*>   StorageMapType;
    typedef std::map<F,FileMapValue>  FileMapType;
 
    FileListInfo finfo;
    StorageMapType  m_storage;
    FileMapType  fileMap;
    H& handler;
    bool gmode;
    bool checksums;
    std::string fid;


 public:

    DataGetHandlerMF(H& in_handler, const FileListInfo& in_filelist,
                     const std::string& filetype_id, 
                     const std::string& header_tag, 
                     bool global_mode=true, bool use_checksums=false);

    ~DataGetHandlerMF() {clearData(); fileMap.clear();}

    const FileListInfo& getFileListInfo() const {return finfo;}

    bool isGlobal() const {return gmode;}

    bool isLocal() const {return !gmode;}


    bool queryData(const F& fkey, const R& rkey);

    bool queryFile(const F& fkey);


    const D& getData(const F& fkey, const R& rkey);

    void removeData(const F& fkey, const R& rkey);

    void removeData(const F& fkey);

    void clearData();


    void getFileMap(XMLHandler& xmlout) const;

    std::map<int,F> getSuffixMap() const;

    std::set<F> getFileKeys() const;

    std::set<R> getKeys(const F& fkey);

    void outputKeys(XMLHandler& xmlout);

    void removeDataMEM(const F& fkey, const R& rkey);

    void removeDataMEM(const F& fkey);

    void clearDataMEM();
    
    void closeAll();  // files closed, but data retained



 private:

    void fail(const F& fkey, const R& rkey);
    
    void fail(const F& fkey);

    void fail(const std::string& msg);
    
    IOMap<R,D>* get_file_ptr(const F& fkey);
    
    void open(FileMapValue& fmv);
    
    void close(FileMapValue& fmv);

          // disallow copies
    DataGetHandlerMF(const DataGetHandlerMF& in);
    DataGetHandlerMF& operator=(const DataGetHandlerMF& in);

};



   // Constructor checks that the information in the headers
   // of all existing files is consistent, then sets up the
   // file map.  Each file is opened (one by one), the header string
   // is read, and then the file is closed.

template <typename H, typename F, typename R, typename D>
DataGetHandlerMF<H,F,R,D>::DataGetHandlerMF(H& in_handler,
                                            const FileListInfo& in_filelist,
                                            const std::string& filetype_id,
                                            const std::string& header_tag,
                                            bool global_mode, bool use_checksums)
  :  finfo(in_filelist), handler(in_handler), gmode(global_mode), 
     checksums(use_checksums), fid(tidyString(filetype_id)) 
{
 for (int suffix=finfo.getMinFileNumber();
          suffix<=finfo.getMaxFileNumber();suffix++){

    std::string filename=finfo.getFileName(suffix);

         // open all existing files and check consistency of headers
    std::string headerxml;
    bool exists;
    {IOMap<R,D> iom(gmode);
     exists=iom.peekHeader(headerxml,filename,fid);}
    if (!exists) continue;

    XMLHandler xmlr; xmlr.set_from_string(headerxml);
    if (!handler.checkHeader(xmlr,suffix)){
       fail("Header string in file is\n"+headerxml+"header info in file "
                  +filename+" does not match info in current Handler\n\n"
                  +"...execution aborted...\n");}

            // extract the file key from this file
    try{
       XMLHandler xmlf(xmlr,header_tag);
       F fkey(xmlf);
       typename FileMapType::iterator it=fileMap.find(fkey);
       if (it!=fileMap.end()){
          fail(std::string("duplicate keys in fileMap in current Handler\n")
               +" ... too confusing to continue\n file suffix "
               +int_to_string(suffix)+" and suffix "+int_to_string((it->second).suffix)
               +" have same file key\n");}
       fileMap.insert(std::make_pair(fkey, FileMapValue(suffix)));}
    catch(const std::exception& xp){
       fail("Could not extract FileKey from file "+filename+"\n");}
    }
}


template <typename H, typename F, typename R, typename D>
IOMap<R,D>* DataGetHandlerMF<H,F,R,D>::get_file_ptr(const F& fkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 if (it==fileMap.end()) return 0;
 if (it->second.fptr==0) open(it->second);
 return it->second.fptr;
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::open(FileMapValue& fmv)
{
 std::string filename=finfo.getFileName(fmv.suffix);
 try {
    fmv.fptr=new IOMap<R,D>(gmode);
    fmv.fptr->openReadOnly(filename,fid,checksums);}
 catch(const std::exception& xp){
    fail("failure opening file "+filename+" in DataGetHandlerMF");}
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::close(FileMapValue& fmv)
{
 delete fmv.fptr; 
 fmv.fptr=0;
}


template <typename H, typename F, typename R, typename D>
const D& DataGetHandlerMF<H,F,R,D>::getData(const F& fkey, const R& rkey)
{
 StorageKey skey(fkey,rkey);
 typename StorageMapType::const_iterator dt=m_storage.find(skey);
 if (dt!=m_storage.end()) return *(dt->second);
 IOMap<R,D> *fptr=get_file_ptr(fkey);
 if (fptr==0) fail(fkey);
 D *result(new D);
 m_storage[skey]=result;
 try {fptr->get(rkey,*result);}
 catch(const std::exception& xp){fail(fkey,rkey);}
 return *result;
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::fail(const F& fkey, const R& rkey)
{
 if ((isPrimaryRank())||(!gmode)){
    std::cout << "DataGetHandlerMF could not find requested record:"<<std::endl;
    std::cout << " File stub: "<< finfo.getFileStub()<<std::endl;
    XMLHandler xmlout;
    xmlout.set_root("FileRecordKey");
    XMLHandler xmlh;
    fkey.output(xmlh);
    xmlout.put_child(xmlh);
    rkey.output(xmlh);
    xmlout.put_child(xmlh);
    std::cout << xmlout.str()<<std::endl;
    std::cout << "...execution aborted..."<<std::endl;}
 clearData(); fileMap.clear();
 errorLaph("aborting");
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::fail(const F& fkey)
{
 if ((isPrimaryRank())||(!gmode)){
    std::cout << "DataGetHandlerMF could not find requested file key:"<<std::endl;
    std::cout << " File stub: "<< finfo.getFileStub()<<std::endl;
    XMLHandler xmlout;
    fkey.output(xmlout);
    std::cout << xmlout.str()<<std::endl;
    std::cout << "...execution aborted..."<<std::endl;}
 clearData(); fileMap.clear();
 errorLaph("aborting");
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::fail(const std::string& msg)
{
 if ((isPrimaryRank())||(!gmode)){
    std::cerr << "DataGetHandlerMF error: "<<msg<<std::endl;}
 clearData(); fileMap.clear();
 errorLaph("aborting");
}


template <typename H, typename F, typename R, typename D>
bool DataGetHandlerMF<H,F,R,D>::queryData(const F& fkey, const R& rkey)
{
 StorageKey skey(fkey,rkey);
 typename StorageMapType::const_iterator dt=m_storage.find(skey);
 if (dt!=m_storage.end()) return true;
 IOMap<R,D> *fptr=get_file_ptr(fkey);
 if (fptr==0) return false;
 return fptr->exist(rkey);
}


template <typename H, typename F, typename R, typename D>
bool DataGetHandlerMF<H,F,R,D>::queryFile(const F& fkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 return (it!=fileMap.end());
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::clearData()
{
 for (typename StorageMapType::iterator
      it=m_storage.begin();it!=m_storage.end();it++) delete it->second;
 m_storage.clear();
 for (typename FileMapType::iterator 
      it=fileMap.begin();it!=fileMap.end();it++) close(it->second);
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::closeAll()
{
 for (typename FileMapType::iterator 
      it=fileMap.begin();it!=fileMap.end();it++){
    if (it->second.fptr!=0){
       it->second.fptr->close();}}
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::removeData(const F& fkey, const R& rkey)
{
 StorageKey skey(fkey,rkey);
 typename StorageMapType::iterator dt=m_storage.find(skey);
 if (dt!=m_storage.end()){
    delete dt->second;
    m_storage.erase(dt);}
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::removeData(const F& fkey)
{
 for (typename StorageMapType::iterator
      it=m_storage.begin();it!=m_storage.end();){
    if (it->first.fkey==fkey){
        delete it->second;
        typename StorageMapType::iterator dt=it;
        it++; m_storage.erase(dt);}
    else
        it++;
    }
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::getFileMap(XMLHandler& xmlout) const
{
 xmlout.set_root("FileMap");
 for (typename FileMapType::const_iterator it=fileMap.begin();
      it!=fileMap.end();it++){
    XMLHandler xmlt("Entry");
    XMLHandler xmlp;
    it->first.output(xmlp); xmlt.put_child(xmlp);
    xmlt.put_child("Suffix",make_string((it->second).suffix));
    xmlout.put_child(xmlt);}
}


template <typename H, typename F, typename R, typename D>
std::map<int,F> DataGetHandlerMF<H,F,R,D>::getSuffixMap() const
{
 std::map<int,F> filekeys;
 for (typename FileMapType::const_iterator it=fileMap.begin();
      it!=fileMap.end();it++)
    filekeys.insert(std::make_pair((it->second).suffix,it->first));
 return filekeys;
} 


template <typename H, typename F, typename R, typename D>
std::set<F> DataGetHandlerMF<H,F,R,D>::getFileKeys() const
{
 std::set<F> filekeys;
 for (typename FileMapType::const_iterator it=fileMap.begin();
      it!=fileMap.end();it++)
    filekeys.insert(it->first);
 return filekeys;
} 


template <typename H, typename F, typename R, typename D>
std::set<R> DataGetHandlerMF<H,F,R,D>::getKeys(const F& fkey)
{
 std::set<R> keys;
 IOMap<R,D>* fptr=get_file_ptr(fkey);
 if (fptr!=0) fptr->getKeys(keys);
 return keys;
}
 


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::outputKeys(XMLHandler& xmlout)
{
 xmlout.set_root("AvailableKeys");
 for (typename FileMapType::iterator it=fileMap.begin();
      it!=fileMap.end();it++){
    std::set<R> keys;
    if (it->second.fptr!=0) 
       it->second.fptr->getKeys(keys);
    else{
       open(it->second);
       it->second.fptr->getKeys(keys);
       close(it->second);} 
    for (typename std::set<R>::const_iterator 
            kt=keys.begin();kt!=keys.end();kt++){
       XMLHandler xmli,xmlk;
       it->first.output(xmli);
       kt->output(xmlk);
       XMLHandler xmlt("Key");
       xmlt.put_child(xmli);
       xmlt.put_child(xmlk);
       xmlout.put_child(xmlt);}}
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::clearDataMEM()
{
 clearData();
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::removeDataMEM(const F& fkey, const R& rkey)
{
 removeData(fkey,rkey);
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMF<H,F,R,D>::removeDataMEM(const F& fkey)
{
 removeData(fkey);
}


   // **************************************************************
   // *                                                            *
   // *                      DataPutHandlerMF                      *
   // *                                                            *
   // **************************************************************

   // "fileMap" is a map associating a file key to a FileMapValue, 
   // which contains a suffix and an IOMap pointer.  Upon construction, 
   // "fileMap" is assigned by opening each file one by one, reading the 
   // header string, and extracting the file key.  None of the files is
   // left open.  "current" is an iterator that keeps track of the
   // current file available for data insertion.
   
   // To insert data, you must first call "open" with a file key.  If the
   // file key refers to a file that is already open, then all that happens
   // is that the current file pointer is changed. If this file is not opened,
   // then it is opened for writing.  If no file is associated with the
   // file key, then a new file is created and opened.   After the "open"
   // command, you then use "putData" to insert records into that one file.
   // Subsequently calling "open" with another file key changes the current
   // file pointer.  The previous file is left open.
   
   // "queryData" and "queryFile" opens the files, and they are left open.
   
   // "flush" and "flushAll" can be used to flush the data in one or all 
   // files.  Flushing means writing out the current IOMap record map to
   // file.  If a program crashing, an IOMap that has not been flushed is
   // corrupted.  If flushed, it will be okay for subsequent programs.
   // Note that flushing does not close the file.  Use "close" or "closeAll"
   // to actually close files.  Keep in mind that while a file is open,
   // its IOMap must maintain the map of record keys in memory.


template <typename H, typename F, typename R, typename D>
class DataPutHandlerMF : public DataPutHandlerBaseMF<H,F,R,D>
{

    struct FileMapValue
    {
      int suffix;
      IOMap<R,D> *fptr;
      FileMapValue(int in_suff) : suffix(in_suff), fptr(0) {}
      ~FileMapValue() {delete fptr;}
    };

    typedef std::map<F,FileMapValue>  FileMapType;


    FileListInfo finfo;
    FileMapType fileMap;
    H& handler;
    bool gmode;
    bool checksums;
    std::string fid;
    typename FileMapType::iterator current;
    int strpfactor, strpunit;


 public:

    DataPutHandlerMF(H& in_handler, const FileListInfo& in_filelist,
                     const std::string& filetype_id,  
                     const std::string& header_tag, 
                     bool global_mode=true, bool use_checksums=false,
                     int striping_factor=1, int striping_unit=0);

    ~DataPutHandlerMF() {fileMap.clear();}

    void setOverWrite() {finfo.setOverWrite();}

    void setNoOverWrite() {finfo.setNoOverWrite();}

    const FileListInfo& getFileListInfo() const {return finfo;}

    bool isGlobal() const {return gmode;}

    bool isLocal() const {return !gmode;}


    void open(const F& fkey);

    void putData(const R& rkey, const D& data);  // insert into current file

    void flush();   // flush current file

    void close();   // close current file

             // calls openFile(fkey) first, then inserts
             
    void putData(const F& fkey, const R& rkey, const D& data);

    void flush(const F& fkey);

    void close(const F& fkey);
    
    void flushAll();

    void closeAll();



    bool queryData(const F& fkey, const R& rkey);

    bool queryData(const R& rkey);  // query in current open file

    bool queryFile(const F& fkey);



    void merge(const FileListInfo& infiles,        // adds data from the "infiles"
               const std::string& header_tag);     // into this file list
                                                   // use in serial code ONLY
                                              
    void merge(const std::vector<FileListInfo>& infiles,  
               const std::string& header_tag);     

    void getFileMap(XMLHandler& xmlout) const;

    std::map<int,F> getSuffixMap() const;

    std::set<F> getFileKeys() const;


 private:


    void fail(const std::string& msg);
    
    void fail(const F& fkey, const R& rkey);

    typename FileMapType::iterator get_file_ptr(const F& fkey);
    
    void openUpdate(FileMapValue& fmv);
    
    void openNew(const F& fkey, FileMapValue& fmv);

    void close(typename FileMapType::iterator it);

    void flush(typename FileMapType::iterator it);


          // disallow copies
    DataPutHandlerMF(const DataPutHandlerMF& in);
    DataPutHandlerMF& operator=(const DataPutHandlerMF& in);

};



   // constructor checks that the information in the headers
   // of all existing files is consistent, then sets up the
   // file map

template <typename H, typename F, typename R, typename D>
DataPutHandlerMF<H,F,R,D>::DataPutHandlerMF(H& in_handler,
                                            const FileListInfo& in_filelist,
                                            const std::string& filetype_id,
                                            const std::string& header_tag,
                                            bool global_mode, bool use_checksums,
                                            int striping_factor, int striping_unit)
  :  finfo(in_filelist), handler(in_handler), gmode(global_mode), checksums(use_checksums), 
     fid(tidyString(filetype_id)), current(fileMap.end()),
     strpfactor(striping_factor), strpunit(striping_unit)
{
 for (int suffix=finfo.getMinFileNumber();
          suffix<=finfo.getMaxFileNumber();suffix++){

    std::string filename=finfo.getFileName(suffix);

         // open all existing files and check consistency of headers
    std::string headerxml;
    bool exists;
    {IOMap<R,D> iom(gmode);
     exists=iom.peekHeader(headerxml,filename,fid);}
    if (!exists) continue;

    XMLHandler xmlr; xmlr.set_from_string(headerxml);
    if (!handler.checkHeader(xmlr,suffix)){
       fail("Header string in file is\n"+headerxml+"header info in file "
                  +filename+" does not match info in current Handler\n\n"
                  +"...execution aborted...\n");}

            // extract the file key from this file
    try{
       XMLHandler xmlf(xmlr,header_tag);
       F fkey(xmlf);
       typename FileMapType::iterator it=fileMap.find(fkey);
       if (it!=fileMap.end()){
          fail(std::string("duplicate keys in fileMap in current Handler\n")
               +" ... too confusing to continue\n file suffix "
               +int_to_string(suffix)+" and suffix "+int_to_string((it->second).suffix)
               +" have same file key\n");}
       fileMap.insert(std::make_pair(fkey, FileMapValue(suffix)));}
    catch(const std::exception& xp){
       fail("Could not extract FileKey from file "+filename+"\n");}
    }
}


     //  Get file pointer corresponding to fkey.  If "fkey" is already in the
     //  fileMap, then check if the file is open (fptr assigned).  If not open,
     //  open the file.  If already open, we are done.  If "fkey" is not in the
     //  fileMap, then we find a new suffix and open a new file.

template <typename H, typename F, typename R, typename D>
typename DataPutHandlerMF<H,F,R,D>::FileMapType::iterator 
        DataPutHandlerMF<H,F,R,D>::get_file_ptr(const F& fkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 if (it!=fileMap.end()){
    if (it->second.fptr==0) openUpdate(it->second);}
 else{
    int findex=0;
    try{ findex=finfo.getFirstAvailableSuffix(gmode);}
    catch(const std::exception& xp){
       fail(std::string("could not create file...no more suffix indices")
                  +" available...execution aborted...\n");}
    std::pair< typename FileMapType::iterator, bool > pr;
    pr=fileMap.insert(std::make_pair(fkey,FileMapValue(findex)));
    if (pr.second==true) it=pr.first;
    else{ fail("DataPutHandlerMF insertion failed on index "+int_to_string(findex));}
    openNew(fkey,it->second);}
 return it;
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::openUpdate(FileMapValue& fmv)
{
 std::string filename=finfo.getFileName(fmv.suffix);
 try {
    fmv.fptr=new IOMap<R,D>(gmode);
    std::string header;
    fmv.fptr->openUpdate(filename,fid,header,'L',1,0,checksums,finfo.isModeOverwrite());}
 catch(const std::exception& xp){
    fail("failure opening file "+filename+" in DataPutHandlerMF");}
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::openNew(const F& fkey, FileMapValue& fmv)
{
 std::string filename=finfo.getFileName(fmv.suffix);
 if (doesFileExist(filename)){
    fail(std::string("file collision: another process has created a file")
         +" during execution of this program..best to abort\n");}
 try {
    fmv.fptr=new IOMap<R,D>(gmode);
    XMLHandler headerxml;
    handler.writeHeader(headerxml,fkey,fmv.suffix);  // write header info 
    fmv.fptr->openNew(filename,fid,headerxml.str(),false,'L',strpfactor,strpunit,checksums,
                      finfo.isModeOverwrite());}
 catch(const std::exception& xp){
    fail("failure opening file "+filename+" in DataPutHandlerMF");}
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::close(typename FileMapType::iterator it)
{
 if (it==fileMap.end()) return;
 delete it->second.fptr; 
 it->second.fptr=0;
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::flush(typename FileMapType::iterator it)
{
 if (it==fileMap.end()) return;
 if (it->second.fptr!=0) it->second.fptr->flush();
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::fail(const std::string& msg)
{
 if ((isPrimaryRank())||(!gmode)){
    std::cerr << "DataPutHandlerMF error: "<<msg<<std::endl;}
 fileMap.clear();
 errorLaph("aborting");
}

template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::fail(const F& fkey, const R& rkey)
{
 if ((isPrimaryRank())||(!gmode)){
    std::cout << "DataPutHandlerMF could not insert requested record:"<<std::endl;
    std::cout << " File stub: "<< finfo.getFileStub()<<std::endl;
    XMLHandler xmlout;
    xmlout.set_root("FileRecordKey");
    XMLHandler xmlh;
    fkey.output(xmlh);
    xmlout.put_child(xmlh);
    rkey.output(xmlh);
    xmlout.put_child(xmlh);
    std::cout << xmlout.str()<<std::endl;
    std::cout << "...execution aborted..."<<std::endl;}
 fileMap.clear();
 errorLaph("aborting");
}


    // opens file, resets the current file pointer

template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::open(const F& fkey)
{
 current=get_file_ptr(fkey);
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::putData(const R& rkey, const D& data)
{
 if (current==fileMap.end()) { fail("No current file; cannot insert Data");}
 if (current->second.fptr==0) openUpdate(current->second);
 try{current->second.fptr->put(rkey,data);}
 catch(const std::exception& xp){
    fail(current->first,rkey);}
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::flush()
{
 flush(current);
}

template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::close()
{
 close(current);
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::flushAll()
{
 for (typename FileMapType::iterator it=fileMap.begin();it!=fileMap.end();it++)
    flush(it);
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::closeAll()
{
 for (typename FileMapType::iterator it=fileMap.begin();it!=fileMap.end();it++)
    close(it);
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::putData(
                 const F& fkey, const R& rkey, const D& data)
{
 open(fkey);
 try{ current->second.fptr->put(rkey,data);}
 catch(const std::exception& xp){
    fail(fkey,rkey);}
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::flush(const F& fkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 flush(it);
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::close(const F& fkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 close(it);
}


template <typename H, typename F, typename R, typename D>
bool DataPutHandlerMF<H,F,R,D>::queryData(
                 const F& fkey, const R& rkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 if (it==fileMap.end()) return false;
 if (it->second.fptr==0) openUpdate(it->second);
 return (it->second.fptr)->exist(rkey);
}


template <typename H, typename F, typename R, typename D>
bool DataPutHandlerMF<H,F,R,D>::queryData(const R& rkey)
{
 if (current==fileMap.end()) { fail("No current file; cannot query Data");}
 if (current->second.fptr==0) openUpdate(current->second);
 return current->second.fptr->exist(rkey);
}


template <typename H, typename F, typename R, typename D>
bool DataPutHandlerMF<H,F,R,D>::queryFile(const F& fkey)
{
 typename FileMapType::iterator it=fileMap.find(fkey);
 return (it!=fileMap.end());
}



template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::getFileMap(XMLHandler& xmlout) const
{
 xmlout.set_root("FileMap");
 for (typename FileMapType::const_iterator it=fileMap.begin();
      it!=fileMap.end();it++){
    XMLHandler xmlt("Entry");
    XMLHandler xmlp;
    it->first.output(xmlp); xmlt.put_child(xmlp);
    xmlt.put_child("Suffix",make_string((it->second).suffix));
    xmlout.put_child(xmlt);}
} 

template <typename H, typename F, typename R, typename D>
std::map<int,F> DataPutHandlerMF<H,F,R,D>::getSuffixMap() const
{
 std::map<int,F> filekeys;
 for (typename FileMapType::const_iterator it=fileMap.begin();
      it!=fileMap.end();it++)
    filekeys.insert(std::make_pair((it->second).suffix,it->first));
 return filekeys;
} 


template <typename H, typename F, typename R, typename D>
std::set<F> DataPutHandlerMF<H,F,R,D>::getFileKeys() const
{
 std::set<F> filekeys;
 for (typename FileMapType::const_iterator it=fileMap.begin();
      it!=fileMap.end();it++)
    filekeys.insert(it->first);
 return filekeys;
} 



template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::merge(const FileListInfo& infiles,
                                      const std::string& header_tag)
{
 DataGetHandlerMF<H,F,R,D> inmerge(handler,infiles,fid,header_tag,gmode);
 std::set<F> fkeys=inmerge.getFileKeys();
 printLaph(make_strf("  number of file keys = %d\n",fkeys.size()));
 for (typename std::set<F>::const_iterator ft=fkeys.begin();ft!=fkeys.end();ft++){
    open(*ft);
    std::set<R> rkeys=inmerge.getKeys(*ft);
    printLaph(make_strf("  file suffix %d: number of records to be added = %d",
              current->second.suffix,rkeys.size()));
    for (typename std::set<R>::const_iterator rt=rkeys.begin();rt!=rkeys.end();rt++){
       const D& data=inmerge.getData(*ft,*rt);
       putData(*rt,data);
       inmerge.removeData(*ft,*rt);}
//    close(*ft);
    inmerge.clearData();
    }
} 

   //  The merge subroutine above causes lustre lock up on ranger when several
   //  jobs are running simultaneously.  Hence, it has been revised below in such
   //  a way to slow down directory searches.

template <typename H, typename F, typename R, typename D>
void DataPutHandlerMF<H,F,R,D>::merge(const std::vector<FileListInfo>& infiles,
                                      const std::string& header_tag)
{
 for (unsigned int k=0;k<infiles.size();k++){
    printLaph(make_strf("merging file list k = %\n",k));
    std::string instub=infiles[k].getFileStub();
    int nfkeys=0;
    for (int suffix=infiles[k].getMinFileNumber();suffix<=infiles[k].getMaxFileNumber();suffix++){
       FileListInfo inf(instub,suffix,suffix);
       DataGetHandlerMF<H,F,R,D> inmerge(handler,inf,fid,header_tag,gmode);
       std::set<F> fkeys=inmerge.getFileKeys();
       nfkeys+=fkeys.size();
       for (typename std::set<F>::const_iterator ft=fkeys.begin();ft!=fkeys.end();ft++){
          open(*ft);
          std::set<R> rkeys=inmerge.getKeys(*ft);
          printLaph(make_strf(" file list %d:  from suffix %d into merged suffix %d: number of records added = %d\n",
                     k,suffix,current->second.suffix,rkeys.size()));
          for (typename std::set<R>::const_iterator rt=rkeys.begin();rt!=rkeys.end();rt++){
             const D& data=inmerge.getData(*ft,*rt);
             putData(*rt,data);
             inmerge.removeData(*ft,*rt);}
//        close(*ft);
          }}
    printLaph(make_strf("  number of file keys merged in list %d was %d\n",k,nfkeys)); 
    }
}

   // **************************************************************
   // *                                                            *
   // *                    DataGetHandlerMFNOM                     *
   // *                                                            *
   // **************************************************************

   // This class mirrors DataPutHanderMF but puts into the NamedObjMap,
   // not into files.  No checksums, and only global mode.
   
   // "clearData", "removeData" only removes the knowledge of the
   // data from this handler.  It does NOT remove the data from the
   // NamedObjMap.  Use "removeDataMEM" to do this.


template <typename H, typename F, typename R, typename D>
class DataGetHandlerMFNOM : public DataGetHandlerBaseMF<H,F,R,D>
{

    struct StorageKey
    {
      F fkey;
      R rkey;
      StorageKey(const F& in_fkey, const R& in_rkey) : fkey(in_fkey), rkey(in_rkey) {}
      StorageKey(const StorageKey& in)  : fkey(in.fkey), rkey(in.rkey) {}
      StorageKey& operator=(const StorageKey& in)
        {fkey=in.fkey; rkey=in.rkey; return *this;}
      bool operator<(const StorageKey& rhs) const
        {return ((fkey<rhs.fkey) || ((fkey==rhs.fkey)&&(rkey<rhs.rkey)));}
    };

    struct MemMapValue
    {
      int suffix;
      std::map<R,D> *mptr;
      MemMapValue(int in_suff) : suffix(in_suff), mptr(0) {}
      MemMapValue(int in_suff, std::map<R,D> *inptr) : suffix(in_suff), mptr(inptr) {}
      ~MemMapValue() {}
    };

    typedef std::map<StorageKey,D*>   StorageMapType;
    typedef std::map<F,MemMapValue>   MemMapType;

    FileListInfo finfo;
    StorageMapType  m_storage;
    MemMapType  memMap;
    H& handler;
    std::string fid;

 public:

    DataGetHandlerMFNOM(H& in_handler, const FileListInfo& in_filelist,
                        const std::string& filetype_id, 
                        const std::string& header_tag);

    ~DataGetHandlerMFNOM() {clearData(); memMap.clear();}

    const FileListInfo& getFileListInfo() const {return finfo;}

    bool isGlobal() const {return true;}

    bool isLocal() const {return false;}


    bool queryData(const F& fkey, const R& rkey);

    bool queryFile(const F& fkey);


    const D& getData(const F& fkey, const R& rkey);

    void removeData(const F& fkey, const R& rkey);

    void removeData(const F& fkey);

    void clearData();


    void getFileMap(XMLHandler& xmlout) const;

    std::map<int,F> getSuffixMap() const;

    std::set<F> getFileKeys() const;

    std::set<R> getKeys(const F& fkey);

    void outputKeys(XMLHandler& xmlout);


    void removeDataMEM(const F& fkey, const R& rkey);

    void removeDataMEM(const F& fkey);
    
    void clearDataMEM();

    void closeAll() {}


 private:

    void fail(const F& fkey, const R& rkey);
    
    void fail(const F& fkey);

    void fail(const std::string& msg);
    
    std::map<R,D>* get_mem_ptr(const F& fkey);
    
    void open(MemMapValue& fmv);
    
    void close(MemMapValue& fmv);

          // disallow copies
    DataGetHandlerMFNOM(const DataGetHandlerMFNOM& in);
    DataGetHandlerMFNOM& operator=(const DataGetHandlerMFNOM& in);

};



   // Constructor checks that the information in the headers
   // of all existing files is consistent, then sets up the
   // file map.  Each file is opened (one by one), the header string
   // is read, and then the file is closed.

template <typename H, typename F, typename R, typename D>
DataGetHandlerMFNOM<H,F,R,D>::DataGetHandlerMFNOM(H& in_handler,
                                            const FileListInfo& in_filelist,
                                            const std::string& filetype_id,
                                            const std::string& header_tag)
  :  finfo(in_filelist), handler(in_handler), fid(tidyString(filetype_id)) 
{
 for (int suffix=finfo.getMinFileNumber();
          suffix<=finfo.getMaxFileNumber();suffix++){

    std::string filename=finfo.getFileName(suffix);
       
         // check consistency of headers of all existing data in NamedObjMap,
         // and check for duplicate file keys

    if (NamedObjMap::query(filename)){
       XMLHandler xmlr;
       NamedObjMap::getXMLInfo(filename,xmlr);
       std::string ftype;
       xmlread(xmlr,"FileType",ftype,"DataGetHandlerMFNOM");
       if (ftype!=filetype_id)
          fail("Wrong filetype id for this object in NamedObjMap\n");
               // extract the file key from this file
       try{
          XMLHandler xmlf(xmlr,header_tag);
          F fkey(xmlf);
          typename MemMapType::iterator it=memMap.find(fkey);
          if (it!=memMap.end()){
             fail(std::string("duplicate keys in memMap in current Handler\n")
                  +" ... too confusing to continue\n file suffix "
                  +int_to_string(suffix)+" and suffix "+int_to_string((it->second).suffix)
                  +" have same file key\n");}
          memMap.insert(std::make_pair(fkey, 
              MemMapValue(suffix,&(NamedObjMap::getData<std::map<R,D>>(filename)))));}
       catch(const std::exception& xp){
          fail("Could not extract FileKey from file "+filename+"\n");}}}
}



template <typename H, typename F, typename R, typename D>
std::map<R,D>* DataGetHandlerMFNOM<H,F,R,D>::get_mem_ptr(const F& fkey)
{
 typename MemMapType::iterator it=memMap.find(fkey);
 if (it==memMap.end()) return 0;
 if (it->second.mptr==0) open(it->second);
 return it->second.mptr;
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMFNOM<H,F,R,D>::open(MemMapValue& fmv)
{
 std::string filename=finfo.getFileName(fmv.suffix);
 try {
    if (!NamedObjMap::query(filename)){
        fail(std::string("Could not open ")+filename);}
    fmv.mptr=&(NamedObjMap::getData<std::map<R,D>>(filename));}
 catch(const std::exception& xp){
    fail("failure opening file "+filename+" in DataGetHandlerMFNOM");}
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMFNOM<H,F,R,D>::close(MemMapValue& fmv)
{
 fmv.mptr=0;
}


template <typename H, typename F, typename R, typename D>
const D& DataGetHandlerMFNOM<H,F,R,D>::getData(const F& fkey, const R& rkey)
{
 StorageKey skey(fkey,rkey);
 typename StorageMapType::const_iterator dt=m_storage.find(skey);
 if (dt!=m_storage.end()) return *(dt->second);
 std::map<R,D> *mptr=get_mem_ptr(fkey);
 if (mptr==0) fail(fkey);
 typename std::map<R,D>::iterator rt=mptr->find(rkey);
 if (rt==mptr->end()){fail(fkey,rkey);}
 m_storage[skey]=&(rt->second);
 return rt->second;
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMFNOM<H,F,R,D>::fail(const F& fkey, const R& rkey)
{
 std::cout << "DataGetHandlerMFNOM could not find requested record:"<<std::endl;
 std::cout << " File stub: "<< finfo.getFileStub()<<std::endl;
 XMLHandler xmlout;
 xmlout.set_root("FileRecordKey");
 XMLHandler xmlh;
 fkey.output(xmlh);
 xmlout.put_child(xmlh);
 rkey.output(xmlh);
 xmlout.put_child(xmlh);
 std::cout << xmlout.str()<<std::endl;
 std::cout << "...execution aborted..."<<std::endl;
 clearData(); memMap.clear();
 errorLaph("aborting");
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMFNOM<H,F,R,D>::fail(const F& fkey)
{
 std::cout << "DataGetHandlerMFNOM could not find requested file key:"<<std::endl;
 std::cout << " File stub: "<< finfo.getFileStub()<<std::endl;
 XMLHandler xmlout;
 fkey.output(xmlout);
 std::cout << xmlout.str()<<std::endl;
 std::cout << "...execution aborted..."<<std::endl;
 clearData(); memMap.clear();
 errorLaph("aborting");
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMFNOM<H,F,R,D>::fail(const std::string& msg)
{
 printLaph(make_strf("DataGetHandlerMFNOM error: %s\n",msg));
 clearData(); memMap.clear();
 errorLaph("aborting");
}


template <typename H, typename F, typename R, typename D>
bool DataGetHandlerMFNOM<H,F,R,D>::queryData(const F& fkey, const R& rkey)
{
 StorageKey skey(fkey,rkey);
 typename StorageMapType::const_iterator dt=m_storage.find(skey);
 if (dt!=m_storage.end()) return true;
 std::map<R,D> *mptr=get_mem_ptr(fkey);
 if (mptr==0) return false;
 typename std::map<R,D>::const_iterator rt=mptr->find(rkey);
 return (rt!=mptr->end());
}



template <typename H, typename F, typename R, typename D>
bool DataGetHandlerMFNOM<H,F,R,D>::queryFile(const F& fkey)
{
 typename MemMapType::const_iterator dt=memMap.find(fkey);
 return (dt!=memMap.end()); 
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMFNOM<H,F,R,D>::clearData()
{
 m_storage.clear();
 for (typename MemMapType::iterator 
      it=memMap.begin();it!=memMap.end();it++) close(it->second);
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMFNOM<H,F,R,D>::removeData(const F& fkey, const R& rkey)
{
 StorageKey skey(fkey,rkey);
 typename StorageMapType::iterator dt=m_storage.find(skey);
 if (dt!=m_storage.end()){
    m_storage.erase(dt);}
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMFNOM<H,F,R,D>::removeData(const F& fkey)
{
 for (typename StorageMapType::iterator
      it=m_storage.begin();it!=m_storage.end();){
    if (it->first.fkey==fkey){
        typename StorageMapType::iterator dt=it;
        it++; m_storage.erase(dt);}
    else
        it++;
    }
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMFNOM<H,F,R,D>::removeDataMEM(const F& fkey, const R& rkey)
{
 StorageKey skey(fkey,rkey);
 typename StorageMapType::const_iterator dt=m_storage.find(skey);
 if (dt!=m_storage.end()){
    m_storage.erase(dt);}
 typename MemMapType::iterator it=memMap.find(fkey);
 if (it==memMap.end()) return;
 if (it->second.mptr==0) open(it->second);
 std::map<R,D> *mptr=it->second.mptr;
 if (mptr==0) return;
 mptr->erase(rkey);  // remove from NamedObjMap
 if (mptr->empty()){
    std::string filename=finfo.getFileName(it->second.suffix);
    NamedObjMap::erase(filename);
    memMap.erase(it);}
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMFNOM<H,F,R,D>::removeDataMEM(const F& fkey)
{
 for (typename StorageMapType::iterator
      it=m_storage.begin();it!=m_storage.end();){
    if (it->first.fkey==fkey){
        typename StorageMapType::iterator dt=it;
        it++; m_storage.erase(dt);}
    else
        it++;
    }
 typename MemMapType::iterator it=memMap.find(fkey);
 if (it==memMap.end()) return;
 if (it->second.mptr==0) open(it->second);
 std::map<R,D> *mptr=it->second.mptr;
 if (mptr==0) return;
 mptr->clear();  // remove from NamedObjMap
 std::string filename=finfo.getFileName(it->second.suffix);
 NamedObjMap::erase(filename);
 memMap.erase(it);
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMFNOM<H,F,R,D>::clearDataMEM()
{
 m_storage.clear();
 for (typename MemMapType::iterator 
      it=memMap.begin();it!=memMap.end();it++){
         std::string filename=finfo.getFileName(it->second.suffix);
         NamedObjMap::erase(filename);
         memMap.erase(it);}
}


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMFNOM<H,F,R,D>::getFileMap(XMLHandler& xmlout) const
{
 xmlout.set_root("FileMap");
 for (typename MemMapType::const_iterator it=memMap.begin();
      it!=memMap.end();it++){
    XMLHandler xmlt("Entry");
    XMLHandler xmlp;
    it->first.output(xmlp); xmlt.put_child(xmlp);
    xmlt.put_child("Suffix",make_string((it->second).suffix));
    xmlout.put_child(xmlt);}
}


template <typename H, typename F, typename R, typename D>
std::map<int,F> DataGetHandlerMFNOM<H,F,R,D>::getSuffixMap() const
{
 std::map<int,F> filekeys;
 for (typename MemMapType::const_iterator it=memMap.begin();
      it!=memMap.end();it++)
    filekeys.insert(std::make_pair((it->second).suffix,it->first));
 return filekeys;
} 


template <typename H, typename F, typename R, typename D>
std::set<F> DataGetHandlerMFNOM<H,F,R,D>::getFileKeys() const
{
 std::set<F> filekeys;
 for (typename MemMapType::const_iterator it=memMap.begin();
      it!=memMap.end();it++)
    filekeys.insert(it->first);
 return filekeys;
} 


template <typename H, typename F, typename R, typename D>
std::set<R> DataGetHandlerMFNOM<H,F,R,D>::getKeys(const F& fkey)
{
 std::set<R> keys;
 std::map<R,D> *mptr=get_mem_ptr(fkey);
 if (mptr!=0){
    for (typename std::map<R,D>::const_iterator it=mptr->begin();
       it!=mptr->end();it++)
          keys.insert(it->first);}
 return keys;
}
 


template <typename H, typename F, typename R, typename D>
void DataGetHandlerMFNOM<H,F,R,D>::outputKeys(XMLHandler& xmlout)
{
 xmlout.set_root("AvailableKeys");
 for (typename MemMapType::iterator it=memMap.begin();
      it!=memMap.end();it++){
    std::set<R> keys;
    keys=getKeys(it->first);
    for (typename std::set<R>::const_iterator 
            kt=keys.begin();kt!=keys.end();kt++){
       XMLHandler xmli,xmlk;
       it->first.output(xmli);
       kt->output(xmlk);
       XMLHandler xmlt("Key");
       xmlt.put_child(xmli);
       xmlt.put_child(xmlk);
       xmlout.put_child(xmlt);}}
}


   // **************************************************************
   // *                                                            *
   // *                    DataPutHandlerMFNOM                     *
   // *                                                            *
   // **************************************************************

   // This class mirrors DataPutHanderMF but puts into the NamedObjMap,
   // not into files.  No checksums, and only global mode.


template <typename H, typename F, typename R, typename D>
class DataPutHandlerMFNOM : public DataPutHandlerBaseMF<H,F,R,D>
{

    struct MemMapValue
    {
      int suffix;
      std::map<R,D> *mptr;
      MemMapValue(int in_suff) : suffix(in_suff), mptr(0) {}
      MemMapValue(int in_suff, std::map<R,D> *inptr) : suffix(in_suff), mptr(inptr) {}
      ~MemMapValue() {}
    };

    typedef std::map<F,MemMapValue>   MemMapType;

    FileListInfo finfo;
    MemMapType  memMap;
    H& handler;
    std::string fid;
    typename MemMapType::iterator current;


 public:

    DataPutHandlerMFNOM(H& in_handler, const FileListInfo& in_filelist,
                        const std::string& filetype_id,  
                        const std::string& header_tag);

    ~DataPutHandlerMFNOM() {memMap.clear();}

    void setOverWrite() {finfo.setOverWrite();}

    void setNoOverWrite() {finfo.setNoOverWrite();}

    const FileListInfo& getFileListInfo() const {return finfo;}

    bool isGlobal() const {return true;}

    bool isLocal() const {return false;}


    void open(const F& fkey);

    void putData(const R& rkey, const D& data);  // insert into current file

    void flush(){}   // flush current file

    void close();   // close current file

             // calls openFile(fkey) first, then inserts
             
    void putData(const F& fkey, const R& rkey, const D& data);

    void flush(const F& fkey){}

    void close(const F& fkey);
    
    void flushAll() {}

    void closeAll();



    bool queryData(const F& fkey, const R& rkey);

    bool queryData(const R& rkey);  // query in current open file

    bool queryFile(const F& fkey);


    void getFileMap(XMLHandler& xmlout) const;

    std::map<int,F> getSuffixMap() const;

    std::set<F> getFileKeys() const;


 private:


    void fail(const std::string& msg);
    
    void fail(const F& fkey, const R& rkey);

    typename MemMapType::iterator get_mem_ptr(const F& fkey);
    
    void openUpdate(MemMapValue& fmv);
    
    void openNew(const F& fkey, MemMapValue& fmv);

    void close(typename MemMapType::iterator it);


          // disallow copies
    DataPutHandlerMFNOM(const DataPutHandlerMFNOM& in);
    DataPutHandlerMFNOM& operator=(const DataPutHandlerMFNOM& in);

};



   // constructor checks that the information in the headers
   // of all existing files is consistent, then sets up the
   // file map

template <typename H, typename F, typename R, typename D>
DataPutHandlerMFNOM<H,F,R,D>::DataPutHandlerMFNOM(H& in_handler,
                                            const FileListInfo& in_filelist,
                                            const std::string& filetype_id,
                                            const std::string& header_tag)
  :  finfo(in_filelist), handler(in_handler), fid(tidyString(filetype_id)), 
     current(memMap.end())
{
 for (int suffix=finfo.getMinFileNumber();
          suffix<=finfo.getMaxFileNumber();suffix++){

    std::string filename=finfo.getFileName(suffix);
       
         // check consistency of headers of all existing data in NamedObjMap,
         // and check for duplicate file keys

    if (NamedObjMap::query(filename)){
       XMLHandler xmlr;
       NamedObjMap::getXMLInfo(filename,xmlr);
       std::string ftype;
       xmlread(xmlr,"FileType",ftype,"DataPutHandlerMFNOM");
       if (ftype!=filetype_id)
          fail("Wrong filetype id for this object in NamedObjMap\n");
               // extract the file key from this file
       try{
          XMLHandler xmlf(xmlr,header_tag);
          F fkey(xmlf);
          typename MemMapType::iterator it=memMap.find(fkey);
          if (it!=memMap.end()){
             fail(std::string("duplicate keys in memMap in current Handler\n")
                  +" ... too confusing to continue\n file suffix "
                  +int_to_string(suffix)+" and suffix "+int_to_string((it->second).suffix)
                  +" have same file key\n");}
          memMap.insert(std::make_pair(fkey, 
              MemMapValue(suffix,&(NamedObjMap::getData<std::map<R,D>>(filename)))));}
       catch(const std::exception& xp){
          fail("Could not extract FileKey from file "+filename+"\n");}}}
}


     //  Get file pointer corresponding to fkey.  If "fkey" is already in the
     //  memMap, then check if the file is open (mptr assigned).  If not open,
     //  open the file.  If already open, we are done.  If "fkey" is not in the
     //  memMap, then we find a new suffix and open a new file.

template <typename H, typename F, typename R, typename D>
typename DataPutHandlerMFNOM<H,F,R,D>::MemMapType::iterator 
        DataPutHandlerMFNOM<H,F,R,D>::get_mem_ptr(const F& fkey)
{
 typename MemMapType::iterator it=memMap.find(fkey);
 if (it!=memMap.end()){
    if (it->second.mptr==0) openUpdate(it->second);}
 else{
    int findex=0;
    try{ findex=finfo.getFirstAvailableSuffixNOM();}
    catch(const std::exception& xp){
       fail(std::string("could not create file...no more suffix indices")
                  +" available...execution aborted...\n");}
    std::pair< typename MemMapType::iterator, bool > pr;
    pr=memMap.insert(std::make_pair(fkey,MemMapValue(findex)));
    if (pr.second==true) it=pr.first;
    else{ fail("DataPutHandlerMFNOM insertion failed on index "+int_to_string(findex));}
    openNew(fkey,it->second);}
 return it;
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMFNOM<H,F,R,D>::openUpdate(MemMapValue& fmv)
{
 std::string filename=finfo.getFileName(fmv.suffix);
 try {
    if (!NamedObjMap::query(filename)){
        fail(std::string("Could not openUpdate ")+filename);}
    fmv.mptr=&(NamedObjMap::getData<std::map<R,D>>(filename));}
 catch(const std::exception& xp){
    fail("failure opening file "+filename+" in DataPutHandlerMFNOM");}
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMFNOM<H,F,R,D>::openNew(const F& fkey, MemMapValue& fmv)
{
 std::string filename=finfo.getFileName(fmv.suffix);
 if (NamedObjMap::query(filename)){
    fail(std::string("namedobj collision: another process has created an obj")
         +" during execution of this program..best to abort\n");}
 try {
    fmv.mptr=&(NamedObjMap::insert<std::map<R,D>>(filename));
    XMLHandler xmlinfo;
    handler.writeHeader(xmlinfo,fkey,fmv.suffix);  // write header info 
    xmlinfo.put_child("FileType",fid);
    NamedObjMap::setXMLInfo(filename,xmlinfo);}
 catch(const std::exception& xp){
    fail("failure opening file "+filename+" in DataPutHandlerMFNOM");}
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMFNOM<H,F,R,D>::close(typename MemMapType::iterator it)
{
 if (it==memMap.end()) return;
 it->second.mptr=0;
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMFNOM<H,F,R,D>::fail(const std::string& msg)
{
 printLaph(make_strf("DataPutHandlerMFNOM error: %s\n",msg));
 memMap.clear();
 errorLaph("aborting");
}

template <typename H, typename F, typename R, typename D>
void DataPutHandlerMFNOM<H,F,R,D>::fail(const F& fkey, const R& rkey)
{
 std::cout << "DataPutHandlerMFNOM could not insert requested record:"<<std::endl;
 std::cout << " File stub: "<< finfo.getFileStub()<<std::endl;
 XMLHandler xmlout;
 xmlout.set_root("FileRecordKey");
 XMLHandler xmlh;
 fkey.output(xmlh);
 xmlout.put_child(xmlh);
 rkey.output(xmlh);
 xmlout.put_child(xmlh);
 std::cout << xmlout.str()<<std::endl;
 std::cout << "...execution aborted..."<<std::endl;
 memMap.clear();
 errorLaph("aborting");
}

    // opens file, resets the current file pointer

template <typename H, typename F, typename R, typename D>
void DataPutHandlerMFNOM<H,F,R,D>::open(const F& fkey)
{
 current=get_mem_ptr(fkey);
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMFNOM<H,F,R,D>::putData(const R& rkey, const D& data)
{
 if (current==memMap.end()) { fail("No current file; cannot insert Data");}
 if (current->second.mptr==0) openUpdate(current->second);
 try{
    if (finfo.isModeOverwrite())
       current->second.mptr->operator[](rkey)=data;
    else if (current->second.mptr->find(rkey)==current->second.mptr->end())
       current->second.mptr->insert(std::make_pair(rkey,data));
    else{ 
       fail("Cannot insert Data: record already exists and overwrite not specified");}}
 catch(const std::exception& xp){
    fail(current->first,rkey);}
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMFNOM<H,F,R,D>::close()
{
 close(current);
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMFNOM<H,F,R,D>::closeAll()
{
 for (typename MemMapType::iterator it=memMap.begin();it!=memMap.end();it++)
    close(it);
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMFNOM<H,F,R,D>::putData(
                 const F& fkey, const R& rkey, const D& data)
{
 open(fkey);
 try{
    if (finfo.isModeOverwrite())
       current->second.mptr->operator[](rkey)=data;
    else if (current->second.mptr->find(rkey)==current->second.mptr->end())
       current->second.mptr->insert(std::make_pair(rkey,data));
    else{ 
       fail("Cannot insert Data: record already exists and overwrite not specified");}}
 catch(const std::exception& xp){
    fail(fkey,rkey);}
}


template <typename H, typename F, typename R, typename D>
void DataPutHandlerMFNOM<H,F,R,D>::close(const F& fkey)
{
 typename MemMapType::iterator it=memMap.find(fkey);
 close(it);
}


template <typename H, typename F, typename R, typename D>
bool DataPutHandlerMFNOM<H,F,R,D>::queryData(
                 const F& fkey, const R& rkey)
{
 typename MemMapType::iterator it=memMap.find(fkey);
 if (it==memMap.end()) return false;
 if (it->second.mptr==0) openUpdate(it->second);
 return ((it->second.mptr)->find(rkey)!=(it->second.mptr)->end());
}


template <typename H, typename F, typename R, typename D>
bool DataPutHandlerMFNOM<H,F,R,D>::queryData(const R& rkey)
{
 if (current==memMap.end()) { fail("No current file; cannot query Data");}
 if (current->second.mptr==0) openUpdate(current->second);
 return (current->second.mptr->find(rkey)!=current->second.mptr->end());
}


template <typename H, typename F, typename R, typename D>
bool DataPutHandlerMFNOM<H,F,R,D>::queryFile(const F& fkey)
{
 typename MemMapType::iterator it=memMap.find(fkey);
 return (it!=memMap.end());
}



template <typename H, typename F, typename R, typename D>
void DataPutHandlerMFNOM<H,F,R,D>::getFileMap(XMLHandler& xmlout) const
{
 xmlout.set_root("FileMap");
 for (typename MemMapType::const_iterator it=memMap.begin();
      it!=memMap.end();it++){
    XMLHandler xmlt("Entry");
    XMLHandler xmlp;
    it->first.output(xmlp); xmlt.put_child(xmlp);
    xmlt.put_child("Suffix",make_string((it->second).suffix));
    xmlout.put_child(xmlt);}
} 

template <typename H, typename F, typename R, typename D>
std::map<int,F> DataPutHandlerMFNOM<H,F,R,D>::getSuffixMap() const
{
 std::map<int,F> filekeys;
 for (typename MemMapType::const_iterator it=memMap.begin();
      it!=memMap.end();it++)
    filekeys.insert(std::make_pair((it->second).suffix,it->first));
 return filekeys;
} 


template <typename H, typename F, typename R, typename D>
std::set<F> DataPutHandlerMFNOM<H,F,R,D>::getFileKeys() const
{
 std::set<F> filekeys;
 for (typename MemMapType::const_iterator it=memMap.begin();
      it!=memMap.end();it++)
    filekeys.insert(it->first);
 return filekeys;
} 



   // **************************************************************
   // *                                                            *
   // *                      DataGetHandlerSF                      *
   // *                                                            *
   // **************************************************************


template <typename H, typename R, typename D>
class DataGetHandlerSF : public DataGetHandlerBaseSF<H,R,D>
{

    std::map<R,D*> m_storage;
    IOMap<R,D> *iomptr;
    H& handler;

 public:

    DataGetHandlerSF(H& in_handler, const std::string& file_name, 
                     const std::string& filetype_id, bool global_mode=true,
                     bool use_checksums=false);

    ~DataGetHandlerSF() {clearData(); delete iomptr;}

    std::string getFileName() const {return iomptr->getFileName();}

    bool isGlobal() const {return iomptr->isGlobal();}

    bool isLocal() const {return iomptr->isLocal();}


    bool queryData(const R& rkey);

    const D& getData(const R& rkey);

    void removeData(const R& rkey);

    void clearData();


    std::set<R> getKeys();

    void outputKeys(XMLHandler& xmlout);

    void removeDataMEM(const R& rkey);

    void clearDataMEM();


 private:

    void fail(const R& rkey);

    void fail(const std::string& msg);

          // disallow copies
    DataGetHandlerSF(const DataGetHandlerSF& in);
    DataGetHandlerSF& operator=(const DataGetHandlerSF& in);

};



   // constructor checks that the information in the header
   // is consistent

template <typename H, typename R, typename D>
DataGetHandlerSF<H,R,D>::DataGetHandlerSF(H& in_handler, const std::string& filename,
                                          const std::string& filetype_id, 
                                          bool global_mode, bool use_checksums)
                       :  handler(in_handler)
{
 std::string headerxml;
 try{
    iomptr=new IOMap<R,D>(global_mode);
    iomptr->openReadOnly(filename,filetype_id,headerxml,use_checksums);}
 catch(const std::exception& xp) {
    fail("could not open file "+filename+" for reading");}

 XMLHandler xmlr; xmlr.set_from_string(headerxml);
 if (!handler.checkHeader(xmlr)){
    fail("Header string in file is \n"+headerxml+"\n header info in file "+filename
         +" does not match info in current Handler\n ...execution aborted...\n");}
}


template <typename H, typename R, typename D>
const D& DataGetHandlerSF<H,R,D>::getData(const R& rkey)
{
 typename std::map<R,D*>::const_iterator dt=m_storage.find(rkey);
 if (dt!=m_storage.end()) return *(dt->second);
 D *result(new D);
 m_storage[rkey]=result;
 try {iomptr->get(rkey,*result);}
 catch(const std::exception& xp){delete result; fail(rkey);}
 return *result;
}


template <typename H, typename R, typename D>
bool DataGetHandlerSF<H,R,D>::queryData(const R& rkey)
{
 return iomptr->exist(rkey);
}


template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::clearData()
{
 for (typename std::map<R,D*>::iterator
      it=m_storage.begin();it!=m_storage.end();it++) delete it->second;
 m_storage.clear();
}


template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::removeData(const R& rkey)
{
 typename std::map<R,D*>::iterator dt=m_storage.find(rkey);
 if (dt!=m_storage.end()){
    delete dt->second;
    m_storage.erase(dt);}
}



template <typename H, typename R, typename D>
std::set<R> DataGetHandlerSF<H,R,D>::getKeys()
{
 std::set<R> keys;
 iomptr->getKeys(keys);
 return keys;
} 


template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::outputKeys(XMLHandler& xmlout)
{
 xmlout.set_root("AvailableKeys");
 std::set<R> keys;
 iomptr->getKeys(keys);
 for (typename std::set<R>::const_iterator 
            kt=keys.begin();kt!=keys.end();kt++){
    XMLHandler xmlk;
    kt->output(xmlk);
    XMLHandler xmlt("Key");
    xmlt.put_child(xmlk);
    xmlout.put_child(xmlt);}
}


template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::fail(const R& rkey)
{
 if ((isPrimaryRank())||(iomptr->isLocal())){
    std::cout << "DataGetHandlerSF could not find requested record:"<<std::endl;
    std::cout << " File name: "<< iomptr->getFileName() <<std::endl;
    XMLHandler xmlout;
    xmlout.set_root("RecordKey");
    XMLHandler xmlh;
    rkey.output(xmlh);
    xmlout.put_child(xmlh);
    std::cout << xmlout.str()<<std::endl;
    std::cout << "...execution aborted..."<<std::endl;}
 clearData(); delete iomptr;
 errorLaph("aborting");
}



template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::fail(const std::string& msg)
{
 if ((isPrimaryRank())||(iomptr->isLocal())){
    std::cerr << "DataGetHandlerSF error: "<<msg<<std::endl;}
 clearData(); delete iomptr;
 errorLaph("aborting");
}


template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::clearDataMEM()
{
 clearData();
}


template <typename H, typename R, typename D>
void DataGetHandlerSF<H,R,D>::removeDataMEM(const R& rkey)
{
 removeData(rkey);
}


   // **************************************************************
   // *                                                            *
   // *                      DataPutHandlerSF                      *
   // *                                                            *
   // **************************************************************


template <typename H, typename R, typename D>
class DataPutHandlerSF : public DataPutHandlerBaseSF<H,R,D>
{

    H& handler;
    IOMap<R,D> *iomptr;


 public:

    DataPutHandlerSF(H& inptr, const std::string& file_name, 
                     const std::string& filetype_id,
                     bool overwrite=false, bool global_mode=true,
                     bool use_checksums=false, int striping_factor=1,
                     int striping_unit=0);

    ~DataPutHandlerSF() {delete iomptr;}

    std::string getFileName() const {return iomptr->getFileName();}

    bool isGlobal() const {return iomptr->isGlobal();}

    bool isLocal() const {return iomptr->isLocal();}


    void putData(const R& rkey, const D& data);

    void flush();

    bool queryData(const R& rkey);


 private:
 
    void fail(const std::string& msg);

    void fail(const R& rkey);


          // disallow copies
    DataPutHandlerSF(const DataPutHandlerSF& in);
    DataPutHandlerSF& operator=(const DataPutHandlerSF& in);

};



   // constructor checks that the information in the header
   // is consistent

template <typename H, typename R, typename D>
DataPutHandlerSF<H,R,D>::DataPutHandlerSF(H& in_handler,
                                          const std::string& file_name,
                                          const std::string& filetype_id,
                                          bool overwrite, bool global_mode,
                                          bool use_checksums, int striping_factor,
                                          int striping_unit)
                     :  handler(in_handler)
{
 XMLHandler headerxml;
 handler.writeHeader(headerxml);  // write header info 
 try{
   iomptr=new IOMap<R,D>(global_mode);
   std::string header(headerxml.str());
   iomptr->openUpdate(file_name,filetype_id,header,'L',
                      striping_factor,striping_unit,use_checksums,overwrite);
   if (!(iomptr->isNewFile())){
      XMLHandler xmlr; xmlr.set_from_string(header);
      if (!handler.checkHeader(xmlr)){
         fail("Header string in file is \n"+header+"\n header info in file "
              +file_name+" does not match info in current Handler\n "
              +"...execution aborted...\n");}}}
 catch(const std::exception& xp) {
    fail("could not open file "+file_name+" for writing"); }
}


template <typename H, typename R, typename D>
void DataPutHandlerSF<H,R,D>::putData(const R& rkey, const D& data)
{
 try{
    iomptr->put(rkey,data);}
 catch(const std::exception& xp){
    fail(rkey);}
}


template <typename H, typename R, typename D>
void DataPutHandlerSF<H,R,D>::flush()
{
 iomptr->flush();
}


template <typename H, typename R, typename D>
bool DataPutHandlerSF<H,R,D>::queryData(const R& rkey)
{
 return iomptr->exist(rkey);
}



template <typename H, typename R, typename D>
void DataPutHandlerSF<H,R,D>::fail(const std::string& msg)
{
 if ((isPrimaryRank())||(iomptr->isLocal())){
    std::cerr << "DataPutHandlerSF error: "<<msg<<std::endl;}
 delete iomptr;
 errorLaph("aborting");
}


template <typename H, typename R, typename D>
void DataPutHandlerSF<H,R,D>::fail(const R& rkey)
{
 if ((isPrimaryRank())||(iomptr->isLocal())){
    std::cout << "DataPutHandlerSF could not insert requested record:"<<std::endl;
    std::cout << " File name: "<< iomptr->getFileName() <<std::endl;
    XMLHandler xmlout;
    xmlout.set_root("RecordKey");
    XMLHandler xmlh;
    rkey.output(xmlh);
    xmlout.put_child(xmlh);
    std::cout << xmlout.str()<<std::endl;
    std::cout << "...execution aborted..."<<std::endl;}
 delete iomptr;
 errorLaph("aborting");
}


   // **************************************************************
   // *                                                            *
   // *                    DataGetHandlerSFNOM                     *
   // *                                                            *
   // **************************************************************

   // This class mirrors DataPutHanderSF but puts into NamedObjMap,
   // not into files.  No checksums, and only global mode.


template <typename H, typename R, typename D>
class DataGetHandlerSFNOM : public DataGetHandlerBaseSF<H,R,D>
{

    H& handler;
    std::map<R,D*> m_storage;
    std::map<R,D>* iomptr;
    std::string m_objname;

 public:

    DataGetHandlerSFNOM(H& in_handler, const std::string& file_name, 
                        const std::string& filetype_id);

    ~DataGetHandlerSFNOM() {clearData();}

    std::string getFileName() const {return m_objname;}

    bool isGlobal() const {return true;}

    bool isLocal() const {return false;}


    bool queryData(const R& rkey);

    const D& getData(const R& rkey);

    void removeData(const R& rkey);

    void clearData();

    void removeDataMEM(const R& rkey);

    void clearDataMEM();


    std::set<R> getKeys();

    void outputKeys(XMLHandler& xmlout);


 private:

    void fail(const R& rkey);

    void fail(const std::string& msg);

          // disallow copies
    DataGetHandlerSFNOM(const DataGetHandlerSFNOM& in);
    DataGetHandlerSFNOM& operator=(const DataGetHandlerSFNOM& in);

};



   // constructor checks that the information in the header
   // is consistent

template <typename H, typename R, typename D>
DataGetHandlerSFNOM<H,R,D>::DataGetHandlerSFNOM(H& in_handler, const std::string& filename,
                                                const std::string& filetype_id)
                       :  handler(in_handler), m_objname(tidyString(filename))
{
 if (m_objname.empty())
    fail("empty object name in DataGetHandlerSFNOM");
 if (!NamedObjMap::query(m_objname)){
    fail("could not find object named "+m_objname+" in NamedObjMap");}
 XMLHandler xmlf;
 NamedObjMap::getXMLInfo(m_objname,xmlf);
 std::string ftype;
 xmlread(xmlf,"FileType",ftype,"DataGetHandlerSFNOM");
 if (ftype!=tidyString(filetype_id))
    fail("Wrong filetype id for this object in NamedObjMap\n");
               // extract the file key from this file
 XMLHandler xmlr;
 NamedObjMap::getXMLInfo(m_objname,xmlr);
 if (!handler.checkHeader(xmlr)){
    fail("Header string is \n"+xmlr.str()+"\n for object "+m_objname
         +" does not match info in current Handler\n ...execution aborted...\n");}
 iomptr=&(NamedObjMap::getData<std::map<R,D>>(m_objname));
}



template <typename H, typename R, typename D>
const D& DataGetHandlerSFNOM<H,R,D>::getData(const R& rkey)
{
 typename std::map<R,D*>::const_iterator dt=m_storage.find(rkey);
 if (dt!=m_storage.end()) return *(dt->second);
 typename std::map<R,D>::iterator rt=iomptr->find(rkey);
 if (rt==iomptr->end()){fail(rkey);}
 D *result=&(rt->second);
 m_storage[rkey]=result;
 return *result;
}


template <typename H, typename R, typename D>
bool DataGetHandlerSFNOM<H,R,D>::queryData(const R& rkey)
{
 typename std::map<R,D*>::const_iterator dt=m_storage.find(rkey);
 if (dt!=m_storage.end()) return true;
 typename std::map<R,D>::const_iterator rt=iomptr->find(rkey);
 return (rt!=iomptr->end());
}


template <typename H, typename R, typename D>
void DataGetHandlerSFNOM<H,R,D>::clearData()
{
 m_storage.clear();
}


template <typename H, typename R, typename D>
void DataGetHandlerSFNOM<H,R,D>::removeData(const R& rkey)
{
 typename std::map<R,D*>::iterator dt=m_storage.find(rkey);
 if (dt!=m_storage.end()){
    m_storage.erase(dt);}
}



template <typename H, typename R, typename D>
void DataGetHandlerSFNOM<H,R,D>::removeDataMEM(const R& rkey)
{
 typename std::map<R,D*>::iterator dt=m_storage.find(rkey);
 if (dt!=m_storage.end()){
    m_storage.erase(dt);}
 iomptr->erase(rkey);
 if (iomptr->empty()){
    NamedObjMap::erase(m_objname);}
}



template <typename H, typename R, typename D>
void DataGetHandlerSFNOM<H,R,D>::clearDataMEM()
{
 m_storage.clear();
 NamedObjMap::erase(m_objname);
}


template <typename H, typename R, typename D>
std::set<R> DataGetHandlerSFNOM<H,R,D>::getKeys()
{
 std::set<R> keys;
 for (typename std::map<R,D>::const_iterator it=iomptr->begin();
    it!=iomptr->end();it++)
    keys.insert(it->first);
 return keys;
} 


template <typename H, typename R, typename D>
void DataGetHandlerSFNOM<H,R,D>::outputKeys(XMLHandler& xmlout)
{
 xmlout.set_root("AvailableKeys");
 std::set<R> keys(getKeys());
 for (typename std::set<R>::const_iterator 
            kt=keys.begin();kt!=keys.end();kt++){
    XMLHandler xmlk;
    kt->output(xmlk);
    XMLHandler xmlt("Key");
    xmlt.put_child(xmlk);
    xmlout.put_child(xmlt);}
}


template <typename H, typename R, typename D>
void DataGetHandlerSFNOM<H,R,D>::fail(const R& rkey)
{
 std::cout << "DataGetHandlerSFNOM could not find requested record:"<<std::endl;
 std::cout << " File name: "<< getFileName() <<std::endl;
 XMLHandler xmlout;
 xmlout.set_root("RecordKey");
 XMLHandler xmlh;
 rkey.output(xmlh);
 xmlout.put_child(xmlh);
 std::cout << xmlout.str()<<std::endl;
 std::cout << "...execution aborted..."<<std::endl;
 clearData();
 errorLaph("aborting");
}



template <typename H, typename R, typename D>
void DataGetHandlerSFNOM<H,R,D>::fail(const std::string& msg)
{
 printLaph(make_strf("DataGetHandlerSFNOM error: %s\n",msg));
 clearData(); 
 errorLaph("aborting");
}

   // **************************************************************
   // *                                                            *
   // *                    DataPutHandlerSFNOM                     *
   // *                                                            *
   // **************************************************************


template <typename H, typename R, typename D>
class DataPutHandlerSFNOM : public DataPutHandlerBaseSF<H,R,D>
{

    H& handler;
    std::map<R,D>* iomptr;
    std::string m_objname;
    bool m_overwrite;


 public:

    DataPutHandlerSFNOM(H& inptr, const std::string& file_name, 
                        const std::string& filetype_id,
                        bool overwrite=false);

    ~DataPutHandlerSFNOM() {}

    std::string getFileName() const {return m_objname;}

    bool isGlobal() const {return true;}

    bool isLocal() const {return false;}


    void putData(const R& rkey, const D& data);

    void flush();

    bool queryData(const R& rkey);


 private:
 
    void fail(const std::string& msg);

    void fail(const R& rkey);


          // disallow copies
    DataPutHandlerSFNOM(const DataPutHandlerSFNOM& in);
    DataPutHandlerSFNOM& operator=(const DataPutHandlerSFNOM& in);

};



   // constructor checks that the information in the header
   // is consistent

template <typename H, typename R, typename D>
DataPutHandlerSFNOM<H,R,D>::DataPutHandlerSFNOM(H& in_handler,
                                          const std::string& file_name,
                                          const std::string& filetype_id,
                                          bool overwrite)
                     :  handler(in_handler), m_objname(tidyString(file_name)),
                        m_overwrite(overwrite)
{
 if (m_objname.empty())
    fail("empty object name in DataGetHandlerSFNOM");
 if (NamedObjMap::query(m_objname)){
    XMLHandler xmlf;
    NamedObjMap::getXMLInfo(m_objname,xmlf);
    std::string ftype;
    xmlread(xmlf,"FileType",ftype,"DataPutHandlerSFNOM");
    if (ftype!=tidyString(filetype_id))
       fail("Wrong filetype id for this object in NamedObjMap\n");
    XMLHandler xmlr;
    NamedObjMap::getXMLInfo(m_objname,xmlr);
    if (!handler.checkHeader(xmlr)){
       fail("Header string is \n"+xmlr.str()+"\n for object "+m_objname
            +" does not match info in current Handler\n ...execution aborted...\n");}
    iomptr=&(NamedObjMap::getData<std::map<R,D>>(m_objname));}
 else{
 try {
    iomptr=&(NamedObjMap::insert<std::map<R,D>>(m_objname));
    XMLHandler xmlinfo;
    handler.writeHeader(xmlinfo);  // write header info 
    xmlinfo.put_child("FileType",tidyString(filetype_id));
    NamedObjMap::setXMLInfo(m_objname,xmlinfo);}
 catch(const std::exception& xp){
    fail("failure opening file "+m_objname+" in DataPutHandlerMFNOM");}}
}



template <typename H, typename R, typename D>
void DataPutHandlerSFNOM<H,R,D>::putData(const R& rkey, const D& data)
{
 if (m_overwrite)
    iomptr->operator[](rkey)=data;
 else if (iomptr->find(rkey)==iomptr->end())
    iomptr->insert(std::make_pair(rkey,data));
 else{ 
    fail(rkey);}
}


template <typename H, typename R, typename D>
void DataPutHandlerSFNOM<H,R,D>::flush()
{}


template <typename H, typename R, typename D>
bool DataPutHandlerSFNOM<H,R,D>::queryData(const R& rkey)
{
 return ((iomptr)->find(rkey)!=(iomptr)->end());
}



template <typename H, typename R, typename D>
void DataPutHandlerSFNOM<H,R,D>::fail(const std::string& msg)
{
 printLaph(make_strf("DataPutHandlerSFNOM error: %s\n",msg));
 errorLaph("aborting");
}


template <typename H, typename R, typename D>
void DataPutHandlerSFNOM<H,R,D>::fail(const R& rkey)
{
 std::cout << "DataPutHandlerSFNOM could not insert requested record:"<<std::endl;
 std::cout << " File name: "<< getFileName() <<std::endl;
 XMLHandler xmlout;
 xmlout.set_root("RecordKey");
 XMLHandler xmlh;
 rkey.output(xmlh);
 xmlout.put_child(xmlh);
 std::cout << xmlout.str()<<std::endl;
 std::cout << "...execution aborted..."<<std::endl;
 errorLaph("aborting");
}


// **************************************************************



   // If the stub of "in_filelist" begins with "NOM_", results are retrieved from
   // "NamedObjMap", otherwise, they are retrieved from files.  


template <typename H, typename F, typename R, typename D>
class DataGetHandlerMFO
{
    DataGetHandlerBaseMF<H,F,R,D> *base_ptr;
 public:
    DataGetHandlerMFO(H& in_handler, const FileListInfo& in_filelist,
                      const std::string& filetype_id, 
                      const std::string& header_tag, 
                      bool global_mode=true, bool use_checksums=false);
    ~DataGetHandlerMFO() {delete base_ptr;}
    const FileListInfo& getFileListInfo() const {return base_ptr->getFileListInfo();}
    bool isGlobal() const {return base_ptr->isGlobal();}
    bool isLocal() const {return base_ptr->isLocal();}
    bool queryData(const F& fkey, const R& rkey) {return base_ptr->queryData(fkey,rkey);}
    bool queryFile(const F& fkey) {return base_ptr->queryFile(fkey);}
    const D& getData(const F& fkey, const R& rkey) {return base_ptr->getData(fkey,rkey);}
    void removeData(const F& fkey, const R& rkey) {base_ptr->removeData(fkey,rkey);}
    void removeData(const F& fkey) {base_ptr->removeData(fkey);}
    void removeDataMEM(const F& fkey, const R& rkey) {base_ptr->removeDataMEM(fkey,rkey);}
    void removeDataMEM(const F& fkey) {base_ptr->removeDataMEM(fkey);}
    void clearData() {base_ptr->clearData();}
    void closeAll() {base_ptr->closeAll();}
    void clearDataMEM() {base_ptr->clearDataMEM();}
    void getFileMap(XMLHandler& xmlout) const {base_ptr->getFileMap(xmlout);}
    std::map<int,F> getSuffixMap() const {return base_ptr->getSuffixMap();}
    std::set<F> getFileKeys() const {return base_ptr->getFileKeys();}
    std::set<R> getKeys(const F& fkey) {return base_ptr->getKeys(fkey);}
    void outputKeys(XMLHandler& xmlout) {base_ptr->outputKeys(xmlout);}
};

   
template <typename H, typename F, typename R, typename D>
DataGetHandlerMFO<H,F,R,D>::DataGetHandlerMFO(H& in_handler,
                                              const FileListInfo& in_filelist,
                                              const std::string& filetype_id,
                                              const std::string& header_tag,
                                              bool global_mode, bool use_checksums)
{
 std::string stub=in_filelist.getFileStub();
 if (stub.find("NOM_")!=0){
    base_ptr=new DataGetHandlerMF<H,F,R,D>(in_handler,in_filelist,filetype_id,
                                           header_tag,global_mode,use_checksums);}
 else{
    if (!global_mode){
        errorLaph("Data to/from NamedObjMap in local mode not currently supported");}
    stub.erase(0,4);
    base_ptr=new DataGetHandlerMFNOM<H,F,R,D>(in_handler,
              FileListInfo(stub,in_filelist.getMinFileNumber(),in_filelist.getMaxFileNumber(),
                           in_filelist.isModeOverwrite()),
              filetype_id,header_tag);}
}
    


   // If the stub of "in_filelist" begins with "NOM_", results are inserted into
   // "NamedObjMap", otherwise, they are written to files.  


template <typename H, typename F, typename R, typename D>
class DataPutHandlerMFO
{
    DataPutHandlerBaseMF<H,F,R,D> *base_ptr;
 public:
    DataPutHandlerMFO(H& in_handler, const FileListInfo& in_filelist,
                      const std::string& filetype_id,  
                      const std::string& header_tag, 
                      bool global_mode=true, bool use_checksums=false,
                      int striping_factor=1, int striping_unit=0);
    ~DataPutHandlerMFO() {delete base_ptr;}
    void setOverWrite() {base_ptr->setOverWrite();}
    void setNoOverWrite() {base_ptr->setNoOverWrite();}
    const FileListInfo& getFileListInfo() const {return base_ptr->getFileListInfo();}
    bool isGlobal() const {return base_ptr->isGlobal();}
    bool isLocal() const {return base_ptr->isLocal();}
    void open(const F& fkey) {base_ptr->open(fkey);}
    void putData(const R& rkey, const D& data) {base_ptr->putData(rkey,data);}  // insert into current file
    void flush() {base_ptr->flush();}   // flush current file
    void close() {base_ptr->close();}   // close current file
    void putData(const F& fkey, const R& rkey, const D& data) {base_ptr->putData(fkey,rkey,data);}
    void flush(const F& fkey) {base_ptr->flush(fkey);}
    void close(const F& fkey) {base_ptr->close(fkey);}
    void flushAll() {base_ptr->flushAll();}
    void closeAll() {base_ptr->closeAll();}
    bool queryData(const F& fkey, const R& rkey) {return base_ptr->queryData(fkey,rkey);}
    bool queryData(const R& rkey) {return base_ptr->queryData(rkey);}  // query in current open file
    bool queryFile(const F& fkey) {return base_ptr->queryFile(fkey);}
    void getFileMap(XMLHandler& xmlout) const {base_ptr->getFileMap(xmlout);}
    std::map<int,F> getSuffixMap() const {return base_ptr->getSuffixMap();}
    std::set<F> getFileKeys() const {return base_ptr->getFileKeys();}
};

 
template <typename H, typename F, typename R, typename D>
DataPutHandlerMFO<H,F,R,D>::DataPutHandlerMFO(H& in_handler,
                                              const FileListInfo& in_filelist,
                                              const std::string& filetype_id,
                                              const std::string& header_tag,
                                              bool global_mode, bool use_checksums,
                                              int striping_factor, int striping_unit)
{
 std::string stub=in_filelist.getFileStub();
 if (stub.find("NOM_")!=0){
    base_ptr=new DataPutHandlerMF<H,F,R,D>(in_handler,in_filelist,filetype_id,
                                           header_tag,global_mode,use_checksums,
                                           striping_factor,striping_unit);}
 else{
    if (!global_mode){
        errorLaph("Data to/from NamedObjMap in local mode not currently supported");}
    stub.erase(0,4);
    base_ptr=new DataPutHandlerMFNOM<H,F,R,D>(in_handler,
              FileListInfo(stub,in_filelist.getMinFileNumber(),in_filelist.getMaxFileNumber(),
                           in_filelist.isModeOverwrite()),
              filetype_id,header_tag);}
}
    

 // *************************************************************************************
 

template <typename H, typename R, typename D>
class DataGetHandlerSFO
{
    DataGetHandlerBaseSF<H,R,D> *base_ptr;
 public:
    DataGetHandlerSFO(H& in_handler, const std::string& file_name, 
                      const std::string& filetype_id, bool global_mode=true,
                      bool use_checksums=false);
    ~DataGetHandlerSFO() {delete base_ptr;}
    std::string getFileName() const {return base_ptr->getFileName();}
    bool isGlobal() const {return base_ptr->isGlobal();}
    bool isLocal() const {return base_ptr->isLocal();}
    bool queryData(const R& rkey) {return base_ptr->queryData(rkey);}
    const D& getData(const R& rkey) {return base_ptr->getData(rkey);}
    void removeData(const R& rkey){base_ptr->removeData(rkey);}
    void removeDataMEM(const R& rkey){base_ptr->removeDataMEM(rkey);}
    void clearData() {base_ptr->clearData();}
    void clearDataMEM() {base_ptr->clearDataMEM();}
    std::set<R> getKeys() {return base_ptr->getKeys();}
    void outputKeys(XMLHandler& xmlout) {base_ptr->outputKeys(xmlout);}
};


template <typename H, typename R, typename D>
DataGetHandlerSFO<H,R,D>::DataGetHandlerSFO(H& in_handler, const std::string& filename,
                                            const std::string& filetype_id, 
                                            bool global_mode, bool use_checksums)
{
 std::string infile(tidyString(filename));
 if (infile.find("NOM_")!=0){
    base_ptr=new DataGetHandlerSF<H,R,D>(in_handler,infile,filetype_id,
                                         global_mode,use_checksums);}
 else{
    if (!global_mode){
        errorLaph("Data to/from NamedObjMap in local mode not currently supported");}
    infile.erase(0,4);
    base_ptr=new DataGetHandlerSFNOM<H,R,D>(in_handler,infile,filetype_id);}
}
    

template <typename H, typename R, typename D>
class DataPutHandlerSFO
{
    DataPutHandlerBaseSF<H,R,D> *base_ptr;
 public:
    DataPutHandlerSFO(H& inptr, const std::string& file_name, 
                      const std::string& filetype_id,
                      bool overwrite=false, bool global_mode=true,
                      bool use_checksums=false, int striping_factor=1,
                      int striping_unit=0);
    ~DataPutHandlerSFO() {delete base_ptr;}
    std::string getFileName() const {return base_ptr->getFileName();}
    bool isGlobal() const {return base_ptr->isGlobal();}
    bool isLocal() const {return base_ptr->isLocal();}
    void putData(const R& rkey, const D& data) {base_ptr->putData(rkey,data);}
    void flush() {base_ptr->flush();}
    bool queryData(const R& rkey) {return base_ptr->queryData(rkey);}  
};


template <typename H, typename R, typename D>
DataPutHandlerSFO<H,R,D>::DataPutHandlerSFO(H& in_handler,
                                            const std::string& file_name,
                                            const std::string& filetype_id,
                                            bool overwrite, bool global_mode,
                                            bool use_checksums, int striping_factor,
                                            int striping_unit)
{
 std::string infile(tidyString(file_name));
 if (infile.find("NOM_")!=0){
    base_ptr=new DataPutHandlerSF<H,R,D>(in_handler,infile,filetype_id,
                                         overwrite,global_mode,use_checksums,
                                         striping_factor,striping_unit);}
 else{
    if (!global_mode){
        errorLaph("Data to/from NamedObjMap in local mode not currently supported");}
    infile.erase(0,4);
    base_ptr=new DataPutHandlerSFNOM<H,R,D>(in_handler,infile,filetype_id,overwrite);}
}
    



// **************************************************************
}
#endif
