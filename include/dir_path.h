#ifndef DIR_PATH_H
#define DIR_PATH_H

namespace LaphEnv {

class DirPath {
  unsigned long store; // every 3 bits is a direction starting from right
                       // maximum of 10 directions (since 32 bits)

public:
  DirPath() : store(0) {}

  DirPath(const DirPath &in_path) : store(in_path.store) {}

  DirPath &operator=(const DirPath &in_path) {
    store = in_path.store;
    return *this;
  }

  DirPath &operator=(int in_dir) {
    store = 0;
    add_dir(in_dir);
    return *this;
  }

  DirPath(int dir) : store(0) { add_dir(dir); }

  DirPath(int dir1, int dir2) : store(0) {
    add_dir(dir1);
    add_dir(dir2);
  }

  DirPath(int dir1, int dir2, int dir3) : store(0) {
    add_dir(dir1);
    add_dir(dir2);
    add_dir(dir3);
  }

  DirPath(int dir1, int dir2, int dir3, int dir4) : store(0) {
    add_dir(dir1);
    add_dir(dir2);
    add_dir(dir3);
    add_dir(dir4);
  }

  DirPath(int dir1, int dir2, int dir3, int dir4, int dir5) : store(0) {
    add_dir(dir1);
    add_dir(dir2);
    add_dir(dir3);
    add_dir(dir4);
    add_dir(dir5);
  }

  DirPath(int dir1, int dir2, int dir3, int dir4, int dir5, int dir6)
      : store(0) {
    add_dir(dir1);
    add_dir(dir2);
    add_dir(dir3);
    add_dir(dir4);
    add_dir(dir5);
    add_dir(dir6);
  }

  DirPath(int dir1, int dir2, int dir3, int dir4, int dir5, int dir6, int dir7)
      : store(0) {
    add_dir(dir1);
    add_dir(dir2);
    add_dir(dir3);
    add_dir(dir4);
    add_dir(dir5);
    add_dir(dir6);
    add_dir(dir7);
  }

  DirPath(int dir1, int dir2, int dir3, int dir4, int dir5, int dir6, int dir7,
          int dir8)
      : store(0) {
    add_dir(dir1);
    add_dir(dir2);
    add_dir(dir3);
    add_dir(dir4);
    add_dir(dir5);
    add_dir(dir6);
    add_dir(dir7);
    add_dir(dir8);
  }

  DirPath(const std::vector<int> &dirs) : store(0) {
    for (int k = 0; k < int(dirs.size()); ++k)
      add_dir(dirs[k]);
  }

  DirPath(const XMLHandler &xml_in) : store(0) {
    XMLHandler xmlr(xml_in, "DirPath");
    std::vector<int> dirs;
    try {
      xmlread(xmlr, "Directions", dirs, "DirPath");
    } catch (const std::exception &xp) {
      errorLaph("could not read DirPath");
    }
    for (int k = 0; k < int(dirs.size()); ++k)
      add_dir(dirs[k]);
  }

  DirPath &addDir(int dir) {
    add_dir(dir);
    return *this;
  }

  int Length() const {
    int count = 0;
    unsigned long temp = store;
    while (temp != 0) {
      temp >>= 3;
      ++count;
    }
    return count;
  }

  bool operator<(const DirPath &rhs) const { return (store < rhs.store); }

  bool operator==(const DirPath &rhs) const { return (store == rhs.store); }

  bool operator!=(const DirPath &rhs) const { return (store != rhs.store); }

  class const_iterator {
    static const unsigned long rightmask = 0x7ul;
    unsigned long m_store;
    const_iterator(unsigned long in_store) {
      unsigned long temp = in_store;
      m_store = 0;
      while (temp != 0) {
        m_store <<= 3;
        m_store |= (temp & rightmask);
        temp >>= 3;
      }
    }

  public:
    const_iterator() {}
    const_iterator &operator=(const const_iterator &rhs) {
      m_store = rhs.m_store;
      return *this;
    }
    int operator*() const { return int(m_store & rightmask) - 4; }
    const_iterator &operator++() {
      m_store >>= 3;
      return *this;
    } // prefix
    const_iterator &operator++(int) {
      m_store >>= 3;
      return *this;
    } // postfix
    bool operator==(const const_iterator &rhs) const {
      return (m_store == rhs.m_store);
    }
    bool operator!=(const const_iterator &rhs) const {
      return (m_store != rhs.m_store);
    }
    friend class DirPath;
  };

  class const_reverse_iterator {
    static const unsigned long rightmask = 0x7ul;
    unsigned long m_store;
    const_reverse_iterator(unsigned long in_store) : m_store(in_store) {}

  public:
    const_reverse_iterator() {}
    const_reverse_iterator &operator=(const const_reverse_iterator &rhs) {
      m_store = rhs.m_store;
      return *this;
    }
    int operator*() const { return int(m_store & rightmask) - 4; }
    const_reverse_iterator &operator++() {
      m_store >>= 3;
      return *this;
    } // prefix
    const_reverse_iterator &operator++(int) {
      m_store >>= 3;
      return *this;
    } // postfix
    bool operator==(const const_reverse_iterator &rhs) const {
      return (m_store == rhs.m_store);
    }
    bool operator!=(const const_reverse_iterator &rhs) const {
      return (m_store != rhs.m_store);
    }
    friend class DirPath;
  };

  const_iterator begin() const { return const_iterator(store); }

  const_iterator end() const { return const_iterator(0); }

  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(store);
  }

  const_reverse_iterator rend() const { return const_reverse_iterator(0); }

  std::string output() const {
    std::ostringstream oss;
    oss << "[";
    for (DirPath::const_iterator it = begin(); it != end(); it++) {
      oss << " " << *it;
    }
    oss << "]";
    return oss.str();
  }

  void output(XMLHandler &xmlout) const {
    xmlout.set_root("DirPath");
    std::string dirs;
    for (DirPath::const_iterator it = begin(); it != end(); it++) {
      dirs += make_string(*it) + " ";
    }
    xmlout.put_child("Directions", dirs);
  }

private:
  void add_dir(const int dir) {
    static const unsigned long rightmask = 0x7ul;
    static const unsigned long full = 0x8000000ul;
    static const unsigned long zero = 0x8ul;
    if ((dir < -3) || (dir > 3)) {
      errorLaph("bad direction in dirPath");
    }
    if (dir == 0)
      return;
    unsigned long udir = (unsigned long)(dir + 4);
    if (store == 0) {
      store = udir;
      return;
    }
    if ((store & rightmask) + udir == zero) {
      store >>= 3;
      return;
    }
    if (store > full) {
      errorLaph("cannot add more directions in dirPath");
    } else {
      store <<= 3;
      store |= udir;
    }
  }
};

//  This takes a set and assigns an integer index ordering
//  to the elements so it "looks" like a vector.  Usage:
//
//     set<Complex> a_set;
//     SetIndexer<Complex> a_vec(a_set);
//
//     Complex z = a_vec[0];
//
//  Warning: do NOT change the underlying set after creating
//  an object of this type!!

template <typename T> class SetIndexer {
  std::vector<typename std::set<T>::const_iterator> m_data;

  // disable copy, default constructor
  SetIndexer();
  SetIndexer(const SetIndexer<T> &rhs);
  SetIndexer &operator=(const SetIndexer<T> &rhs);

public:
  SetIndexer(const std::set<T> &data) {
    m_data.resize(data.size());
    int index = 0;
    for (typename std::set<T>::const_iterator it = data.begin();
         it != data.end(); it++, index++)
      m_data[index] = it;
  }

  ~SetIndexer() {}

  const T &operator[](int k) const { return *m_data[k]; }

  int size() const { return m_data.size(); }
};

} // namespace LaphEnv
#endif
