#ifndef STOP_WATCH_H
#define STOP_WATCH_H

#include "timer.h"

namespace LaphEnv {

// ***********************************************************************
// *                                                                     *
// *   Objects of this class can be used for getting timings on the      *
// *   host.  Typical usage is                                           *
// *                                                                     *
// *      StopWatch rolex;                                               *
// *      rolex.start();   // starts the timer                           *
// *      ....do something...                                            *
// *      rolex.stop();   // halts the timer                             *
// *      cout << "total time is "<<rolex.getTimeInSeconds()<<endl;      *
// *                                                                     *
// *   The member "getTimeInSeconds()" returns the total accumulated     *
// *   time.  You can call a "stop", then a subsequent "start", and      *
// *   the time will be accumulated.  Use "reset" to set the timer       *
// *   back to zero.  "getLastIntervalInSeconds()" returns the           *
// *   time reported between the last start and stop.  Attempting        *
// *   to reset a running timer results in a fatal execution error.      *
// *                                                                     *
// ***********************************************************************

class StopWatch {

  quda::Timer<false> m_timer;

  StopWatch(const StopWatch &intimer) = delete; // no copy constructor

  StopWatch &operator=(const StopWatch &in) = delete;

public:
  StopWatch() {}

  ~StopWatch() {}

  void start(); // warning if running

  void stop(); // warning if stopped

  void reset(); // error if running (must be stopped)

  double getTimeInSeconds() const; // error if running (must be stopped)

  double getLastIntervalInSeconds() const; // ok if running

  bool isRunning() const;
};

// **********************************************************************
} // namespace LaphEnv
#endif
