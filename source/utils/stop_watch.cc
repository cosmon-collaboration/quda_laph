#include "stop_watch.h"
#include "laph_stdio.h"

namespace LaphEnv {

void StopWatch::start() {
  if (m_timer.running) {
    printLaph("Warning: calling start() on running StopWatch");
  } else {
    m_timer.start();
  }
}

void StopWatch::stop() {
  if (m_timer.running) {
    m_timer.stop();
  } else {
    printLaph("Warning: calling stop() on stopped StopWatch");
  }
}

void StopWatch::reset() { m_timer.reset(nullptr, nullptr, 0); }

double StopWatch::getTimeInSeconds() const {
  if (m_timer.running) {
    errorLaph("Can only query total time from a stopped StopWatch");
  }
  return m_timer.time;
}

double StopWatch::getLastIntervalInSeconds() const {
  return m_timer.last_interval;
}

bool StopWatch::isRunning() const { return m_timer.running; }

} // namespace LaphEnv
