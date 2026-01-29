#pragma once

#include <iomanip>
#include <sys/time.h>

/**
 * \brief Time stamp class
 *
 * This class is used to store time stamps and compute elapsed times.
 * Simplified version of the FEAT3 class of the same name.
 *
 * \author Dirk Ribbrock
 * \author Peter Zajac
 */
class TimeStamp
{
private:
  /// Our time-stamp.
  timeval _time;

public:
  /// Constructor.
  TimeStamp()
  {
    stamp();
  }

  /**
   * \brief Stamps the current time-stamp.
   *
   * This function updates the time-stamp to the current time.
   *
   * \returns \c *this
   */
  TimeStamp& stamp()
  {
    gettimeofday(&_time, 0);
    return *this;
  }

  /**
   * \brief Calculates the time elapsed between two time stamps.
   *
   * \param[in] before
   * A time stamp that represents a previous moment.
   *
   * \returns
   * The time elapsed between the time stamps \p before and \c this in seconds.
   */
  double elapsed(const TimeStamp& before) const
  {
    return double(_time.tv_sec - before._time.tv_sec) + 1E-6 * double(_time.tv_usec - before._time.tv_usec);
  }

  /**
   * \brief Calculates the time elapsed between the time stamp and now.
   *
   * \note
   * This function does \b not update (stamp) this time stamp.
   *
   * \returns
   * The time elapsed between this time stamp and now in seconds.
   */
  double elapsed_now() const
  {
    return TimeStamp().elapsed(*this);
  }

  /**
   * \brief Calculate the time elapsed between two time stamps in microseconds.
   *
   * \param[in] before
   * A time stamp that represents a previous moment.
   *
   * \returns
   * The time elapsed between the time stamps \p before and \c this in microseconds.
   */
  long long elapsed_micros(const TimeStamp& before) const
  {
    return 1000000ll * (long long)(_time.tv_sec - before._time.tv_sec) +
           (long long)(_time.tv_usec - before._time.tv_usec);
  }

  /**
   * \brief Calculates the time elapsed between the time stamp and now in microseconds.
   *
   * \note
   * This function does \b not update (stamp) this time stamp.
   *
   * \returns
   * The time elapsed between this time stamp and now in microseconds.
   */
  long long elapsed_micros_now() const
  {
    return TimeStamp().elapsed_micros(*this);
  }

  /**
   * \brief Formats an elapsed time in microseconds as a string.
   *
   * This function formats the given elapsed time as a string.
   *
   * \param[in] micros
   * The elapsed time to be formatted.
   *
   * \returns
   * The time elapsed as a formatted string, for example <code>h:m:ss.mmm</code>.
   */
  static std::string format_micros(long long micros)
  {
    std::ostringstream oss;
    oss << std::setfill('0');

    // check whether the time is negative
    if(micros < 0)
    {
      oss << '-';
      micros = -micros;
    }

    long long hours = micros / 3600000000ll;
    long long minutes = (micros / 60000000ll) % 60ll;
    long long seconds = (micros / 1000000ll) % 60ll;
    long long milliseconds = (micros / 1000ll) % 1000ll;

    if(hours > 0)
    {
      oss << std::setw(2) << hours << ":";
    }

    if(hours > 0 || minutes > 0)
    {
      oss << std::setw(2) << minutes << ":";
    }

    if(hours > 0 || minutes > 0 || seconds > 0)
    {
      oss << std::setw(2) << seconds << "." << std::setw(3) << milliseconds;
    }

    // return formatted string
    return oss.str();
  }

  /**
   * \brief Return the time elapsed between two time stamps as a string.
   *
   * See TimeStamp::format_micros() for more information about the formatting options.
   *
   * \param[in] before
   * A time stamp that represents a previous moment.
   *
   * \param[in] format
   * Specifies string formatting to be used.
   *
   * \returns
   * The time elapsed between the time stamps \p before and \c this as a formatted string,
   * for example <code>h:mm:ss.mm</code>.
   */
  std::string elapsed_string(const TimeStamp& before) const
  {
    return format_micros(elapsed_micros(before));
  }

  /**
   * \brief Calculates the time elapsed between the time stamp and now as a string.
   *
   * \note
   * This function does \b not update (stamp) this time stamp.
   *
   * See TimeStamp::format_micros() for more information about the formatting options.
   *
   * \param[in] format
   * Specifies string formatting to be used.
   *
   * \returns
   * The time elapsed between this time stamp and now as a formatted string,
   * for example <code>h:mm:ss.mm</code>.
   */
  std::string elapsed_string_now() const
  {
    return TimeStamp().elapsed_string(*this);
  }

  /**
   * \brief Comparison operator.
   *
   * \param[in] other
   * Another time-stamp.
   *
   * \returns
   * \c true, if \c this time-stamp has been taken earlier than \p other, otherwise \c false.
   */
  bool operator<(const TimeStamp& other) const
  {
    return (_time.tv_sec < other._time.tv_sec) ||
           ((_time.tv_sec == other._time.tv_sec) && (_time.tv_usec < other._time.tv_usec));
  }
}; // class TimeStamp

/**
 * \brief Stop-Watch class
 *
 * This class implements an incremental stop-watch.
 *
 * \author Peter Zajac
 */
class StopWatch
{
private:
  TimeStamp _start_time;
  TimeStamp _stop_time;
  bool _running;
  long long _micros;

public:
  StopWatch() : _running(false), _micros(0ll)
  {
  }

  /// Resets the elapsed time.
  void reset()
  {
    _running = false;
    _micros = 0ll;
  }

  /// Starts the stop-watch.
  void start()
  {
    if(!_running)
      _start_time.stamp();
    _running = true;
  }

  /// Stops the stop-watch and increments elapsed time.
  void stop()
  {
    if(_running)
    {
      _stop_time.stamp();

      // update elapsed time
      _micros += _stop_time.elapsed_micros(_start_time);
    }
    _running = false;
  }

  /// Returns \c true if the stop-watch is currently running.
  bool running() const
  {
    return _running;
  }

  /// Returns the total elapsed time in seconds.
  double elapsed() const
  {
    return 1E-6 * double(elapsed_micros());
  }

  /// Returns the total elapsed time in micro-seconds.
  long long elapsed_micros() const
  {
    if(!_running)
      return _micros;
    else
    {
      // stop-watch is currently running, so compute current time
      return _micros + TimeStamp().elapsed_micros(_start_time);
    }
  }

  /**
   * \brief Return the time elapsed in the stop-watch.
   *
   * See TimeStamp::format_micros() for more information about the formatting options.
   *
   * \param[in] format
   * Specifies the output string format to be used.
   *
   * \returns
   * The time elapsed between in the stop-watch as a formatted string <code>h:mm:ss.nnn</code>.
   */
  std::string elapsed_string() const
  {
    return TimeStamp::format_micros(elapsed_micros());
  }

  /**
   * \brief Returns the formatted percentage that this stop watch elapsed time relative to another total runtime stop
   * watch
   *
   * \param[in] total_watch
   * The stop watch that represents the total runtime
   *
   * \param[in] prefix
   * The prefix string
   *
   * \param[in] postfix
   * The postfix string
   *
   * \returns The formatted percentage of this stop watch elapsed time relative to watch_total
   */
  std::string
  percent_of(const StopWatch& total_watch, const std::string& prefix = " [", const std::string& postfix = "% ]") const
  {
    return prefix + std::to_string(100.0 * this->elapsed() / total_watch.elapsed()) + postfix;
  }
}; // class StopWatch
