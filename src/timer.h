/**
 * @file timer.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-09-15
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#include<chrono>
#include<ctime>
#include<iostream>

using namespace std::chrono;

class ChronoTimer{

  system_clock::time_point _start;
  system_clock::time_point _end;
  std::time_t _t;

public:

  ChronoTimer()
  {
    _t = system_clock::to_time_t(system_clock::now());
  }

  void
  start()
  {
    _start = system_clock::now();
  }

  void 
  end()
  {
    _end = system_clock::now();
  }

  double
  elapsed_msec() const
  {
    return duration_cast<milliseconds>(_end-_start).count();
  }

  double
  elapsed_sec() const
  {
    return duration_cast<seconds>(_end-_start).count();
  }

  double
  elapsed_microsec() const
  {
    return duration_cast<microseconds>(_end-_start).count();
  }

  void
  stamp_time() const
  {
    std::cout << std::ctime(&_t) <<std::endl;
  }

};

