#ifndef MARVEL_CONFIG_H
#define MARVEL_CONFIG_H
#include <chrono>
#include <iostream>
namespace marvel{

//TIMING
#define __TIME_BEGIN__\
  {auto _start = chrono::system_clock::now();

#define __TIME_END__(_str)\
  auto _end = chrono::system_clock::now();\
  auto _duration = chrono::duration_cast<chrono::microseconds>(_end - _start); \
  cout << _str << " cost " <<\
  double(_duration.count()) *chrono::microseconds::period::num / chrono::microseconds::period::den \
  << "seconds" << endl;}



}
#endif
