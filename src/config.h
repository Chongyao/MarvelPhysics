#ifndef MARVEL_CONFIG_H
#define MARVEL_CONFIG_H
#include <chrono>
#include <iostream>
#include <list>
#include <string>
#include <assert.h> 

namespace marvel{

//TIMING
class TIMING{
 public :
  static void begin(){
    starts_.push_back(std::chrono::system_clock::now());
  }
  static void end(const std::string& info){
    assert(!starts_.empty());
      
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - starts_.back());
    std::cout << info << " cost "
              << double(duration.count()) *std::chrono::microseconds::period::num / std::chrono::microseconds::period::den
         << "seoncds" << std::endl;
    starts_.pop_back();
  }
 private:
  static std::list<std::chrono::system_clock::time_point> starts_;
  
};

#define __TIME_BEGIN__\
  TIMING::begin();

#define __TIME_END__(_str)\
  TIMING::end(_str);



}
#endif
