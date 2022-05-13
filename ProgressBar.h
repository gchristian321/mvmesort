#ifndef ProgressBar_HEADER_
#define ProgressBar_HEADER_
#include <chrono>
#include <string>

class ProgressBar {
public:
  ProgressBar(long long max_entries, int bar_width = 40):
    progress(0.),
    barWidth(bar_width),
    entryCount(max_entries),
    t1(std::chrono::system_clock::now())
  { }
  void operator()(long long entry)
  {
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
      if(!i) std::cout << "[";
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
    }
    auto t2 = std::chrono::system_clock::now();
    double elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
    double total_time = elapsed_time / progress;
    double time_left = total_time - elapsed_time;
    time_left/=1e3; // --> seconds
    auto time_remaining_message = ConvertTime(time_left);
    if(time_left > 100*3600 || std::isnan(time_left) || std::isinf(time_left))
      { time_remaining_message = "00:00:00"; }
    std::cout << "] " << int(progress * 100.) << " % [" <<
      time_remaining_message << "]\r";
    std::cout.flush();
    progress = (double)entry/entryCount;
  }
  void Reset()
  {
    progress = 0.;
    t1 = std::chrono::system_clock::now();
  }
  std::string ConvertTime(double time_left)
  {
    if(time_left < 1) return "00:00:00";
    int seconds=0,minutes=0,hours=0;
    if(time_left >= 3600) {
      hours = int(time_left/3600);
      time_left -= hours*3600;
    }
    if(time_left >= 60) {
      minutes = int(time_left/60);
      time_left -= minutes*60;
    }
    seconds = time_left;
    std::string time_remaining_message =
      Form("%.2i:%.2i:%.2i",hours,minutes,seconds);
    return time_remaining_message;
  }
private:
  float progress;
  int barWidth;
  long long entryCount;
  std::chrono::time_point<std::chrono::system_clock, std::chrono::duration<long, std::ratio<1, 1000000000> > > t1;
};

#endif
