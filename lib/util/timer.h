#ifndef ABSMC_TIMER_H
#define ABSMC_TIMER_H

#include <string>
#include <chrono>
#include <sstream>

namespace util {

class Timer {
public:

    /// constructor: create timer and start it
    Timer() {
        start();
    }

    /// print time
    static std::string printTime() {
        auto temp = std::chrono::system_clock::to_time_t( std::chrono::system_clock::now() );
        std::stringstream out;
        out << std::put_time(std::localtime(&temp), "%F %T");
        return out.str();
    }

    void start() {
        t1 = std::chrono::system_clock::now();
    }

    void stop() {
        t2 = std::chrono::system_clock::now();
    }

    /// print elapsed time (seconds)
    std::string printElapsedTime() {
        auto t_elapsed = t2 - t1;
        std::stringstream out;
        out << std::chrono::duration_cast<std::chrono::seconds>( t_elapsed ).count();
        return out.str() + " s";
    }

private:
    std::chrono::time_point<std::chrono::system_clock> t1;
    std::chrono::time_point<std::chrono::system_clock> t2;
};

} // end namespace util

#endif
