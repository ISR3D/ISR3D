#ifndef ABSMC_LOGGER_H
#define ABSMC_LOGGER_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "timer.h"

namespace util {

class Logger : public std::ostream {
public:
    friend Logger& logger();

    /// constructor: create config object and load configuration from ascii file
    Logger(std::string const& fileName ="") :
        std::ostream(&buffer), buffer(logFile) {
        if (fileName != "") {
            logFile.open(fileName);
        }
    }

    /// set log file name
    void setFile(std::string const& fileName) {
        logFile.close();
        logFile.open(fileName);
        if (!logFile.is_open()) {
            throw std::runtime_error("Unable to open log file");
        }
    }

    void closeFile() {
        logFile.close();
    }

    /// print log line
    void print(std::string const& inputString);

    /// close file on exit
    ~Logger() {
        logFile.close();
    }

private:
    class LogStreamBuf: public std::stringbuf
    {
        std::ofstream& logFile;
    public:
        LogStreamBuf(std::ofstream& str)
            :logFile(str)
        {}

        ~LogStreamBuf() {
            if (pbase() != pptr()) {
                putOutput();
            }
        }

        virtual int sync() {
            putOutput();
            return 0;
        }

        void putOutput() {
            // Called by destructor.
            // destructor can not call virtual methods.
            std::string tstamp = Timer::printTime() + ": ";
            std::cout << tstamp << str();
            if (logFile.is_open()) {
                logFile << tstamp << str();
                logFile.flush();
            }
            str("");
        }
    };

    LogStreamBuf buffer;

    std::string fileName;
    std::ofstream logFile;
};

/// This makes the logger executable-specific
inline Logger& logger() {
    static Logger instance;
    return instance;
}

void Logger::print(std::string const& inputString){
    std::cout << inputString << std::endl;
    if (logFile.is_open()) {
        logFile << inputString << std::endl;
    }
}

} // end namespace util

#endif
