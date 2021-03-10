#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "config.h"

namespace util {

//////////////////// Class Config: implementation ///////////////////////////

void Config::readFile(std::string const& fileName)
{

    std::ifstream configFile(fileName.c_str());
    if (!configFile.is_open()) {
        error("readFile()","Could not open config file "+fileName, true);
        return;
    }
    parseInputStream(configFile);
    configFile.close();
}

void Config::readString(std::string const& inputString)
{
    std::istringstream inputStrStream(inputString);
    parseInputStream(inputStrStream);
}

template <typename T>
void Config::add(std::string const& key, T const& val)
{
    std::ostringstream oss;
    oss << val;
    cfg_map.insert(std::make_pair(key,oss.str()));
}

void Config::printAll(std::ostream & ostr) const
{
    std::map<std::string,std::string>::const_iterator it;
    for (it=cfg_map.begin(); it!=cfg_map.end(); it++) {
        ostr << (*it).first << "=" << (*it).second << "\n";
    }
    return;
}

bool Config::hasValue(std::string const& cfg_key) const
{
    std::map<std::string,std::string>::const_iterator it = cfg_map.find(cfg_key);
    return (it!=cfg_map.end());
}

template<typename T>
T Config::getValue(std::string const& cfg_key) const
{

    std::map<std::string,std::string>::const_iterator it = cfg_map.find(cfg_key);
    if (it==cfg_map.end()) {
        error("getValue()","key \""+cfg_key+"\" not found.", true);
    }

    std::string cfg_val = (*it).second;
    std::istringstream iss(cfg_val);
    T val;
    iss >> val;
    return val;
}

void Config::parseInputStream(std::istream& inputStream)
{
    std::string line;
    std::string::size_type loc;

    while (!inputStream.eof() ) {
        // read one line at a time
        getline(inputStream,line);
        // strip off comments (starting by "#")
        loc = line.find("#", 0);
        if (loc!=std::string::npos) line.erase(loc);

        // ignore empty lines
        loc = line.find_first_not_of(" ", 0);
        if (loc==std::string::npos) continue;
    // ignore lines with {

        // split line at first blank
        loc=line.find(" ", 0);
        if (loc==std::string::npos || loc==0 || loc+1==line.length()) {
            error("parseInputStream()","ignoring line: "+line, false);
            continue;
        }
        std::string key = line.substr(0,loc);
        std::string val = line.substr(loc+1);
        add(key,val);
    }
    return;
}

void Config::error(std::string const& methodName, std::string const& message, bool doExit) const
{
    std::cerr << "Config::" << methodName << ":" << message << std::endl;
    if (doExit) exit(EXIT_FAILURE);
}

//////////////////// Class Config: template instantiation ///////////////////

template
void Config::add(std::string const& key, bool const& val);

template
void Config::add(std::string const& key, char const& val);
template
void Config::add(std::string const& key, unsigned char const& val);

template
void Config::add(std::string const& key, int const& val);
template
void Config::add(std::string const& key, unsigned int const& val);

template
void Config::add(std::string const& key, long int const& val);
template
void Config::add(std::string const& key, unsigned long int const& val);

template
void Config::add(std::string const& key, short int const& val);
template
void Config::add(std::string const& key, unsigned short int const& val);

template
void Config::add(std::string const& key, float const& val);

template
void Config::add(std::string const& key, double const& val);



template
bool  Config::getValue(std::string const& cfg_key) const;

template
char  Config::getValue(std::string const& cfg_key) const;
template
unsigned char Config::getValue(std::string const& cfg_key) const;

template
int  Config::getValue(std::string const& cfg_key) const;
template
unsigned int Config::getValue(std::string const& cfg_key) const;

template
long int  Config::getValue(std::string const& cfg_key) const;
template
unsigned long int Config::getValue(std::string const& cfg_key) const;

template
short int  Config::getValue(std::string const& cfg_key) const;
template
unsigned short int Config::getValue(std::string const& cfg_key) const;

template
float  Config::getValue(std::string const& cfg_key) const;
template
double Config::getValue(std::string const& cfg_key) const;

template
std::string Config::getValue(std::string const& cfg_key) const;

} // end namespace util

