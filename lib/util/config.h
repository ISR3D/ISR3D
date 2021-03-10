#ifndef ABSMC_CONFIG_H
#define ABSMC_CONFIG_H

#include <string>
#include <iostream>
#include <map>

namespace util {

/// A quick-and dirty ascii config file parser.
/** Config file syntax is as follows:
 *    - empty lines are ignored
 *    - parts of a line starting with "#" are ignored (regarded as comments)
 *       ("#" need not be at beginning)
 *    - other lines are regarded as statements taking the form "key value"
 *    - strings must not contain spaces; these will be stripped off
 *    - supported types are all built-in types; for other types to work
 *      correctly, the getValue() template must be specialized
 */
class Config {
    public:
        /// constructor: create config object and load configuration from ascii file
        Config(std::string const& fileName="")
        {
            if (fileName!="") readFile(fileName);
        }
        /// load configuration from text file
        void readFile(std::string const& fileName);
        /// load configuration from text file
        void readString(std::string const& inputString);
        /// add a key/value pair
        template <typename T>
        void add(std::string const& key, T const& val);
        /// print all (key,value) pairs
        void printAll(std::ostream& ostr) const;
        /// get value of parameter indexed by cfg_key
        template <typename T> T getValue(std::string const& cfg_key) const;
        /// test if a key exists
        bool hasValue (std::string const& cfg_key) const;
private:
    std::map<std::string,std::string> cfg_map;
    void parseInputStream(std::istream& inputStream) ;
    void error(std::string const& methodName, std::string const& message, bool doExit) const;
private:
    friend Config& config();
};

inline Config& config()
{
    static Config instance;
    return instance;
}

} // end namespace util

#endif
