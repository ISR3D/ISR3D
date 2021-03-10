#ifndef ABSMC_MISC_UTIL_H
#define ABSMC_MISC_UTIL_H

#include <cstddef>
#include <vector>
#include <sstream>

namespace util  {

void tokenizeString(std::string const& str, std::vector<std::string>& tokens, char delimiter)
{
    std::istringstream iss(str);
    std::string token;
    while (std::getline(iss, token, delimiter)) {
        if (!token.empty() ) {tokens.push_back(token); }
    }
}

template<typename T>
void tokenizeAndConvert(std::string const& str, std::vector<T>& tokens, char delimiter)
{
    std::istringstream iss(str);
    std::string tokenAsString;
    T tokenAsT;
    while (std::getline(iss, tokenAsString, delimiter)) {
        if (!tokenAsString.empty() ) {
            std::istringstream tokenAsStringStream(tokenAsString);
            tokenAsStringStream >> tokenAsT;
            tokens.push_back(tokenAsT);
        }
    }
}

} // end namespace util

#endif
