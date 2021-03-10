#ifndef ABSMC_FILE_UTIL_H
#define ABSMC_FILE_UTIL_H

#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cassert>

namespace util  {

std::string createFileName(std::string const& path, std::string const& baseName, std::string const& extension="", int tag=0, int tagWidth=0)
{
    assert(tag >= 0);
    assert(tagWidth >= 0);

    static const char pathSep = '/';

    std::stringstream fNameStream;
    if (path.length() > 0) { fNameStream << path << pathSep; }
    fNameStream << baseName;
    if (tagWidth>0) {  fNameStream << "_" << std::setfill('0') << std::setw(tagWidth) << tag; }
    if (extension.length() > 0) { fNameStream << '.' << extension; }
    return fNameStream.str();
}


struct FileNameParts {
    FileNameParts(std::string const& directory_, std::string const& baseName_, std::string const& extension_, std::string const& fileName_)
        : directory(directory_), baseName(baseName_), extension(extension_), fileName(fileName_) { }

    std::string directory;
    std::string baseName;
    std::string extension;
    std::string fileName;
};


FileNameParts parseFileName(std::string const& fullName)
{
    static const char dirSep = '/';
    static const char extSep = '.';

    size_t lastDirSepPos = fullName.find_last_of(dirSep);
    const std::string directory = fullName.substr(0, lastDirSepPos+1);
    const std::string fileName = fullName.substr(lastDirSepPos+1);

    size_t lastExtSepPos = fileName.find_last_of(extSep);
    const std::string baseName = fileName.substr(0, lastExtSepPos);
    const std::string extension = fileName.substr(lastExtSepPos+1);

    return FileNameParts(directory, baseName, extension, fileName);
}

} // end namespace util

#endif
