#pragma once
#include "logger.h"
#include <sys/stat.h>  
#include <string>
#include <sstream>
#include <vector>

namespace utils{
     
    template<typename T, typename... Args>
    std::unique_ptr<T> make_unique(Args&&... args) {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }


    inline void checkFolder(const std::string& path) {
        struct stat info;
        if (stat(path.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR)) {
            try {
                system(("mkdir -p " + path).c_str());
                Logger::info("Folder created");
            } catch (const std::runtime_error& e) {
                Logger::error(e.what());
                exit(EXIT_FAILURE);
            } 
        }
    }
    
    inline std::string to_string(std::vector<double> A, std::string sep= ", "){ 
        std::ostringstream oss;
        for (const double& val : A){
            oss << val << sep;
        } 
        return oss.str();
    }
}