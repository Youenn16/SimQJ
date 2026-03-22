#pragma once
#include <string>
#include <iostream>
#include <iomanip>
#include <type_traits>
#include "sim_math.h"
#include "version.h"

#define RED     "\033[31m"
#define GREEN     "\033[32m"
#define BLUE     "\033[34m"
#define MAGENTA     "\033[35m"
#define YELLOW     "\033[33m"
#define RESET   "\033[0m" 

class Logger {
public:

    template<typename... Args>
    static void log(const Args&... args) {
        std::cout << YELLOW << "[LOG] " << RESET ;
        using expander = int[];
        (void)expander{0, (print_item(args), 0)...};
        
        std::cout << std::endl;
    }


    template<typename... Args>
    static void debug(const Args&... args) {
        std::cout << MAGENTA << "[DEBUG] " << RESET ;
        using expander = int[];
        (void)expander{0, (print_item(args), 0)...};
        
        std::cout << std::endl;
    }

    template<typename... Args>
    static void error(const Args&... args) {
        std::cerr << RED << "[ERROR] " << RESET ;
        using expander = int[];
        (void)expander{0, (print_item(args), 0)...};
        
        std::cout << std::endl;
    }

     template<typename... Args>
    static void info(const Args&... args) {
        std::cout << GREEN << "[INFO] " << RESET ;
        using expander = int[];
        (void)expander{0, (print_item(args), 0)...};
        
        std::cout << std::endl;
    }


    template<typename T>
    static void appendArg(std::ostringstream& oss, const T& arg) {
        typedef typename std::decay<T>::type DecayedT; 
        if (std::is_same<DecayedT, std::string>::value || std::is_same<DecayedT, const char*>::value || std::is_same<DecayedT, char*>::value) { 
            oss  << arg;
        } else { 
            oss << arg<< ",    " ;
        }
    }
 

    template<typename... Args>
    static void info_big(std::string title,std::string color,  const Args&... args) {
        int width = 70;
        int padding = (width + title.size()) / 2.0;
        std::cout << color << std::string(width, '-')<< RESET << std::endl;
        std::cout << color << std::setw(padding) << title << RESET << std::endl;

        std::ostringstream oss; 
        (void)std::initializer_list<int>{ (appendArg(oss, args), 0)... }; 

         padding = (width + oss.str().size()) / 2.0;
        std::cout  << std::setw(padding) << oss.str() << RESET << std::endl;
       
       
        std::cout << color << std::string(width, '-')<< RESET << std::endl; 
        std::cout  << std::endl; 
    }

    static int visualWidth(const std::string& s) {
        int count = 0;
        for (size_t i = 0; i < s.size(); ) {
            unsigned char c = s[i];
            if      (c < 0x80) i += 1;  // ASCII
            else if (c < 0xE0) i += 2;  // 2-byte sequence
            else if (c < 0xF0) i += 3;  // 3-byte sequence (your box chars)
            else               i += 4;  // 4-byte sequence
            count++;
        }
        return count;
    }
 

     static void SimQJStarter() {
        std::string args = "   ██████╗ ██╗███╗   ███╗       ██████╗      ██╗";
        int width = 70;
        int padding    = (width - visualWidth(args)) / 2.0;
        std::cout << BLUE << std::string(width, '=')<< RESET << std::endl;
        std::cout << BLUE << std::string(padding, ' ')  << "   ██████╗ ██╗███╗   ███╗       ██████╗      ██╗" << RESET << std::endl;
        std::cout << BLUE << std::string(padding, ' ')  << "  ██╔════╝ ██║████╗ ████║      ██╔═══██╗     ██║" << RESET << std::endl;
        std::cout << BLUE << std::string(padding, ' ')  << "  ███████╗ ██║██╔████╔██║      ██║   ██║     ██║" << RESET << std::endl;
        std::cout << BLUE << std::string(padding, ' ')  << "  ╚════██║ ██║██║╚██╔╝██║      ██║▄▄ ██║     ██║" << RESET << std::endl;
        std::cout << BLUE << std::string(padding, ' ')  << "  ███████║ ██║██║ ╚═╝ ██║      ╚██████╔╝ ██╗ ██║" << RESET << std::endl;
        std::cout << BLUE << std::string(padding, ' ')  << "  ╚══════╝ ╚═╝╚═╝     ╚═╝       ╚══██═╝  ╚████" << RESET << std::endl;
        
        args = std::string("Version    ") + VERSION_STRING;
        padding = (width + args.size()) / 2.0;
        std::cout << BLUE << std::setw(padding) << args << RESET << std::endl;


        std::cout << BLUE << std::string(width, '=')<< RESET << std::endl; 
        std::cout  << std::endl; 
 
         
    }








       
    static void spacer() {
        std::cout << std::endl;
    }
 

    static void loadingBar(double progress, int total, double time=-1, int barWidth = 70) {
        float ratio = progress / total;
        int filled = (int)(ratio * barWidth);

        if (progress > 0 ) {
            std::cout << "\033[3A"; // Move up 3 lines (Empty line, Bar, Empty line)
        }
 
        std::cout << "\033[K\n"; // Clear line and move down
        std::cout << "\r\033[K ["; // Reset to start of line, Clear it, then print
        for (int i = 0; i < barWidth; i++) {
            if (i < filled) std::cout << "=";
            else if (i == filled) std::cout << ">";
            else std::cout << " ";
        } 
        std::cout << "] " << (int)(ratio * 100) << "%,     t : "  << (time > 0 ? std::to_string(time) + " s" : "") << std::endl;
        std::cout << "\033[K\n" << std::flush;// Clear line and move down


    }

private:
    // C++17 version 

    // template<typename T>
    // static void print_item(const T& item) {
    //     if constexpr (std::is_same<T, sim_math::Matrix>::value) {
    //         item.print();
    //     } else {
    //         std::cout << item;
    //     }
    //     std::cout << " "; 
    // }

    // C++ 11 version


    static void print_item(const sim_math::Matrix& item) {
        item.print();
        std::cout << " ";
    }

    static void print_item(const lapack_complex_double item) { 
        std::cout << sim_math::creal(item)<<" + "<< sim_math::cimag(item) << "j    ";
    }


    template<typename T>
    static void print_item(const T& item) {
        std::cout << item << " ";
    }
};