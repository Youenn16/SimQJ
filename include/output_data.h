#pragma once
#include <iostream> // string
#include <iomanip> // precision
#include <fstream> // ofstream
#include <vector> // vectors
#include "simulation.h"
#include "config.h"
#include "logger.h"




class Simulation; // forward declaration

class OutputData {

    // std::unique_ptr<std::vector<double>> energy_history;



    // std::vector<double>* ent_array_list;
    // std::vector<double>* info_qj_ent;
    // std::vector<double>* info_nh_ent;
    // std::vector<double>* info_qj;
    // std::vector<double>* delta_ent_qj;
    // std::vector<double>* delta_ent_nh;
    

    // std::unique_ptr<std::vector<double> > time_array;
    
    // std::vector<double> norm;
    // std::vector<double> norm;
    // std::vector<double> delta;
    // std::vector<double> obs;
    // std::vector<double> proba_jump;
    // std::vector<double> proba_jump_site;
    // std::vector<double> jump_time;
    // std::vector<double> jump_pb_value;
    // std::vector<int> jump_index;
    // std::vector<int> error_index;

    public:
       

        typedef void (*DataChannelFunc)(Simulation& simu, DataChannel &var);
        
        struct ActiveDataChannels {
            DataChannelFunc func; 
            DataChannel* var;
        }; 
 
         
        OutputData(std::string path, std::ofstream& file ): file_(file), path_(path) {}


        

        void setupDataChannels(Config& config);
        void collectData(Simulation& simu);
        void collectProbas(Simulation& simu, double r, double r_site);
        void collectEntanglementStatistics(Simulation& simu);
      
        void saveData(Simulation& simu); 
        void reinit(); 
        int spacing_ = 12; // For .dat files saving


        std::ofstream&  file() {
            return file_;
        }

        std::string path() {
            return path_;
        }

        void set_path(std::string new_path) {
            path_ = new_path;
        }
        

    private: 
        /**
         * @brief Writes a formatted field to a file with a fixed column width.
         *        Appends a comma separator unless it is the last field in the row.
         *
         * @tparam T      Type of the value to write (e.g. double, int, std::string)
         * @param file    The output file stream to write to
         * @param isLast  True if this is the last field in the row (no comma appended)
         * @param value   The value to write
         */
        template<typename T>
        void writeField(std::ofstream& file, bool isLast, const T& value) {
            file << std::setw(spacing_) << value;
            if (!isLast) file << ",";
        }


        void save_vector(std::ofstream& file_, std::string path,std::vector<double>& data, int prec=6);    
        void save_vector(std::ofstream& file_, std::string path,std::vector<int>& data);
        void save_csv(std::ofstream &file_1, std::string header, std::string path, std::vector<int> &output_indices, std::vector<ActiveDataChannels> &activeDataChannels);
        void save_ent_statistics(std::ofstream &file_1,   std::string path_qj, std::string path_nh, std::vector<int> &output_indices, std::vector<ActiveDataChannels> &activeDataChannels);

        std::ofstream& file_;
        std::string path_;

        std::vector<ActiveDataChannels> activeDataChannels;

        static void entanglementChannel(Simulation& simu, DataChannel &var ); 
        static void normChannel(Simulation& simu, DataChannel &var);
        static void timeChannel(Simulation& simu, DataChannel &var); 
        static void densityChannel(Simulation& simu, DataChannel &var);  
        static void entanglementStatisticsChannel(Simulation &simu, DataChannel &var);
        static void probasChannel(Simulation &simu, DataChannel &var);

       
        
};
 