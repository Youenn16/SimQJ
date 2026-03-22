#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <set>

#define NO_CLICK_KEY "no_click"
#define HIGH_ORDER_KEY "high_order"
#define SIZE_KEY "L"
#define TIME_STEP_KEY "dt"
#define TIME_MAX_KEY "time_max"
#define INIT_STATE_KEY "init_state"
#define NB_PART_KEY "nb_part"


#define RECORD_PROBA_KEY "record_probas"
#define RECORD_ENTANGLEMENT_STATISTIC_KEY "record_entantanglement_statistics"
#define ENTANGLEMENT_KEY "entanglement"
#define ENTANGLEMENT_STATISTIC_KEY "entanglement_statistics"
#define NORM_KEY "norm"
#define TIME_KEY "time"
#define DENSITY_KEY "density"
#define PROBA_KEY "proba"


template<typename T> struct Tag {};


struct DataChannel{
        std::string name;
        std::vector<std::vector<double>> data;
        int target_size;
        // std::vector<double>* targets[];
        int param;
};

struct Value {
    enum Type { DOUBLE, BOOL, STRING } type;
    double  d;
    bool    b;
    std::string s;
};

struct Config {
    
    std::string model_name_;
    std::map<std::string, Value> parameters; // store 
    std::vector<DataChannel> monitoring;

    bool checkParam(std::string param_name){
        return parameters.find(param_name) != parameters.end();
    }



    // Unified entry point so the call site stays clean
    template<typename T>
    T getParam(std::string param_name) {
        return getParam(param_name, Tag<T>{});
    }

    double getParam(std::string param_name, Tag<double>) {
        return parameters.find(param_name)->second.d;
    }

    bool getParam(std::string param_name, Tag<bool>) {
        return parameters.find(param_name)->second.b;
    }

    std::string getParam(std::string param_name, Tag<std::string>) {
        return parameters.find(param_name)->second.s;
    } 
     

    void setParam ( const std::string& key, double val) {
        parameters[key] = {Value::DOUBLE, val, false, ""};
    }
    void setParam( const std::string& key, bool val) {
        parameters[key] = {Value::BOOL, 0.0, val, ""};
    }
    void setParam( const std::string& key, std::string val) {
        parameters[key] = {Value::STRING, 0.0, false, val};
    }
 
};

// Helper to trim whitespace
inline std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}


// auto contains = [](std::initializer_list<std::string> list, const std::string& key) {
//     return std::find(list.begin(), list.end(), key) != list.end();
// };

inline DataChannel generateDataChannel(std::string name_, int param_, int target_size_, int RESERVE_MEM){
    DataChannel var;
    var.name = name_;
    var.param = param_;
    var.target_size = target_size_;
    var.data.resize(target_size_);
    for (size_t i = 0; i < target_size_; i++)
    {
        var.data[i].reserve(RESERVE_MEM); 
    } 

    return var;
}

inline Config parseConfig(const std::string& filename ) {
    Config config;
    std::ifstream file(filename);
    std::string line, currentSection;
    int RESERVE_MEM=5000; 
    
    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
 
        if (line.back() == ':') {
            currentSection = line.substr(0, line.size() - 1);
            continue;
        }

        if(config.checkParam("dt") && config.checkParam("time_max") ){
            RESERVE_MEM = 3 * config.getParam<double>("time_max")/config.getParam<double>("dt");
        }
        
        if (currentSection == "ModelParameters" || currentSection == "SimulationParameters") {

            size_t sep = line.find('=');
            if (sep != std::string::npos) {
                std::string key = trim(line.substr(0, sep));
                std::string valStr = trim(line.substr(sep + 1));
                std::string lower = valStr;
                std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower); 

                // When adding a config parameters, if it is of type bool or string, it must be added in the following
                std::set<std::string> boolKeys = {"no_click","high_order", "high_order","record_probas"};
                std::set<std::string> stringKeys = {"init_state"};
                
                if(key == "name"){
                    config.model_name_ = lower;
                }else if(stringKeys.count(key) ){
                    config.setParam( key,lower);
                }else if(boolKeys.count(key) ){
                    config.setParam(key,lower=="true"?true: false);
                }else{
                    config.setParam(key,std::stod(valStr));
                }
            } 
        } 
        else if (currentSection == "DataMonitoring") {
            size_t startBracket = line.find('(');
            size_t endBracket = line.find(')');
            std::string tmp_name =trim(line.substr(0, startBracket));

            int tmp_target_size = 1;
            if(tmp_name == ENTANGLEMENT_STATISTIC_KEY){ 
                config.setParam(RECORD_ENTANGLEMENT_STATISTIC_KEY, true);
                continue;

            }
           
            if (startBracket != std::string::npos && endBracket != std::string::npos) {// case where channels have a param

                std::string args = line.substr(startBracket + 1, endBracket - startBracket - 1);
                std::stringstream ss(args);
                std::string item;
                // Here you build the different channels where data will be saved, typically if you have entanglement(2,4,8),
                // in the config file, three channels will be built with always name="entanglement", param=2,  param=4 or param=8
                // and a target_size=1 meaning only one list per channel. When histogramms are built target_size=4 allowing many lists. 
                while (std::getline(ss, item, ',')) {
 
                    config.monitoring.push_back(generateDataChannel(tmp_name, std::stoi(trim(item)), tmp_target_size,RESERVE_MEM));
 
                }
            } else { // case where no param are passed
                
                config.monitoring.push_back(generateDataChannel(tmp_name, -1, tmp_target_size,RESERVE_MEM)); 
            }

           
        }
    }


    if(config.checkParam(RECORD_PROBA_KEY) && config.getParam<bool>(RECORD_PROBA_KEY)){
        config.monitoring.push_back(generateDataChannel("proba", -1, 2,RESERVE_MEM));  
    }

    // Handling the entanglement statistic case 
    if(config.checkParam(RECORD_ENTANGLEMENT_STATISTIC_KEY) ){

        for (auto &var : config.monitoring)
        {
            if(var.name == ENTANGLEMENT_KEY){
                config.monitoring.push_back(generateDataChannel(ENTANGLEMENT_STATISTIC_KEY, var.param, 6,RESERVE_MEM));
            }
        } 
        
    }


    return config;
}
