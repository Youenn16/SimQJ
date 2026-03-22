#include "sim_math.h"
#include "models/ising_model.h"
#include "logger.h"
#include "simulation.h"
#include "output_data.h"
#include "config.h"
#include "utils.h"
#include <iostream>
#include <filesystem>

#include <sys/stat.h>
#include <stdexcept>
#include <getopt.h>
#include <filesystem>


static const std::string srcDir = std::string(__FILE__).substr(0, std::string(__FILE__).find_last_of("/\\"));




#define VERSION 0.7
#define VERSION_ID 1000


static bool unique_run = true;
static bool verbose = false;
static int batch_id = -1;
static int loop_size = -1;
static std::string output_path = srcDir + "/../build/res/"; 
static std::string config_path = srcDir + "/../config/config.param"; 
 


void showHelp(const char* programName) {
    std::cout << "Usage: " << programName << " [OPTIONS]\n\n"
                          << "  --version              Version "<< VERSION << "\n" 
                          << "  -v, --verbose          Enable verbose output\n" 
                          << "  -o, --output <folder>  Output folder\n" 
                          << "  -b, --batch_id         Batch Id \n"
                          << "  -l, --loop_number      Number of iterations\n"
                          << "  -c, --config_path      Config path\n"
                          << "  -h, --help             Show this help\n"; 
}


int main(int argc, char* argv[]) {
    int opt;  
 

    static struct option long_options[] = {
        {"version", no_argument,       nullptr, VERSION_ID},  // --version
        {"verbose", no_argument,       nullptr, 'v'},  // --verbose  
        {"output",  required_argument, nullptr, 'o'},  // --output <folder>
        {"batch_id",   required_argument, nullptr, 'b'},  // --batch_id <n>
        {"loop_number",   required_argument, nullptr, 'l'},  // --loop_number <n>
        {"config_path",   required_argument, nullptr, 'c'},  // --config_path <file>
        {"help",    no_argument,       nullptr, 'h'},  // --help
        {nullptr,   0,                 nullptr,  0 }   // sentinel, always end with this
    };
    // The string "vo:" means:
    //   'v' = flag (no argument)
    //   'o' = option WITH argument (the ':' means it requires a value)
    while ((opt = getopt_long(argc, argv, "vfo:n:hb:l:c:", long_options, nullptr)) != -1) {
        switch (opt) {
            case VERSION_ID:
                std::cout << "Version : " << VERSION << std::endl;
                return 1;
            case 'v':
                verbose = true;
                break;
            case 'o':
                output_path = srcDir+ "/../" + std::string(optarg);   // optarg points to the option's value
                break;
            case 'b':
                batch_id = std::stoi(optarg); 
                break;
            case 'l':
                loop_size = std::stoi(optarg); 
                break;
            case 'c':
                config_path = srcDir+  "/../" + std::string(optarg); 
                break;
            case 'h': 
                showHelp(argv[0]);
                return 1;
                break;
            case '?':    // unknown option or missing argument
                std::cerr << "Usage: " << argv[0] << " [-v] [-o output]\n";
                return 1;
        }
    } 

    
    Logger::SimQJStarter();

    srand(51301); 
    Config config = parseConfig(config_path); // Configuration file
    
    try {
        utils::checkFolder(output_path);
    } catch (const std::runtime_error& e) {
        Logger::error(e.what());
        exit(EXIT_FAILURE);
    }
    std::ofstream file_data; 
    OutputData output(output_path, file_data); // Results saving
    output.setupDataChannels(config); // Preparing the relevant quantities to save

    // Model Mapping
    std::unique_ptr<ModelBase> model;  
    if (config.model_name_ == "isingchain_density"){
        model = utils::make_unique<IsingModel>(config);
    }else{
        Logger::error("Model not supported (or mistake in the model name config)" );
        exit(EXIT_FAILURE);
    } 


     if(loop_size < 0 || batch_id < 0){
        loop_size = 1;
        batch_id = -1;
    }
    
    // Simulation initialisation
    Simulation simu(output, config );
    simu.init();
    model->init(simu); // always initialized model after simulation



    simu.createLogs(batch_id);
   

    for (size_t k = 0; k < loop_size; k++)
    {  
        if(batch_id>-1){ 
            output.set_path(output_path+"traj_"+std::to_string(k)+ "_");
        }else{
            output.set_path(output_path+"traj_");
        }


        model->simulate(simu); 
        output.saveData(simu);
        simu.reinit();
        output.reinit();

    }
    simu.clean();



    
    return 0;
}