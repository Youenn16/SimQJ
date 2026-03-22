#pragma once
#include "../model.h"
#include "simulation.h"
#include "sim_math.h"
#include "config.h"


class IsingModel : public Model<IsingModel> {

public:
    static constexpr double DEFAULT_J = 1.0;
    IsingModel(Config& config);
    IsingModel(double gamma,  double h, double kappa, double J = DEFAULT_J );



    // IsingModel create(Config& config);

    // void simulate(data::OutputData& output, double time, double dt, int seed) const;
    // void train_impl();
    void simulate(Simulation& simu) ;

    void init(Simulation& simu) ;

    void createLogs() ;

private:
    double gamma_;
    double J_;
    double h_;
    double kappa_;
 
    std::unique_ptr<sim_math::Matrix> H;

    



};

