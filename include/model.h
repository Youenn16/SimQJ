#pragma once
#include <string>
#include <iostream>
#include "simulation.h"

class ModelBase {
public:
    virtual ~ModelBase() = default;  
    virtual void simulate(Simulation& simulation) = 0;     
    virtual void init(Simulation& simulation) = 0;     
};
 
template<typename Derived>
class Model : public ModelBase {
public:
    explicit Model(std::string name, std::string jump_name, bool particle_conservation)
        : name_(std::move(name)), jump_name_(std::move(jump_name)), particle_conservation_(particle_conservation) 
         {}

    // // Non-virtual public API
    // void simulate(data::OutputData& output, double time, double dt, int seed) const {
    //     return derived().simulate(output, time, dt, seed);
    // }

    void simulate(Simulation& simulation) override {
        return derived().simulate(simulation);
    }

    void init(Simulation& simulation) {
        return derived().init(simulation);
    }

    void createLogs(){ 
        Logger::info(" Model Creation Logs"); 
    }

    // void train(int epochs) {
    //     epochs_ = epochs;
    //     trained_ = true;
    //     derived().train_impl();
    // }

    const std::string& name() const {
        return name_;
    }

protected:
    // Access to derived safely for const objects
    const Derived& derived() const {
        return static_cast<const Derived&>(*this);
    }
    // Access to derived safely for mutable objects
    Derived& derived() {
        return static_cast<Derived&>(*this);
    }

protected:
    std::string name_;
    std::string jump_name_; 
    bool particle_conservation_;

};

// name 
// jump type 
// particle conservation 
// simulation method 
// first order / high order 
