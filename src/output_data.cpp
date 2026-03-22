#pragma once
#include "output_data.h"
#include "simulation.h"

void OutputData::entanglementChannel(Simulation &simu, DataChannel &var)
{
    // If entanglement statistic is recorded, entanglement is saved 
    var.data[0].push_back(simu.get_entanglement(var.param));
}
 

void OutputData::normChannel(Simulation &simu, DataChannel &var)
{
    var.data[0].push_back(simu.norm());
}

void OutputData::densityChannel(Simulation &simu, DataChannel &var)
{
    var.data[0].push_back(simu.get_density(var.param));
}

void OutputData::timeChannel(Simulation &simu, DataChannel &var)
{
    var.data[0].push_back(simu.time());
}

void OutputData::entanglementStatisticsChannel(Simulation &simu, DataChannel &var)
{ 
}

void OutputData::probasChannel(Simulation &simu, DataChannel &var)
{ 
}
 


void OutputData::setupDataChannels(Config &config)
{
    std::map<std::string, DataChannelFunc> registry;
    registry[ENTANGLEMENT_KEY] = entanglementChannel; 
    registry[NORM_KEY] = normChannel;
    registry[TIME_KEY] = timeChannel;
    registry[DENSITY_KEY] = densityChannel; 
    registry[ENTANGLEMENT_STATISTIC_KEY] = entanglementStatisticsChannel;
    registry[PROBA_KEY] = probasChannel; 

    for (auto &var : config.monitoring)
    {
        auto it = registry.find(var.name); // search channel saving action corresponding to monitoring requirement
        if (it != registry.end())
        {
            ActiveDataChannels af; // create the channel composed of the function to save, and the saving list + params
            af.func = it->second;
            af.var = &var;
            activeDataChannels.push_back(af); // list of all the channels
        }
        else
        {
            Logger::error("Unknown monitoring variable: ",var.name ) ;
        }
    }
    Logger::info("[Outputs] Channels configured ! ");
}


void OutputData::collectData(Simulation &simu)
{
    for (size_t i = 0; i < activeDataChannels.size(); ++i)
    {
        if (activeDataChannels[i].var->name != PROBA_KEY && activeDataChannels[i].var->name != ENTANGLEMENT_STATISTIC_KEY )
        {
            activeDataChannels[i].func(simu, *activeDataChannels[i].var);
        }
    }
}

void OutputData::collectProbas(Simulation &simu,double r, double r_site)
{
    if(simu.record_probas_){
        for (size_t i = 0; i < activeDataChannels.size(); ++i)
        {
            if (activeDataChannels[i].var->name == PROBA_KEY)
            {
                activeDataChannels[i].var->data[0].push_back(r);
                activeDataChannels[i].var->data[1].push_back(r_site); 
            }
        }
    }

}


void OutputData::collectEntanglementStatistics(Simulation &simu)
{
    if(simu.record_entanglement_statistics_){
        for (size_t i = 0; i < activeDataChannels.size(); ++i)
            {
                if (activeDataChannels[i].var->name == ENTANGLEMENT_STATISTIC_KEY)
                {
                    int l = activeDataChannels[i].var->param;
                    if(simu.jump_case() && simu.before_jump()){
                        activeDataChannels[i].var->data[0].push_back(simu.get_current_entanglement(l,0)-simu.get_current_entanglement(l,1)); // delta_ent_nh
                        activeDataChannels[i].var->data[1].push_back(simu.get_current_entanglement(l,1)); // ent_nh info
                    }else if(simu.jump_case() && !simu.before_jump()){
                        activeDataChannels[i].var->data[2].push_back(simu.get_current_entanglement(l,1)-simu.get_current_entanglement(l,0)); // delta_ent_qj
                        activeDataChannels[i].var->data[3].push_back(simu.get_current_entanglement(l,0)); // ent_qj info
                        activeDataChannels[i].var->data[4].push_back(simu.get_index_jump()); // ent_qj index
                        activeDataChannels[i].var->data[5].push_back(simu.time()); // ent_qj index
                    }
                     
                }
            }
    }

}

void OutputData::saveData(Simulation &simu)
{
    std::string filename = "";

    std::vector<int> entanglement_indices;
    std::vector<int> entanglement_statistics_indices;
    std::vector<int> density_indices;
    
    for (int i = 0; i < activeDataChannels.size(); i++)
    {
        if (activeDataChannels[i].var->name == ENTANGLEMENT_KEY)
        {
            entanglement_indices.push_back(i);
        }
        else if (activeDataChannels[i].var->name == ENTANGLEMENT_STATISTIC_KEY)
        {
            entanglement_statistics_indices.push_back(i);
        }
        else if (activeDataChannels[i].var->name == DENSITY_KEY)
        {
            density_indices.push_back(i);
        }else if(activeDataChannels[i].var->name == PROBA_KEY){
            filename = path() + "proba_jump.dat";
            save_vector(file(), filename, activeDataChannels[i].var->data[0]);
            filename = path() + "proba_jump_site.dat";
            save_vector(file(), filename, activeDataChannels[i].var->data[1]);
 
        }
        else // 1 channel case
        {
            filename = path() +   activeDataChannels[i].var->name + ".dat";
            save_vector(file(), filename, activeDataChannels[i].var->data[0]);
             
        }
    }

    if (!entanglement_indices.empty())
    {
        filename = path() + "entanglement.dat";
        save_csv(file(), "S_l=",filename, entanglement_indices, activeDataChannels);
    }

    if (!entanglement_statistics_indices.empty())
    {
        save_ent_statistics(file(), path() + "entanglement_stat_qj.dat",path() + "entanglement_stat_nh.dat", entanglement_statistics_indices, activeDataChannels);
    }

     if (!density_indices.empty())
    {
        filename = path() + "densities.dat";
        save_csv(file(), "n_i=",filename, density_indices, activeDataChannels);
    }
 

    Logger::info("[Outputs] Saved ! ");
    Logger::spacer();
}

void OutputData::reinit()
{
    for (size_t i = 0; i < activeDataChannels.size(); ++i)
    {
        for (auto &inner : (activeDataChannels[i].var)->data)
        {
            inner.clear();
        }
    }
}

void OutputData::save_vector(std::ofstream &file, std::string path, std::vector<double> &data, int prec)
{


    file.open(path);
    if (file.is_open())
    {
        for (double d : data)
        {
            file << std::fixed << std::setprecision(prec) << d << '\n';
        }
        file.close();
    }
}

void OutputData::save_vector(std::ofstream &file, std::string path, std::vector<int> &data)
{
    file.open(path);
    if (file.is_open())
    {
        for (int d : data)
        {
            file << d << '\n';
        }
        file.close();
    }
}

 



// Save data under the format 
// header param(output indice 1) header param(output indice 2) ...  header param(output indice n)
void OutputData::save_csv(std::ofstream &file_1, std::string header, std::string path, std::vector<int> &output_indices, std::vector<ActiveDataChannels> &activeDataChannels)
{

    file_1.open(path);
    if (file_1.is_open())
    {

        file_1 << std::fixed << std::setprecision(6);

        
        int index_ent = output_indices[0];
        int nb_point = (activeDataChannels[index_ent].var->data[0]).size();

        // HEADER WRITING
        for (size_t i = 0; i < output_indices.size(); ++i)
        {
            index_ent = output_indices[i];
            std::string value = header + std::to_string(activeDataChannels[index_ent].var->param);
            writeField<std::string>(file_1, i == output_indices.size() - 1, value);
        }

        file_1 << '\n';

        // DATA WRITING
        for (size_t i = 0; i < nb_point; ++i)
        {

            for (size_t k = 0; k < output_indices.size(); ++k)
            {
                index_ent = output_indices[k];
                std::string value = std::to_string( (activeDataChannels[index_ent].var->data)[0][i]);
                writeField<std::string>(file_1, k == output_indices.size() - 1, value); 
            }

            file_1 << '\n';
        }

        file_1.close();
    }
}



void OutputData::save_ent_statistics(std::ofstream &file_1,   std::string path_qj, std::string path_nh, std::vector<int> &output_indices, std::vector<ActiveDataChannels> &activeDataChannels)
{
    std::string value = "";
    // DELTA ENT QJ
    file_1.open(path_qj);
    if (file_1.is_open())
    {

        file_1 << std::fixed << std::setprecision(6);

        int index_ent = output_indices[0];
        int nb_point = (activeDataChannels[index_ent].var->data[2]).size(); //delta ent_qj size

        // HEADER WRITING
        file_1 << std::setw(spacing_) << "t";
        file_1 << std::setw(spacing_) << "site";
        for (size_t i = 0; i < output_indices.size(); ++i)
        {
            index_ent = output_indices[i];
            value = "S_l=" + std::to_string(activeDataChannels[index_ent].var->param);
            writeField<std::string>(file_1, false, value);
        }
        for (size_t i = 0; i < output_indices.size(); ++i)
        {
            index_ent = output_indices[i];
            value = "DS_l=" + std::to_string(activeDataChannels[index_ent].var->param);
            writeField<std::string>(file_1, i == output_indices.size() - 1, value);
        }


        file_1 << '\n';
        // DATA WRITING
        for (size_t i = 0; i < nb_point; ++i)
        {
            index_ent = output_indices[0];
            value = std::to_string( (activeDataChannels[index_ent].var->data)[5][i]);
            writeField<std::string>(file_1, false, value); 
            value = std::to_string( (activeDataChannels[index_ent].var->data)[4][i]) ;
            writeField<std::string>(file_1, false, value); 

            for (size_t k = 0; k < output_indices.size(); ++k)
            {
                index_ent = output_indices[k];
                value = std::to_string( (activeDataChannels[index_ent].var->data)[3][i]); // ent before qj
                writeField<std::string>(file_1, false, value);  
            }

            for (size_t k = 0; k < output_indices.size(); ++k)
            {
                index_ent = output_indices[k];
                value = std::to_string( (activeDataChannels[index_ent].var->data)[2][i]); // delta_ent_qj
                writeField<std::string>(file_1, k == output_indices.size() - 1, value); 
            }

            file_1 << '\n';
        }

        file_1.close();
    }


    // DELTA ENT NH
    file_1.open(path_nh);
    if (file_1.is_open())
    {

        file_1 << std::fixed << std::setprecision(6);
 
        int index_ent = output_indices[0];
        int nb_point = (activeDataChannels[index_ent].var->data[2]).size(); //delta ent_qj size

        // HEADER WRITING
        for (size_t i = 0; i < output_indices.size(); ++i)
        {
            index_ent = output_indices[i];
            value = "S_l=" + std::to_string(activeDataChannels[index_ent].var->param);
            writeField<std::string>(file_1, false, value);
        }
        for (size_t i = 0; i < output_indices.size(); ++i)
        {
            index_ent = output_indices[i];
            value = "DS_l=" + std::to_string(activeDataChannels[index_ent].var->param);
            writeField<std::string>(file_1, i == output_indices.size() - 1, value);
        }


        file_1 << '\n';
        // DATA WRITING
        for (size_t i = 0; i < nb_point; ++i)
        { 

            for (size_t k = 0; k < output_indices.size(); ++k)
            {
                index_ent = output_indices[k];
                value = std::to_string( (activeDataChannels[index_ent].var->data)[1][i]); // ent before qj
                writeField<std::string>(file_1, false, value);  
            }

            for (size_t k = 0; k < output_indices.size(); ++k)
            {
                index_ent = output_indices[k];
                value = std::to_string( (activeDataChannels[index_ent].var->data)[0][i]); // delta_ent_qj
                writeField<std::string>(file_1, k == output_indices.size() - 1, value); 
            }

            file_1 << '\n';
        }

        file_1.close();
    }
}