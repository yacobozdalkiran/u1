#ifndef FIELD_H
#define FIELD_H

#include <random>
#include <vector>

#include "../geometry/geometry.h"

class GaugeField{
    private:
    int L;
    int V;
    std::vector<double> links; //On stocke uniquement la phase
    
    public:
        GaugeField(const Geometry &geo);
        //Start
        void hot_start(std::mt19937_64& rng);
        void cold_start(); 
        //Access
        double get_link(int site, int mu) const;
        void set_link(int site, int mu, double new_value);
        void add_to_link(int site, int mu, double value);
        //Plaquette/Staples
        double plaquette(const Geometry& geo, int site) const;
};

#endif
