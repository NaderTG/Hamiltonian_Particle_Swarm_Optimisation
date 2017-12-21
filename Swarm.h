/*
 * This program is an implementation of particle swarm optimisation algorithm for finding the global extremum of a function using stochastic Hamiltonian formulation.
 * Copyright (C) 2017 Nader Ganaba.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Swarm.h
 *
 *  Created by Nader on 19/03/2017.
 */

#ifndef Swarm_h
#define Swarm_h

#include "Particle.h"
#include "Params.h"
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <gsl/gsl_rng.h>

template<int N>
class Swarm{
private:
    int num_population_;
    std::vector<Particle<N>*> population_;
    double* p_best_;
    double f_g_best_; //The cost at the global best
public:
    
    Swarm(int size){
        
        p_best_ = (double*) calloc (N,sizeof(double));
        num_population_ = size;
        f_g_best_ = 100000.0;
        
    }
    template<class System>
    void showGlobalBest( System cost_func){
        std::cout << "Best position = [";
        for(int i = 0; i < N-1; i++){
            std::cout <<p_best_[i]  << ", ";
        }
        std::cout <<p_best_[N-1] << "]\n";
        double f_val = 0.0;
        cost_func(p_best_, f_val );
        
        std::cout << "The cost functional is then = " << f_val << std::endl;
        
    }
    
    int getPopulationSize(){return num_population_;}
    template<class System>
    void initPopulation(Params& params_val, System cost_func){
        //
        double* position;
        double* velocity;
        double f_val = 0.0;
        position = (double*) calloc (N,sizeof(double));
        velocity = (double*) calloc (N,sizeof(double));
        
        const gsl_rng_type * T;
        gsl_rng * r;
        
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        for(int i = 0; i < params_val.size_pop_; i++){
            
            for(int k = 0; k < params_val.dims_; k++){
                // position[k] = random use params_val.xmin_ and stuff
                // velocity[k] =
                // temp = (double) i / (double )(params_val.size_pop_ -1) ;
                position[k] = (params_val.xmax_ - params_val.xmin_)*(gsl_rng_uniform (r)) + params_val.xmin_;
                velocity[k] = (params_val.vmax_ - params_val.vmin_)*(gsl_rng_uniform (r)) + params_val.vmin_;
            }
            
            //cells_domain.push_back(new Cell(_position, idx, 0, 0));
            population_.push_back(new Particle<N>(i, position, velocity));
            population_[i]->initPositionBest();
            //Compute f_val
            cost_func(position, f_val );
            population_[i]->setCostVal( f_val );
            //                //set global best
            if(f_val <= f_g_best_){
                f_g_best_ = f_val;
                for(int k = 0; k <params_val.dims_ ; k++){
                    p_best_[k] =position[k];
                }
            }
            
            
        }
        
        gsl_rng_free (r);
    }
    
    void printPopulation(){
        for(int i = 0; i < num_population_; i++){
            population_.at(i)->printInfo();
        }
    }
    template<class System>
    void updateParticles(Params& params_val, System cost_func){
        double c1, c2;
        c1 = params_val.c1_;
        c2 = params_val.c2_;
        double f_val = 0.0;
        
        int counter = 0;
        
        //while(error > params_val._tolerance && counter < params_val.max_iter_){
        while(counter < params_val.max_iter_){
            for(int i = 0; i < num_population_; i++){
                
                population_[i]->updateVelocity(p_best_, params_val.inertia_, params_val.c1_, params_val.c2_);
                population_[i]->updatePosition();
                cost_func(population_[i]->getPosition(), f_val);
                population_[i]->setCostVal(f_val);
                population_[i]->updatePositionBest();
                
                //update the global best
                if( f_val< f_g_best_){
                    f_g_best_ = f_val;
                    population_[i]->copyPosition(p_best_);
                    
                }
            }
            //compute error
            counter++;
        }
    }
    
    //                void eval_cost( System cost){
    //                   //
    //                    double f_val;
    //                    for(int i = 0; i < num_population_; i++){
    //
    //                        f_val = cost(population_[i]->getPosition());
    //                        population_[i]->setCostVal(f_val);
    //                        population_[i]->updatePositionBest();
    //                    }
    //
    //                }
    
    
    friend std::ostream& operator<<(std::ostream& os,  Swarm& _sw){
        os << "Needs to be done";
    }
    
    
};




#endif /* Swarm_h */
