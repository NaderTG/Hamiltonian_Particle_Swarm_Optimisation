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
 *  Particle.h
 *
 *  Created by Nader on 19/03/2017.
 */

#ifndef Particle_h
#define Particle_h

#include <stdlib.h>
#include <iostream>
#include <gsl/gsl_rng.h>

template<int N>
class Particle
{
    double* position_;
    double* velocity_;
    double* p_best_; //Position of personal best
    double* g_best_; //Global best
    
    int part_id_; //population id
    double f_val_, f_val_old_; //value of functional at personal best
    int part_rank_; //Ranking among the population
    
    gsl_rng * r;
    
public:
    Particle(int id_val){
        position_ = (double*) calloc (N,sizeof(double));
        velocity_ = (double*) calloc (N,sizeof(double));
        p_best_    = (double*) calloc (N,sizeof(double));
        g_best_    = (double*) calloc (N,sizeof(double));
        f_val_    = 0.0;
        f_val_old_ = 0.0;
        part_id_   = id_val;
        
        gsl_rng_env_setup();
        r = gsl_rng_alloc (gsl_rng_default);
        gsl_rng_set(r, (time(NULL)*(id_val + 1)));
        
    }
    
    Particle(int id_val, double* pos){
        position_ = (double*) calloc (N,sizeof(double));
        velocity_ = (double*) calloc (N,sizeof(double));
        p_best_    = (double*) calloc (N,sizeof(double));
        g_best_    = (double*) calloc (N,sizeof(double));
        f_val_    = 0.0;
        f_val_old_ = 0.0;
        
        part_id_ = id_val;
        
        gsl_rng_env_setup();
        r = gsl_rng_alloc (gsl_rng_default);
        gsl_rng_set(r, (time(NULL)*(id_val + 1)));
        
        for(int i = 0; i < N; i++){
            position_[i] = pos[i];
        }
    }
    
    Particle(int id_val, double* pos, double* vel){
        position_ = (double*) calloc (N,sizeof(double));
        velocity_ = (double*) calloc (N,sizeof(double));
        p_best_    = (double*) calloc (N,sizeof(double));
        g_best_    = (double*) calloc (N,sizeof(double));
        f_val_    = 0.0;
        f_val_old_ = 0.0;
        
        part_id_ = id_val;
        
        gsl_rng_env_setup();
        r = gsl_rng_alloc (gsl_rng_default);
        gsl_rng_set(r, (time(NULL)*(id_val + 1)));
        
        for(int i = 0; i < N; i++){
            position_[i] = pos[i];
            velocity_[i] = vel[i];
        }
    }
    
    void printInfo(){
        std::cout <<"Particle [" << part_id_<< "] \t";
        std::cout << "Position [ ";
        for(int i = 0; i < N-1; i++){
            std::cout <<  position_[i] <<", ";
        }
        std::cout  << position_[N-1] << "] \t";
        
        std::cout << "Velocity [ ";
        for(int i = 0; i < N-1; i++){
            std::cout  <<  velocity_[i] <<", ";
        }
        std::cout  << velocity_[N-1] << "]\t";
        std::cout << "F_val = " << f_val_ << std::endl;
        
    }
    
    void testRnd(){
        double u = gsl_rng_uniform (r);
        double u2 = gsl_rng_uniform (r);
        std::cout << "r = " << u << ",\t" << u2 << std::endl;
    }
    
    friend std::ostream& operator<<(std::ostream& output,  Particle<N>& part){
        
        output << "Particle [" << part.getID() << "] \t";
        output << "Position [ ";
        for(int i = 0; i < N-1; i++){
            output << part.position_[i] <<", ";
        }
        output <<part.position_[N-1] << "] \t";
        
        output << "Velocity [ ";
        for(int i = 0; i < N-1; i++){
            output << part.velocity_[i] <<", ";
        }
        output <<part.velocity_[N-1] << "]\n";
        
        
        return output;
    }
    
    
    double operator()(int i) {return position_[i];}
    
    
    
    void updateGlobalBest(double *g_best_vec){
        //Put error thingy
        for(int i = 0 ; i < N; i++){
            g_best_[i] = g_best_vec[i];
        }
    }
    
    void initPositionBest(){
        
        for(int i = 0 ; i < N; i++){
            p_best_[i] = position_[i];
        }
        
    }
    
    double getPositionBest(int i) {return p_best_[i];}
    
    void updatePositionBest(){
        if( f_val_  < f_val_old_){
            for(int i = 0 ; i < N; i++){
                p_best_[i] = position_[i];
            }
        }
    }
    
    
    
    double* getPosition(){return position_;}
    double getPosition(int i) {return position_[i];}
    void copyPosition(double* a){
        for(int k = 0; k < N; k++){
            a[k] = position_[k];
        }
    }
    void setPosition(int i, double val){position_[i] = val;}
    void setPosition(double* val){
        for(int i = 0; i < N; i++){
            position_[i] = val[i];
        }
    }
    double* getVelocity(){return velocity_;}
    double getVelocity(int i) {return velocity_[i];}
    void setVelocity(int i, double val){velocity_[i] = val;}
    
    
    
    void updateVelocity(double* g_best, double weight, double c1, double c2){
        double r1 = 0.0, r2 = 0.0;
        
        
        for(int i  = 0; i < N; i++){
            r1 = c1*gsl_rng_uniform (r);
            r2 = c2*gsl_rng_uniform (r);
            velocity_[i] = weight*velocity_[i] + r1*(p_best_[i] -position_[i]) + r2*(g_best[i] - position_[i]);
            
        }
        
    }
    
    void updatePosition(){
        
        for(int i  = 0; i < N; i++){
            position_[i] += velocity_[i];
        }
        
    }
    
    void updateRank(int rank){
        part_rank_ = rank;
    }
    void setCostVal(double val) {f_val_old_ = f_val_; f_val_ = val;}
    double getCostVal(){
        return f_val_;
    }
    
    int getID(){return part_id_;}
    
};




#endif /* Particle_h */
