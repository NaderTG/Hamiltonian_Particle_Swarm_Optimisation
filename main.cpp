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
 *  main.cpp
 *
 *  Created by Nader on 19/03/2017.
 */

#include <iostream>
#include <vector>
#include "Particle.h"
#include "Params.h"
#include "Swarm.h"
#include "system_func.h"
#include <gsl/gsl_rng.h>

void testSwarm(){
    Params pam("Rastrigin");
    const int dims = 2;
    Swarm<dims> parts(10);
    parts.initPopulation(pam, Rastrigin<dims>);
    // parts.printPopulation();
    parts.updateParticles(pam, Rastrigin<dims>);
    std::cout << "\n";
    parts.showGlobalBest(Rastrigin<dims>);
    
    
    
}



int main(int argc, const char * argv[]) {

    std::cout << "Finding the global minimum for Rastrigin function using Particle Swarm Optimisation!\n";
    testSwarm(); //Not finished

    return 0;
}
