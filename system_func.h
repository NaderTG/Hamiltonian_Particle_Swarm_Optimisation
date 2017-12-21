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
 *  system_func.h
 *  Hamiltonian_PSO
 *
 *  Created by Nader on 19/03/2017.
 */

#ifndef system_func_h
#define system_func_h

#define _USE_MATH_DEFINES // for C++
#include <cmath>

#define POW2(x) (x*x)
template<int N>
void Rastrigin( const double* x , double &y )
{
    double temp1 = 0.0, temp2 = 0.0;
    double A = 10.0;
    for(int i = 0; i < N; i++){
        temp1 += x[i]*x[i];
        temp2 += A*cos(2.0*M_PI*x[i]);
    }

    y = A*N +temp1 - temp2;
    
}

template<int N>
void Rosenbrock( const double* x , double &y )
{
    double temp1 = 0.0, temp2 = 0.0;
    
    for(int i = 0; i < N-1; i++){
        temp1 += (x[i]-1.0)*(x[i]-1.0);
        temp2 += 100.0*(x[i+1] - x[i]*x[i])*(x[i+1] - x[i]*x[i]) ;
    }
    
    y =  temp1 + temp2;
}

template<int N>
void Goldstein( const double* x , double &y ){
    double temp1 = 0.0, temp2 = 0.0;
    
    temp1 = 1 + (x[0] + x[1] + 1)*(x[0] + x[1] + 1)*(19 - 14*x[0] + 3*x[0]*x[0] - 14*x[1] + 6*x[0]*x[1] + 3*x[1]*x[1]);
    temp2 = 30 + (2*x[0] - 3*x[1])*(2*x[0] - 3*x[1])*(18 - 32*x[0] + 12*x[0]*x[0] + 48*x[1] - 36*x[0]*x[1] + 27*x[1]*x[1]);
    y =  temp1*temp2;
}

template<int N>
void calcState(const double* u , double *x , double *y){
    //The dynamical constraint is
    // \dot{x_1} = x_2      x_1(0) = 0
    // \dot{x_2} = u        x_2(0) = 0
    
    //The norm \| [x_1^{n+1} - x_1^n - \Delta t x^n_2 ,x^{n+1}_2 - x^n_2  - u^n ] \|^2 is given explicitly by
    //
    double dt = 1.5/15.0;
    x[0] = M_PI; y[0] = 1;
    for(int i = 1; i < N; i++){
        
        x[i] = x[i-1] +dt*y[i-1];
        y[i] = y[i-1] +dt*u[i] - dt*x[i] + dt*y[i-1]*(1.0 - x[i]*x[i]);
    }
    
    
}

template<int N>
void controlCost( const double* u , double &y ){
    
    double* x1 = (double*) calloc (N,sizeof(double));
    double* x2 =  (double*) calloc (N,sizeof(double));
    
    calcState<N>(u, x1, x2);
    
    double sum = 0.0;
    double dt = 1.5/15.0;
    // sum = 0.5*u[0]*u[0];
    double temp1, temp2;
    for(int i = 1; i < N-1; i++){
        temp1 = x1[i] - x1[i-1] -dt*x2[i-1];
        temp2 = x2[i] - x2[i-1] -dt*u[i]  + dt*x1[i] - dt*x2[i-1]*(1.0 - x1[i]*x1[i]);
        sum += dt*0.5*u[i]*u[i] + 0.5*dt*temp1*temp1 + 0.5*dt*temp2*temp2;
        
    }
    sum += 0.5*x1[N]*x1[N];
    delete(x1);
    delete(x2);
    
    y = sum;
}
#endif /* system_func_h */
