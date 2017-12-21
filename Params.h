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
 *  Params.h
 *
 *  Created by Nader on 19/03/2017.
 */

#ifndef Params_h
#define Params_h

#include <cmath>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream> // ifstream
#include <sstream> // stringstream


namespace params_config {
    
    struct data: std::map <std::string, std::string>
    {
        
        bool iskey( const std::string& s ) const
        {
            return count( s ) != 0;
        }
    };
    //---------------------------------------------------------------------------
    // The extraction operator reads configuration::data until EOF.
    // Invalid data is ignored.
    //
    std::istream& operator >> ( std::istream& ins, data& _data )
    {
        std::string s, key, value;
        
        // For each (key, value) pair in the file
        while (std::getline( ins, s ))
        {
            std::string::size_type begin = s.find_first_not_of( " \f\t\v" );
            
            // Skip blank lines
            if (begin == std::string::npos) continue;
            
            // Skip commentary
            if (std::string( "#;" ).find( s[ begin ] ) != std::string::npos) continue;
            
            // Extract the key value
            std::string::size_type end = s.find( '=', begin );
            key = s.substr( begin, end - begin );
            
            // (No leading or trailing whitespace allowed)
            key.erase( key.find_last_not_of( " \f\t\v" ) + 1 );
            
            // No blank keys allowed
            if (key.empty()) continue;
            
            // Extract the value (no leading or trailing whitespace allowed)
            begin = s.find_first_not_of( " \f\n\r\t\v", end + 1 );
            end   = s.find_last_not_of(  " \f\n\r\t\v" ) + 1;
            
            value = s.substr( begin, end - begin );
            
            // Insert the properly extracted (key, value) pair into the map
            _data[ key ] = value;
        }
        
        return ins;
    }
    
    
    std::ostream& operator << ( std::ostream& outs, const data& _data )
    {
        data::const_iterator iter;
        for (iter = _data.begin(); iter != _data.end(); iter++)
            outs << iter->first << " = " << iter->second << std::endl;
        return outs;
    }
    
    
    
}

class Params{
public:
    double xmin_, xmax_, vmin_, vmax_;
    double c1_, c2_;
    int size_pop_;  //Size of population
    int dims_;      //Dimensions
    // int _epochs;
    int max_iter_;
    double inertia_;
    double tolerance_;
    
    std::string filename;
    Params(std::string _filename);
    
    friend std::ostream& operator<<(std::ostream& os, const Params& pm);
    
};

Params::Params(std::string _filename){
    
    
    
    params_config::data myconfig;
    params_config::data::const_iterator iter;
    filename = _filename;
    std::string filename_1 = filename + ".ini";
    std::ifstream f(filename_1.c_str());
    
    f >> myconfig;
    f.close();
    
    //xmin
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("x min" == iter->first){
            xmin_  = std::stod (iter->second);
        }
    }
    //xmax
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("x max" == iter->first){
            xmax_  = std::stod (iter->second);
        }
    }
    
    //vmin
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("v min" == iter->first){
            vmin_  = std::stod (iter->second);
        }
    }
    
    //v max
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("v max" == iter->first){
            vmax_ = std::stod (iter->second);
        }
    }
    
    //c1
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("c1" == iter->first){
            c1_ = std::stod (iter->second);
        }
    }
    
    //c2
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("c2" == iter->first){
            c2_ = std::stod (iter->second);
        }
    }
    
    //tolerance
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("tolerance" == iter->first){
            tolerance_ = std::stod (iter->second);
        }
    }
    
    //inertia
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("inertia" == iter->first){
            inertia_ = std::stod (iter->second);
        }
    }
    
    //size pop
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("particle size" == iter->first){
            size_pop_  = std::stoi (iter->second);
        }
    }
    
    //Max Epoch
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("max epoch" == iter->first){
            max_iter_  = std::stoi (iter->second);
        }
    }
    
    //Dimensions
    for (iter = myconfig.begin(); iter != myconfig.end(); iter++){
        if("dimensions" == iter->first){
            dims_  = std::stoi (iter->second);
        }
    }
    
}

std::ostream& operator<<(std::ostream& os, const Params& pm)
{
    
    os << "Parameters:\n";
    os << "Position limits = [" <<pm.xmin_ << ", " << pm.xmax_ << "]\n";
    os << "Velocity limits = [" << pm.vmin_ << ", " << pm.vmax_ << "]\n";
    os << "Population size = " << pm.dims_ << std::endl;
    os << "Population dimension = " << pm.size_pop_ << std::endl;
    os << "Iteration params = [ " << pm.max_iter_ << ", " << pm.inertia_ << ", " << pm.tolerance_ << " ]\n";
    return os;
}


#endif /* Params_h */
