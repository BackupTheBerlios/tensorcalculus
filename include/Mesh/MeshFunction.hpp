/*
 * Copyright (C) 2010 Philipp Wähnert
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __MESHFUNCTION_HPP
#define __MESHFUNCTION_HPP

#include <istream>
#include <sstream>
#include <vector>

#include "Utilities/Utilities.hpp"

namespace TensorCalculus {

  template<typename T>
  class MeshFunction {
  private:
    int num_functions;
    std::vector< std::vector<T> > data;
    
    void read_data_from_stream(std::istream& stream);
    
  public:
    MeshFunction(std::istream& stream) {
      read_data_from_stream(stream);
    }
    
    const std::vector<T>& getValues(int vertex) const {
      return data[vertex];
    }    
    
    const int countValues() const { return data.size(); }
    const int countFunctions() const { return num_functions; }
  };
  
  template<typename T>
  void MeshFunction<T>::read_data_from_stream(std::istream& stream) {
    std::string line;
    int num_vertices;
    
    if (!getline(stream, line)) throw StreamFormatError();
    {
      std::stringstream sstream(line);
      if (!(sstream >> num_vertices)) throw StreamFormatError();
      if (!(sstream >> num_functions)) throw StreamFormatError();
      data.reserve(num_vertices);
    }
  
    for (int vertex = 0; vertex < num_vertices; ++vertex) {
      std::vector<T> new_data;
      new_data.reserve(num_functions);
      
      if (!getline(stream, line)) throw StreamFormatError();
      std::stringstream sstream(line);
      for (int i = 0; i < num_functions; ++i) {
        T new_value;
        if (!(sstream >> new_value)) throw StreamFormatError();
        new_data.push_back(new_value);
      }

      data.push_back(new_data);
    }
  }

} // namespace TensorCalculus

#endif // __MESHFUNCTION_HPP
