/*
 * Copyright (C) 2011 Henry Auer
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

#ifndef __FULLTENSOR_HPP
#define __FULLTENSOR_HPP

#include <iostream>
#include "Representation/TensorRepresentation.hpp"
#include "Matrix/MatrixOperators.hpp"
#include "Utilities/Index.hpp"

namespace TensorCalculus{
  
 template <typename T> class FullTensor{
  
  protected:
   std::vector <T> v;
   std::vector <int> componentDimensions;
   int  d;
   long vSize;
  
  private:
   void init (const std::vector<int> &componentDimensions, const std::vector<T> &v){
      (*this).v = v;
      (*this).componentDimensions = componentDimensions;
      (*this).d = componentDimensions.size();
      (*this).vSize = v.size();
   }
   void init (const char *filename){
      (*this).read(filename);
   }
   void init(){
      d = 0;
      vSize = 0;
      v.resize(0);
      componentDimensions.resize(0);
   }
  
  public:
   FullTensor (const char *filename){
    init(&filename[0]);
   }
   FullTensor (const std::vector<int> &componentDimensions, const std::vector<T> &v){
       init(componentDimensions, v);
   }

   FullTensor (){
      init();
   }
   
  std::vector <T> getV () const {
   return v;
  }
  std::vector <int> getComponentDimensions () const {
   return componentDimensions;
  }

  int getComponentDimension (int mu) const {
     return componentDimensions[mu];
  }

  T getVofIndex (int mu) const {
     return v[mu];
  }

  void read(const char *filename){
     std::ifstream from(&filename[0]);
     if(from.fail()){
        throw std::invalid_argument("no input file ");
     }
      char t[64];

      from.setf(std::ios::scientific, std::ios::floatfield);
      from >> t >> t >> t >> t;
      from >> d;
      componentDimensions.resize(d);
      vSize = 1;

      for (int mu=0; mu<d; mu++){
       from >> t;
       from >> t;
       from >> componentDimensions[mu];
       vSize *= componentDimensions[mu];
      }
      v.resize(vSize);

      for (int i=0; i<vSize; i++){
       for(int mu=0; mu<d; mu++){ //index in Datei �berspringen
        from >> t;
       }
       from >> v[i];
      }

      std::cout<<"read file: "<<filename<<std::endl;
      std::cout<<"d = "<<d<<std::endl;
      std::cout<<"componentDimensions ="<<std::endl;
      using namespace VectorOperators;
      std::cout<<componentDimensions<<std::endl;
      //std::cout<<v <<std::endl;
  }

  void write2disk(const char *filename) {
     std::ofstream fout(&filename[0]);

     fout.setf(std::ios::scientific, std::ios::floatfield);
     fout<<"FullTensor in d = "<< d <<std::endl;
     int prodOfComponent = 1;
     //std::vector<int> partialProdOfComponent(d+1);
     //partialProdOfComponent[0] = 1;
     for(int mu=0; mu<d; mu++){
        fout<<"n["<<mu<<"] = "<<componentDimensions[mu]<<std::endl;
        //partialProdOfComponent[mu+1] = partialProdOfComponent[mu] * componentDimension[mu];
        prodOfComponent *= componentDimensions[mu];
     }

     Index index(componentDimensions);
     index.begin();

     for(int i=0; i<prodOfComponent; i++){
        using namespace VectorOperators;
        fout<<index.getCurrent()<<"  "<<v[i]<<std::endl;
        ++index;
     }
   } //end write2disk
  void setVofIndex (const std::vector<int> &index, T value){
     //todo: abfragen, ob index passt
     int position   = 0;
     int stepLenght = 1;
     for(int mu=0; mu<d; mu++){
        position += index[mu] * stepLenght;
        stepLenght *= componentDimensions[mu];
     }
     v[position] = value;
  }


 }; // end class FullTensor
} //end namespace TensorCalculus

#endif /* __FULLTENSOR_HPP */
