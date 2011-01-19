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

#ifndef SPINADAPT_HPP_
#define SPINADAPT_HPP_

#include"vector"

namespace TensorCalculus{

// v_mnop = <mn||op> = <mn|op> - <mn|po>
void readAOSpinAdapted (std::vector<double> &fullAdaptAO);
void spinAdaptedAO (const int aoCount, const std::vector<double> &fullAO, std::vector<double> &fullAdaptAO);
void readAOSpinAdaptedchemical (std::vector<double> &fullAdaptAO);
void spinAdaptedAOchemical (const int aoCount, const std::vector<double> &fullAO, std::vector<double> &fullAdaptAO);

}; //end namespace TensorCalculus

#endif /* SPINADAPT_HPP_ */
