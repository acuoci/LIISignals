/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2018 Alberto Cuoci                                       |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_RungeKutta4thOrder_H
#define OpenSMOKE_RungeKutta4thOrder_H

#include <vector>

namespace OpenSMOKE
{
	//!  A class for solving ODE systems beased on the 4th order Runge-Kutta method
	/*!
	This class provides the tools for solving ODE systems beased on the 4th order Runge-Kutta method
	*/

	class RungeKutta4thOrder
	{

	public:

		/**
		*@brief Default constructor
		*/
		RungeKutta4thOrder();

		/**
		*@brief Sets the initial conditions (together with memory allocation)
		*@param	n number of equations
		*@param	tInitial initial time
		*@param	uInitial initial values of unknowns
		*/
		void SetInitialConditions(const unsigned int n, const double tInitial, const double* uInitial);

		/**
		*@brief Sets the constant time step
		*@param	dt time step
		*/
		void SetTimeStep(const double dt);

		/**
		*@brief Sets the final time
		*@param	tFinal the final time
		*/
		void SetFinalTime(const double tFinal);

		/**
		*@brief Sets the ODE function
		*@param	*f pointer to the ODE function
		*/
		void SetOdeSystem(void(*f)(const double t, const double* u, double* dudt))
		{
			f_ = f;
		}

		/**
		*@brief Solves the ODE system
		*/
		void Solve();

		/**
		*@brief		Returns the solution
		*@return	the solution matrix (rows=number of steps, colums=time|unknown1|unknown2|...)
		*/
		const std::vector< std::vector<double> >& solution() const { return solution_;  }

	private:

		/**
		*@brief Allocates memory
		*/
		void MemoryAllocation();

		/**
		*@brief Pointer to the ODE function 
		*/
		void(*f_)(const double t, const double* u, double* dudt);


	private:

		std::vector<double> uInitial_;	//!< initial values of unknowns

		std::vector<double> u0_;		//!< vector for internal use
		std::vector<double> u1_;		//!< vector for internal use
		std::vector<double> u2_;		//!< vector for internal use
		std::vector<double> u3_;		//!< vector for internal use

		std::vector<double> f0_;		//!< vector for internal use
		std::vector<double> f1_;		//!< vector for internal use
		std::vector<double> f2_;		//!< vector for internal use
		std::vector<double> f3_;		//!< vector for internal use

		std::vector< std::vector<double> > solution_;	//!< the solution matrix(rows = number of steps, colums = time | unknown1 | unknown2 | ...)


		unsigned int n_;		//!< number of equations/unknowns
		double tInitial_;		//!< initial time
		double tFinal_;			//!< final time
		double dt_;				//!< time step

	};

}

#include "RungeKutta4thOrder.hpp"

#endif /* OpenSMOKE_RungeKutta4thOrder_H */
