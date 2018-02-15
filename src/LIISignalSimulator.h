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

#ifndef OpenSMOKE_LIISignalSimulator_H
#define OpenSMOKE_LIISignalSimulator_H

#include "LIISignalModel.h"

namespace OpenSMOKE
{
	//!  A class for simulating the temporal evolution of the LII Signal
	/*!
	This class provides the tools for simulating the temporal evolution of the LII Signal
	*/

	class LIISignalSimulator
	{
	
	public:

		/**
		*@brief Default constructor
		*@param		lii reference to a LIISignalModel object
		*/
		LIISignalSimulator(LIISignalModel& lii);

		/**
		*@brief Simulates the temporal evolution of temperature and mass of soot particles
		*/
		void Solve();

		/**
		*@brief Differential equations to be solved
		*@param	t current time
		*@param	u current values of unknowns (Tp and mp)
		*@param	dudt current values of derivatives of unknowns (Tp and mp) with respect to time
		*/
		void Equations(const double t, const double* u, double* dudt) const;

		/**
		*@brief Prints the results on a file
		*@param	filename file name
		*/
		void Print(const std::string filename) const;

		/**
		*@brief Sets the initial diameter of soot particles
		*@param	dp0 initial diameter (in m)
		*/
		void SetInitialDiameter(const double dp0);

		/**
		*@brief Sets the initial temperature of soot particles
		*@param	Tp0 initial temperature (in K)
		*/
		void SetInitialSootParticlesTemperature(const double Tp0);

		/**
		*@brief Sets the gas mixture temperature
		*@param	Tg gas mixture temperature (in K)
		*/
		void SetGasMixtureTemperature(const double Tg);

		/**
		*@brief Sets the gas mixture pressure
		*@param	p gas mixture pressure (in Pa)
		*/
		void SetGasMixturePressure(const double p);

		/**
		*@brief Sets the final time of integration
		*@param	tf final time of integration (in s)
		*/
		void SetFinalTime(const double tf);

		/**
		*@brief Sets the time step of integration
		*@param	dt time step of integration (in s)
		*/
		void SetTimeStep(const double dt);

	
	private:

		LIISignalModel&	lii_;			//!< lii signal model

		std::vector<double> t_;			//!< solution vector: time (in s)
		std::vector<double> Tp_;		//!< solution vector: soot particle temperature (in K)
		std::vector<double> mp_;		//!< solution vector: soot particle mass (in kg)
		std::vector<double> dp_;		//!< solution vector: soot particle diameter (in m)
		std::vector<double> J_;			//!< solution vector: evaporation mass flow rate (in kg/s)
		std::vector<double> SLII_;		//!< solution vector: LII signal (in W/m)
		std::vector<double> Qabs_;		//!< solution vector: absorbed heat (in W)
		std::vector<double> Qcon_;		//!< solution vector: conduction heat (in W)
		std::vector<double> Qeva_;		//!< solution vector: evaporation heat (in W)
		std::vector<double> Qrad_;		//!< solution vector: radiation heat (in W)
		std::vector<double> Qtot_;		//!< solution vector: total heat (in W)

		double Tp0_;					//!< initial temperature of soot particles (in K)
		double dp0_;					//!< initial diameter of soot particles (in m)
		double p_;						//!< pressure (in Pa)
		double Tg_;						//!< gas mixture temperature (in K)
		double tf_;						//!< final time of integration (in s)
		double dt_;						//!< time step of integration (in s)

	private:

		static const double pi_;		//!< pi

	};

}

#include "LIISignalSimulator.hpp"

#endif /* OpenSMOKE_LIISignalSimulator_H */
