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

#include "LogNormalDistribution.h"
#include "UserDefinedDistribution.h"
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

		enum ParticleSizeDistributionFunction { DISTRIBUTION_MONODISPERSED, DISTRIBUTION_LOGNORMAL, DISTRIBUTION_USER_DEFINED };

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

		/**
		*@brief Sets the time window for averaging
		*@param	dt_window the time window for averaging (in s)
		*/
		void SetTimeWindow(const double dt_window);

		/**
		*@brief Sets the log-normal particle size function
		*@param	log_normal_pdf the log normal particle size distribution function
		*/
		void SetLogNormalParticleSizeDistributionFunction(OpenSMOKE::LogNormalDistribution& log_normal_pdf);

		/**
		*@brief Sets a user defined distribution function
		*@param	user_defined_pdf the user defined particle size distribution function
		*/
		void SetUserDefinedParticleSizeDistributionFunction(OpenSMOKE::UserDefinedDistribution& user_defined_pdf);

		/**
		*@brief Returns the LII signal (in W/m)
		*@return the LII signal (in W/m)
		*/
		const std::vector<double>& SLII() const { return SLII_; }

		/**
		*@brief Returns the normalized signal
		*@return the normalized signal
		*/
		const std::vector<double>& nSLII() const { return nSLII_; }

		/**
		*@brief Returns the particle temperature (in K)
		*@return the particle temperature (in K)
		*/
		const std::vector<double>& Tp() const { return Tp_; }

		/**
		*@brief Returns the LII signal (in W/m) averaged over the user defined window time
		*@return the averaged LII signal (in W/m)
		*/
		const std::vector<double>& SLII_averaged() const { return SLII_averaged_; }

		/**
		*@brief Returns the normalized signal averaged over the user defined window time
		*@return the averaged normalized signal
		*/
		const std::vector<double>& nSLII_averaged() const { return nSLII_averaged_; }

		/**
		*@brief Returns the particle temperature (in K) averaged over the user defined window time
		*@return the averaged particle temperature (in K)
		*/
		const std::vector<double>& Tp_averaged() const { return Tp_averaged_; }


	private:

		/**
		*@brief Simulates the temporal evolution of temperature and mass of soot particles according to a monodispersed distribution
		*/
		void SolveMonodispersedDistribution();

		/**
		*@brief Simulates the temporal evolution of temperature and mass of soot particles according to a log-normal distribution
		*/
		void SolveLogNormalDistribution();

		/**
		*@brief Simulates the temporal evolution of temperature and mass of soot particles according to a user-defined distribution
		*/
		void SolveUserDefinedDistribution();

		/**
		*@brief Reconstructs the normalized signal at the end of the simulations
		*/
		void NormalizedSignal();

		/**
		*@brief Averages the quantities of interest over the user-specified window
		*/
		void Averaging();


	private:

		LIISignalModel & lii_;			//!< lii signal model

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
		std::vector<double> nSLII_;		//!< solution vector: normalized LII signal

		std::vector<double> Tp_averaged_;			//!< solution vector: averaged soot particle temperature (in K)
		std::vector<double> SLII_averaged_;			//!< solution vector: averaged LII signal (in W/m)
		std::vector<double> nSLII_averaged_;		//!< solution vector: averaged normalized LII signal


		double Tp0_;					//!< initial temperature of soot particles (in K)
		double dp0_;					//!< initial diameter of soot particles (in m)
		double p_;						//!< pressure (in Pa)
		double Tg_;						//!< gas mixture temperature (in K)
		double tf_;						//!< final time of integration (in s)
		double dt_;						//!< time step of integration (in s)
		double dt_window_;				//!< window time for averaging (in s)

		ParticleSizeDistributionFunction distribution_;			//!< particle size diameter distribution function
		OpenSMOKE::LogNormalDistribution* log_normal_pdf_;		//!< pointer to a log-normal distribution
		OpenSMOKE::UserDefinedDistribution* user_defined_pdf_;	//!< pointer to a user-defined distribution

	private:

		static const double pi_;		//!< pi

	};

}

#include "LIISignalSimulator.hpp"

#endif /* OpenSMOKE_LIISignalSimulator_H */
