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

#ifndef OpenSMOKE_LIISignalModel_H
#define OpenSMOKE_LIISignalModel_H

#include "GasMixture.h"
#include "SootParticles.h"

namespace OpenSMOKE
{
	//!  A class for modeling the LII Signal
	/*!
	This class provides the tools for modeling the LII Signal
	*/

	class LIISignalModel
	{

	public:

		enum LaserTemporalIntensityModel { LASER_MODEL_GAUSSIAN };

		enum HeatCondutionModel {	HEAT_CONDUCTION_FREE_MOLECULAR, 
									HEAT_CONDUCTION_CONTINUUM,
									HEAT_CONDUCTION_TRANSITION_MCCOY_CHA,
									HEAT_CONDUCTION_TRANSITION_FUCHS };
	
	public:

		/**
		*@brief Default constructor
		*@param		gas reference to a GasMixture object
		*@param		gas reference to a SootParticles object
		*/
		LIISignalModel(GasMixture& gas, SootParticles& soot);

		/**
		*@brief Calculation of radiative power
		*@param		Tp	particle temperature (in K)
		*@param		Tg	gas temperature (in K)
		*@param		dp	particle diameter (in m)
		*@return	the radiative power (in W)
		*/
		double QRadiation(const double Tp, const double Tg, const double dp) const;

		/**
		*@brief Calculation of absorption power
		*@param		t the current time (in s)
		*@param		dp	particle diameter (in m)
		*@return	the absorption power (in W)
		*/
		double QAbsorption(const double t, const double dp) const;

		/**
		*@brief Calculation of thermal conduction power
		*@param		Tp	particle temperature (in K)
		*@param		Tg	gas temperature (in K)
		*@param		p	pressure (in Pa)
		*@param		dp	particle diameter (in m)
		*@return	the thermal conduction power (in W)
		*/
		double QConduction(const double Tp, const double Tg, const double p, const double dp) const;

		/**
		*@brief Calculation of evaporation power
		*@param		Tp	particle temperature (in K)
		*@param		J   evaporation mass flow rate (in kg/s)
		*@return	the evaporation power (in W)
		*/
		double QEvaporation(const double Tp, const double J) const;

		/**
		*@brief Calculation of evaporation mass flow rate
		*@param		Tp	particle temperature (in K)
		*@param		Tg	gas temperature (in K)
		*@param		p	pressure (in Pa)
		*@param		dp	particle diameter (in m)
		*@return	the evaporation mass flow rate (in kg/s)
		*/
		double JEvaporation(const double Tp, const double Tg, const double p, const double dp) const;

		/**
		*@brief Calculation of LII signal
		*@param		dp	particle diameter (in m)
		*@param		Tp	particle temperature (in K)
		*@return	the LII signal (in W/m)
		*/
		double LIISignal(const double dp, const double Tp) const;

		/**
		*@brief Returns a reference to the gas mixture model
		*@return	the gas mixture model
		*/
		const GasMixture& gas() const { return gas_; }

		/**
		*@brief Returns a reference to the soot particle model
		*@return	the soot particle model
		*/
		const SootParticles& soot() const { return soot_; }

	public:

		/**
		*@brief Sets the total emission coefficient
		*@param theta the total emission coefficient (for black bodies is equal to 1)
		*/
		void SetTotalEmissionCoefficient(const double theta);

		/**
		*@brief Sets the soot absorption function
		*@param theta the soot absorption function
		*/
		void SetSootAbsorptionFunction(const double Em);

		/**
		*@brief Sets the laser excitation wavelength
		*@param lambda the laser excitation wavelength (in m)
		*/
		void SetLaserExcitationWavelength(const double lambda);

		/**
		*@brief Sets the thermal accomodation factor
		*@param alpha the thermal accomodation factor
		*/
		void SetThermalAccomodationCoefficient(const double alpha);

		/**
		*@brief Sets the evaporation coefficient
		*@param beta the evaporation coefficient
		*/
		void SetEvaporationCoefficient(const double beta);

		/**
		*@brief Sets the laser fluence
		*@param phi the laser fluence (in J/m2)
		*/
		void SetLaserFluence(const double phi);

		/**
		*@brief Sets the laser mean time
		*@param phi the laser mean time (in s)
		*/
		void SetLaserMeanTime(const double mu);

		/**
		*@brief Sets the laser FHMW
		*@param phi the laser FHMW (in s)
		*/
		void SetLaserFHMW(const double FHMW);

		/**
		*@brief Sets the heat conduction model
		*@param flag the heat conduction model
		*/
		void SetHeatCondutionModel(const HeatCondutionModel flag);

	private:

		/**
		*@brief Calculation of laser temporal intensity
		*@param		t the current time (s)
		*@return	the laser temporal intensity (in J/m2/s)
		*/
		double LaserTemporalIntensity(const double t) const;

		/**
		*@brief Calculation of adsorption cross section
		*@param		dp	particle diameter (in m)
		*@return	the adsorption cross section (in m2)
		*/
		double AdsorptionCrossSection(const double dp) const;

		/**
		*@brief Calculation of flux in free molecular regime
		*@param		Tp	particle temperature (in K)
		*@param		dp	particle diameter (in m)
		*@return	the flux (in 1/m2/s)
		*/
		double FluxFreeMolecularRegime(const double Tp, const double dp) const;

		/**
		*@brief Calculation of flux in continuum regime
		*@param		Tp	particle temperature (in K)
		*@param		Tg	gas temperature (in K)
		*@param		p	pressure (in Pa)
		*@param		dp	particle diameter (in m)
		*@return	the flux (in 1/m2/s)
		*/
		double FluxContinuumRegime(const double Tp, const double Tg, const double p, const double dp) const;

		/**
		*@brief Calculation of conductive power in free molecular regime for Fuchs model
		*@param		Tp		particle temperature (in K)
		*@param		Tdelta	boundary layer temperature (in K)
		*@param		p		pressure (in Pa)
		*@param		dp		particle diameter (in m)
		*@return	the	conductive power (in W)
		*/
		double QConductionTransitionFuchsFreeMolecular(const double Tp, const double Tdelta, const double p, const double dp) const;
		
		/**
		*@brief Calculation of conductive power in continuum regime for Fuchs model
		*@param		Tdelta	boundary layer temperature (in K)
		*@param		Tg		gas temperature (in K)
		*@param		p		pressure (in Pa)
		*@param		dp		particle diameter (in m)
		*@return	the	conductive power (in W)
		*/
		double QConductionTransitionFuchsContinuum(const double Tdelta, const double Tg, const double p, const double dp) const;

		/**
		*@brief Calculation of soot spectral emissivity
		*@param		dp	particle diameter (in m)
		*/
		double SootSpectralEmissivity(const double dp) const;

		/**
		*@brief Calculation of spectral response function of the detection system 
		*@param		lambda	wave length (in m)
		*@return	the spectral response function of the detection system (in 1/m)
		*/
		double Omega(const double lambda) const;

	
	private:

		GasMixture&		gas_;		//!< gas mixture model
		SootParticles&	soot_;		//!< soot particles model

		double theta_;				//!< total emission coefficient (for black-body equal to 1)

		double Em_;					//!< Soot absorption function
		double lambda_ex_;			//!< Laser excitation wavelength (m)

		double alpha_;							//!< thermal accomodation coefficient	
		HeatCondutionModel conduction_model_;	//!< heat conduction model

		double beta_;				//!< evaporation coefficient

		double phi_;				//!< laser fluence (J/m2)
		double mu_;					//!< mean time of Gaussian signal (in s)
		double FHMW_;				//!< FWHM of temporal laser beam profile (in s)

		LaserTemporalIntensityModel laser_model_;	//!< model for temporal laser intensity

		double lambda_min_;			//!< Bandpass wavelength min (in m)
		double lambda_max_;			//!< Bandpass wavelength max (in m)


	private:

		static const double sigmaSB_;	//!< the Stefan-Boltzmann constant (in W/m2/K4)
		static const double kB_;		//!< the Boltzmann constant (in m2.kg/s2/K)
		static const double R_;			//!< the constant of ideal gases (in J/kmol/K)
		static const double pi_;		//!< pi
		static const double Nav_;		//!< Avogadro's number (in 1/kmol)
		static const double h_;			//!< Planck's constant (in m2.kg/s)
		static const double c_;			//!< speed of light (in m/s)

	};

}

#include "LIISignalModel.hpp"

#endif /* OpenSMOKE_LIISignalModel_H */
