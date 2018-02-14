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

#ifndef OpenSMOKE_SootParticles_H
#define OpenSMOKE_SootParticles_H

namespace OpenSMOKE
{

	//!  A class describing the properties of soot particles
	/*!
	This class provides the correlations for describing the properties of soot particles
	*/

	class SootParticles
	{

	public:

		enum CorrelationVaporPressure		{ PV_HOFMAN_2007, PV_GOULDER_2002, PV_MICHELSEN_2008 };
		enum CorrelationVaporizationHeat	{ HV_HOFMAN_2007, HV_MICHELSEN_2008 };
		enum CorrelationMolecularWeight		{ MV_HOFMAN_2007, MV_MICHELSEN_2008 };
	
	public:

		/**
		*@brief Default constructor
		*/
		SootParticles();

		/**
		*@brief		Calculation of density of soot
		*@return	density (in kg/m3)
		*/
		double Density();

		/**
		*@brief Calculation of constant pressure specific heat (mass units)
		*@param		T	temperature (in K)
		*@return		the Cp (in J/kg/K)
		*/
		double MassSpecificHeatConstantPressure(const double T);

		/**
		*@brief Calculation of derivative of constant pressure specific heat
				with respect to temperature (mass units)
		*@param		T	temperature (in K)
		*@return		the derivative dCp/dT (in J/kg/K2)
		*/
		double DerivativeMassSpecificHeatConstantPressure(const double T);

		/**
		*@brief Calculation of vapor pressure 
		*@param		T	temperature (in K)
		*@return		vapor pressure (in Pa)
		*/
		double VaporPressure(const double T);

		/**
		*@brief Calculation of molecular weight of soot vapors
		*@param		T	temperature (in K)
		*@return		molecular weight (in kg/kmol)
		*/
		double MolecularWeight(const double T);

		/**
		*@brief Calculation of molecular cross section of soot vapors
		*@param		T	temperature (in K)
		*@return		cross section (in m2)
		*/
		double CrossSection(const double T);

		/**
		*@brief Calculation of vaporization heat
		*@param		T	temperature (in K)
		*@return		vaporization heat (in J/kmol)
		*/
		double VaporizationHeat(const double T);

		/**
		*@brief Calculation of diffusion coefficient of soot vapors
		*@param		T		temperature (in K)
		*@param		p		pressure (in Pa)
		*@param		gamma	heat capacity ratio of gaseous mixture
		*@return	the diffusion coefficient (in m2/s)
		*/
		double DiffusionCoefficient(const double T, const double p, const double gamma);


	public:

		/**
		*@brief Sets the correlation to be adopted for vapor pressure calculation
		*@param	flag the correlation 
		*/
		void SetCorrelationVaporPressure(const CorrelationVaporPressure flag);

		/**
		*@brief Sets the correlation to be adopted for vaporization heat calculation
		*@param	flag the correlation
		*/
		void SetCorrelationVaporizationHeat(const CorrelationVaporizationHeat flag);

		/**
		*@brief Sets the correlation to be adopted for molecular weight calculation
		*@param	flag the correlation
		*/
		void SetCorrelationMolecularWeight(const CorrelationMolecularWeight flag);


	private:

		CorrelationMolecularWeight	mv_correlation_;		//!< the correlation for evaluating the molecular weight
		CorrelationVaporPressure	pv_correlation_;		//!< the correlation for evaluating the vapor pressure
		CorrelationVaporizationHeat hv_correlation_;		//!< the correlation for evaluating the vaporization heat

		static const double kB_;	//!< the Boltzmann constant (in m2.kg/s2/K)
		static const double R_;		//!< the constant of ideal gases (in J/kmol/K)
		static const double pi_;	//!< pi

	};

}

#include "SootParticles.hpp"

#endif /* OpenSMOKE_SootParticles_H */
