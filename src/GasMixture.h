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

#ifndef OpenSMOKE_GasMixture_H
#define OpenSMOKE_GasMixture_H

namespace OpenSMOKE
{
	//!  A class describing the properties of a gaseous mixture
	/*!
	This class provides the correlations for describing the properties of a gaseous mixture
	*/

	class GasMixture
	{
	
	public:

		/**
		*@brief Default constructor
		*/
		GasMixture();

		/**
		*@brief Calculation of constant pressure specific heat (mass units)
		*@param		T	temperature (in K)
		*@return		the Cp (in J/kg/K)
		*/
		double MassSpecificHeatConstantPressure(const double T);

		/**
		*@brief Calculation of constant volume specific heat (mass units)
		*@param		T	temperature (in K)
		*@return		the Cv (in J/kg/K)
		*/
		double MassSpecificHeatConstantVolume(const double T);

		/**
		*@brief Calculation of constant pressure specific heat (mole units)
		*@param		T	temperature (in K)
		*@return		the Cp (in J/kmol/K)
		*/
		double MoleSpecificHeatConstantPressure(const double T);

		/**
		*@brief Calculation of constant volume specific heat (mole units)
		*@param		T	temperature (in K)
		*@return		the Cv (in J/kmol/K)
		*/
		double MoleSpecificHeatConstantVolume(const double T);

		/**
		*@brief Calculation of heat capacity ratio gamma=Cp/Cv
		*@param		T	temperature (in K)
		*@return		the heat capacity ratio (dimensionless)
		*/
		double Gamma(const double T);

		/**
		*@brief Calculation of thermal conductivity
		*@param		T	temperature (in K)
		*@return		the thermal conductivity (in W/m/s)
		*/
		double ThermalConductivity(const double T);

		/**
		*@brief Sets the molecular weight of the mixture
		*@param	MW	molecular weight (in kg/kmol)
		*/
		void SetMolecularWeight(const double MW);

		/**
		*@brief		Returns the molecular weight
		*@return	molecular weight (in kg/kmol)
		*/
		double M() const { return mw_; }

	
	private:

		double mw_;					//!< the gas molecular weight (in kg/kmol)

	private:

		static const double R_;		//!< the constant of ideal gases (in J/kmol/K)

	};

}

#include "GasMixture.hpp"

#endif /* OpenSMOKE_GasMixture_H */
