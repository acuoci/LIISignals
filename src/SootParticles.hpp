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

namespace OpenSMOKE
{
	const double SootParticles::kB_ = 1.38064852e-23;	// [m2.kg/s2/K]
	const double SootParticles::R_ = 8314.4598;			// [J/kmol/K]
	const double SootParticles::pi_ = 3.14159265359;	// [-]

	SootParticles::SootParticles()
	{
		pv_correlation_ = PV_HOFMAN_2007;
		mv_correlation_ = MV_HOFMAN_2007;
		hv_correlation_ = HV_HOFMAN_2007;
		density_ = 1860.;	
	}

	void SootParticles::SetCorrelationVaporPressure(const CorrelationVaporPressure flag)
	{
		pv_correlation_ = flag;
	}

	void SootParticles::SetCorrelationVaporizationHeat(const CorrelationVaporizationHeat flag)
	{
		hv_correlation_ = flag;
	}

	void SootParticles::SetCorrelationMolecularWeight(const CorrelationMolecularWeight flag)
	{
		mv_correlation_ = flag;
	}

	void SootParticles::SetDensity(const double density)
	{
		density_ = density;
	}

	double SootParticles::Density() const
	{
		return density_;
	}

	double SootParticles::MassSpecificHeatConstantPressure(const double T) const
	{
		const double a = 1878.;			// [J/kg/K]
		const double b = 0.1082;		// [J/kg/K2]
		const double c = -1.5149e8;		// [JK/kg]

		return a + b * T + c / T / T;	// [J/kg/K]
	}

	double SootParticles::DerivativeMassSpecificHeatConstantPressure(const double T) const
	{
		const double b = 0.1082;		// [J/kg/K2]
		const double c = -1.5149e8;		// [JK/kg]

		return b - 2.*c / T / T / T;	// [J/kg/K2]
	}

	double SootParticles::VaporPressure(const double T) const
	{
		if (pv_correlation_ == PV_HOFMAN_2007)
		{
			return std::exp(-111.4 + T*(0.0906 +T*(-2.764e-5 + T*(4.175e-9 - T*2.488e-13))));   // [Pa]
		}

		else if (pv_correlation_ == PV_GOULDER_2002)
		{
			return std::exp(-122.96 + T*(0.0906 + T*(-2.764e-5 + T*(4.175e-9 - T*2.488e-13)))) * 101325;   // [Pa]
		}

		else if (pv_correlation_ == PV_MICHELSEN_2008)
		{
			return std::exp(-198.11 + T*(0.20732 +T*(-9.6899e-05 + T*(2.3826e-08 + T*(-2.9423e-12 + T*1.4363e-16))))) * 1e5;   //[Pa]
		}

		return 0;
	}

	double SootParticles::MolecularWeight(const double T) const
	{
		if (mv_correlation_ == MV_HOFMAN_2007)
		{
			return (0.01718 + T*(6.865e-7 + T*(2.996e-9 +T*(-8.595e-13 + T*1.049e-16)))) * 1000;   // [kg/kmol]
		}

		else if (mv_correlation_ == MV_MICHELSEN_2008)
		{
			return (82.299 + T*(-0.19406 + T*(1.9146e-4 + T*(-8.4923e-8 + T*(1.9306e-11 +T*(-2.2032e-15 + T*1.0052e-19))))));   //[kg/kmol]
		}

		return 0;
	}

	double SootParticles::CrossSection(const double T) const
	{
		const double M = MolecularWeight(T) / 1000.;	// [kg/mol]
		return (1.8e-19 +M*(-1.857e-17 + M*(1.404e-15 +M*(-2.593e-14 + M*(2.075e-13 - M*6.667e-13)))));   //[m2]
	}

	double SootParticles::VaporizationHeat(const double T) const
	{
		if (hv_correlation_ == HV_HOFMAN_2007)
		{
			return 1e3*(205398 + T*(736.6 + T*(-0.4071 + T*(1.199e-4 + T*(-1.795e-8+T*1.072e-12)))));   // [J/kmol]
		}

		else if (hv_correlation_ == HV_MICHELSEN_2008)	// Michelsen (2008)
		{
			return 1e3*(9.2892e5 + T*(-609.8 + T*(0.63252 + T*(-2.9484e-4 + T*(7.0101e-8 + T*(-8.4111e-12 + T*4.088e-16))))));   //[J/kmol]
		}

		return 0;
	}

	double SootParticles::DiffusionCoefficient(const double T, const double p, const double gamma) const
	{
		const double Mv = MolecularWeight(T);		// [kg/kmol]
		const double sigma = CrossSection(T);		// [m2]
		const double f = (9.*gamma - 5.) / 4.;		// [-]

		return f * kB_*T / (4.*sigma*p) * std::sqrt(R_*T / pi_ / Mv);
	}
}
