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
	const double GasMixture::R_ = 8314.4598;				// [J/kmol/K]

	GasMixture::GasMixture()
	{
		mw_ = 28.;
	}

	void GasMixture::SetMolecularWeight(const double MW)
	{
		mw_ = MW;
	}

	double GasMixture::MassSpecificHeatConstantPressure(const double T) const
	{
		// N2 thermodynamic coefficients from CRECK Modeling

		const double Ti = 1050.;	// [K]
		const double Tm = 3800.;	// [K]

		if (T <= Ti)		// Low-temperature region
		{
			return (1.143621e+03+T*(-7.243422e-01+T*(1.588338e-03+T*(-1.114793e-06+T* 2.738496e-10))));		// [J/kg/K]
		}
		else if (T <= Tm)	// High-temperature region
		{
			return (8.051735e+02 + T * (5.649814e-01 + T * (-2.535527e-04 + T * (5.466132e-08 - T*4.591919e-12))));		// [J/kg/K]
		}
		else
		{
			return (8.051735e+02 + Tm * (5.649814e-01 + Tm * (-2.535527e-04 + Tm * (5.466132e-08 - Tm * 4.591919e-12))));		// [J/kg/K]

		}
	}

	double GasMixture::MassSpecificHeatConstantVolume(const double T) const
	{
		const double Cp = MoleSpecificHeatConstantPressure(T);	// [J/kmol/K]

		return (Cp - R_) / mw_;									// [J/kg/K]
	}

	double GasMixture::MoleSpecificHeatConstantPressure(const double T) const
	{
		return MassSpecificHeatConstantPressure(T)*mw_;
	}

	double GasMixture::MoleSpecificHeatConstantVolume(const double T) const
	{
		return MassSpecificHeatConstantVolume(T)*mw_;
	}

	double GasMixture::Gamma(const double T) const
	{
		const double Cp = MoleSpecificHeatConstantPressure(T);	// [J/kmol/K]
		const double Cv = Cp-R_;								// [J/kmol/K]

		return Cp / Cv;
	}

	double GasMixture::ThermalConductivity(const double T) const
	{
		// N2 coefficients from CRECK Modeling
		const double logT = std::log(T);
		return std::exp(1.188512e-02 + logT*(-2.896464e+00 + logT*(5.514232e-01 - logT*2.729252e-02)));  // [W/m/K]
	}
}
