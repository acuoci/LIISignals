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
	const double LIISignalModel::kB_ = 1.38064852e-23;		// [m2.kg/s2/K]
	const double LIISignalModel::sigmaSB_ = 5.670367e-8;	// [W/m2/K4]
	const double LIISignalModel::R_ = 8314.4598;			// [J/kmol/K]
	const double LIISignalModel::pi_ = 3.14159265359;		// [-]
	const double LIISignalModel::Nav_ = 6.022e26;			// [1/kmol]


	LIISignalModel::LIISignalModel(GasMixture& gas, SootParticles& soot) :
		gas_(gas),
		soot_(soot)
	{
		// Radiation
		{
			theta_ = 1.;	// total emission coefficient (for black-body equal to 1)
		}

		// Absorption
		{
			Em_ = 0.4;					// soot absorption function
			lambda_ex_ = 1064e-9;		// laser excitation wavelength (m)
		}

		// Heat conduction
		{
			conduction_model_ = HEAT_CONDUCTION_FREE_MOLECULAR;		// heat condution model
			alpha_ = 0.3;											// thermal accomodation coefficient	
		}

		// Evaporation
		{
			beta_ = 1;	// evaporation coefficient
		}

		// Laser temporal intensity signal
		{
			laser_model_ = LASER_MODEL_GAUSSIAN;	// model
			phi_ = 1000.;							// laser fluence (in J/m2)
			mu_ = 60e-9;							// mean time of Gaussian signal (in s)
			FHMW_ = 7e-9;							// FWHM of temporal laser beam profile (s)
		}
	}

	void LIISignalModel::SetHeatCondutionModel(const HeatCondutionModel flag)
	{
		conduction_model_ = flag;
	}

	void LIISignalModel::SetTotalEmissionCoefficient(const double theta)
	{
		theta_ = theta;
	}
	
	void LIISignalModel::SetSootAbsorptionFunction(const double Em)
	{
		Em_ = Em;
	}

	void LIISignalModel::SetLaserExcitationWavelength(const double lambda)
	{
		lambda_ex_ = lambda;
	}

	void LIISignalModel::SetThermalAccomodationCoefficient(const double alpha)
	{
		alpha_ = alpha;
	}

	void LIISignalModel::SetEvaporationCoefficient(const double beta)
	{
		beta_ = beta;
	}

	void LIISignalModel::SetLaserFluence(const double phi)
	{
		phi_ = phi;
	}

	void LIISignalModel::SetLaserMeanTime(const double mu)
	{
		mu_ = mu;
	}

	void LIISignalModel::SetLaserFHMW(const double FHMW)
	{
		FHMW_ = FHMW;
	}

	double LIISignalModel::QRadiation(const double Tp, const double Tg, const double dp)
	{
		return pi_ * (dp*dp)*theta_*sigmaSB_*(std::pow(Tp, 4.) - std::pow(Tg, 4.));
	}

	double LIISignalModel::QAbsorption(const double t, const double dp)
	{
		const double F = LaserTemporalIntensity(t);
		const double Cabs = AdsorptionCrossSection(dp);

		return Cabs*F;
	}

	double LIISignalModel::QConduction(const double Tp, const double Tg, const double p, const double dp)
	{
		if (conduction_model_ == HEAT_CONDUCTION_FREE_MOLECULAR)
		{
			const double fa = 1./(gas_.Gamma(Tg)-1.);
			const double fb = 1./(gas_.Gamma(Tp) - 1.);
			const double fc = 1./(gas_.Gamma((Tg+Tp)/2.)-1.);
			const double I = (fa + 4.*fc + fb) / 6.;
			const double gammaStar = 1.+1./I;

			return alpha_ * pi_* (dp*dp) * p / 8. * std::sqrt(8.*R_*Tg / pi_ / gas_.M()) *
								(gammaStar + 1.) / (gammaStar - 1.)*(Tp / Tg - 1.);
		}

		else if (conduction_model_ == HEAT_CONDUCTION_CONTINUUM)
		{
			const double fa = gas_.ThermalConductivity(Tg);
			const double fc = gas_.ThermalConductivity((Tg+Tp)/2.);
			const double fb = gas_.ThermalConductivity(Tp);
			const double I = (Tp - Tg) / 6.*(fa + 4.*fc + fb);

			return 2.*pi_*dp*I;
		}

		else if (conduction_model_ == HEAT_CONDUCTION_TRANSITION_MCCOY_CHA)
		{
			const double kg = gas_.ThermalConductivity(Tg);
			const double Mg = gas_.M();
			const double gamma = gas_.Gamma(Tg);
			const double f = (9 * gamma - 5.) / 4.;
			const double G = 8.*f / alpha_ / (gamma + 1.);
			const double lambda = kg / f / p * (gamma - 1.)*std::sqrt(pi_*Mg*Tg / 2. / R_);


			return 2.*kg*pi_*dp*dp*(Tp-Tg)/(dp+G*lambda);
		}

		return 0;
	}

	double LIISignalModel::QEvaporation(const double Tp, const double J)
	{
		const double deltaHvs = soot_.VaporizationHeat(Tp);	// [J/kmol]
		const double Mvs = soot_.MolecularWeight(Tp);       // [kg/kmol]

		return deltaHvs/Mvs*J;								// [W]

	}

	double LIISignalModel::JEvaporation(const double Tp, const double Tg, const double p, const double dp)
	{
		const double Mvs = soot_.MolecularWeight(Tp);			// [kg/kmol]
		const double NFMs = FluxFreeMolecularRegime(Tp, dp);	// [1/m2/s]
		const double NCs  = FluxContinuumRegime(Tp, Tg, p, dp);	// [1/m2/s]

		const double Nvs = 1. / (1. / NFMs + 1. / NCs);			// [1/m2/s]

		return pi_ * (dp*dp) * Nvs*Mvs / Nav_;					// [kg/s]
	}

	double LIISignalModel::LaserTemporalIntensity(const double t)
	{
		if (laser_model_ == LASER_MODEL_GAUSSIAN)
		{
			const double sigma = FHMW_ / (2. * std::sqrt(2. * log(2.)));
			return phi_ / sigma / std::sqrt(2.*pi_)*std::exp(-( (t - mu_)*(t - mu_))/(2.*sigma*sigma));
		}

		return 0;
	}

	double LIISignalModel::AdsorptionCrossSection(const double dp)
	{
		return pi_*pi_/lambda_ex_* (dp*dp*dp) * Em_;
	}

	double LIISignalModel::FluxFreeMolecularRegime(const double Tp, const double dp)
	{
		const double pvs = soot_.VaporPressure(Tp);				// [Pa]
		const double Mvs = soot_.MolecularWeight(Tp);			// [kg/kmol]

		return beta_*pvs/kB_/Tp * std::sqrt(R_*Tp/(2.*pi_*Mvs));   // [1/m2/s]

	}

	double LIISignalModel::FluxContinuumRegime(const double Tp, const double Tg, const double p, const double dp)
	{
		const double pvs = soot_.VaporPressure(Tp);									// [Pa]
		const double Gammavs = soot_.DiffusionCoefficient(Tp, p, gas_.Gamma(Tg));	// [m2/s]

		return 2. * pvs / kB_ / Tp * Gammavs / dp;   // [1/m2/s]
	}

	double LIISignalModel::SootSpectralEmissivity(const double dp)
	{
		return 4.*AdsorptionCrossSection(dp)/pi_/(dp*dp);
	}
}
