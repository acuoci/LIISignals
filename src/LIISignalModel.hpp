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
	const double LIISignalModel::h_ = 6.62607004e-34;		// [m2.kg/s]
	const double LIISignalModel::c_ = 299792458;			// [m/s]

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

		// Detection window
		{
			lambda_min_ = 545e-9;
			lambda_max_ = 555e-9;
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

	void LIISignalModel::SetDetectionWindow(const double lambda_min, const double lambda_max)
	{
		lambda_min_ = lambda_min;
		lambda_max_ = lambda_max;
	}

	double LIISignalModel::QRadiation(const double Tp, const double Tg, const double dp) const
	{
		return pi_ * (dp*dp)*theta_*sigmaSB_*(std::pow(Tp, 4.) - std::pow(Tg, 4.));
	}

	double LIISignalModel::QAbsorption(const double t, const double dp) const
	{
		const double F = LaserTemporalIntensity(t);
		const double Cabs = AdsorptionCrossSection(dp);

		return Cabs * F;
	}

	double LIISignalModel::QConduction(const double Tp, const double Tg, const double p, const double dp) const
	{
		if (conduction_model_ == HEAT_CONDUCTION_FREE_MOLECULAR)
		{
			const double fa = 1. / (gas_.Gamma(Tg) - 1.);
			const double fb = 1. / (gas_.Gamma(Tp) - 1.);
			const double fc = 1. / (gas_.Gamma((Tg + Tp) / 2.) - 1.);
			const double I = (fa + 4.*fc + fb) / 6.;
			const double gammaStar = 1. + 1. / I;

			return alpha_ * pi_* (dp*dp) * p / 8. * std::sqrt(8.*R_*Tg / pi_ / gas_.M()) *
				(gammaStar + 1.) / (gammaStar - 1.)*(Tp / Tg - 1.);
		}

		else if (conduction_model_ == HEAT_CONDUCTION_CONTINUUM)
		{
			const double fa = gas_.ThermalConductivity(Tg);
			const double fc = gas_.ThermalConductivity((Tg + Tp) / 2.);
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


			return 2.*kg*pi_*dp*dp*(Tp - Tg) / (dp + G * lambda);
		}

		else if (conduction_model_ == HEAT_CONDUCTION_TRANSITION_FUCHS)
		{
			const double eps = 1e-6;
			double Tmin = (1. - eps)*Tg + eps * Tp;
			double Tmax = (1. - eps)*Tp + eps * Tg;
			double Tdelta = (Tg + Tp) / 2.;

			double fmin = QConductionTransitionFuchsContinuum(Tmin, Tg, p, dp) -
				QConductionTransitionFuchsFreeMolecular(Tp, Tmin, p, dp);

			double fmax = QConductionTransitionFuchsContinuum(Tmax, Tg, p, dp) -
				QConductionTransitionFuchsFreeMolecular(Tp, Tmax, p, dp);

			double fdelta = QConductionTransitionFuchsContinuum(Tdelta, Tg, p, dp) -
				QConductionTransitionFuchsFreeMolecular(Tp, Tdelta, p, dp);

			for (unsigned int k = 0; k < 12; k++)
			{
				if (fmin*fdelta < 0.)
				{
					Tmax = Tdelta;
					fmax = fdelta;

					Tdelta = (Tmin + Tdelta) / 2.;
					fdelta = QConductionTransitionFuchsContinuum(Tdelta, Tg, p, dp) -
						QConductionTransitionFuchsFreeMolecular(Tp, Tdelta, p, dp);
				}
				else
				{
					Tmin = Tdelta;
					fmin = fdelta;

					Tdelta = (Tmax + Tdelta) / 2.;
					fdelta = QConductionTransitionFuchsContinuum(Tdelta, Tg, p, dp) -
						QConductionTransitionFuchsFreeMolecular(Tp, Tdelta, p, dp);
				}
			}

			return QConductionTransitionFuchsFreeMolecular(Tp, Tdelta, p, dp);;
		}

		return 0;
	}

	double LIISignalModel::QEvaporation(const double Tp, const double J) const
	{
		const double deltaHvs = soot_.VaporizationHeat(Tp);	// [J/kmol]
		const double Mvs = soot_.MolecularWeight(Tp);       // [kg/kmol]

		return deltaHvs / Mvs * J;								// [W]

	}

	double LIISignalModel::JEvaporation(const double Tp, const double Tg, const double p, const double dp) const
	{
		const double Mvs = soot_.MolecularWeight(Tp);			// [kg/kmol]
		const double NFMs = FluxFreeMolecularRegime(Tp, dp);	// [1/m2/s]
		const double NCs = FluxContinuumRegime(Tp, Tg, p, dp);	// [1/m2/s]

		const double Nvs = 1. / (1. / NFMs + 1. / NCs);			// [1/m2/s]

		return pi_ * (dp*dp) * Nvs*Mvs / Nav_;					// [kg/s]
	}

	double LIISignalModel::LaserTemporalIntensity(const double t) const
	{
		if (laser_model_ == LASER_MODEL_GAUSSIAN)
		{
			const double sigma = FHMW_ / (2. * std::sqrt(2. * log(2.)));
			return phi_ / sigma / std::sqrt(2.*pi_)*std::exp(-((t - mu_)*(t - mu_)) / (2.*sigma*sigma));
		}

		return 0;
	}

	double LIISignalModel::AdsorptionCrossSection(const double dp) const
	{
		return pi_ * pi_ / lambda_ex_ * (dp*dp*dp) * Em_;
	}

	double LIISignalModel::FluxFreeMolecularRegime(const double Tp, const double dp) const
	{
		const double pvs = soot_.VaporPressure(Tp);				// [Pa]
		const double Mvs = soot_.MolecularWeight(Tp);			// [kg/kmol]

		return beta_ * pvs / kB_ / Tp * std::sqrt(R_*Tp / (2.*pi_*Mvs));   // [1/m2/s]

	}

	double LIISignalModel::FluxContinuumRegime(const double Tp, const double Tg, const double p, const double dp) const
	{
		const double pvs = soot_.VaporPressure(Tp);									// [Pa]
		const double Gammavs = soot_.DiffusionCoefficient(Tp, p, gas_.Gamma(Tg));	// [m2/s]

		return 2. * pvs / kB_ / Tp * Gammavs / dp;   // [1/m2/s]
	}

	double LIISignalModel::QConductionTransitionFuchsFreeMolecular(const double Tp, const double Tdelta, const double p, const double dp) const
	{
		const double fa = 1. / (gas_.Gamma(Tdelta) - 1.);
		const double fb = 1. / (gas_.Gamma(Tp) - 1.);
		const double fc = 1. / (gas_.Gamma((Tdelta + Tp) / 2.) - 1.);
		const double I = (fa + 4.*fc + fb) / 6.;
		const double gammaStar = 1. + 1. / I;

		return	alpha_ * pi_* (dp*dp) * p / 8. * std::sqrt(8.*R_*Tdelta / pi_ / gas_.M()) *
			(gammaStar + 1.) / (gammaStar - 1.)*(Tp / Tdelta - 1.);
	}

	double LIISignalModel::QConductionTransitionFuchsContinuum(const double Tdelta, const double Tg, const double p, const double dp) const
	{
		const double fa = gas_.ThermalConductivity(Tg);
		const double fc = gas_.ThermalConductivity((Tg + Tdelta) / 2.);
		const double fb = gas_.ThermalConductivity(Tdelta);
		const double I = (Tdelta - Tg) / 6.*(fa + 4.*fc + fb);

		const double kg = gas_.ThermalConductivity(Tdelta);
		const double Mg = gas_.M();
		const double gamma = gas_.Gamma(Tdelta);
		const double f = (9 * gamma - 5.) / 4.;

		const double lambda = kg / f / p * (gamma - 1.)*std::sqrt(pi_*Mg*Tdelta / 2. / R_);

		const double Lambda1 = 1. + (2.*lambda / dp);
		const double Lambda2 = 1. + (2.*lambda / dp)*(2.*lambda / dp);

		const double coefficient = 0.2*std::pow(Lambda1, 5.) - 1. / 3.*Lambda2*std::pow(Lambda1, 3.) + 2. / 15.*std::pow(Lambda2, 2.5);
		const double delta = std::pow(dp / 2., 3.) / lambda / lambda * coefficient - dp / 2.;

		return 4.*pi_*(delta + dp / 2.)*I;
	}


	double LIISignalModel::SootSpectralEmissivity(const double dp) const
	{
		return 4.*AdsorptionCrossSection(dp) / pi_ / (dp*dp);
	}

	double LIISignalModel::Omega(const double lambda) const
	{
		return 1.;	// TODO
	}


	double LIISignalModel::LIISignal(const double dp, const double Tp) const
	{
		// Integral calculation
		const unsigned int n = 20;
		const double dlambda = (lambda_max_ - lambda_min_) / double(n - 1);
		const double epsilon = SootSpectralEmissivity(dp);

		double I = 0.;
		for (unsigned int i = 0; i < n; i++)
		{
			const double lambda_a = lambda_min_ + i * dlambda;
			const double lambda_b = lambda_a + dlambda;
			const double omega_a = Omega(lambda_a);
			const double omega_b = Omega(lambda_b);
			const double f_a = omega_a * epsilon / std::pow(lambda_a, 5.) / (std::exp(h_*c_ / lambda_a / kB_ / Tp) - 1.);
			const double f_b = omega_b * epsilon / std::pow(lambda_b, 5.) / (std::exp(h_*c_ / lambda_b / kB_ / Tp) - 1.);

			I += dlambda * 0.50*(f_a + f_b);
		}

		return 2 * (pi_*pi_)*h_*(c_*c_)*(dp*dp)*I;
	}
}
