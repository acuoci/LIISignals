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

#include "RungeKutta4thOrder.h"

namespace OpenSMOKE
{
	LIISignalSimulator* ptSimulator;

	void ptEquations(const double t, const double* u, double* dudt)
	{
		ptSimulator->Equations(t, u, dudt);
	}

	const double LIISignalSimulator::pi_ = 3.14159265359;		// [-]

	LIISignalSimulator::LIISignalSimulator(LIISignalModel& lii) :
		lii_(lii)
	{
		ptSimulator = this;

		Tp0_ = 1700.10;			// initial temperature of soot particles (in K)
		dp0_ = 20e-9;			// initial diameter of soot particles (in m)

		Tg_ = 1700.;			// temperature of gas mixture (in K)
		p_ = 1e5;				// pressure of gas mixture (in Pa)

		tf_ = 2e-6;				// final time of integration (in s)
		dt_ = 2e-10;			// time step of integration (in s)
	}


	void LIISignalSimulator::SetInitialDiameter(const double dp0)
	{
		dp0_ = dp0;
	}

	void LIISignalSimulator::SetInitialSootParticlesTemperature(const double Tp0)
	{
		Tp0_ = Tp0;
	}

	void LIISignalSimulator::SetGasMixtureTemperature(const double Tg)
	{
		Tg_ = Tg;
	}

	void LIISignalSimulator::SetGasMixturePressure(const double p)
	{
		p_ = p;
	}

	void LIISignalSimulator::SetFinalTime(const double tf)
	{
		tf_ = tf;
	}

	void LIISignalSimulator::SetTimeStep(const double dt)
	{
		dt_ = dt;
	}

	void LIISignalSimulator::Solve()
	{
		const double mp0 = pi_ / 6.*std::pow(dp0_, 3.)*lii_.soot().Density();

		OpenSMOKE::RungeKutta4thOrder rk4;

		std::vector<double> y(2);
		y[0] = Tp0_;
		y[1] = mp0;

		rk4.SetInitialConditions(2, 0., y.data());
		rk4.SetFinalTime(tf_);
		rk4.SetTimeStep(dt_);
		rk4.SetOdeSystem(ptEquations);
		rk4.Solve();

		// Recovering the solution
		unsigned int n = rk4.solution()[0].size();
		t_ = rk4.solution()[0];
		Tp_ = rk4.solution()[1];
		mp_ = rk4.solution()[2];
		dp_.resize(n);
		J_.resize(n);
		SLII_.resize(n);
		Qabs_.resize(n);
		Qcon_.resize(n);
		Qeva_.resize(n);
		Qrad_.resize(n);
		Qtot_.resize(n);

		
		for (unsigned int k = 0; k < n; k++)
		{
			dp_[k] = std::pow(6.*mp_[k] / pi_ / lii_.soot().Density(), 1. / 3.);

			J_[k] = lii_.JEvaporation(Tp_[k], Tg_, p_, dp_[k]);
			SLII_[k] = lii_.LIISignal(dp_[k], Tp_[k]);

			Qabs_[k] = lii_.QAbsorption(t_[k], dp_[k]);
			Qcon_[k] = lii_.QConduction(Tp_[k], Tg_, p_, dp_[k]);
			Qeva_[k] = lii_.QEvaporation(Tp_[k], J_[k]);
			Qrad_[k] = lii_.QRadiation(Tp_[k], Tg_, dp_[k]);
			Qtot_[k] = Qabs_[k] - (Qcon_[k] + Qeva_[k] + Qrad_[k]);
			
		}
	}

	void LIISignalSimulator::Equations(const double t, const double* u, double* dudt) const
	{
		// Recover main unknowns
		const double Tp = u[0];
		const double mp = u[1];

		// Reconstruct diameter
		const double dp = std::pow(6.*mp / pi_ / lii_.soot().Density(), 1. / 3.);

		// Vaporization mass flux
		const double J = lii_.JEvaporation(Tp, Tg_, p_, dp);

		// Heats
		const double Qabs = lii_.QAbsorption(t, dp);
		const double Qcon = lii_.QConduction(Tp, Tg_, p_, dp);
		const double Qeva = lii_.QEvaporation(Tp, J);
		const double Qrad = lii_.QRadiation(Tp, Tg_, dp);
		const double Qtot = Qabs - (Qcon + Qeva + Qrad);

		// Heat balance equation
		const bool include_derivative_of_cp_ = false;
		double coefficient = lii_.soot().MassSpecificHeatConstantPressure(Tp);		
		if (include_derivative_of_cp_ == true)
			coefficient += Tp*lii_.soot().DerivativeMassSpecificHeatConstantPressure(Tp);
		const double dTpdt = Qtot / mp / coefficient;

		// Mass balance equation
		const double dmpdt = -J;

		// Recover ODE rhs
		dudt[0] = dTpdt;
		dudt[1] = dmpdt;
	}

	void LIISignalSimulator::Print(const std::string filename) const
	{
		std::ofstream fOut(filename.c_str(), std::ios::out);
		fOut.setf(std::ios::scientific);

		fOut << std::left << std::setw(15) << "t[s]";
		fOut << std::left << std::setw(15) << "LII[W/m]";
		fOut << std::left << std::setw(15) << "Tp[K]";
		fOut << std::left << std::setw(15) << "dp[m]";
		fOut << std::left << std::setw(15) << "Qabs[W]";
		fOut << std::left << std::setw(15) << "Qcon[W]";
		fOut << std::left << std::setw(15) << "Qeva[W]";
		fOut << std::left << std::setw(15) << "Qrad[W]";
		fOut << std::left << std::setw(15) << "Qtot[W]";
		fOut << std::left << std::setw(15) << "J[kg/s]";
		fOut << std::endl;

		for (unsigned int k = 0; k < t_.size(); k++)
		{
			fOut << std::left << std::setw(15) << t_[k];
			fOut << std::left << std::setw(15) << SLII_[k];		// TODO
			fOut << std::left << std::setw(15) << Tp_[k];
			fOut << std::left << std::setw(15) << dp_[k];
			fOut << std::left << std::setw(15) << Qabs_[k];
			fOut << std::left << std::setw(15) << Qcon_[k];
			fOut << std::left << std::setw(15) << Qeva_[k];
			fOut << std::left << std::setw(15) << Qrad_[k];
			fOut << std::left << std::setw(15) << Qtot_[k];
			fOut << std::left << std::setw(15) << J_[k];
			fOut << std::endl;
		}

		fOut.close();
	}
}
