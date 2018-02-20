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

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

#include "SootParticles.h"
#include "GasMixture.h"
#include "LIISignalModel.h"
#include "RungeKutta4thOrder.h"
#include "LogNormalDistribution.h"
#include "LIISignalSimulator.h"

int main()
{
	const double tf = 2e-6;			// total time of simulation (in s)
	const double dt = 2e-10;		// time step of integration (in s)

	std::ifstream fInput("input", std::ios::in);
	std::string dummy;

	std::getline(fInput, dummy);	// distribution type: monodispersed vs polydispersed
	std::string distribution;
	fInput >> distribution;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// heat conduction model: free-molecular, continuum, transition-mccoy-cha, transition-fuchs
	std::string heat_conduction;
	fInput >> heat_conduction;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// rectangular detection bandpass: bandpass wavelength min. (in nm)	
	double lambda_min; fInput >> lambda_min;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// Rectangular detection bandpass: bandpass wavelength max. (in nm)	
	double lambda_max; fInput >> lambda_max;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// laser wavelength (in nm)
	double lambda; fInput >> lambda;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// Full width at half maximum of temporal laser beam profile (in ns)
	double fwhm; fInput >> fwhm;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// Soot absorption function E(m)
	double Em; fInput >> Em;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// Laser fluence (in J/cm2)
	double J; fInput >> J;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// gas temperature (in K)
	double Tg; fInput >> Tg;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// pressure (in Pa)
	double p; fInput >> p;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// molecular weight of gas mixture (in kg/kmol)	
	double mw; fInput >> mw;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// thermal accommodation coefficient alpha		
	double alpha; fInput >> alpha;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// mass accommodation coefficient beta	
	double beta; fInput >> beta;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// particle density (in kg/m3)
	double density; fInput >> density;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// initial particle temperature (in K)
	double Tp0; fInput >> Tp0;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// initial mean diameter (in nm)
	double dp_mean; fInput >> dp_mean;
	std::getline(fInput, dummy);

	OpenSMOKE::SootParticles soot;
	{
		soot.SetDensity(density);
	}
	OpenSMOKE::GasMixture gas;
	{
		gas.SetMolecularWeight(mw);
	}

	OpenSMOKE::LIISignalModel lii(gas, soot);
	{
		lii.SetLaserFluence(J*1.e4);
		lii.SetLaserExcitationWavelength(lambda/1.e9);
		lii.SetLaserFHMW(fwhm/1.e9);
		lii.SetSootAbsorptionFunction(Em);
		lii.SetThermalAccomodationCoefficient(alpha);
		lii.SetEvaporationCoefficient(beta);
		lii.SetDetectionWindow(lambda_min/1.e9, lambda_max/1.e9);
		
		if (heat_conduction == "free-molecular")
			lii.SetHeatCondutionModel(OpenSMOKE::LIISignalModel::HEAT_CONDUCTION_FREE_MOLECULAR);
		else if (heat_conduction == "continuum")
			lii.SetHeatCondutionModel(OpenSMOKE::LIISignalModel::HEAT_CONDUCTION_CONTINUUM);
		else if (heat_conduction == "transition-mccoy-cha")
			lii.SetHeatCondutionModel(OpenSMOKE::LIISignalModel::HEAT_CONDUCTION_TRANSITION_MCCOY_CHA);
		else if (heat_conduction == "transition-fuchs")
			lii.SetHeatCondutionModel(OpenSMOKE::LIISignalModel::HEAT_CONDUCTION_TRANSITION_FUCHS);
		else
		{
			std::cout << "Wrong heat conduction model. Available options: free-molecular, continuum, transition-mccoy-cha, transition-fuchs" << std::endl;
			std::cout << "Press enter to continue..." << std::endl;
			getchar();
			exit(-1);
		}
	}

	OpenSMOKE::LIISignalSimulator simulator(lii);
	{
		simulator.SetInitialSootParticlesTemperature(Tp0);
		simulator.SetGasMixtureTemperature(Tg);
		simulator.SetGasMixturePressure(p);
		simulator.SetFinalTime(tf);
		simulator.SetTimeStep(dt);
		simulator.SetInitialDiameter(dp_mean/1.e9);
	}

	if (distribution == "monodispersed")
	{
		simulator.Solve();
		simulator.Print("SolutionMonodispersed.out");

	}
	else if (distribution == "polydispersed")
	{
		std::getline(fInput, dummy);						// std. deviation (in nm)
		double sigmap; fInput >> sigmap;
		std::getline(fInput, dummy);

		std::getline(fInput, dummy);						// number of quadrature points
		unsigned int n; fInput >> n;
		std::getline(fInput, dummy);

		// From mean diameter dp and std. deviation sigmap
		const double sigma = std::exp(std::sqrt(std::log(sigmap*sigmap / std::exp(2.*std::log(dp_mean)) + 1.)));
		const double d_cmd = std::exp(std::log(dp_mean) - std::log(sigma)*std::log(sigma) / 2.);

		// Print distribution of file
		OpenSMOKE::LogNormalDistribution log_normal_pdf(d_cmd, sigma, n);
		log_normal_pdf.Print("LogNormalDistribution.out");

		simulator.SetLogNormalParticleSizeDistributionFunction(log_normal_pdf);
		simulator.Solve();
		simulator.Print("SolutionPolydispersed.out");
	}

	std::cout << "Successfully done! Press enter to continue..." << std::endl;
	getchar();
	return 0;
}
