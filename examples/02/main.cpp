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
#include <sstream>

#include "SootParticles.h"
#include "GasMixture.h"
#include "LIISignalModel.h"
#include "RungeKutta4thOrder.h"
#include "UserDefinedDistribution.h"
#include "LIISignalSimulator.h"

void parse_file(const std::string file_name, int& number_of_lines, int& number_of_columns);

int main()
{
	std::ifstream fInput("input", std::ios::in);
	std::string dummy;

	std::getline(fInput, dummy);	// file containing the list of operating conditions
	std::string input_file_name;
	fInput >> input_file_name;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// file containing the list of operating conditions
	std::string output_file_name;
	fInput >> output_file_name;
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

	std::getline(fInput, dummy);	// thermal accommodation coefficient alpha		
	double alpha; fInput >> alpha;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// mass accommodation coefficient beta	
	double beta; fInput >> beta;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// particle density (in kg/m3)
	double density; fInput >> density;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// total time of simulation (in ns)
	double tf; fInput >> tf;	tf /= 1e9;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// time step (in ns)
	double dt; fInput >> dt;	dt /= 1e9;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// list output times (in ns)
	int n_times; fInput >> n_times;
	std::vector<double> out_times(n_times);
	std::vector<int>	out_steps(n_times);
	for (int i = 0; i < n_times; i++)
	{
		fInput >> out_times[i];	out_times[i] /= 1e9;
		out_steps[i] = static_cast<int>(out_times[i] / dt);
	}
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// temperature threshold
	double threshold_temperature; fInput >> threshold_temperature;
	std::getline(fInput, dummy);

	std::getline(fInput, dummy);	// mole fraction threshold
	double threshold_mole_fractions; fInput >> threshold_mole_fractions;
	std::getline(fInput, dummy);

	fInput.close();


	int nlines = 0;
	int ncolumns = 0;
	parse_file(input_file_name, nlines, ncolumns);
	int nabscissas = (ncolumns - 3)/3;

	std::ifstream fData(input_file_name, std::ios::in);
	std::ofstream fOutput_SLII(output_file_name+".SLII", std::ios::out);
	std::ofstream fOutput_nSLII(output_file_name+".nSLII", std::ios::out);
	std::ofstream fOutput_Tp(output_file_name+".Tp", std::ios::out);
	fOutput_SLII.setf(std::ios::scientific);
	fOutput_nSLII.setf(std::ios::scientific);
	fOutput_Tp.setf(std::ios::scientific);

	for (unsigned int j = 0; j < out_steps.size(); j++)
	{
		std::stringstream number; number << out_times[j] * 1e9;
		std::string label = number.str() + "(ns)";
		fOutput_SLII << std::left << std::setw(16) << std::fixed << std::setprecision(0) << label;
		fOutput_nSLII << std::left << std::setw(16) << std::fixed << std::setprecision(0) << label;
		fOutput_Tp << std::left << std::setw(16) << std::fixed << std::setprecision(0) << label;
	}
	fOutput_SLII << std::endl;
	fOutput_nSLII << std::endl;
	fOutput_Tp << std::endl;

	std::string line;
	std::getline(fData, line);	// first line is header
	const double start = std::clock();
	for (int k=1;k<nlines;k++)
	{
		double Tg;	fData >> Tg;
		double mwg;	fData >> mwg;
		double Pg;	fData >> Pg;

		std::vector<bool> presence(nabscissas);
		std::vector<double> dp;
		for (int j = 0; j < nabscissas; j++)
		{
			double dummy;
			fData >> dummy;
			if (dummy > 0.)
			{
				dp.push_back(dummy);
				presence[j] = true;
			}
			else
			{
				presence[j] = false;
			}
		}

		std::vector<double> np;
		for (int j = 0; j < nabscissas; j++)
		{
			double dummy;
			fData >> dummy;
			if (presence[j] == true)
				np.push_back(dummy);
		}

		std::vector<double> mole_fractions;
		double sum_mole_fractions = 0.;
		for (int j = 0; j < nabscissas; j++)
		{
			double dummy;
			fData >> dummy;
			if (presence[j] == true)
			{
				mole_fractions.push_back(dummy);
				sum_mole_fractions += dummy;
			}
		}

		
		if (sum_mole_fractions >= threshold_mole_fractions && Tg > threshold_temperature)
		{
			std::cout << "Point: " << k  << " over " << nlines << std::endl;

			std::vector<double> y(dp.size());
			for (unsigned int j = 0; j < dp.size(); j++)
				y[j] = np[j] * mole_fractions[j];

			OpenSMOKE::SootParticles soot;
			{
				soot.SetDensity(density);
			}
			OpenSMOKE::GasMixture gas;
			{
				gas.SetMolecularWeight(mwg);
			}

			OpenSMOKE::LIISignalModel lii(gas, soot);
			{
				lii.SetLaserFluence(J*1.e4);
				lii.SetLaserExcitationWavelength(lambda / 1.e9);
				lii.SetLaserFHMW(fwhm / 1.e9);
				lii.SetSootAbsorptionFunction(Em);
				lii.SetThermalAccomodationCoefficient(alpha);
				lii.SetEvaporationCoefficient(beta);
				lii.SetDetectionWindow(lambda_min / 1.e9, lambda_max / 1.e9);

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

			// User-defined distribution function
			OpenSMOKE::UserDefinedDistribution user_defined_pdf(dp.data(), y.data(), dp.size());
			user_defined_pdf.Print("UserDefinedDistribution.out");

			// LII Signal Simulator
			OpenSMOKE::LIISignalSimulator simulator(lii);
			{
				simulator.SetInitialSootParticlesTemperature(Tg + 0.10);
				simulator.SetGasMixtureTemperature(Tg);
				simulator.SetGasMixturePressure(Pg);
				simulator.SetFinalTime(tf);
				simulator.SetTimeStep(dt);
				simulator.SetInitialDiameter(user_defined_pdf.mu() / 1.e9);
			}

			// Simulation
			simulator.SetUserDefinedParticleSizeDistributionFunction(user_defined_pdf);
			simulator.Solve();
			simulator.Print("SolutionUserDefined.out");

			// Write
			for (unsigned int j = 0; j < out_steps.size(); j++)
			{
				fOutput_SLII << std::left << std::setw(16) << std::scientific << std::setprecision(6) << simulator.SLII()[out_steps[j]];
				fOutput_nSLII << std::left << std::setw(16) << std::scientific << std::setprecision(6) << simulator.nSLII()[out_steps[j]];
				fOutput_Tp << std::left << std::setw(16) << std::scientific << std::setprecision(6) << simulator.Tp()[out_steps[j]];
			}
		}
		else
		{
			std::cout << "Point: " << k << " over " << nlines << " (skipped)" << std::endl;

			// Write
			for (unsigned int j = 0; j < out_steps.size(); j++)
			{
				fOutput_SLII << std::left << std::setw(16) << std::scientific << std::setprecision(6) << 0.;
				fOutput_nSLII << std::left << std::setw(16) << std::scientific << std::setprecision(6) << 0.;
				fOutput_Tp << std::left << std::setw(16) << std::scientific << std::setprecision(6) << 0.;
			}
		}

		fOutput_SLII << std::endl;
		fOutput_nSLII << std::endl;
		fOutput_Tp << std::endl;

		const double total_duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		const double time_per_cell = total_duration / static_cast<double>(k);

		std::cout << "Elapsed CPU time (s): " << total_duration << std::endl;
		std::cout << "Estimated remaining time (s): " << (nlines-k)*time_per_cell << std::endl;
	}

	fOutput_SLII.close();
	fOutput_nSLII.close();
	fOutput_Tp.close();

	std::cout << "Successfully done! Press enter to continue..." << std::endl;
	getchar();
	return 0;
}

void parse_file(const std::string file_name, int& number_of_lines, int& number_of_columns)
{
	std::ifstream myfile;
	
	myfile.open(file_name, std::ios::in);

	std::string line;
	number_of_lines = 0;
	number_of_columns = 0;
	while (std::getline(myfile, line))
	{
		++number_of_lines;

		if (number_of_lines == 2)
		{
			double number;
			std::istringstream ss(line);
			while (ss >> number)
				number_of_columns++;
		}
	}

	myfile.close();

	std::cout << "Number of lines in data file: " << number_of_lines << std::endl;
	std::cout << "Number of cols  in data file: " << number_of_columns << std::endl;
}
