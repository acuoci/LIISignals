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

#include <algorithm>

namespace OpenSMOKE
{

	UserDefinedDistribution::UserDefinedDistribution(const double* x, const double* y, const unsigned int n)
	{
		// Memory allocation
		n_ = n;
		x_.resize(n_);
		y_.resize(n_);
		p_.resize(n_);

		// Fill the vectors
		double sum = 1.e-32;
		for (unsigned int i = 0; i < n_; i++)
		{
			x_[i] = x[i];
			y_[i] = y[i];
			sum += y_[i];
		}

		// Probabilities
		for (unsigned int i = 0; i < n_; i++)
			p_[i] = y_[i] / sum;

		// Mean
		mu_ = 0.;
		for (unsigned int i = 0; i < n_; i++)
			mu_ += p_[i] * x_[i];

		// Variance
		sigmap2_ = 0.;
		for (unsigned int i = 0; i < n_; i++)
			sigmap2_ += p_[i] * x_[i] * x_[i];
		sigmap2_ -= mu_ * mu_;

		// Summary on the screen
		std::cout << "-----------------------------------------------" << std::endl;
		std::cout << "User defined Distribution function: summary    " << std::endl;
		std::cout << "-----------------------------------------------" << std::endl;
		std::cout << " * mean:          " << mu_ << std::endl;
		std::cout << " * std deviation: " << std::sqrt(sigmap2_) << std::endl;
		std::cout << " * min element:   " << *std::min_element(x_.begin(), x_.end()) << std::endl;
		std::cout << " * max element:   " << *std::max_element(x_.begin(), x_.end()) << std::endl;
		std::cout << "-----------------------------------------------" << std::endl;
	}

	void UserDefinedDistribution::Print(const std::string filename) const
	{
		std::ofstream fOut(filename.c_str(), std::ios::out);
		fOut.setf(std::ios::scientific);

		fOut << std::left << std::setw(15) << "x(udm)";
		fOut << std::left << std::setw(15) << "y(a.u.)";
		fOut << std::left << std::setw(15) << "p";
		fOut << std::endl;

		for (unsigned int k = 0; k < n_; k++)
		{
			fOut << std::left << std::setw(15) << x_[k];
			fOut << std::left << std::setw(15) << y_[k];
			fOut << std::left << std::setw(15) << p_[k];
			fOut << std::endl;
		}

		fOut.close();
	}
}
