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
	const double LogNormalDistribution::pi_ = 3.14159265359;	// [-]

	LogNormalDistribution::LogNormalDistribution(const double cmd, const double sigma, const unsigned int n)
	{
		n_ = n;
		cmd_ = cmd;
		sigma_ = sigma;
		ln_cmd_ = std::log(cmd_);
		ln_sigma_ = std::log(sigma_);

		// Internal constants
		c1_ = 1. / std::sqrt(2.*pi_) / ln_sigma_;
		c2_ = 2.*ln_sigma_*ln_sigma_;

		// Maximum value
		max_x_ = std::exp(-ln_sigma_ * ln_sigma_ + ln_cmd_);
		max_ = (*this)(max_x_);

		// Boundaries for numerical integration
		const double x_tilde = 1e-6;
		const double ln_sigma_2_ = ln_sigma_ * ln_sigma_;
		max_bc_ = std::exp(ln_sigma_2_ + ln_cmd_ +
			std::sqrt(-2.*ln_sigma_2_*ln_cmd_ + 2.*std::log(cmd_ / x_tilde)*ln_sigma_2_));
		min_bc_ = std::exp(ln_sigma_2_ + ln_cmd_ -
			std::sqrt(-2.*ln_sigma_2_*ln_cmd_ + 2.*std::log(cmd_ / x_tilde)*ln_sigma_2_));

		// Moments
		mu_ = std::exp(ln_cmd_ + ln_sigma_2_ / 2.);
		sigmap2_ = (std::exp(ln_sigma_2_) - 1.)*std::exp(2.*ln_cmd_ + ln_sigma_2_);
		gamma1_ = std::sqrt((std::exp(ln_sigma_2_) - 1.))*(2. + std::exp(ln_sigma_2_));
		gamma2_ = std::exp(4.*ln_sigma_2_) + 2.*std::exp(3.*ln_sigma_2_) + 3.*std::exp(2.*ln_sigma_2_) - 6.;

		// Summary on the screen
		std::cout << "--------------------------------------------" << std::endl;
		std::cout << "Log-Normal Distribution function: summary   " << std::endl;
		std::cout << "--------------------------------------------" << std::endl;
		std::cout << " * CMD:              " << cmd_ << std::endl;
		std::cout << " * sigma:            " << sigma_ << std::endl;
		std::cout << " * ln(CMD):          " << ln_cmd_ << std::endl;
		std::cout << " * ln(sigma):        " << ln_sigma_ << std::endl;
		std::cout << " * xmin:             " << min_bc_ << std::endl;
		std::cout << " * xmax:             " << max_bc_ << std::endl;
		std::cout << " * mean dp:          " << mu_ << std::endl;
		std::cout << " * std.dev. sigmap:  " << std::sqrt(sigmap2_) << std::endl;
		std::cout << " * variance sigma2p: " << sigmap2_ << std::endl;
		std::cout << " * skewness gamma1:  " << gamma1_ << std::endl;
		std::cout << " * kurtosis gamma2:  " << gamma2_ << std::endl;
		std::cout << "--------------------------------------------" << std::endl;
		std::cout << std::endl;

		// Calculated Pdf and probabilities
		CalculatedPdf();
	}

	double LogNormalDistribution::operator() (const double x) const
	{
		return c1_ / x * std::exp(-std::pow(std::log(x) - ln_cmd_, 2.) / c2_);
	}

	void LogNormalDistribution::SetNumberofIntervals(const unsigned int n)
	{
		n_ = n;
		CalculatedPdf();
	}

	void LogNormalDistribution::CalculatedPdf()
	{
		x_.resize(n_);
		pdf_.resize(n_);
		p_.resize(n_);

		// Calculate Pdf (uniform distribution)
		double cdf = 0.;
		{
			const double dx = (max_bc_ - min_bc_) / double(n_);

			for (unsigned int i = 0; i < n_; i++)
			{
				x_[i] = min_bc_ + dx / 2. + i * dx;
				pdf_[i] = (*this)(x_[i]);
				cdf += dx * pdf_[i];
			}

			for (unsigned int i = 0; i < n_; i++)
				p_[i] = pdf_[i] * dx / cdf;
		}

		// Check accuracy: sum
		{
			const double e = (cdf - 1.) * 100.;

			std::cout << "---------------------------------------------" << std::endl;
			std::cout << "Log-Normal Distribution function: accuracy   " << std::endl;
			std::cout << "---------------------------------------------" << std::endl;
			std::cout << " * cdf (analytical):  " << 1. << std::endl;
			std::cout << " * cdf (numerical):   " << cdf << std::endl;
			std::cout << " * error(%):          " << e << std::endl;
			std::cout << "---------------------------------------------" << std::endl;

			if (std::fabs(e) > 1.)
			{
				std::cout << "The numerical accuracy in describing the log-normal distribution is too low." << std::endl;
				std::cout << "Please increase the number of intervals for numerical integration." << std::endl;
			}

			std::cout << std::endl;
		}

		// Check accuracy: mean
		{
			double mu = 0.;
			for (unsigned int i = 0; i < n_; i++)
				mu += p_[i] * x_[i];

			const double e = (mu - mu_) / mu_ * 100.;

			std::cout << "---------------------------------------------" << std::endl;
			std::cout << "Log-Normal Distribution function: accuracy   " << std::endl;
			std::cout << "---------------------------------------------" << std::endl;
			std::cout << " * mean (analytical): " << mu_ << std::endl;
			std::cout << " * mean (numerical):  " << mu << std::endl;
			std::cout << " * error(%):          " << e << std::endl;
			std::cout << "---------------------------------------------" << std::endl;

			if (std::fabs(e) > 1.)
			{
				std::cout << "The numerical accuracy in describing the log-normal distribution is too low." << std::endl;
				std::cout << "Please increase the number of intervals for numerical integration." << std::endl;
			}

			std::cout << std::endl;
		}

		// Check accuracy: variance
		{
			double sigmap2 = 0.;
			for (unsigned int i = 0; i < n_; i++)
				sigmap2 += p_[i] * x_[i] * x_[i];
			sigmap2 -= mu_ * mu_;

			const double e = (sigmap2 - sigmap2_) / sigmap2_ * 100.;

			std::cout << "---------------------------------------------" << std::endl;
			std::cout << "Log-Normal Distribution function: accuracy   " << std::endl;
			std::cout << "---------------------------------------------" << std::endl;
			std::cout << " * variance (analytical): " << sigmap2_ << std::endl;
			std::cout << " * variance (numerical):  " << sigmap2 << std::endl;
			std::cout << " * error(%):              " << e << std::endl;
			std::cout << "---------------------------------------------" << std::endl;

			if (std::fabs(e) > 1.)
			{
				std::cout << "The numerical accuracy in describing the log-normal distribution is too low." << std::endl;
				std::cout << "Please increase the number of intervals for numerical integration." << std::endl;
			}

			std::cout << std::endl;
		}
	}

	void LogNormalDistribution::Print(const std::string filename) const
	{
		std::ofstream fOut(filename.c_str(), std::ios::out);
		fOut.setf(std::ios::scientific);

		fOut << std::left << std::setw(15) << "x(udm)";
		fOut << std::left << std::setw(15) << "pdf(1/udm)";
		fOut << std::left << std::setw(15) << "p";
		fOut << std::endl;

		for (unsigned int k = 0; k < n_; k++)
		{
			fOut << std::left << std::setw(15) << x_[k];
			fOut << std::left << std::setw(15) << pdf_[k];
			fOut << std::left << std::setw(15) << p_[k];
			fOut << std::endl;
		}

		fOut.close();
	}
}
