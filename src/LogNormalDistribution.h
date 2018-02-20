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

#ifndef OpenSMOKE_LogNormalDistribution_H
#define OpenSMOKE_LogNormalDistribution_H

namespace OpenSMOKE
{
	//!  A class describing the properties of a gaseous mixture
	/*!
	This class provides the correlations for describing the properties of a gaseous mixture
	*/

	class LogNormalDistribution
	{
	
	public:

		/**
		*@brief Default constructor
		*@param	cmd		median
		*@param	sigma	geometric width
		*/
		LogNormalDistribution(const double cmd, const double sigma, const unsigned int n);

		/**
		*@brief Returns the probability distribution value at point x
		*@param		x	point at which the pdf is required
		*@return		the pdf
		*/
		double operator() (const double x) const;

		/**
		*@brief Returns the maximum value of the pdf
		*@return		the maximum value
		*/
		double max() const { return max_; }

		/**
		*@brief Returns the coordinate of maximum value of the pdf
		*@return		the coordinate of maximum value
		*/
		double max_x() const { return max_x_; }

		/**
		*@brief Returns the minimum coordinate from which the pdf is larger than a minimum threshold
		*@return		the minimum coordinate
		*/
		double min_bc() const { return min_bc_; }

		/**
		*@brief Returns the maximum coordinate up to which the pdf is larger than a minimum threshold
		*@return		the maximum coordinate
		*/
		double max_bc() const { return max_bc_; }

		/**
		*@brief Returns number of discretization intervals
		*@return		the number of discretization intervals
		*/
		unsigned int n() const { return n_; }

		/**
		*@brief Returns the coordinates of quadrature points
		*@return		the coordinates of quadrature points
		*/
		const std::vector<double>& x() const { return x_; }

		/**
		*@brief Returns the pdf values at the quadrature points
		*@return		the pdf values at the quadrature points
		*/
		const std::vector<double>& pdf() const { return pdf_; }

		/**
		*@brief Returns the probabilities at the quadrature points
		*@return		the probabilities at the quadrature points
		*/
		const std::vector<double>& p() const { return p_; };

		/**
		*@brief Sets the number of intervals (i.e. quadrature points) for numerical integration
		*@param	n number of quadrature points
		*/
		void SetNumberofIntervals(const unsigned int n);

		/**
		*@brief Prints the pdf and the probabilities at the quadrature points
		*/
		void Print(const std::string filename) const;

	private:

		/**
		*@brief Calculates the pdfs and the probabilities at the quadrature points
		*/
		void CalculatedPdf();

	private:

		double sigma_;				//!< geometric width
		double cmd_;				//!< median
		double ln_sigma_;			//!< logarithm of geometric width
		double ln_cmd_;				//!< logarithm of median

		double max_;				//!< max value of pdf
		double max_x_;				//!< position of maximum value of pdf
		double min_bc_;				//!< minimum value of interval along independent coordinat for numerical integration
		double max_bc_;				//!< maximum value of interval along independent coordinat for numerical integration

		double c1_;					//!< internal coefficient
		double c2_;					//!< internal coefficient

		double mu_;					//!< mean value
		double sigmap2_;			//!< variance
		double gamma1_;				//!< skewness
		double gamma2_;				//!< kurtosis
			
		unsigned int n_;			//!< number of quadrature points

		std::vector<double> x_;		//!< positions of quadrature points
		std::vector<double> p_;		//!< probabilities at quadrature points
		std::vector<double> pdf_;	//!< pdfs at quadrature points

	private:

		static const double pi_;	//!< pi

	};

}

#include "LogNormalDistribution.hpp"

#endif /* OpenSMOKE_LogNormalDistribution_H */
