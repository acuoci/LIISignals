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

#ifndef OpenSMOKE_UserDefinedDistribution_H
#define OpenSMOKE_UserDefinedDistribution_H

namespace OpenSMOKE
{
	//!  A class describing the properties of a gaseous mixture
	/*!
	This class provides the correlations for describing the properties of a gaseous mixture
	*/

	class UserDefinedDistribution
	{

	public:

		/**
		*@brief Default constructor: Sets a user defined distribution function
		*@param	x abscissas
		*@param y ordinates (not necessarily normalized)
		*@param	n number of points
		*/
		UserDefinedDistribution(const double* x, const double* y, const unsigned int n);

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
		*@brief Returns the probabilities at the quadrature points
		*@return		the probabilities at the quadrature points
		*/
		const std::vector<double>& p() const { return p_; };

		/**
		*@brief Returns the mean value
		*@return		the mean value
		*/
		double mu() const { return mu_; };

		/**
		*@brief Prints the pdf and the probabilities at the quadrature points
		*/
		void Print(const std::string filename) const;


	private:

		double mu_;					//!< mean value
		double sigmap2_;			//!< variance

		unsigned int n_;			//!< number of quadrature points

		std::vector<double> x_;		//!< positions of quadrature points
		std::vector<double> p_;		//!< probabilities at quadrature points
		std::vector<double> y_;		//!< un-normalized probabilities
	};

}

#include "UserDefinedDistribution.hpp"

#endif /* OpenSMOKE_UserDefinedDistribution_H */
