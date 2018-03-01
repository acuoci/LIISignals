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
	RungeKutta4thOrder::RungeKutta4thOrder()
	{

	}

	void RungeKutta4thOrder::SetInitialConditions(const unsigned int n, const double tInitial, const double* uInitial)
	{
		tInitial_ = tInitial;
		n_ = n;

		MemoryAllocation();

		uInitial_.resize(n_);
		for (unsigned int i = 0; i < n_; i++)
			uInitial_[i] = uInitial[i];
	}

	void RungeKutta4thOrder::SetTimeStep(const double dt)
	{
		dt_ = dt;
	}

	void RungeKutta4thOrder::SetFinalTime(const double tFinal)
	{
		tFinal_ = tFinal;
	}

	void RungeKutta4thOrder::MemoryAllocation()
	{
		u0_.resize(n_);
		u1_.resize(n_);
		u2_.resize(n_);
		u3_.resize(n_);

		f0_.resize(n_);
		f1_.resize(n_);
		f2_.resize(n_);
		f3_.resize(n_);

		solution_.resize(n_ + 1);
	}

	void RungeKutta4thOrder::Solve()
	{
		// Number of steps
		const unsigned int nsteps = (unsigned int)((tFinal_ - tInitial_) / dt_);

		// Resize solution matrix
		for (unsigned int i = 0; i < n_ + 1; i++)
			solution_[i].resize(nsteps + 1);

		// Initial conditions
		u0_ = uInitial_;
		solution_[0][0] = tInitial_;
		for (unsigned int i = 0; i < n_; i++)
			solution_[i + 1][0] = u0_[i];

		// Loop
		for (unsigned int k = 0; k < nsteps; k++)
		{
			// Current time
			const double t0 = tInitial_ + k * dt_;

			// Evalution at t0
			f_(t0, u0_.data(), f0_.data());

			const double t1 = t0 + dt_ / 2.;
			for (unsigned int i = 0; i < n_; i++)
				u1_[i] = u0_[i] + dt_ * f0_[i] / 2.;
			f_(t1, u1_.data(), f1_.data());

			const double t2 = t0 + dt_ / 2.;
			for (unsigned int i = 0; i < n_; i++)
				u2_[i] = u0_[i] + dt_ * f1_[i] / 2.;
			f_(t2, u2_.data(), f2_.data());

			const double t3 = t0 + dt_;
			for (unsigned int i = 0; i < n_; i++)
				u3_[i] = u0_[i] + dt_ * f2_[i];
			f_(t3, u3_.data(), f3_.data());

			for (unsigned int i = 0; i < n_; i++)
				u0_[i] += dt_ * (f0_[i] + 2.*f1_[i] + 2.*f2_[i] + f3_[i]) / 6.;

			solution_[0][k + 1] = t3;
			for (unsigned int i = 0; i < n_; i++)
				solution_[i + 1][k + 1] = u0_[i];
		}
	}
}
