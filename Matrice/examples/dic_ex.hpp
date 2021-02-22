/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#pragma once
#include <io/io.hpp>
#include <algs/correlation.hpp>
#include <algs/graphics.hpp>
#include <algs/geometry.hpp>

DGE_MATRICE_BEGIN
namespace example {

namespace fs = dgelom::fs;
using corr_optim_t = dgelom::correlation_optimizer::icgn_bic_1<float>;
using raw_image_t = corr_optim_t::matrix_type;
using smooth_image_t = corr_optim_t::smooth_image_t;

/// <summary>
/// \brief Example of correlation optimizer for image matching.
/// </summary>
/// <param name="'apath'">Path to current location of .</param>
/// <param name="'dfolder'">Folder where the images are stored.</param>
/// <returns>Make no sense.</returns>
int corr_optimer_eng(fs::path&& apath, std::string&& dfolder) try
{
	///\brief Attach image file info to a data loader (tiff herein).
	decltype(auto) path = apath.append(dfolder);
	auto image_loader = dgelom::make_loader(path,
		dgelom::io::tiff<raw_image_t::value_t>());

	///\brief Get reference image info.
	const auto ref_image = image_loader.forward().front();
	const auto rows = ref_image.rows(), cols = ref_image.cols();

	///\brief Set optimizer options and mesh computing domain.
	const auto options = corr_optim_t::options_type{ 15 };
	const auto x_space = dgelom::linspace<int>::_(
		options.radius() << 1, cols - options.radius() << 1, 40);
	const auto y_space = dgelom::linspace<int>::_(
		options.radius() << 1, rows - options.radius() << 1, 40);

	auto disp_x = corr_optim_t::matrix_type(x_space.size() * y_space.size(),
		image_loader.depth(), 0.);

	smooth_image_t f{ ref_image };
	for (; !image_loader.end(); ) {
		const auto cur_image = image_loader.forward().front();
		decltype(f) g{ cur_image };

		auto disp_col = disp_x.cview(image_loader.pos() - 1);
		corr_optim_t solver(f, g, options);
		size_t npoint = 0;
		for (const auto y : y_space) {
			for (const auto x : x_space) {
				solver.init(x, y);
				auto p = decltype(solver)::param_type();
				std::cout << "-IT-" << "-----DISP----"
					<< "----QUAD LOSS---" << "----RHO LOSS----"
					<< "---||dp||----" << "ID: " << npoint << "\n";
				for (auto nit = 0; nit < options.maxiters(); ++nit) {
					auto [squared_loss, rho_loss, nodp] = solver.robust_sol(p);
					std::cout << " "
						<< std::setiosflags(std::ios::left)
						<< dgelom::str(nit, 2).append("  ")
						<< std::setprecision(7)
						<< std::setw(15)
						<< p[0]
						<< std::setw(15)
						<< squared_loss
						<< std::setw(15)
						<< rho_loss
						<< std::setw(15)
						<< nodp << "\n";
					if (nodp < options.tol(1))
						break;
					if (squared_loss < options.tol(100))
						break;
				}
				disp_col[npoint] = p[0];
				++npoint;
			}
		}
	}
	dgelom::IO::CSV csv(path.parent_path().string() + "\\res_noise_10\\Welsch_S0.01.csv");
	csv.open(dgelom::IO::app);
	for (auto row = disp_x.rwbegin(); row != disp_x.rwend(); ++row) {
		csv.append(row.begin(), row.end());
	}
	csv.close();

	std::cout << " >> [Matrice message] finish!\n";

	return 0;
}
catch (dgelom::exception::error& e) {
	std::cout << "Err: " << e.what()
		<< "in function: " << e.location()._func
		<< ". (See line " << e.location()._line
		<< " in file '" << e.location()._file << "')\n";
}
catch (std::exception& e) {
	std::cout << "Unknown exception with info: [" << e.what() << "]" << std::endl;
}
catch (...) {
	std::cout << "Unknown exception." << std::endl;
}
}
DGE_MATRICE_END