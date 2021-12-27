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
#include <array>

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
/// <param name="'apath'">Path to current location.</param>
/// <param name="'dfolder'">Folder where the images are stored.</param>
/// <returns>Make no sense.</returns>
int corr_optimer_eng(fs::path&& apath, auto dfolder, auto val)
try {
	///\brief Attach image file info to a data loader (tiff herein).
	const auto path = apath.append(dfolder);
	auto image_loader = dgelom::make_loader(path,
		dgelom::io::tiff<raw_image_t::value_t>());

	///\brief Get reference image info.
	const auto ref_image = image_loader.forward().front();
	const auto rows = ref_image.rows(), cols = ref_image.cols();
	smooth_image_t f{ ref_image };

	///\brief Set optimizer options and mesh computing domain.
	const auto options = corr_optim_t::options_type{ 20 };
	const auto subsize = options.radius() << 1;
	/*const */auto x_space = dgelom::make_linspace(
		subsize, cols - subsize, 2);
	/*const */auto y_space = dgelom::make_linspace(
		subsize, rows - subsize, 2);
	x_space(0) = 933, x_space(1) = 1523;
	y_space(0) = 915, y_space(1) = 921;
	auto disp_x = corr_optim_t::matrix_type(x_space.size()*y_space.size()*2, image_loader.depth(), 0.);
	auto error = dgelom::Matrix_<float, ::dynamic, 2>(image_loader.depth());

	///\brief Move loader back one step to start from reference.
	image_loader.shift(-1);

	auto loss_scale = val;
	for (; !image_loader.end(); ) {
		const auto cur_image = image_loader.forward().front();
		decltype(f) g{ cur_image };

		const auto idx = image_loader.pos() - 1;
		auto disp_col = disp_x.cview(idx);
		corr_optim_t solver(f, g, options);
		solver.set_loss_scale(loss_scale);

		size_t npoint = 0;
		for (const auto y : y_space) {
			for (const auto x : x_space) {
				solver.init(x, y);
				auto p = decltype(solver)::param_type();
				std::cout << "-IT-" << "-----DISP----"
					<< "----QUAD LOSS---" << "----RHO LOSS----"
					<< "---||dp||----" << "ID: " << npoint << "\n";
				for (auto nit = 0; nit < min(10,options.maxiters()); ++nit) {
					auto [squared_loss, rho_loss, nodp] = solver.robust_sol(p);
					//auto [squared_loss, nodp] = solver(p);
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
					if (nodp < options.tol(10))
						break;
					if (squared_loss < options.tol(100))
						break;
				}
				disp_col[(npoint)<<1] = p[0];
				disp_col[(npoint)<<1|1] = p[3];
				++npoint;
			}
		}
		
		//const auto u = disp_x.cbegin(idx);
		//const auto [mean, stdv] = corr::eval_perf(u.begin(), u.end(), /*Groundtruth=*/idx*0.05);
		//error.view<0>(idx) = { float(mean), float(stdv) };
	}

	dgelom::IO::CSV csv;
	std::array<std::string, 2> labels{ "Mean", "Stdv" };

	//csv.reset(path.parent_path().append("res\\error_uniform.csv"));
	//csv.open(dgelom::IO::app);
	//csv.append("Impulse-Welsch-"+dgelom::str(loss_scale));
	////csv.append("Impulse-Gaussian");
	//for (auto col = error.cwbegin(); col!=error.cwend(); ++col) {
	//	csv.append(std::move(labels[col.pos()]), 
	//		col.begin(), col.end());
	//}
	//csv.close();

	csv.reset(path.parent_path().append("res\\disp_4pts_0.0005.csv"));
	csv.open(dgelom::IO::app);
	for (auto row = disp_x.rwbegin(); row != disp_x.rwend(); ++row) {
		csv.append(row.begin(), row.end());
	}
	csv.close();

	std::cout << " >> [Matrice message] finish for scale " + dgelom::str(loss_scale) + "\n";

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

void multi_points_optimizer(std::vector<std::array<corr_optim_t::value_type, 4>>&& nodes, fs::path&& path)
try {
	using value_t = corr_optim_t::value_type;
	const auto refimg = io::imread(path.string() + "\\tension_00.tif").matrix<value_t>();
	const auto curimg = io::imread(path.string() + "\\tension_01.tif").matrix<value_t>();

	const auto opts = corr_optim_t::options_type{ 20 };
	auto disp = corr_optim_t::matrix_type(nodes.size(), 2, 0.);

	smooth_image_t f(refimg), g(curimg);
	auto optim = corr_optim_t{f, g, opts};
	auto rview = disp.rwbegin();
	for (decltype(auto) node : nodes) {
		optim.init(node[0], node[1]);
		auto p = decltype(optim)::param_type{};
		p[0] = node[2] - node[0];
		p[3] = node[3] - node[1];
		for (auto i = 0; i < min(10, opts.maxiters()); ++i) {
			const auto [coef, error] = optim(p);
			if (error < opts.tol()) {
				break;
			}
		}
		rview = { p[0], p[3] };
		++rview;
	}
}
catch (dgelom::exception::error& e) {
	std::cout << "Err: " << e.what()
		<< "in function: " << e.location()._func
		<< ". (See line " << e.location()._line
		<< " in file '" << e.location()._file << "')\n";
}
}
DGE_MATRICE_END