/**************************************************************
	This example shows how to refine image patch matching
	using Matrice correlation module.
 **************************************************************/
#pragma once
#include <io/io.hpp>
#include <algs/correlation.hpp>
#include <algs/graphics.hpp>
#include <algs/geometry.hpp>

namespace fs = dgelom::fs;
using corr_optim_t = dgelom::correlation_optimizer::icgn_bic_1<float>;
using raw_image_t = corr_optim_t::matrix_type;
using smooth_image_t = corr_optim_t::smooth_image_t;

int corr_optimer_eng() try 
{
	const auto path = dgelom::fs::current_path().append("noise_10");
	auto image_loader = dgelom::make_loader(path, dgelom::io::tiff<raw_image_t::value_t>());

	const auto ref_image = image_loader.forward().front();
	const auto rows = ref_image.rows(), cols = ref_image.cols();

	const auto options = corr_optim_t::options_type{ 15 };
	const auto x_space = dgelom::linspace<int>::_(options._Radius * 2, cols - options._Radius * 2, 40);
	const auto y_space = dgelom::linspace<int>::_(options._Radius * 2, rows - options._Radius * 2, 40);

	auto disp_x = corr_optim_t::matrix_type(x_space.size() * y_space.size(), image_loader.depth(), 0.);

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
					<< "----QUAD LOSS---" << "----RHO LOSS----" << "---||dp||----" << "ID: " << npoint << "\n";
				for (auto nit = 0; nit < options._Maxits; ++nit) {
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
					if (nodp < 1.0e-6)
						break;
					if (squared_loss < 1.e-2)
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
	std::cout << e.what() << "in function: " << e.location()._func
		<< ". \n    See line " << e.location()._line
		<< " in file: " << e.location()._file << std::endl;
}
catch (std::exception& e) {
	std::cout << "Unknown exception with info: [" << e.what() << "]" << std::endl;
}
catch (...) {
	std::cout << "Unknown exception." << std::endl;
}