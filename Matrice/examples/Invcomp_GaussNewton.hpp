/**************************************************************
	This example shows how to refine image patch matching
	using Matrice correlation module.
 **************************************************************/
#pragma once
#include <iostream>
#include "../include/Matrice/io/io.hpp"
#include "../include/Matrice/algs/correlation/_optim.h"

int main() try {
	// image path "...\\images"
	const std::string _Impath = "...\\images";

	// create image loader, which loads images from path: _Impath
	dgelom::io::data_loader_uint8 _Loader(
		dgelom::io::directory{ _Impath, {} },
		[](auto&& _)->auto{
		auto _Img = cv::imread(_, cv::IMREAD_GRAYSCALE);
		return dgelom::move(data);
	});

	// using 1-st order IC-GN with bicubic spline interpolation
	using optimzer_t = dgelom::correlation_optimizer::icgn_bic<float, 1>;
	using options_t = dgelom::correlation_optimizer::options;
	using frame_t = optimzer_t::interp_type;

	// load reference image
	const auto _Ref = _Loader.forward().front();
	frame_t f(_Ref);

	// load current image
	const auto _Cur = _Loader.forward().front();
	frame_t g(_Cur);

	// create optimizer with 31-by-31 window
	optimzer_t optim(f, g, options_t{ 15 });

	// given any initial estimation {x,y} = {50, 50}
	optimzer_t::point_type p{ 50, 50 };
	// initialize optimizer
	optim.init(p);
	// set initial parameter vector to zero
	optimzer_t::param_type par{ 0 };
	// solving
	for (const auto _Its : dgelom::range(1, optim.options()._Maxits)) {
		auto[_Coef, _Tol] = optim(par);
		if (_Tol < options_t::_Mytol<decltype(_Tol)>){
			break;
		}
	}
	
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