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
#include <valarray>
#include <private/autograd/_ad_exps.h>

DGE_MATRICE_BEGIN
namespace example {
void autodiff_test() {
	// Make a variable and a function expression...
	auto x = ade::make_variable(0.);
	const auto f = 2.*x + ade::sin(x);

	// Eval function value and derivative...
	auto value = f();
	auto df = f.deriv();

	// Update variable value.
	x = 3.1415926/2;

	// Eval function value and derivative...
	value = f();
	df = f.deriv();
}
}
DGE_MATRICE_END