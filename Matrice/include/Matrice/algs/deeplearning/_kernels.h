/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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
***********************************************************************/
#pragma once
#include "../../core/matrix.h"

MATRICE_NAMESPACE_BEGIN_
namespace dl { namespace kernel {
// \alias template for dynamic vector type
template<typename _Ty> using vector_type = Matrix<_Ty>;

template<typename _Ty>
MATRICE_HOST_FINL _Ty _Gaussian(const vector_type<_Ty>& x, const vector_type<_Ty>& x_pri, const _Ty _sigma) {
	vector_type<_Ty> _Diff = x - x_pri;
	auto _Prod = (_Diff*_Diff).sum();
	return (std::exp(-0.5 *_Prod / (_sigma * _sigma)));
}


} }
_MATRICE_NAMESPACE_END