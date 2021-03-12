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
*********************************************************************/
#include "core/matrix.h"
#include "algs/spline/_spline.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ty>
_Bspline<_Ty, bicubic_tag>::_Bspline(const matrix_type& _Data) {
	_Mybase::_Mycoef = make_shared_matrix(matrix_type(_Data.shape()));
	_Precompute(_Data);
}
template<typename _Ty>
void _Bspline<_Ty, bicubic_tag>::_Precompute(const matrix_type& _Data) {
	matrix_type _Buff(_Data.shape());
}

template class _Bspline<float, bicubic_tag>;
template class _Bspline<double, bicubic_tag>;
_DETAIL_END
DGE_MATRICE_END