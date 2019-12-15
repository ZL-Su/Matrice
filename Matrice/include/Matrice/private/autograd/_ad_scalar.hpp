/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

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
#include "_ad_ops.h"

DGE_MATRICE_BEGIN
_INTERNAL_BEGIN
template<typename _Dty, bool _Enable> struct _Auto_diff_spec_op;
_INTERNAL_END
_DETAIL_BEGIN
/**
 * \class _Auto_diff_scalar<>
 * \brief A scalar type with auto-differentiation capability
 * \param <_Dty> the matrix type used to hold the derivatives. 
				The value type and the number of derivatives to compute
				are determined from this type. _Dty can be: 
				\c Matrix_<f,3,1> for 3 derivatives, or 
				\c Matrix_<f,::dynamic,1> if the num of ders is deduced
				at runtime.
 * This class represents a scalar value while tracking its derivaties using Matrice's expression.
 * It supports almost all primitive math functions.
 * It also can be used as the scalar type of a Matrix_<> object, in this case, the expression mechanism only occurs at the top Matrix_ level, while the derivatives are computed right away.
 */
template<typename _Dty> 
class _Auto_diff_scalar
	: public internal::_Auto_diff_spec_op<_Dty, is_matrix_v<remove_all_t<_Dty>>>
{
	using _Mybase = internal::_Auto_diff_spec_op<_Dty, is_matrix_v<remove_all_t<_Dty>>>;
public:
	using dervt_type = remove_all_t<_Dty>;
	using value_type = typename dervt_type::value_type;
};
_DETAIL_END

DGE_MATRICE_END