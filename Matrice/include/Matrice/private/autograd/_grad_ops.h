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
#include "_grad_utils.h"

DGE_MATRICE_BEGIN _DETAIL_BEGIN
struct _Grad_tag {};
struct _Grad_op_tag : _Grad_tag {};
struct _Constant_tag : _Grad_op_tag {};
struct _Constant_tag : _Grad_op_tag {};
struct _Linear_tag : _Grad_op_tag {};
struct _Square_tag : _Grad_op_tag {};

namespace grad {
template<typename _Tag> struct _Op{};

template<> struct _Op<_Constant_tag> {
	template<typename _Ty, typename _Gradty>
	MATRICE_GLOBAL_INL static auto& eval(const _Ty& _Val, _Gradty& _Grad) {
		return (_Grad = 0*_Grad);
	}
};
template<> struct _Op<_Linear_tag> {
	template<typename _Ty, typename _Gradty>
	MATRICE_GLOBAL_INL static auto& eval(const _Ty& _Coef, _Gradty& _Grad) {
		return (_Grad = _Coef*_Grad);
	}
};
template<> struct _Op<_Square_tag> {
	template<typename _Ty, typename _Gradty>
	MATRICE_GLOBAL_INL static auto& eval(const _Ty& _Val, _Gradty& _Grad) {
		return (_Grad = 2 * _Val * _Grad);
	}
};
}
_DETAIL_END DGE_MATRICE_END
