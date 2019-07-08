/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
#include "../functions.h"

DGE_MATRICE_BEGIN
namespace dnn {
_DETAIL_BEGIN
template<typename _Ty>
void _conv2d_impl(_TAG device_tag::cpu, const Matrix<_Ty>& _In, Matrix<_Ty>& _Out) {

}

template<typename _Ty>
void _conv2d_impl(_TAG device_tag::gpu, const Matrix<_Ty>& _In, Matrix<_Ty>& _Out) {

}
_DETAIL_END

template<typename _Tag, typename _Ty>
MATRICE_HOST_INL void _conv2d(const Matrix<_Ty>& _In, Matrix<_Ty>& _Out) {
	static_assert(is_floating_point_v<_Ty>, "Only float or double type is allowed.");
	detail::_conv2d_impl(_Tag(), _In, _Out);
}
}
DGE_MATRICE_END