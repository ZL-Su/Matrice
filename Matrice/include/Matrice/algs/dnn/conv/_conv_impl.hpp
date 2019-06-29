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
template<typename _Ty, typename = enable_if_t<is_floating_point_v<_Ty>>>
void _conv2d_impl(const Matrix<_Ty>& _In, Matrix<_Ty>& _Out, _TAG device_tag::cpu) {

}

template<typename _Ty, typename = enable_if_t<is_floating_point_v<_Ty>>>
void _conv2d_impl(const Matrix<_Ty>& _In, Matrix<_Ty>& _Out, _TAG device_tag::gpu) {

}
_DETAIL_END

template<typename _Tag, typename _Ty>
void _conv2d() {
	detail::_conv2d_impl(Matrix<_Ty>(), Matrix<_Ty>(), _Tag());
}
}
DGE_MATRICE_END