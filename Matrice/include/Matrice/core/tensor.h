/**************************************************************************
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
**************************************************************************/
#pragma once

#include "../private/tensor/_tensor.hpp"

DGE_MATRICE_BEGIN

template<typename _Ty, size_t _Depth = 0>
using Tensor = detail::_Tensor<_Ty, _Depth>;

//template<typename _Ty, size_t... _Shape>
//using tensor_ = detail::_Tensor_impl<_Ty, _Shape...>;

//template<typename _Ty, int... _Shape>
//using multi_matrix = detail::_Multi_matrix<_Ty, _Shape...>;
DGE_MATRICE_END
