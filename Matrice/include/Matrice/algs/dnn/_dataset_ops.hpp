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

#include "forward.hpp"
#include "util/_macros.h"

MATRICE_ALG_BEGIN(dnn)
template<typename _Ty,  typename _Uy = _Ty> MATRICE_HOST_INL 
auto data_iter(size_t nbatch, const Matrix<_Ty>& feats, const Matrix<_Uy>& labels) {
	const auto num_examples = min(feats.size(), labels.size());

}
MATRICE_ALG_END(dnn)