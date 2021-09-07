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
#pragma once

#include "_type_traits.h"

#ifdef MATRICE_ENABLE_CXX20
DGE_MATRICE_BEGIN
template<typename _Ty>
concept real_type = is_integral_v<_Ty> || is_floating_point_v<_Ty>;
template<typename _My>
concept matrix_type = is_matrix_v<_My>;
DGE_MATRICE_END
#endif