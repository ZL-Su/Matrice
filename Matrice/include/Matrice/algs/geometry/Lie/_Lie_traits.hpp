/*  ************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
*	************************************************************************/
#pragma once

#include "../private/_type_traits.h"
#include "forward.hpp"
#include "_Lie_fwd.hpp"

DGE_MATRICE_BEGIN
template<size_t _Dim, typename _Ty>
struct traits<detail::_SO<_Dim, _Ty>> {
	using value_type = _Ty;
	using group_type = detail::_SO<_Dim, _Ty>;
	using properties = internal::_Lie_group_prop<group_type>;
	using vector_type = detail::Matrix_<value_type, _Dim, 1>;
	using jacobian_type = detail::Matrix_<value_type, _Dim, _Dim>;
};
DGE_MATRICE_END