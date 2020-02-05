/**************************************************************************
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
**************************************************************************/
#pragma once
#include "_Lie_group_base.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
/**
 *\brief Class template for Special Orthogonal (Rotation) group.
 *\param <_Dim> dimension of SO(_Dim) group.
 */
template<size_t _Dim, typename _Ty = default_type> 
class _SO : public _Lie_group_base<_SO<_Dim, _Ty>> {
	using _Myt = _SO;
	using _Mybase = _Lie_group_base<_Myt>;
public:

};
_DETAIL_END
DGE_MATRICE_END