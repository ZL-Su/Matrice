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
#include "_Lie_base.hpp"
#include "../quaternion.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
/**
 *\brief Class template for Special Orthogonal (Rotation) group.
 *\param <_Dim> dimension of SO(_Dim) group.
 */
template<typename _Ty = default_type> 
class _SO2 : public _Lie_group_base<_SO2<_Ty>> {
	using _Myt = _SO2;
	using _Mybase = _Lie_group_base<_SO2<_Ty>>;
public:

};
_DETAIL_END
DGE_MATRICE_END