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
#define MATRICE_MAKE_LIE_GROUP_BEGIN(NAME) \
template<typename _Ty = default_type> \
class NAME : public _Lie_group_base<NAME<_Ty>> { \
    using _Myt = NAME; \
    using _Mybase = _Lie_group_base<NAME<_Ty>>; \
public:
#define MATRICE_MAKE_LIE_GROUP_END(NAME) }

/**
 *\brief Class template for 2D Special Orthogonal (Rotation) group.
 */
MATRICE_MAKE_LIE_GROUP_BEGIN(_SO2)

MATRICE_MAKE_LIE_GROUP_END(_SO2);
/**
 *\brief Class template for 3D Special Orthogonal (Rotation) group.
 */
MATRICE_MAKE_LIE_GROUP_BEGIN(_SO3)

MATRICE_MAKE_LIE_GROUP_END(_SO3);

/**
 *\brief Class template for 2D Special Euclidean (Rigid motion) group.
 */
MATRICE_MAKE_LIE_GROUP_BEGIN(_SE2)

MATRICE_MAKE_LIE_GROUP_END(_SE2);
/**
 *\brief Class template for 3D Special Euclidean (Rigid motion) group.
 */
MATRICE_MAKE_LIE_GROUP_BEGIN(_SE3)

MATRICE_MAKE_LIE_GROUP_END(_SE3);

#undef MATRICE_MAKE_LIE_GROUP_BEGIN
#undef MATRICE_MAKE_LIE_GROUP_END
_DETAIL_END
DGE_MATRICE_END