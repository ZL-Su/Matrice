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

#include "../util/utils.h"
#include "../forward.hpp"

DGE_MATRICE_BEGIN
struct bicubic_tag {};
struct biquintic_tag {};
struct biseptic_tag {};

_DETAIL_BEGIN
// Forward declarations
template<typename _Derived> class _Spline;
template<typename _Ty, typename _Cat, typename... _Args> class _Bspline {};

// Traits for CRTP
template<typename _Ty, typename _Cat, typename... _Args>
struct traits<_Bspline<_Ty, _Cat, _Args...>> {
	using value_type = _Ty;
	using category = _Cat;
	using type = detail::_Bspline<value_type, category, _Args...>;
};

_DETAIL_END

template<typename _Ty>
using bicubic_bspline_t = detail::_Bspline<_Ty, bicubic_tag>;

DGE_MATRICE_END

#include "spline/_spline.hpp"