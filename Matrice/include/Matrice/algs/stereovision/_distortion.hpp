/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2022, Zhilong(Dgelom) Su, all rights reserved.

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
#include "core.hpp"
#include "internal/expr_base.hpp"

MATRICE_NAMESPACE_BEGIN(vision)

struct fwd { enum { forward = 0 }; };
struct bwd { enum { backward = 0 }; };

template<typename _Ty, class _Tag = void>
requires is_floating_point_v<_Ty>
class distortion : public xpr::__xpr__ {
	using _Myt = distortion;
	using _Mydt = Vec4_<_Ty>;
public:
	using value_t = _Mydt::value_type;

	explicit distortion(_Mydt::const_initlist factors) noexcept
		:_Mydata(factors) {
	}

	/// <summary>
	/// \brief Apply distortion to ideal image coordiantes (x, y).
	/// </summary>
	/// <returns> Distorted coordinates in the normalized image domain.</returns> 
	MATRICE_GLOBAL_FINL auto apply(value_t x, value_t y) const noexcept {
		return _Apply(x, y);
	}

private:
	_Mydt _Mydata;
	MATRICE_GLOBAL_FINL auto _Apply(value_t x, value_t y) const noexcept;
};
template<typename _Ty, class _Tag> requires is_floating_point_v<_Ty>
MATRICE_GLOBAL_FINL auto distortion<_Ty, _Tag>::_Apply(value_t x, value_t y)const noexcept
{
	const auto r_2 = sqsum(x, y);
	const auto [k1, k2, p1, p2] = _Mydata.unbind();
	const auto tmp = 1 + k1*r_2 + k2*sq(r_2) + 2*(p1*x + p2*y);
	const auto x_d = tmp * x + p2 * r_2;
	const auto y_d = tmp * y + p1 * r_2;
	MATRICE_USE_STD(make_tuple);
	return make_tuple(x_d, y_d);
}

template<typename _Ty>
	requires is_floating_point_v<_Ty>
class distortion<_Ty, bwd> : public distortion<_Ty, void> {
	using _Myt = distortion;
	using _Mybase = distortion<_Ty, void>;
	using _Mydt = Vec4_<_Ty>;
public:
	using value_t = _Mydt::value_type;

};

template<typename _Ty>
	requires is_floating_point_v<_Ty>
class distortion<_Ty, fwd> : public distortion<_Ty, void> {
	using _Myt = distortion;
	using _Mybase = distortion<_Ty, void>;
	using _Mydt = Vec4_<_Ty>;
public:
	using value_t = _Mydt::value_type;

};

template<typename _Ty> using fwd_dmodel_t = distortion<_Ty, fwd>;
template<typename _Ty> using bwd_dmodel_t = distortion<_Ty, bwd>;

MATRICE_NAMESPACE_END(vision)
