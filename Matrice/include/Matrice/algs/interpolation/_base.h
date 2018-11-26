/***********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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

#include "../../core/matrix.h"
#include "../../core/vector.h"
#include "../../core/solver.h"
#include "../../private/_range.h"
#include "../../private/_tag_defs.h"

MATRICE_ALGS_BEGIN
enum {
	INTERP     = 18,

	BSPLINE    = 2307,
	BILINEAR   = 2152,
	BICUBIC    = 2156,
	BIQUINTIC  = 2160,
	BISEPTIC   = 2164,
	_BICBSPL = INTERP | BICUBIC | BSPLINE,
	_BIQNSPL = INTERP | BIQUINTIC | BSPLINE,
	_BISPSPL = INTERP | BISEPTIC | BSPLINE,
};

// \Forward declaration
//template<typename _Ty, typename _Tag> class _Spline_interpolation;
template<typename _Ty, typename _Tag> class _Spline_interpolation;

// \Interpolation traits definition
template<typename _Ty> struct interpolation_traits {};
template<typename _Ty, typename _Tag>
struct interpolation_traits<_Spline_interpolation<_Ty, _Tag>> {
	using value_type = _Ty;
	using matrix_type = Matrix<value_type>;
	using category = _Tag;
	using type = _Spline_interpolation<value_type, category>;
	static constexpr auto option = INTERP|BSPLINE;
};
template<typename _Ty>
using interpolation_traits_t = typename interpolation_traits<_Ty>::type;

// \Interpolation auto dispatching
template<typename _Ty, typename _Tag> struct auto_interp_dispatcher {
	using type = _Spline_interpolation<_Ty, _Tag>;
};
template<typename _Ty = float, typename _Tag = _TAG bicspl_tag>
using auto_interp_dispatcher_t = typename auto_interp_dispatcher<_Ty, _Tag>::type;

template<typename _Derived> class _Interpolation_base {
	using _Myt = _Interpolation_base;
	using _Mydt = _Derived;
	using _Mytraits = interpolation_traits<_Mydt>;
public:
	using category = category_type_t<_Mytraits>;
	using value_type = typename _Mytraits::value_type;
	using matrix_type = typename _Mytraits::matrix_type;
	using point_type = Vec2_<value_type>;
	static constexpr auto option = _Mytraits::option;

	_Interpolation_base(const matrix_type& _Data) noexcept
		: _Mydata(_Data) { static_cast<_Mydt*>(this)->_Coeff_impl(); }
	_Interpolation_base(const _Myt& _Other) noexcept
		: _Mydata(_Other._Mydata), _Mycoeff(_Other._Mycoeff) {}

	/**
	 * \get interpolation coeff. matrix.
	 */
	MATRICE_HOST_INL auto& operator()() const { return (_Mycoeff); }

	/**
	 * \get the interpolated value at _Pos. 
	 */
	MATRICE_HOST_INL auto operator()(const point_type& _Pos) const {
		return (this)->_Value_at(_Pos);
	}
	MATRICE_HOST_INL auto operator()(value_type _X, value_type _Y) const {
		return (this)->_Value_at(point_type(_X, _Y));
	}

	/**
	 * \get the interpolated gradient value at _Pos.
	 */
	MATRICE_HOST_INL auto grad(const point_type& _Pos) const {
		return std::make_tuple((this)->_Gradx_at(_Pos), (this)->_Grady_at(_Pos));
	}
	template<axis _Axis, typename = std::enable_if_t<_Axis==axis::x||_Axis==axis::y>>
	MATRICE_HOST_INL auto grad(const point_type& _Pos) const {
		if constexpr (_Axis == axis::x) return (this)->_Gradx_at(_Pos);
		if constexpr (_Axis == axis::y) return (this)->_Grady_at(_Pos);
	}

private:
	MATRICE_HOST_INL auto _Value_at(const point_type& _Pos) const;
	MATRICE_HOST_INL auto _Gradx_at(const point_type& _Pos) const;
	MATRICE_HOST_INL auto _Grady_at(const point_type& _Pos) const;

	std::add_pointer_t<_Mydt> _Mydt_this = static_cast<_Mydt*>(this);

protected:
	const matrix_type& _Mydata;
	const value_type _Myeps = value_type(1.0e-7);
	matrix_type _Mycoeff;
};

MATRICE_ALGS_END
#include "_base.inl"