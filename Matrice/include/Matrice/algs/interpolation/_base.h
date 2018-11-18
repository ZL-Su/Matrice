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

#pragma region <!-- Codes in the region will be deprecated -->

template<typename _Ty, size_t _Opt> class BilinearInterp;
template<typename _Ty, size_t _Opt> class BicubicSplineInterp;
template<typename _Ty, size_t _Opt> class BiquinticSplineInterp;
template<typename _Ty, size_t _Opt> class BisepticSplineInterp;

template<typename _Ty, size_t _Opt = INTERP|BILINEAR>
struct interpolation_traits
{ using type = BilinearInterp<_Ty, _Opt>; };
template<typename _Ty> 
struct interpolation_traits<_Ty, INTERP|BICUBIC|BSPLINE>
{ using type = BicubicSplineInterp<_Ty, INTERP|BICUBIC|BSPLINE>; };
template<typename _Ty>
struct interpolation_traits<_Ty, INTERP|BIQUINTIC|BSPLINE>
{ using type = BiquinticSplineInterp<_Ty, INTERP|BIQUINTIC|BSPLINE>; };
template<typename _Ty>
struct interpolation_traits<_Ty, INTERP|BISEPTIC|BSPLINE>
{ using type = BisepticSplineInterp<_Ty, INTERP|BISEPTIC|BSPLINE>; };

template<typename _Ty, size_t _Opt>
using interpolation_traits_t = typename interpolation_traits<_Ty, _Opt>::type;

template<typename _Ty, typename _Derived> class InterpBase_
{
	using derived_t = _Derived;
public:
	using value_t = _Ty;
	using value_type = value_t;
	using matrix_t = Matrix<value_type>;

	template<typename... _Args>
	MATRICE_GLOBAL_FINL InterpBase_(const _Args&... args);

	MATRICE_GLOBAL_FINL auto& operator() () const { return m_coeff; }

protected:
	template<typename... _Args>
	MATRICE_GLOBAL_FINL void _Bspline_coeff(const _Args& ...args);

	const value_type m_eps = value_type(1.0e-7);
	matrix_t m_coeff;
};

#pragma endregion

template<typename _Ty, std::size_t _Opt> class _Spline_interpolation;

template<typename _Ty, std::size_t _Opt>
struct interpolation_traits<_Spline_interpolation<_Ty, _Opt>> {
	using value_type = _Ty;
	using matrix_type = Matrix<value_type>;
	using type = _Spline_interpolation<value_type, _Opt>;
	static constexpr auto option = _Opt;
};
template<typename _Ty, std::size_t _Opt> struct auto_interp_dispatcher {
	using type = _Spline_interpolation<_Ty, _Opt>;
};
template<typename _Ty, std::size_t _Opt>
using auto_interp_dispatcher_t = typename auto_interp_dispatcher<_Ty, _Opt>::type;

template<typename _Derived> class _Interpolation_base {
	using _Myt = _Interpolation_base;
	using _Mydt = _Derived;
	using _Mytraits = interpolation_traits<_Mydt>;
public:
	using value_type = typename _Mytraits::value_type;
	using matrix_type = typename _Mytraits::matrix_type;
	using point_type = Vec2_<value_type>;
	static constexpr auto option = _Mytraits::option;

	_Interpolation_base(const matrix_type& _Data) : _Mydata(_Data) {
		static_cast<_Mydt*>(this)->_Coeff_impl();
	}

	/**
	 * \get interpolation coeff. matrix.
	 */
	MATRICE_HOST_INL auto& operator()() const {
		return (_Mycoeff);
	}

	/**
	 * \get the interpolated value at _Pos. 
	 */
	MATRICE_HOST_INL auto operator()(const point_type& _Pos) const {
		return static_cast<const _Mydt*>(this)->_Value_at(_Pos);
	}

	/**
	 * \get the interpolated gradient value at _Pos.
	 */
	MATRICE_HOST_INL auto grad(const point_type& _Pos) const {
		return std::make_tuple(
			static_cast<const _Mydt*>(this)->_Gradx_at(_Pos),
			static_cast<const _Mydt*>(this)->_Grady_at(_Pos)
			);
	}
	template<axis _Axis, typename = std::enable_if_t<_Axis==axis::x||_Axis==axis::y>>
	MATRICE_HOST_INL auto grad(const point_type& _Pos) const {
		if constexpr (_Axis == axis::x)
			return static_cast<const _Mydt*>(this)->_Gradx_at(_Pos);
		if constexpr (_Axis == axis::y)
			return static_cast<const _Mydt*>(this)->_Grady_at(_Pos);
	}

protected:
	const matrix_type& _Mydata;
	const value_type _Myeps = value_type(1.0e-7);
	matrix_type _Mycoeff;
};

MATRICE_ALGS_END
#include "_base.inl"