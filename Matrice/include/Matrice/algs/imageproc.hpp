/***************************************************************************
This file is part of Matrice, an effcient and elegant C++ library for SC.
      Copyright(C) 2018, Zhilong (Dgelom) Su (su-zl@seu.edu.cn), 
		                   all rights reserved.

This program is free software : you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.
****************************************************************************/
#pragma once
#include <type_traits>
#include "../private/_range.h"
#include "../algs/interpolation.h"

DGE_MATRICE_BEGIN

template<typename _Pixty, typename = std::enable_if_t<std::is_arithmetic_v<_Pixty>>>
void _Grayscale_stretch(_Pixty* Img, std::size_t rows, std::size_t cols) {
	using pixel_t = _Pixty;
	using iterator = std::add_pointer_t<pixel_t>;

	const std::size_t N = rows * cols;
	iterator _Begin = Img, _End = Img + N;

	pixel_t _Max = pixel_t(0), _Min = pixel_t(255);
	for (; _Begin != _End; ++_Begin) {
		_Max = _Max > *_Begin ? _Max : *_Begin;
		_Min = _Min < *_Begin ? _Min : *_Begin;
	}

	float_t _Scal = 255.f / float_t(_Max - _Min);
	for (_Begin -= N; _Begin != _End; ++_Begin)
		*_Begin = pixel_t(_Scal * (*_Begin - _Min));
}

namespace types { template<typename _Ty, int _M, int _N> class Matrix_; }

namespace detail {
template<typename _Ty, int _M, int _N> class _Multi_matrix;

template<typename _Ty> struct gradient_traits {};
template<std::size_t _Opt> 
struct splineitp_option{static constexpr auto value = _Opt;};
template<std::size_t _Opt>
MATRICE_HOST_INL constexpr auto splineitp_option_v = splineitp_option<_Opt>::value;
template<> struct splineitp_option<BSPL3> { static constexpr auto value = bcspline; };
template<> struct splineitp_option<BSPL5> { static constexpr auto value = bqspline; };
template<> struct splineitp_option<BSPL7> { static constexpr auto value = bsspline; };

template<typename _Ty, std::size_t _Opt> class _Gradient_impl;

template<typename _Ty, std::size_t _Opt> 
struct gradient_traits<_Gradient_impl<_Ty, _Opt>> {
	using value_type = _Ty;
	static constexpr auto option = _Opt;
	using image_type = types::Matrix_<value_type, 0, 0>;
	using matrix_type = types::Matrix_<value_type, 0, 0>;
};

template<std::size_t _Opt> struct _Grad_range_clip {};
template<> struct _Grad_range_clip<BSPL3> {
	template<typename _Ity>
	MATRICE_GLOBAL_INL static auto value(const _Ity _L, const _Ity _U) {
		return range(_L + 1, _U - 3);
	}
};
template<> struct _Grad_range_clip<BSPL5> {
	template<typename _Ity>
	MATRICE_GLOBAL_INL static auto value(const _Ity _L, const _Ity _U) {
		return range(_L + 2, _U - 4);
	}
};
template<> struct _Grad_range_clip<BSPL7> {
	template<typename _Ity>
	MATRICE_GLOBAL_INL static auto value(const _Ity _L, const _Ity _U) {
		return range(_L + 3, _U - 5);
	}
};
/**
 * \Base class for gradient computation in interpolation way.
 */
template<typename _Derived> class _Interpolated_gradient_base {
	using _Mytraits = gradient_traits<_Derived>;
	using _Myitp = interpolation<typename _Mytraits::value_type, splineitp_option_v<_Mytraits::option>>;
public:
	using value_type = typename _Mytraits::value_type;
	using image_type = typename _Mytraits::image_type;
	using matrix_type = typename _Mytraits::matrix_type;
	using point_type = typename _Myitp::type::point_type;
	
	_Interpolated_gradient_base(const image_type& _Image)
		: _Myimg(_Image), _Myop(_Myitp(_Myimg)) {
	}
	_Interpolated_gradient_base(const image_type& _Image, const typename _Myitp::type& _Op)
		: _Myimg(_Image), _Myop(_Op) {
	}

	/**
	 * \Get image gradient: $\fract{\partial I}{\partial _Axis}(_Pos)$
	 */
	MATRICE_HOST_INL auto at(const point_type& _Pos) const {
		return _Myop.grad(_Pos);
	}
	template<axis _Axis>
	MATRICE_HOST_INL auto at(const point_type& _Pos) const { 
		return _Myop.grad<_Axis>(_Pos); 
	}
	/**
	 * \Get image gradient w.r.t. _Axis in rect. range [_L, _R) | [_U, _D)
	 */
	template<axis _Axis> 
	MATRICE_HOST_INL auto at(int _L, int _R, int _U, int _D) const { 
		using _My_range = _Grad_range_clip<_Mytraits::option>;
		matrix_type _Grad(_D - _U, _R - _L);
		for (const auto _Idy : _My_range::value(_U,_D)) {
			auto _Row = _Grad.rbegin(_Idy - _U);
			for (const auto _Idx : _My_range::value(_L, _R)) {
				_Row[_Idx - _L] = at<_Axis>(point_type(_Idx, _Idy));
			}
		}

		return std::forward<matrix_type>(_Grad);
	}

protected:
	const image_type& _Myimg;
	typename _Myitp::type _Myop;
};

template<typename _Ty> class _Gradient_impl<_Ty, SOBEL>  {
	using _Myt = _Gradient_impl;
public:
	using value_type = _Ty;
	using image_type = Matrix<value_type>;
	_Gradient_impl(const image_type& _Image) : _Myimg(_Image) {}

private:
	const image_type& _Myimg;
};

template<typename _Ty> class _Gradient_impl<_Ty, BSPL3>
	: public _Interpolated_gradient_base<_Gradient_impl<_Ty, BSPL3>> {
	using _Myt = _Gradient_impl;
	using _Mybase = _Interpolated_gradient_base<_Myt>;
public:
	using typename _Mybase::image_type;
	using typename _Mybase::value_type;

	_Gradient_impl(const image_type& _Img) : _Mybase(_Img) {}
};
}

/**
 * \TEMPLATE class for image gradient computation.
 * \PARAMS : <_Ty> the data type; <_Alg> the gradient operator.
 */
template<typename _Ty, std::size_t _Alg = BSPL3>
using gradient = detail::_Gradient_impl<_Ty, _Alg>;

DGE_MATRICE_END