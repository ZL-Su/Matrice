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

enum class axis {x, y, z};

template<typename _Pixty, typename = std::enable_if_t<std::is_arithmetic_v<_Pixty>>>
void _Grayscale_stretch(_Pixty* Img, std::size_t rows, std::size_t cols) {
	using pixel_t = _Pixty;
	using iterator = std::add_pointer_t<pixel_t>;

	const std::size_t N = rows * cols;
	iterator _Begin = Img, _End = Img + N;

	pixel_t _Max = pixel_t(0), _Min = pixel_t(255);
	for (; _Begin != _End; ++_Begin)
	{
		_Max = _Max > *_Begin ? _Max : *_Begin;
		_Min = _Min < *_Begin ? _Min : *_Begin;
	}

	float_t _Scal = 255.f / float_t(_Max - _Min);
	for (_Begin -= N; _Begin != _End; ++_Begin)
		*_Begin = pixel_t(_Scal * (*_Begin - _Min));
}

namespace types {
template<typename _Ty, int _M, int _N> class Matrix_;
}

enum {SOBEL = 0, BSPL3 = 3, BSPL5 = 5, BSPL7 = 7};

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
	using image_type = types::Matrix_<unsigned char, 0, 0>;
	using matrix_type = types::Matrix_<value_type, 0, 0>;
};

/**
 * \Base class for gradient computation in interpolation way.
 */
template<typename _Derived> class _Gradient_base {
	using _Mytraits = gradient_traits<_Derived>;
	using _Myitp_type = interpolation<value_type, splineitp_option_v<_Mytraits::option>>;
public:
	using value_type = typename _Mytraits::value_type;
	using image_type = typename _Mytraits::image_type;
	using matrix_type = typename _Mytraits::matrix_type;
	
	_Gradient_base(const image_type& _Image) : _Myimg(_Image) {
		_Myop = _Myitp_type(_Myimg)();
	}

protected:
	const image_type& _Myimg;
	typename _Myitp_type::kernel_type _Myop;
};

template<typename _Ty>
class _Gradient_impl<_Ty, SOBEL>  {
	using _Myt = _Gradient_impl;
	using _Mybase = _Gradient_base<_Myt>;
public:
	using typename _Mybase::image_type;
	_Gradient_impl(const image_type& _Image) : _Mybase(_Image) {}

};
template<typename _Ty>
class _Gradient_impl<_Ty, BSPL3>
	: public _Gradient_base<_Gradient_impl<_Ty, BSPL3>> {
	using _Myt = _Gradient_impl;
	using _Mybase = _Gradient_base<_Myt>;
	using _Myitp_type = interpolation<typename _Mybase::value_type, bcspline>;
public:
	using typename _Mybase::image_type;
	using typename _Mybase::value_type;
	_Gradient_impl(const image_type& _Image) 
		: _Mybase(_Image), _Myop(_Myitp_type(_Mybase::_Myimg)()) {
	}

	template<axis _Axis, typename _Valty>
	MATRICE_HOST_INL auto wrt(_Valty _x, _Valty _y) const {
		if constexpr (_Axis = axis::x) return _Myop.gradx_at();
		if constexpr (_Axis = axis::y);
	}

private:
	typename _Myitp_type::kernel_type _Myop;
};
}

DGE_MATRICE_END