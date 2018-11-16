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
namespace dgelom {

template<typename _Pixty, typename = typename std::enable_if<std::is_arithmetic<_Pixty>::value>::type>
void _Grayscale_stretch(_Pixty* Img, std::size_t rows, std::size_t cols)
{
	typedef _Pixty                         pixel_t;
	typedef pixel_t*                      iterator;

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

template<typename _Ty, std::size_t _Op> class _Gradient_impl;
template<typename _Ty, std::size_t _Op> 
struct gradient_traits<_Gradient_impl<_Ty, _Op>> {
	using value_type = _Ty;
	static constexpr auto op_value = _Op;
	using image_type = types::Matrix_<unsigned char, 0, 0>;
	using matrix_type = types::Matrix_<value_type, 0, 0>;
};

template<typename _Derived> class _Gradient_base {
	using _Mytraits = gradient_traits<_Derived>;
public:
	using value_type = typename _Mytraits::value_type;
	using image_type = typename _Mytraits::image_type;
	using matrix_type = typename _Mytraits::matrix_type;

	_Gradient_base(const image_type& _Image) : _Myimg(_Image) {}


	/**
	 * \get image gradient with respective to _Axis
	 * \Template param _Axis: 0 for x axis, 1 for y axis
	 */
	template<std::size_t _Axis> auto& get() {
		return (_Mygrad[_Axis]);
	}
	template<std::size_t _Axis> const auto& get() const {
		return (_Mygrad[_Axis]);
	}

private:
	const image_type& _Myimg;
	_Multi_matrix<value_type, 0, 0> _Mygrad;
};

template<typename _Ty>
class _Gradient_impl<_Ty, SOBEL> 
	: public _Gradient_base<_Gradient_impl<_Ty, SOBEL>> {
	using _Myt = _Gradient_impl;
	using _Mybase = _Gradient_base<_Myt>;
public:
	using typename _Mybase::image_type;
	_Gradient_impl(const image_type& _Image) : _Mybase(_Image) {}



};

}
