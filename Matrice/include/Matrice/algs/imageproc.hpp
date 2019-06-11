/***************************************************************************
This file is part of Matrice, an effcient and elegant C++ library for SC.
Copyright(C) 2018-2019, Zhilong (Dgelom) Su (su-zl@seu.edu.cn), all rights reserved.

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

template<typename _Pixty, MATRICE_ENABLE_IF(is_arithmetic_v<_Pixty>)>
void _Grayscale_stretch(_Pixty* Img, size_t rows, size_t cols) {
	using pixel_t = _Pixty;
	using iterator = std::add_pointer_t<pixel_t>;

	const auto N = rows * cols;
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

namespace types { 
	template<typename _Ty, int _M, int _N> class Matrix_; 
}

namespace detail {
template<typename _Ty, int _M, int _N> class _Multi_matrix;

template<typename _Ty> struct gradient_traits {};
template<typename _Tag> struct splineitp_option {using tag = _Tag;};
template<typename _Tag>
using splineitp_option_t = typename splineitp_option<_Tag>::tag;
template<> struct splineitp_option<_TAG _Itped_grad_tag::bicspl> { 
	using tag = _TAG _Bspline_itp_tag::bicubic; 
};
template<> struct splineitp_option<_TAG _Itped_grad_tag::biqspl> {
	using tag = _TAG _Bspline_itp_tag::biquintic;
};
template<> struct splineitp_option<_TAG _Itped_grad_tag::bisspl> {
	using tag = _TAG _Bspline_itp_tag::biseptic;
};

template<typename _Ty, typename _Tag> class _Gradient_impl;

template<typename _Ty, typename _Tag>
struct gradient_traits<_Gradient_impl<_Ty, _Tag>> {
	using value_type = _Ty;
	using image_type = types::Matrix_<value_type, 0, 0>;
	using matrix_type = types::Matrix_<value_type, 0, 0>;
	using category = _Tag;
};

template<typename _Tag> struct _Grad_range_clip {};
template<> struct _Grad_range_clip<_TAG _Itped_grad_tag::bicspl> {
	template<typename _Ity>
	MATRICE_GLOBAL_INL static auto _(const _Ity _L, const _Ity _U) {
		return range(_L + 1, _U - 3);
	}
};
template<> struct _Grad_range_clip<_TAG _Itped_grad_tag::biqspl> {
	template<typename _Ity>
	MATRICE_GLOBAL_INL static auto _(const _Ity _L, const _Ity _U) {
		return range(_L + 2, _U - 4);
	}
};
template<> struct _Grad_range_clip<_TAG _Itped_grad_tag::bisspl> {
	template<typename _Ity>
	MATRICE_GLOBAL_INL static auto _(const _Ity _L, const _Ity _U) {
		return range(_L + 3, _U - 5);
	}
};
/**
 * \Base class for gradient computation in interpolation way.
 */
template<typename _Derived> class _Interpolated_gradient_base {
	using _Mytraits = gradient_traits<_Derived>;
	using _Myitp = interpolation<typename _Mytraits::value_type, category_type_t<_Mytraits>>;
public:
	using category = category_type_t<_Mytraits>;
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
	MATRICE_HOST_INL value_type at(const point_type& _Pos) const {
		return _Myop.grad(_Pos);
	}
	template<axis _Axis>
	MATRICE_HOST_INL value_type at(const point_type& _Pos) const {
		return _Myop.grad<_Axis>(_Pos); 
	}
	/**
	 * \Get image gradient w.r.t. _Axis in rect. range [_L, _R) | [_U, _D)
	 */
	template<axis _Axis> 
	MATRICE_HOST_INL decltype(auto)at(int _L, int _R, int _U, int _D) const { 
		using _My_range = _Grad_range_clip<category>;
		const auto _Ry = _My_range::_(_U, _D);
		const auto _Rx = _My_range::_(_L, _R);

		matrix_type _Grad(_D - _U, _R - _L, zero<value_type>);
		for (const auto _Idy : _Ry) {
			auto _Row = _Grad.rbegin(_Idy - _U);
			for (const auto _Idx : _Rx) {
				_Row[_Idx - _L] = at<_Axis>(point_type(_Idx, _Idy));
			}
		}

		return std::forward<matrix_type>(_Grad);
	}
	/**
	 * \Get image interpolator
	 */
	MATRICE_HOST_INL auto& itp() { return (_Myop); }
	MATRICE_HOST_INL const auto& itp() const { return (_Myop); }

protected:
	const image_type& _Myimg;
	typename _Myitp::type _Myop;
};

template<typename _Ty> 
class _Gradient_impl<_Ty, _TAG _Sobel_grad_tag> {
	using _Myt = _Gradient_impl;
public:
	using value_type = _Ty;
	using image_type = Matrix<value_type>;
	_Gradient_impl(const image_type& _Image) : _Myimg(_Image) {}

	MATRICE_HOST_INL auto at(diff_t _x, diff_t _y) const {
		auto gx = zero<value_type>, gy = gx;

		if (_x > 0 && _y > 0 && _x < _Myimg.cols() && _y < _Myimg.rows()) {
			const auto pu = _Myimg[_y - 1], pc = _Myimg[_y], pl = _Myimg[_y + 1];

			const auto I00 = pu[_x - 1], I01 = pu[_x], I02 = pu[_x + 1];
			const auto I10 = pc[_x - 1], I11 = pc[_x], I12 = pc[_x + 1];
			const auto I20 = pl[_x - 1], I21 = pl[_x], I22 = pl[_x + 1];

			gx = (I01 - I00) + two<value_type>*(I12 - I10) + (I22 - I20);
			gy = (I20 - I00) + two<value_type>*(I21 - I01) + (I22 - I02);
		}

		return tuple<value_type, value_type>(gx, gy);
	}
	template<typename _Op>
	MATRICE_HOST_INL auto eval(_Op&& _op) const {
		return _op(_Myimg);
	}

private:
	const image_type& _Myimg;
};

template<typename _Ty> 
class _Gradient_impl<_Ty, _TAG _Itped_grad_tag::bicspl>
	: public _Interpolated_gradient_base<_Gradient_impl<_Ty, _TAG _Itped_grad_tag::bicspl>> {
	using _Myt = _Gradient_impl;
	using _Mybase = _Interpolated_gradient_base<_Myt>;
public:
	using typename _Mybase::image_type;
	using typename _Mybase::value_type;

	_Gradient_impl(const image_type& _Img) : _Mybase(_Img) {}
};

template<typename _Ty>
class _Gradient_impl<_Ty, _TAG _Itped_grad_tag::biqspl>
	: public _Interpolated_gradient_base<_Gradient_impl<_Ty, _TAG _Itped_grad_tag::biqspl>> {
	using _Myt = _Gradient_impl;
	using _Mybase = _Interpolated_gradient_base<_Myt>;
public:
	using typename _Mybase::image_type;
	using typename _Mybase::value_type;

	_Gradient_impl(const image_type& _Img) : _Mybase(_Img) {}
};

template<typename _Ty>
class _Gradient_impl<_Ty, _TAG _Itped_grad_tag::bisspl>
	: public _Interpolated_gradient_base<_Gradient_impl<_Ty, _TAG _Itped_grad_tag::bisspl>> {
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
template<typename _Ty, typename _Tag = _TAG _Itped_grad_tag::bicspl>
using gradient = detail::_Gradient_impl<_Ty, _Tag>;
template<typename _Ty>
using trivial_imgrad = detail::_Gradient_impl<_Ty, _TAG _Sobel_grad_tag>;

DGE_MATRICE_END