/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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

#include "_base.h"

MATRICE_ALGS_BEGIN
_DETAIL_BEGIN

/**
 *\brief N-dimensional grid iterator (nested loop with variable depth).
 *\param <_ND> template argument for grid dimension, viz., _ND = 2 or 3.
 */
template<uint8_t _ND=2> class _Grid_iterator {
	static_assert(_ND <= 4, "Dimensions cannot be greater than 4.");
	using _Myt = _Grid_iterator;
public:
	using idx_array = Vec_<size_t, _ND>;
	explicit _Grid_iterator(const idx_array& _dims)
		: _Outer(_dims), _Pos(0) {
		fill(_Inner, 0);
		_Is_end = (_Inner == _Outer);
	}
	explicit _Grid_iterator(size_t _dim) : _Pos(0) {
		fill(_Outer, _dim), fill(_Inner, 0);
	}

	MATRICE_HOST_INL size_t operator[](size_t i) const { 
		return _Inner[i]; 
	}
	MATRICE_HOST_INL const idx_array& operator*() const { 
		return _Inner; 
	}
	MATRICE_HOST_INL _Myt& operator++() {
		_Is_end = std::true_type::value;
		for (auto k = _ND; k--;) {
			if (++_Inner[k] < _Outer[k]) {
				_Is_end = std::false_type::value;
				break;
			}
			_Inner[k] = 0;
		}
		++_Pos;
		return (*this);
	}
	MATRICE_HOST_INL operator bool() const {
		return (!_Is_end);
	}
	MATRICE_HOST_INL size_t pos() const { 
		return _Pos; 
	}

private:
	bool _Is_end = std::true_type::value;
	idx_array _Outer, _Inner;
	typename idx_array::value_type _Pos;
};

// \define lattice types
enum class _Lattice_tag {approx, sparse, dense};
// \forward declaration of lattice type
template<typename _Ty, uint8_t _ND, _Lattice_tag _Tag> class _Lattice;
// \retrieve lattice traits
template<typename _Ty, uint8_t _ND, _Lattice_tag _Tag>
struct traits<_Lattice<_Ty, _ND, _Tag>> {
	using value_type = _Ty;
	static constexpr auto dimension = _ND;
	_Lattice_tag lattice = _Tag;
};

/**
 *\brief lattice base class
 *\param <_Derived> template type, which is a derived type of the base type 
 */
template<typename _Derived> class _Lattice_base {
	using _Mytraits = traits<_Derived>;
	using _Myt = _Lattice_base;
	using _Mydt = _Derived;
public:
	using value_type = typename _Mytraits::value_type;
	using idx_array = Vec_<size_t, _Mytraits::dimension>;
	using point_type = Vec_<value_type, _Mytraits::dimension>;
	struct options_type {
		point_type _Min, _Max;//min and max points of the given scatter data
		idx_array _Gridesize; //grid size
	};

	/**
	 *\brief Evaluate the interpolated value at point p
	 *\param [p] holds grid coordinates
	 */
	MATRICE_HOST_INL value_type operator()(const point_type& p) const {
		return static_cast<const _Mydt*>(this)->_Eval(p);
	}

	/**
	 *\brief calculate interpolated residual
	 *\param [_cbegin, _cend] iterator range of grid points
	 */
	template<typename _CooIt, typename _ValIt>
	MATRICE_HOST_INL value_type residual(const _CooIt& _cbegin, const _CooIt& _cend, _ValIt& _vbegin)const {
		value_type _Res = zero<value_type>;

		_CooIt p = _cbegin;
		_ValIt v = _vbegin;
		for (; p != _cend; ++p, ++v) {
			(*v) -= (*this)(*p);
			_Res = max(_Res, abs(*v));
		}

		return (_Res);
	}

	/**
	 *\brief helper function to report interpolation status
	 *\param [os] i/o stream.
	 */
	MATRICE_HOST_INL void report(std::iostream& os) const {
		return static_cast<const _Mydt*>(this)->_Report(os);
	}

protected:
	static value_type _Bicubic_kernel(value_type t, size_t k) {
		assert(k < 4);
		switch (k) {
		case 0: return pow(1 - t, 3) / 6; break;
		case 1: return pow(t, 3) / 2 - sqr(t) + 0.6667; break;
		case 2: return (sqr(t) - pow(t, 3) + t) / 2 + 0.16667; break;
		case 3: return pow(t, 3) / 6.; break;
		default: return zero<value_type>;
		}
	}
	static bool _Is_boxed(const idx_array& lt, const idx_array& rb, const idx_array& p) {
		for (auto i = 0; i < _ND; ++i) {
			if (p[i] < lt[i] || p[i] > rb[i])
				return std::false_type::value;
		}
		return std::true_type::value;
	}
};

template<typename _Ty, uint8_t _N>
class _Lattice<_Ty, _N, _Lattice_tag::approx>
	: public _Lattice_base<_Lattice<_Ty, _N, _Lattice_tag::approx>> 
{
	using _Myt = _Lattice;
	using _Mybase = _Lattice_base<_Myt>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::point_type;
	_Lattice(std::function<value_type(const point_type&)>&& _f)
		: _Myfunc(_f) {}

	MATRICE_HOST_INL value_type _Eval(const point_type& p) const {
		return (_Myfunc(p));
	}

	MATRICE_HOST_INL void _Report(std::iostream& io) const {
		os << "[Matrice report] initial approximation";
	}

private:
	std::function<value_type(const point_type&)> _Myfunc;
};

template<typename _Ty, uint8_t _N>
class _Lattice<_Ty, _N, _Lattice_tag::dense>
	: public _Lattice_base<_Lattice<_Ty, _N, _Lattice_tag::dense>>
{
	using _Myt = _Lattice;
	using _Mybase = _Lattice_base<_Myt>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::point_type;
	using typename _Mybase::idx_array;
	using typename _Mybase::options_type;

	template<typename CIt, typename VIt>
	_Lattice(const CIt& cbegin, const CIt& cend, const VIt& vbegin, const options_type& opts)
		:{

	}

	MATRICE_HOST_INL value_type _Eval(const point_type& p) const {
		return (_Myfunc(p));
	}

	MATRICE_HOST_INL void _Report(std::iostream& io) const {
		os << "[Matrice report] initial approximation";
	}

private:
	point_type _Myinv;
	options_type _Myopts;
	Matrix<value_type> _Myphi;
};

template<typename _Ty, uint8_t _ND>
class _Lattice<_Ty, _ND, _Lattice_tag::sparse>
	: public _Lattice_base<_Lattice<_Ty, _ND, _Lattice_tag::sparse>>
{
	using _Myt = _Lattice;
	using _Mybase = _Lattice_base<_Myt>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::point_type;
	_Lattice(std::function<value_type(const point_type&)>&& _f)
		: _Myfunc(_f) {}

	MATRICE_HOST_INL value_type _Eval(const point_type& p) const {
		return (_Myfunc(p));
	}

	MATRICE_HOST_INL void _Report(std::iostream& io) const {
		os << "[Matrice report] initial approximation";
	}

private:
	std::function<value_type(const point_type&)> _Myfunc;
};

_DETAIL_END

/**
 *\brief multilevel bicubic spline interpolation interface
 *\param <_Ty> a scalar template type
 */
template<typename _Ty> 
class _Spline_interpolation<_Ty, _TAG mbicspl_tag> : public
	_Interpolation_base<_Spline_interpolation<_Ty, _TAG mbicspl_tag>>
{
	static_assert(is_scalar_v<_Ty>, "template type _Ty must be a scalar.");
	using _Myt = _Spline_interpolation<_Ty, _TAG mbicspl_tag>;
	using _Mybase = _Interpolation_base<_Myt>;
public:
	using typename _Mybase::category;
	using typename _Mybase::value_type;
	using typename _Mybase::point_type;
	using typename _Mybase::matrix_type;
	using _Mybase::_Interpolation_base;
	static constexpr dim = category::dimension;
	 

	MATRICE_HOST_FINL void _Coeff_impl();
};

MATRICE_ALGS_END