/**************************************************************************
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
**************************************************************************/
#pragma once

#include "_base.h"

MATRICE_ALGS_BEGIN

template<typename _Ty, std::size_t _Opt> class _Spline_interpolation {};

template<typename _Ty>
class _Spline_interpolation<_Ty, _BICBSPL>
	: public _Interpolation_base<_Spline_interpolation<_Ty, _BICBSPL>>
{
	using _Myt = _Spline_interpolation;
	using _Mybase = _Interpolation_base<_Myt>;
public:
	using typename _Mybase::matrix_type;
	using typename _Mybase::value_type;
	using typename _Mybase::point_type;
	static constexpr auto _N1 = 4, _N2 = 3;

	_Spline_interpolation(const matrix_type& _Data) : _Mybase(_Data) {}
	
	MATRICE_HOST_FINL auto _Value_at(const point_type& _Pos) const {
		const auto _Ix = floor<int>(_Pos.x), _Iy = floor<int>(_Pos.y);
		const auto _Dx = _Pos.x - _Ix, _Dy = _Pos.y - _Iy;

		Matrix_<value_type, _N1, 1> _Dxs{ 1., _Dx, _Dx*_Dx, _Dx*_Dx*_Dx };
		Matrix_<value_type, 1, _N1> _Dys{ 1., _Dy, _Dy*_Dy, _Dy*_Dy*_Dy };

		const auto _Coeff = _Mycoeff.block(_Ix-1, _Ix+3, _Iy-1, _Iy+3).eval<4,4>();
		auto _Temp = _L.t().mul(_Coeff.mul(_L).eval()).eval();

		return (_Dys.mul(_Temp.mul(_Dxs).eval()))(0);
	}
	MATRICE_HOST_FINL auto _Gradx_at(const point_type& _Pos) const {
		const auto _Ix = floor<int>(_Pos.x), _Iy = floor<int>(_Pos.y);
		const auto _Dx = _Pos.x - _Ix, _Dy = _Pos.y - _Iy;

		Matrix_<value_type, _N2, 1> _Dxs{ 1., _Dx, _Dx*_Dx };
		Matrix_<value_type, 1, _N1> _Dys{ 1., _Dy, _Dy*_Dy, _Dy*_Dy*_Dy };

		const auto _Coeff = _Mycoeff.block(_Ix-1, _Ix+3, _Iy-1, _Iy+3).eval<4,4>();
		auto _Temp = _L.t().mul(_Coeff.mul(_R).eval()).eval();

		return (_Dys.mul(_Temp.mul(_Dxs).eval()))(0);
	}
	MATRICE_HOST_FINL auto _Grady_at(const point_type& _Pos) const {
		const auto _Ix = floor<int>(_Pos.x), _Iy = floor<int>(_Pos.y);
		const auto _Dx = _Pos.x - _Ix, _Dy = _Pos.y - _Iy;

		Matrix_<value_type, _N1, 1> _Dxs{ 1., _Dx, _Dx*_Dx, _Dx*_Dx*_Dx };
		Matrix_<value_type, 1, _N2> _Dys{ 1., _Dy, _Dy*_Dy };

		const auto _Coeff = _Mycoeff.block(_Ix-1, _Ix+3, _Iy-1, _Iy+3).eval<4,4>();
		auto _Temp = _R.t().mul(_Coeff.mul(_L).eval()).eval();

		return (_Dys.mul(_Temp.mul(_Dxs).eval()))(0);
	}

	MATRICE_HOST_INL void _Coeff_impl();

private:
	using _Mybase::_Mycoeff;
	const Matrix_<value_type, 4, 4> _L{ 1,-3,3,-1,4,0,-6,3,1,3,3,-3,0,0,0,1 };
	const Matrix_<value_type, 4, 3> _R{ -3,6,-3,0,-12,9,3,6,-9,0,0,3 };
};

template<typename _Ty>
class _Spline_interpolation<_Ty, _BIQNSPL>
	: public _Interpolation_base<_Spline_interpolation<_Ty, _BIQNSPL>>
{
	using _Myt = _Spline_interpolation;
	using _Mybase = _Interpolation_base<_Myt>;
public:
	using typename _Mybase::matrix_type;
	using typename _Mybase::value_type;
	using typename _Mybase::point_type;

	_Spline_interpolation(const matrix_type& _Data) : _Mybase(_Data) {}

	MATRICE_HOST_INL void _Coeff_impl();

	MATRICE_HOST_INL auto _Value_at(const point_type& _Pos) {

	}

private:
	using _Mybase::_Mycoeff;
	const Matrix_<value_type, 4, 4> _A{ 1, -3, 3, -1, 4, 0, -6, 3, 1, 3, 3, -3, 0, 0, 0, 1 };
	const Matrix_<value_type, 4, 3> _B{ -3, 6, -3, 0, -12, 9, 3, 6, -9, 0, 0, 3 };
};

template<typename _Ty>
class _Spline_interpolation<_Ty, _BISPSPL>
	: public _Interpolation_base<_Spline_interpolation<_Ty, _BISPSPL>>
{
	using _Myt = _Spline_interpolation;
	using _Mybase = _Interpolation_base<_Myt>;
public:
	using typename _Mybase::matrix_type;
	using typename _Mybase::value_type;
	using typename _Mybase::point_type;

	_Spline_interpolation(const matrix_type& _Data) : _Mybase(_Data) {}

	MATRICE_HOST_INL void _Coeff_impl();

	MATRICE_HOST_FINL auto _Value_at(const point_type& _Pos) const {
		const auto _Ix = static_cast<int>(_Pos.x), _Iy = static_cast<int>(_Pos.y);
		const auto _Dx = _Pos.x - _Ix, _Dy = _Pos.y - _Iy;
	}

private:
	using _Mybase::_Mycoeff;
	const Matrix_<value_type, 4, 4> _A{ 1, -3, 3, -1, 4, 0, -6, 3, 1, 3, 3, -3, 0, 0, 0, 1 };
	const Matrix_<value_type, 4, 3> _B{ -3, 6, -3, 0, -12, 9, 3, 6, -9, 0, 0, 3 };
};

///<brief> Class BicubicSplineInterp will be deprecated </brief>
template<typename _Ty, std::size_t _Opt = _BICBSPL>
class BicubicSplineInterp : public InterpBase_<_Ty, BicubicSplineInterp<_Ty, _Opt>>
{
	using base_t = InterpBase_<_Ty, BicubicSplineInterp<_Ty, _Opt>>;
public:
	enum { option = _Opt };
	using typename base_t::value_t;
	using typename base_t::matrix_t;

	MATRICE_GLOBAL_FINL BicubicSplineInterp(const matrix_t& _Data)
		:base_t(_Data) {}

	MATRICE_HOST_INL void _Bspline_coeff(const matrix_t& _Data);

private:
	using base_t::m_coeff;
	const types::Matrix_<value_t, 4, 4> m_icoef{1, -3, 3, -1, 4, 0, -6, 3, 1, 3, 3, -3, 0, 0, 0, 1};
	const types::Matrix_<value_t, 4, 3> m_gcoef{-3, 6, -3, 0, -12, 9, 3, 6, -9, 0, 0, 3};
};
MATRICE_ALGS_END