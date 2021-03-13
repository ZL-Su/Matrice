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

#include "../spline.h"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Derived> class _Spline {
	using _Mytraits = traits<_Derived>;
	using _Myt = _Spline;
public:
	using value_type = typename _Mytraits::value_type;
	using matrix_type = Matrix<value_type>;

	template<typename... _Ty>
	MATRICE_HOST_INL auto operator()(_Ty&&..._Pos) const noexcept {
		return derived()._Value_at(_Pos...);
	}

	template<typename... _Ty>
	MATRICE_HOST_INL auto grad(_Ty&&..._Pos) const noexcept {
		return derived()._Grad_at(_Pos...);
	}

	MATRICE_GLOBAL_INL const _Derived& derived() const noexcept {
		return *static_cast<_Derived*>(this);
	}
	MATRICE_GLOBAL_INL _Derived& derived() noexcept {
		return *static_cast<_Derived*>(this);
	}

protected:
	shared_matrix_t<value_type> _Mycoef;
};

template<typename _Ty> 
class _Bspline<_Ty, bicubic_tag>
	: public _Spline<_Bspline<_Ty, bicubic_tag>> {
	using _Myt = _Bspline<_Ty, bicubic_tag>;
	using _Mybase = _Spline<_Myt>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::matrix_type;

	_Bspline(const matrix_type& _Data);

private:
	void _Precompute(const matrix_type& _Data);
};
_DETAIL_END
DGE_MATRICE_END