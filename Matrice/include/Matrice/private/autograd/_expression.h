/*********************************************************************
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
*********************************************************************/
#pragma once
#include "../_matrix_exp.hpp"

///<brief> construct computational expression system </brief>
///<method> exp::forward(): for expression evaluation </method>
///<method> exp::backward(): for derivatives computation based on chain rule </method>

DGE_MATRICE_BEGIN namespace exp {
_DETAIL_BEGIN

template<typename _Derived> class _Expression_base
{
	using _Myt = _Expression_base;
	using _Mydt = _Derived;
public:
	_Expression_base() {};
	~_Expression_base() {};

	/**
	 * \return value of an expresstion
	 */
	MATRICE_GLOBAL_INL auto forward() const {
		return (static_cast<const _Mydt*>(this)->operator());
	}
	/**
	 * \return gradient of an expresstion
	 */
	MATRICE_GLOBAL_INL auto backward() const {

	}
private:

};

_DETAIL_END
} DGE_MATRICE_END