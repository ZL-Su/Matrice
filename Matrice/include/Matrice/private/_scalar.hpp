/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library for
3D Vision and Photo-Mechanics.
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

#include "_type_traits.h"

DGE_MATRICE_BEGIN

template<typename _Ty = default_type>
class Scalar {
	static_assert(is_scalar_v<remove_all_t<_Ty>>, 
		"_Ty in Scalar must be a primitive scalar type.");
	using _Myt = Scalar;
public:
	enum { Size = 1, rows_at_compiletime = 1, cols_at_compiletime = 1 };
	using value_type = remove_all_t<_Ty>;
	using reference = value_type&;
	using pointer = value_type*;

	MATRICE_GLOBAL_INL Scalar() noexcept
		: _Myval(0) {
	}
	template<typename _Uy>
	MATRICE_GLOBAL_INL Scalar(const _Uy s) noexcept 
		: _Myval(s) {
	}

	MATRICE_GLOBAL_INL operator reference() const noexcept {
		return (_Myval);
	}
	MATRICE_GLOBAL_INL operator reference() noexcept {
		return (_Myval);
	}

	MATRICE_GLOBAL_INL operator pointer() const noexcept {
		return (&_Myval);
	}
	MATRICE_GLOBAL_INL operator pointer() noexcept {
		return (&_Myval);
	}

	template<typename _Uy>
	MATRICE_GLOBAL_INL _Myt& operator=(const _Uy arg) noexcept {
		if constexpr (is_pointer_v<_Uy>) _Myval = *arg;
		if constexpr (is_expression_v<_Uy>) _Myval = arg;
		return (*this);
	}

	MATRICE_GLOBAL_INL decltype(auto)operator()(auto)const noexcept {
		return (_Myval);
	}
	MATRICE_GLOBAL_INL decltype(auto)operator()(auto) noexcept {
		return (_Myval);
	}

	MATRICE_GLOBAL_INL constexpr decltype(auto)(size)() const noexcept {
		return 1;
	}

	/**
	 *\brief Compute scalar derivative.
	 */
	MATRICE_GLOBAL_INL const _Myt grad()const noexcept {
		return _Myt();
	}
	/**
	 *\brief Compute reciprocal of this scalar.
	 */
	MATRICE_GLOBAL_INL const _Myt inv()const noexcept {
		return _Myt(_Myval == 0 ? 0 : 1 / _Myval);
	}

private:
	value_type _Myval;
};

template<typename _Ty> struct is_scalar<Scalar<_Ty>> {
	static constexpr auto value = true;
};
template<typename _Ty> struct remove_all<Scalar<_Ty>> {
	using type = typename Scalar<_Ty>::value_type;
};
DGE_MATRICE_END