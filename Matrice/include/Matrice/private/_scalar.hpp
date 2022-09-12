/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library for
3D Vision and Photo-Mechanics.
Copyright(C) 2018-2022, Zhilong(Dgelom) Su, all rights reserved.

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
		"_Ty in Scalar<_Ty> must be a primitive scalar type.");
	using _Myt = Scalar;
public:
	enum { Size = 1, rows_at_compiletime = 1, cols_at_compiletime = 1 };
	using value_type = remove_all_t<_Ty>;
	using value_t = value_type;
	using reference = value_type&;
	using pointer = value_type*;

	MATRICE_GLOBAL_INL Scalar() noexcept
		: _Myval(0) {
	}
	MATRICE_GLOBAL_INL Scalar(shape_t<3>) noexcept
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

	MATRICE_GLOBAL_FINL decltype(auto)operator()(size_t)const noexcept {
		return (_Myval);
	}
	MATRICE_GLOBAL_FINL decltype(auto)operator()(size_t) noexcept {
		return (_Myval);
	}

	MATRICE_GLOBAL_INL constexpr decltype(auto)(size)() const noexcept {
		return 1;
	}
	MATRICE_GLOBAL_INL constexpr decltype(auto)(cols)() const noexcept {
		return 1;
	}
	MATRICE_GLOBAL_INL constexpr decltype(auto)(rows)() const noexcept {
		return 1;
	}
	MATRICE_GLOBAL_INL constexpr auto(shape)() const noexcept {
		return shape_t<3>{1, 1, 1};
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

	MATRICE_GLOBAL_INL decltype(auto)sum() const noexcept {
		return (_Myval);
	}
	MATRICE_GLOBAL_INL decltype(auto)reduce() const noexcept {
		return (_Myval);
	}

private:
	value_type _Myval;
};

template<> struct traits<Scalar<float>> {
	using value_type = float;
	using value_t = value_type;
};
template<> struct traits<Scalar<double>> {
	using value_type = double;
	using value_t = value_type;
};

template<typename _Ty> struct is_scalar<Scalar<_Ty>> {
	static constexpr auto value = true;
};
template<typename _Ty> struct remove_all<Scalar<_Ty>> {
	using type = typename Scalar<_Ty>::value_type;
};
DGE_MATRICE_END