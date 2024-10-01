/**********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2022, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#pragma once

#include "util/_macros.h"
#include "private/_type_traits.h"
#include "private/_shape.hpp"

DGE_MATRICE_BEGIN
namespace xpr {
/// <summary>
/// \brief Provide a unified identity for expressions.
/// </summary>
struct __xpr__{};

/// <summary>
/// \breif Provide a unified identity for expreesion operators
/// </summary>
struct __op__ {};

template<class _Derived> 
class Base : public __xpr__ {
	using _Myt = Base;
	using _Mytraits = traits<_Derived>;
public:
	using derived_t = _Derived;
	using traits_t = _Mytraits;

	MATRICE_GLOBAL_FINL auto eval() const noexcept {
		typename _Mytraits::return_type _Res(derived().shape());
		for (auto idx = 0; idx < _Res.size(); ++idx) {
			_Res(idx) = derived()(idx);
		}
		return forward<decltype(_Res)>(_Res);
	}

	MATRICE_GLOBAL_FINL derived_t& derived() noexcept {
		return *((derived_t*)(this));
	}
	MATRICE_GLOBAL_FINL const derived_t& derived() const noexcept {
		return *((derived_t*)(this));
	}
};

template<class _Ty, class _Uy, 
	typename T = conditional_t<is_scalar_v<_Ty>, Scalar<_Ty>, _Ty>,
	typename U = conditional_t<is_scalar_v<_Uy>, Scalar<_Uy>, _Uy>>
static inline constexpr auto max_rows_v = max_integer_v<T::rows_at_compiletime, U::rows_at_compiletime>;
template<class _Ty, class _Uy, 
	typename T = conditional_t<is_scalar_v<_Ty>, Scalar<_Ty>, _Ty>,
	typename U = conditional_t<is_scalar_v<_Uy>, Scalar<_Uy>, _Uy>>
static inline constexpr auto max_cols_v = max_integer_v<T::cols_at_compiletime, U::cols_at_compiletime>;
template<class _Ty, class _Uy, 
	typename T = conditional_t<is_scalar_v<_Ty>, Scalar<_Ty>, _Ty>,
	typename U = conditional_t<is_scalar_v<_Uy>, Scalar<_Uy>, _Uy>>
static inline constexpr auto min_rows_v = min_integer_v<T::rows_at_compiletime, U::rows_at_compiletime>;
template<class _Ty, class _Uy, 
	typename T = conditional_t<is_scalar_v<_Ty>, Scalar<_Ty>, _Ty>,
	typename U = conditional_t<is_scalar_v<_Uy>, Scalar<_Uy>, _Uy>>
static inline constexpr auto min_cols_v = min_integer_v<T::cols_at_compiletime, U::cols_at_compiletime>;
template<typename T, typename U>
static inline constexpr auto common_rows_v = conditional_size_v<is_scalar_v<T> || is_scalar_v<U>,
	max_rows_v<T, U>, min_rows_v<T, U>>;
template<typename T, typename U>
static inline constexpr auto common_cols_v = conditional_size_v<is_scalar_v<T> || is_scalar_v<U>,
	max_cols_v<T, U>, min_cols_v<T, U>>;
}

/// <summary>
/// \brief Expression concept. 
/// _Xpr is a concept iff it's derived from class xpr::__xpr__.
/// </summary>
template<typename _Xpr>
concept Expr = std::is_base_of_v<xpr::__xpr__, _Xpr>;

/// <summary>
/// \brief Expression concept. 
/// _Xpr is a concept iff it's derived from class xpr::__xpr__.
/// </summary>
template<typename _Op>
concept Operator = std::is_base_of_v<xpr::__op__, _Op>;
template<typename _Op> concept Functor = Operator<_Op>;


DGE_MATRICE_END