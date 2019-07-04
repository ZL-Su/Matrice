/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
**************************************************************************/
#pragma once

#include <type_traits>
#include "_macros.h"
#include "_std_wrapper.h"
#include "_exception.h"

DGE_MATRICE_BEGIN

static_assert(sizeof(void *) == 8, "MATRICE supports 64 bit only");

template<typename _Ty = long double>
MATRICE_GLOBAL_INL constexpr _Ty pi{ static_cast<_Ty>(3.14159265358979323846264338327950288419716939937510582097494459) };

namespace detail {
	template<typename _Ty> struct string_to_numval {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return (_Str); }
	};
	template<> struct string_to_numval<int> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stoi(_Str); }
	};
	template<> struct string_to_numval<long> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stol(_Str); }
	};
	template<> struct string_to_numval<float> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stof(_Str); }
	};
	template<> struct string_to_numval<double> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stod(_Str); }
	};
	template<> struct string_to_numval<long double> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stold(_Str); }
	};
	template<> struct string_to_numval<long long> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stoll(_Str); }
	};
	template<> struct string_to_numval<unsigned long> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stoul(_Str); }
	};
	template<> struct string_to_numval<unsigned long long> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stoull(_Str); }
	};
}

/**
 * \brief  cast a string to a user specified type T.
 * \example: 
		const auto val = dgelom::stonv<float>("1.0"); //val = 1.0f;
 */
template<typename T = std::string> MATRICE_HOST_FINL
constexpr T stonv(const std::string& _Str) noexcept {
#if (defined _DEBUG || MATRICE_DEBUG)
	DGELOM_CHECK(!_Str.empty(), "_Str should not be empty.")
#endif
	return detail::string_to_numval<T>::value(_Str); 
}
template<typename T = std::string> MATRICE_HOST_FINL
constexpr T cast_string_to(std::string&& _Str) noexcept {
#if (defined _DEBUG || MATRICE_DEBUG)
	DGELOM_CHECK(!_Str.empty(), "_Str should not be empty.")
#endif
		return detail::string_to_numval<T>::value(_Str);
}

/**
 * \brief  cast T-typed numeric to std::string
 * \example: 
		const auto str = dgelom::cast_to_string(3.14159); //str = "3.14159"
 */
template<typename T> MATRICE_HOST_INL
std::string cast_to_string(T _Val) noexcept {
	return std::to_string(_Val);
}

/**
 * \brief  append a T-typed element into tuple _Tpl
 */
template<typename T, typename... U> MATRICE_HOST_FINL
tuple<U..., T> tuple_append(const tuple<U...>& _Tpl, const T& _Val) {
	return std::tuple_cat(_Tpl, std::make_tuple(_Val));
}

/**
 * \packs the first _N element from _E into a tuple
 */
template<size_t _N> struct tuple_n {
	template<typename U> MATRICE_HOST_FINL static auto _(const U& _E) {
		return tuple_append(tuple_n<_N - 1>::_(_E), _E);
	}
	template<typename U> MATRICE_HOST_FINL static auto _(const U* _E) {
		return tuple_append(tuple_n<_N - 1>::_(_E), _E[_N]);
	}
	template<typename U, typename F> MATRICE_HOST_FINL static auto _(const U* _E, F&& _Op) {
		return tuple_append(tuple_n<_N - 1>::_(_E, _Op), _Op(_E[_N]));
	}
};
template<> struct tuple_n<0> {
	template<typename U> MATRICE_HOST_FINL static auto _(const U& _E) {
		return std::make_tuple(_E);
	}
	template<typename U> MATRICE_HOST_FINL static auto _(const U* _E) {
		return std::make_tuple(_E[0]);
	}
	template<typename U, typename F> MATRICE_HOST_FINL static auto _(const U* _E, F&& _Op) {
		return std::make_tuple(_Op(_E[0]));
	}
};

/**
 *\brief get size of built-in array
 */
template<typename _Ty, size_t _N>
MATRICE_GLOBAL_INL constexpr size_t size(_Ty(&)[_N]) noexcept { 
	return (_N); 
}
/**
 *\brief get size of std::array
 */
template<typename _Ty, size_t _N>
MATRICE_HOST_INL constexpr size_t size(std::array<_Ty, _N>&) noexcept {
	return (_N); 
}
/**
 *\brief get compile-time size of a matrix_<,,>
 */
template<template<typename,int,int>class _My, typename _Ty, int _M, int _N>
MATRICE_HOST_INL constexpr size_t size(const _My<_Ty, _M, _N>&) noexcept {
	return (_M*_N);
}
/**
 *\brief get size of container which must has method size()
 */
template<typename _Cont>
MATRICE_HOST_INL constexpr size_t size(const _Cont& _) noexcept {
	return (_.size()); 
}

/**
 *\brief call func[_Fn] in a lock guarded way.
 *\param [args] variadic argument(s) in function "func".
 */
template<typename _Mtx, typename _Fn, typename... _Args>
MATRICE_HOST_INL decltype(auto)locked_call(_Mtx& mtx, _Fn&& func, _Args...args){
	std::lock_guard<_Mtx> __guard__(mtx);
	return func(args...);
}

/**
 * \unroll a linear index "idx" to 2d indices [y, x].
 */
template<typename _Ity>
MATRICE_GLOBAL_INL decltype(auto) unroll_linear_index(_Ity idx, _Ity width) noexcept {
	const auto y = safe_div(idx, width);
	const auto x = idx - y * width;
	return std::make_tuple(y, x);
}

/**
 * \transform functor definitions
 */
struct transforms {
	template<typename _Ty> struct scale {
		using value_type = _Ty;
		template<typename _Uy>
		MATRICE_GLOBAL_INL scale(const _Uy& _Scale = _Uy(1)) : _Myscale(_Scale) {}
		MATRICE_GLOBAL_INL auto operator()(const value_type& _Val)const { return (_Myscale*_Val); }
		value_type _Myscale = 1.;
	};
	template<typename _Ty> struct clamp {
		using value_type = _Ty;
		template<typename _Uy>
		MATRICE_GLOBAL_INL clamp(const _Uy& _Lower, const _Uy& _Upper) : _Mylower(_Lower),_Myupper(_Upper) {}
		MATRICE_GLOBAL_INL auto operator()(const value_type& _Val)const { return min(max(_Val,_Mylower),_Myupper); }
		value_type _Mylower = std::numeric_limits<value_type>::min();
		value_type _Myupper = std::numeric_limits<value_type>::max();
	};
	template<typename _Ty> struct relu {
		using value_type = _Ty;
		MATRICE_GLOBAL_INL relu() {}
		MATRICE_GLOBAL_INL auto operator()(const value_type& _Val)const { return max(_Myzero,_Val); }
		value_type _Myzero = (value_type)(0);
	};
};
DGE_MATRICE_END
