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
***********************************************************************/
#pragma once
#include <cstdint>
#include <queue>
#include <random>
#include <tuple>
#include <memory>
#include <type_traits>

namespace dgelom {
	using std::uint8_t;
	using std::uint16_t;
	using std::uint32_t;
	using std::uint64_t;
	using std::size_t;
	using std::add_pointer_t;

	using std::shared_ptr;
	using std::unique_ptr;

	using std::tuple;
	using std::priority_queue;

	using std::mt19937;
	using std::uniform_int;
	using std::uniform_real;
	using std::normal_distribution;

	using diff_t = std::ptrdiff_t;

	template<typename _Ty>
	using initlist = std::initializer_list<_Ty>;
	template<typename _Ty> 
	using nested_initlist = initlist<initlist<_Ty>>;

	template<typename _Ty, typename _Uy = _Ty>
	using pair_t = std::pair<_Ty, _Uy>;

#if defined(_HAS_CXX14) || defined(_HAS_CXX17)
	using std::enable_if_t;
	using std::common_type_t;
	using std::is_integral_v;
	using std::is_arithmetic_v;
	//using std::is_scalar_v;
	using std::is_pointer_v;
	using std::is_class_v;
	using std::is_function_v;
#else
#define MATRICE_MAKE_TS_HELPER(TYPE) \
	template<typename _Ty> \
	inline constexpr auto is_##TYPE##_v \
	= std::is_##TYPE##<_Ty>::value;

	template<bool _Test, typename _Ty = void>
	using enable_if_t = typename std::enable_if<_Test, _Ty>::type;
	template<typename... _Ty>
	using common_type_t = typename std::common_type<_Ty...>::type;
	
	MATRICE_MAKE_TS_HELPER(integral);
	MATRICE_MAKE_TS_HELPER(arithmetic);
	//MATRICE_MAKE_TS_HELPER(scalar);
	MATRICE_MAKE_TS_HELPER(pointer);
	MATRICE_MAKE_TS_HELPER(class);
	MATRICE_MAKE_TS_HELPER(function);

#undef MATRICE_MAKE_TS_HELPER
#endif

	using std::get;
	using std::move;
	using std::forward;

#ifndef MATRICE_STD_NUMLIMITS
#define MATRICE_STD_NUMLIMITS std::numeric_limits<_Ty>
#endif
	// Variable Template for std::numeric_limits<_Ty>::min()
	template<typename _Ty>
	constexpr auto min_v = MATRICE_STD_NUMLIMITS::min();
	// Variable Template for std::numeric_limits<_Ty>::max()
	template<typename _Ty>
	constexpr auto max_v = MATRICE_STD_NUMLIMITS::max();
	// Variable Template for std::numeric_limits<_Ty>::epsilon()
	template<typename _Ty>
	constexpr auto epsilon_v = MATRICE_STD_NUMLIMITS::epsilon();
	// Variable Template for std::numeric_limits<_Ty>::infinity()
	template<typename _Ty>
	constexpr auto infinity_v = MATRICE_STD_NUMLIMITS::infinity();
	// Variable Template for std::numeric_limits<_Ty>::quiet_NaN()
	template<typename _Ty>
	constexpr auto quiet_nan_v = MATRICE_STD_NUMLIMITS::quiet_NaN();

#undef MATRICE_STD_NUMLIMITS
}
