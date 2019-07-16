/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include <cstdint>
#include <queue>
#include <random>
#include <tuple>
#include <type_traits>

namespace dgelom {
	using std::uint8_t;
	using std::uint16_t;
	using std::uint32_t;
	using std::uint64_t;
	using std::size_t;
	
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
}
