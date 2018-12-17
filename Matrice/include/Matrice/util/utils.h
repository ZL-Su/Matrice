#pragma once

#include "_macros.h"
#include <type_traits>

#if (defined __enable_cuda__ && !defined __disable_cuda__)
#include <device_functions.h>
#include <thrust\complex.h>
#endif

DGE_MATRICE_BEGIN
template<typename T1, typename T2, typename _Ret = std::common_type_t<T1, T2>, typename = std::enable_if_t<std::is_arithmetic_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret add(const T1& a, const T2& b) { return a + b; }
template<typename T1, typename T2, typename _Ret = std::common_type_t<T1, T2>, typename = std::enable_if_t<std::is_arithmetic_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret sub(const T1& a, const T2& b) { return a - b; }
template<typename T1, typename T2, typename _Ret = std::common_type_t<T1, T2>, typename = std::enable_if_t<std::is_arithmetic_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret mul(const T1& a, const T2& b) { return a * b; }
template<typename T1, typename T2, typename _Ret = std::common_type_t<T1, T2>, typename = std::enable_if_t<std::is_arithmetic_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret div(const T1& a, const T2& b) { return a / b; }
template<typename T1, typename T2, typename _Ret = std::common_type_t<T1, T2>, typename = std::enable_if_t<std::is_arithmetic_v<_Ret>>> 
MATRICE_GLOBAL_FINL constexpr _Ret max(const T1& a, const T2& b) { return a < b ? b : a; }
template<typename T1, typename T2, typename _Ret = std::common_type_t<T1, T2>, typename = std::enable_if_t<std::is_arithmetic_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret min(const T1& a, const T2& b) { return a < b ? a : b; }
template<typename T, typename _Ret = T, typename = std::enable_if_t<std::is_arithmetic_v<T>>> 
MATRICE_GLOBAL_FINL constexpr _Ret sqrt(const T& x) { return std::sqrt(x); }
template<typename T, typename _Ret = T, typename = std::enable_if_t<std::is_arithmetic_v<T>>> 
MATRICE_HOST_FINL constexpr _Ret abs(const T& x) { return std::abs(x); }
template<typename T, typename _Ret = T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
MATRICE_HOST_FINL constexpr _Ret exp(const T& x) { return std::exp(x); }
template<typename T, typename _Ret = T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
MATRICE_HOST_FINL constexpr _Ret log(const T& x) { return std::log(x); }
template<typename T, typename _Ret = T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
MATRICE_HOST_FINL constexpr _Ret log2(const T& x) { return std::log2(x); }
template<typename T, typename _Ret = T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
MATRICE_HOST_FINL constexpr _Ret log10(const T& x) { return std::log10(x); }
template<typename T, typename _Ret = T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
MATRICE_HOST_FINL constexpr _Ret floor(const T& x) { return static_cast<_Ret>(std::floor(x)); }
template<typename T, typename _Ret = T, typename = std::enable_if_t<std::is_arithmetic_v<T>>> 
MATRICE_HOST_FINL constexpr _Ret ceil(const T& x) { return static_cast<_Ret>(std::ceil(x)); }
template<typename T1, typename T2, typename _Ret = std::common_type_t<T1, T2>, typename = std::enable_if_t<std::is_arithmetic_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret pow(const T1& x, const T2& y) { return std::pow(x,y); }

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
 * \Convert a string to a user specified (numerical) value.
 * Example: auto _Ret = dgelom::stonv<float>("1.0");
 */
template<typename T = std::string> MATRICE_HOST_FINL
constexpr T stonv(const std::string& _Str) { return detail::string_to_numval<T>::value(_Str); }

#ifdef _HAS_CXX17
//template<typename... _Args> MATRICE_GLOBAL_FINL
//constexpr auto plus(_Args const&... args) { return (... + args); }
//template<typename... _Args> MATRICE_GLOBAL_FINL
//constexpr auto minus(_Args const&... args) { return (... - args); }
//template<typename... _Args> MATRICE_GLOBAL_FINL
//constexpr auto multiply(_Args const&... args) { return (... * args); }
///**
// * \is _Value satisfies anyone of given _Conditions: x == _Cond1 || ... || x = _Condn
// */
//template<typename T, typename... Ts>
//MATRICE_GLOBAL_INL constexpr bool is_eqto_one_in(T&& _Value, Ts&&... _Conditions) {
//	return ((_Conditions == _Value) || ...);
//}
///**
// * \is _Value satisfies all of given _Conditions: x == _Cond1 && ... && x = _Condn
// */
//template<typename T, typename... Ts>
//MATRICE_GLOBAL_INL constexpr bool is_eqto_all_in(T&& _Value, Ts&&... _Conditions) {
//	return ((_Conditions == _Value) && ...);
//}
#endif

/**
 * \get T-typed zero value
 */
template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
struct zero { static const constexpr T value = T(0);};
template<typename T> MATRICE_GLOBAL_INL constexpr auto zero_v = zero<T>::value;

/**
 * \append a T-typed element into tuple _Tpl
 */
template<typename T, typename... U> MATRICE_HOST_FINL
std::tuple<U..., T> tuple_append(const std::tuple<U...>& _Tpl, const T& _Val) {
	return std::tuple_cat(_Tpl, std::make_tuple(_Val));
}

/**
 * \packs the first _N element from _E into a tuple
 */
template<std::size_t _N> struct tuple_n {
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

DGE_MATRICE_END
