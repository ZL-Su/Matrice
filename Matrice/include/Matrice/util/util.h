#pragma once

#include "_macros.h"
#include <type_traits>

#ifdef __enable_cuda__
#include <device_functions.h>
#include <thrust\complex.h>
#endif

MATRICE_NAMESPACE_BEGIN_

template<typename T1, typename T2, typename _Ret = typename std::decay<decltype(true ? T1() : T2())>::type> MATRICE_GLOBAL_FINL
constexpr _Ret max(const T1& a, const T2& b) { return a < b ? b : a; }
template<typename T1, typename T2, typename _Ret = typename std::decay<decltype(true ? T1() : T2())>::type> MATRICE_GLOBAL_FINL
constexpr _Ret min(const T1& a, const T2& b) { return a < b ? a : b; }

template<typename T, typename _Ret = T, typename = std::enable_if_t<std::is_arithmetic_v<T>>> 
MATRICE_GLOBAL_FINL constexpr _Ret sqrt(const T& x) {
	return std::sqrt(x);
}
template<typename T, typename = std::enable_if_t<std::is_class_v<T>>>
MATRICE_GLOBAL_FINL constexpr auto sqrt(const T& x) {
	return (x.sqrt());
}
template<typename T, typename _Ret = T> MATRICE_HOST_FINL
constexpr _Ret abs(const T& x) { return std::abs(x); }
template<typename T, typename _Ret = T> MATRICE_HOST_FINL
constexpr _Ret log(const T& x) { return std::log(x); }
template<typename T, typename _Ret = T> MATRICE_HOST_FINL
constexpr _Ret floor(const T& x) { return static_cast<_Ret>(std::floor(x)); }
template<typename T, typename _Ret = T> MATRICE_HOST_FINL
constexpr _Ret ceil(const T& x) { return static_cast<_Ret>(std::ceil(x)); }
template<typename T1, typename T2, typename _Ret = typename std::decay<decltype(true ? T1() : T2())>::type> MATRICE_GLOBAL_FINL
constexpr _Ret pow(const T1& x, const T2& y) { return std::pow(x,y); }

namespace details {
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
template<typename T> MATRICE_HOST_FINL
constexpr T stonv(const std::string& _Str) { return details::string_to_numval<T>::value(_Str);
}


#ifdef _HAS_CXX17
template<typename... _Args> MATRICE_GLOBAL_FINL
constexpr auto plus(_Args const&... args) { return (... + args); }
template<typename... _Args> MATRICE_GLOBAL_FINL
constexpr auto minus(_Args const&... args) { return (... - args); }
template<typename... _Args> MATRICE_GLOBAL_FINL
constexpr auto multiply(_Args const&... args) { return (... * args); }
#endif
template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
struct zero { static const constexpr T value = T(0);};

_MATRICE_NAMESPACE_END
