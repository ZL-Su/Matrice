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

template<typename T, typename _Ret = T> MATRICE_GLOBAL_FINL
constexpr _Ret sqrt(const T& x) { return std::sqrt(x); }
template<typename T, typename _Ret = T> MATRICE_GLOBAL_FINL
constexpr _Ret abs(const T& x) { return std::abs(x); }
template<typename T, typename _Ret = T> MATRICE_GLOBAL_FINL
constexpr _Ret log(const T& x) { return std::log(x); }
template<typename T1, typename T2, typename _Ret = typename std::decay<decltype(true ? T1() : T2())>::type> MATRICE_GLOBAL_FINL
constexpr _Ret pow(const T1& x, const T2& y) { return std::pow(x,y); }

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
