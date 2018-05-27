#pragma once
#include <type_traits>
#include "_macros.h"

MATRICE_NAMESPACE_BEGIN_

template<typename T1, typename T2, typename _Ret = typename std::decay<decltype(true ? T1() : T2())>::type>
constexpr _Ret max(const T1& a, const T2& b) { return a < b ? b : a; }
template<typename T1, typename T2, typename _Ret = typename std::decay<decltype(true ? T1() : T2())>::type>
constexpr _Ret min(const T1& a, const T2& b) { return a < b ? a : b; }
template<typename... _Args>
constexpr auto plus(_Args const&... args) { return (... + args); }
template<typename... _Args>
constexpr auto minus(_Args const&... args) { return (... - args); }
template<typename... _Args>
constexpr auto multiply(_Args const&... args) { return (... * args); }

_MATRICE_NAMESPACE_END
