
#pragma once
#include <type_traits>
#include "../private/_expr_type_traits.h"
#ifdef __use_mkl__
#include <mkl.h>
#else
#include <fkl.h>
#endif

namespace dgelom { namespace math {

template<typename T> inline constexpr
T _Det(const T* data, int n, typename std::enable_if<is_float32<T>::value>::type* = 0) 
{ 
	return fblas::_sdetm(data, n); 
}
template<typename T> inline constexpr
T _Det(const T* data, int n, typename std::enable_if<is_float64<T>::value>::type* = 0)
{
	return fblas::_ddetm(data, n);
}

} }
