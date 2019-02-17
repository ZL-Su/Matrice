/***************************************************************************
This file is part of Matrice, an effcient and elegant C++ library for SC.
Copyright(C) 2018, Zhilong (Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.
****************************************************************************/
#pragma once
#include <xutility>
#include<algorithm>
#include "_macros.h"
#include "../arch/ixpacket.h"

MATRICE_NAMESPACE_BEGIN_
template<class _InIt, class _OutIt> MATRICE_HOST_FINL
void transform(const _InIt _First, const _InIt _Last, _OutIt _Dest) {
	static_cast<void>(_First == _Last);
	auto _UFirst = _First; size_t i = 0;
	for (; _UFirst != _Last; ++_UFirst, ++i) _Dest[i] = *_UFirst;
}
template<class _InIt, class _OutIt> MATRICE_HOST_FINL
void transform(const _InIt _First, const _InIt _Last, _OutIt _Dest, size_t _Dstep) {
	using namespace std;
	auto _UFirst = (_First);
	const auto _ULast = (_Last);
	auto _UDest = (_Dest);
	for (; _UFirst != _ULast; ++_UFirst, _UDest += _Dstep)
		*_UDest = *_UFirst;
}
template<class _Fn, class _InIt, class _OutIt> MATRICE_HOST_INL
void transform(_Fn _Func, const _InIt _First, const _InIt _Last, size_t _Stride, _OutIt _Dest, size_t _Invpos = 1) {
	using namespace std;
	auto _UFirst = (_First + _Stride - _Invpos);
	const auto _ULast = (_Last);
	auto _UDest = (_Dest);
	for (; (_UFirst + _Invpos) != _ULast; _UFirst += _Stride, (void)++_UDest)
	{
		*_UDest = _Func(*_UFirst);
	}
	*(_UDest) = _Func(*_UFirst);
}
template<class _Fn, class _InIt, class _OutIt> MATRICE_HOST_INL
void transform(_Fn _Func, const _InIt _First, const _InIt _Last, _OutIt _Dest, size_t _Stride) {
	using namespace std;
	auto _UFirst = (_First);
	const auto _ULast = (_Last);
	auto _UDest = (_Dest);
	for (; _UFirst != _ULast; ++_UFirst, (void)(_UDest+=_Stride))
	{
		*_UDest = _Func(*_UFirst);
	}
}

/**
 *\brief transform a forward container to another
 *\param [_Src, _Dst] containers with forward iterator
 */
template<typename _Fwd1, typename _Fwd2> MATRICE_HOST_INL 
auto transform(const _Fwd1& _Src, _Fwd2& _Dst)->_Fwd2& {
	transform(_Src.begin(), _Src.end(), _Dst.begin());
	return (_Dst);
}

/**
 *\brief Process container entry with a unary functor.
 *\param [_Cont] a container with forward iterator
 *\param [_Func] a unary functor to be applied to every dereferenced iterator
 */
template<typename _Fwdty, typename _Fn, 
	typename = std::enable_if_t<std::is_class_v<_Fwdty>>>
MATRICE_HOST_FINL void for_each(_Fwdty& _Cont, _Fn&& _Func) {
	std::for_each(_Cont.begin(), _Cont.end(), _Func);
}

/**
 *\brief Fill a container
 *\param [_Cont] container with forward iterator
 *\param [_Val] value to be filled into _Cont
 */
template<typename _Fwdty, typename _T, 
	typename = std::enable_if_t<std::is_class_v<_Fwdty>>>
MATRICE_HOST_FINL void fill(_Fwdty& _Cont, const _T& _Val) { 
	std::fill(_Cont.begin(), _Cont.end(), _Val); 
}
/**
 *\brief Fill a range [_First, _Last) with a fixed interval step
 *\param [_First, _Last] the first and last iterator
 *\param [_Stride] interval step
 *\param [_Val] value to be filled into _Cont
 */
template<typename _FwdIt, typename _T, 
	typename = std::enable_if_t<std::is_scalar_v<_T>>>
MATRICE_HOST_FINL void fill(_FwdIt _First, _FwdIt _Last, size_t _Stride, const _T& _Val) {
	auto _UFirst = (_First);
	const auto _ULast = (_Last);
	for (; _UFirst < _ULast; _UFirst += _Stride) *_UFirst = _Val;
}

// \return: sum of the range [_First, _Last)
template<typename _InIt>
MATRICE_GLOBAL_INL auto reduce(_InIt _First, _InIt _Last) {
	static_cast<void>(_First == _Last);
	typename std::pointer_traits<_InIt>::element_type _Ret = 0;
#if MATRICE_SIMD_ARCH==MATRICE_SIMD_AVX
	using packed_t = typename simd::template Packet_<decltype(_Ret)>;
	decltype(auto) _Size = std::distance(_First, _Last);
	decltype(auto) _Step = packed_t::size << 1;
	decltype(auto) _N = _Size / _Step;
	for (auto i = 0, j = 0; i < _N; j = (++i)*_Step) {
		const auto _Pos = _First + j;
		_Ret += (packed_t(_Pos)+packed_t(_Pos+packed_t::size)).reduce();
	}
	for (_First += _N * _Step; _First != _Last; ++_First) _Ret += *_First;
#else
	for (; _First != _Last; ++_First) _Ret += *_First;
#endif
	return (_Ret);
}

// \return: sum of range [_First, _Last) with step := _Stride
template<typename _InIt, 
	typename = std::enable_if_t<std::is_pointer_v<_InIt>>>
MATRICE_GLOBAL_INL auto reduce(_InIt _First, _InIt _Last, index_t _Stride) {
	static_cast<void>(_First == _Last);
	remove_reference_t<decltype(*_First)> _Ret = 0;
	for (; _First < _Last; _First += _Stride) _Ret += *(_First);
	return (_Ret);
}
template<typename _Ty, typename _Op,  typename _InIt = _Ty*>
MATRICE_GLOBAL_INL _Ty reduce(_InIt _First, _InIt _Last, _Op _op)
{
	static_cast<void>(_First == _Last);
	_Ty _Ret = 0;
	for (; _First != _Last; ++_First) _Ret += _op(*_First);
	return (_Ret);
}
template<typename _InIt, typename _Op, typename = std::enable_if_t<std::is_pointer_v<_InIt>&&std::is_function_v<_Op>>>
MATRICE_GLOBAL_INL auto reduce(_InIt _First, _InIt _Last, _Op _op)
{
	static_cast<void>(_First == _Last);
	typename remove_reference<decltype(_First[0])>::type _Ret = 0;
	for (; _First != _Last; ++_First) _Ret += _op(*_First);
	return (_Ret);
}
template<template<typename> class _Op, typename _InIt, typename _Scalar, typename _Func>
MATRICE_GLOBAL_INL auto reduce(_InIt _First, _InIt _Last, _Scalar _Value, _Func _Fn, _Op<_Scalar> _op = _Op<_Scalar>())
{
	static_assert(std::is_arithmetic_v<_Scalar>, "Oops, template parameter '_Scalar' is illegal!");
	static_cast<void>(_First == _Last);
	typename remove_reference<decltype(_First[0])>::type _Ret = 0;
	for (; _First != _Last; ++_First) _Ret += _Fn(_op(*_First, _Value));
	return (_Ret);
}
_DETAIL_BEGIN
template<std::size_t _N>  struct _Reduce_n {
	template<typename _Ty> static
	MATRICE_GLOBAL_INL auto value(const _Ty* _Data[[_N]]) {
		return (_Reduce_n<_N - 1>::value(_Data) + _Data[_N]);
	}
};
template<> struct _Reduce_n<0> {
	template<typename _Ty> static
	MATRICE_GLOBAL_INL auto value(const _Ty* _Data[[]]) { return (_Data[0]); }
};

template<std::size_t _N>  struct _Powers_n {
	template<typename _Ty> static
	MATRICE_GLOBAL_INL auto value(const _Ty _Val) {
		return (_Powers_n<_N - 1>::value(_Val)*_Val);
	}
};
template<> struct _Powers_n<0> {
	template<typename _Ty> static
		MATRICE_GLOBAL_INL auto value(const _Ty _Val) { return _Ty(1); }
};
_DETAIL_END

template<std::size_t _Size = 1> using reduce_n_t = detail::_Reduce_n<_Size - 1>;
template<std::size_t _Expo = 1> using powers_n_t = detail::_Powers_n<_Expo>;
DGE_MATRICE_END
