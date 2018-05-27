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
#include<algorithm>
#include "_macros.h"

MATRICE_NAMESPACE_BEGIN_
using namespace std;

template<class _Fn, class _InIt, class _OutIt> inline
void transform(_Fn _Func, const _InIt _First, const _InIt _Last, size_t _Stride, _OutIt _Dest, size_t _Invpos = 1)
{
	auto _UFirst = _Unchecked(_First + _Stride - _Invpos);
	_DEBUG_RANGE(_UFirst, _Last);
	const auto _ULast = _Unchecked(_Last);
	auto _UDest = _Unchecked_n(_Dest, _Idl_distance<_InIt>(_UFirst, _ULast));
	for (; (_UFirst + _Invpos) != _ULast; _UFirst += _Stride, (void)++_UDest)
	{
		*_UDest = _Func(*_UFirst);
	}
	*(_UDest) = _Func(*_UFirst);
}
template<class _Fn, class _InIt, class _OutIt> inline
void transform(_Fn _Func, const _InIt _First, const _InIt _Last, _OutIt _Dest, size_t _Stride)
{
	_DEBUG_RANGE(_First, _Last);
	auto _UFirst = _Unchecked(_First);
	const auto _ULast = _Unchecked(_Last);
	auto _UDest = _Unchecked_n(_Dest, _Idl_distance<_InIt>(_UFirst, _ULast));
	for (; _UFirst != _ULast; ++_UFirst, (void)(_UDest+=_Stride))
	{
		*_UDest = _Func(*_UFirst);
	}
}
_MATRICE_NAMESPACE_END
