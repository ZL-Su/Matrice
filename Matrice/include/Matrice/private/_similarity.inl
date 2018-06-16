/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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
**************************************************************************/
#pragma once
#include <algorithm>
#include "../algs/similarity.h"
#include "../arch/ixpacket.h"

MATRICE_ALGS_BEGIN
template<typename T> MATRICE_GLOBAL_INL
T Metric_<metric_fn::L1, T>::eval(const pointer _Othr) const
{
	value_t _Score = value_t(0);
#ifdef __AVX__
	using Packet = simd::Packet_<value_t, 4>;
	const auto _N = _Size / Packet::size;
	for (auto i = 0, j = 0; i < _N; j = (++i)*Packet::size) {
		_Score += (Packet(_Data + j) - Packet(_Othr + j)).abs().reduce();
	}
	for (auto i = _N * Packet::size; i < _Size; ++i) {
		_Score += std::abs(_Data[i] - _Othr[i]);
	}
#else
	types::Matrix<value_t> _Left(_Size, 1, _Data), _Right(_Size, 1, _Othr);
	decltype(auto) _Diff = _Left - _Right;
	for (auto i = 0; i < _Size; ++i) {
		_Score += std::abs(_Diff(i));
	}
#endif
	return (_Score);
}

template<typename T> MATRICE_GLOBAL_INL
T Metric_<metric_fn::L2, T>::eval(const pointer _Othr) const
{
	value_t _Score = value_t(0);
#ifdef __AVX__
	using Packet = simd::Packet_<value_t, 4>;
	const auto _N = _Size / Packet::size;
	for (auto i = 0, j = 0; i < _N; j = (++i)*Packet::size) {
		decltype(auto) _Diff = Packet(_Data + j) - Packet(_Othr + j);
		_Score += (_Diff*_Diff).reduce();
	}
	for (auto i = _N * Packet::size; i < _Size; ++i) {
		_Score += (_Data[i] - _Othr[i])*(_Data[i] - _Othr[i]);
	}
#else
	types::Matrix<value_t> _Left(_Size, 1, _Data), _Right(_Size, 1, _Othr);
	auto _Diff = _Left - _Right;
	for (auto i = 0; i < _Size; ++i) {
		decltype(auto) _Val = _Diff(i);
		_Score += _Val*_Val;
	}
#endif
	return (std::sqrt(_Score));
}

template<typename T> MATRICE_GLOBAL_INL
void Metric_<metric_fn::ZNCC, T>::_Init()
{
	_Size = _Radius << 1 | 1; _Size *= _Size;

	auto _Begin = _Data, _End = _Data + _Size;
	_Option[0] = reduce<value_t>(_Begin, _End) / _Size;
	_Option[1] = reduce<value_t>(_Begin, _End, [&](auto It)->value_t { *It -= _Option[0]; return (*It**It); });
}
template<typename T> MATRICE_GLOBAL_INL
T Metric_<metric_fn::ZNCC, T>::eval(const pointer _Othr) const
{
	value_t _Toption[2], _Score = value_t(0);
	_Toption[0] = reduce<value_t>(_Othr, _Othr + _Size) / _Size;
	auto _Op = [](auto _val)->auto{ return _val * _val; };

	auto _Begin = _Data;
	_Toption[1] = reduce<value_t>(_Othr, _Othr + _Size, [&](auto It)->value_t{
		auto _Diff = *It - _Toption[0];
		_Score += _Diff * *(_Begin++);
		return _Diff*_Diff; 
	});

	return (_Score/std::sqrt(_Option[1]*_Toption[1]));
}
MATRICE_ALGS_END
