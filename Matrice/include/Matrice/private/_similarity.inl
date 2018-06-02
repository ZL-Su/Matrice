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
#include "../algs/similarity.h"

MATRICE_ALGS_BEGIN
template<typename T>
inline void Metric_<metric_fn::ZNCC, T>::_Init()
{
	_Size = _Radius << 1 | 1; _Size *= _Size;

	auto _Begin = _Data, _End = _Data + _Size;
	_Option[0] = reduce<value_t>(_Begin, _End) / _Size;
	_Option[1] = reduce<value_t>(_Begin, _End, [&](auto It)->value_t { *It -= _Option[0]; return (*It**It); });
}
template<typename T>
inline T Metric_<metric_fn::ZNCC, T>::eval(const pointer _oth) const
{
	value_t _Toption[2], _Score = value_t(0);
	_Toption[0] = reduce<value_t>(_oth, _oth + _Size) / _Size;
	auto _Op = [](auto _val)->auto{ return _val * _val; };

	auto _Begin = _Data;
	_Toption[1] = reduce<value_t>(_oth, _oth + _Size, [&](auto It)->value_t{
		auto _Diff = *It - _Toption[0];
		_Score += _Diff * *(_Begin++);
		return _Diff*_Diff; 
	});

	return (_Score/std::sqrt(_Option[1]*_Toption[1]));
}
MATRICE_ALGS_END
