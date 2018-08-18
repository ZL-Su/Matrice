/***************************************************************************
This file is part of Matrice, an effcient and elegant C++ library for SC.
      Copyright(C) 2018, Zhilong (Dgelom) Su (su-zl@seu.edu.cn), 
		                   all rights reserved.

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

#include <type_traits>
#include "../core/matrix.h"

MATRICE_ALGS_BEGIN
using std::size_t;
enum class metric_fn { L1, L2, ZNCC };

template<metric_fn Fn, typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type> class Metric_ {};
template<typename T> class Metric_<metric_fn::L1, T> final
{
	using value_t = T;
	using pointer = T * ;
public:
	MATRICE_GLOBAL_INL explicit Metric_(const pointer _data, int _n) noexcept
		: _Data(_data), _Size(_n) {};

	MATRICE_GLOBAL_INL value_t eval(const pointer _oth) const;

private:
	pointer _Data;
	size_t  _Size;
};
template<typename T> class Metric_<metric_fn::L2, T> final
{
	using value_t = T;
	using pointer = T * ;
public:
	MATRICE_GLOBAL_INL explicit Metric_(const pointer _data, int _n) noexcept
		: _Data(_data), _Size(_n) {};

	MATRICE_GLOBAL_INL value_t eval(const pointer _oth) const;
private:
	pointer _Data;
	size_t  _Size;
};
template<typename T> class Metric_<metric_fn::ZNCC, T> final
{
	using value_t = T;
	using pointer = T * ;
public:
	MATRICE_GLOBAL_INL explicit Metric_(const pointer _data, int _radius) noexcept
		: _Data(_data), _Radius(_radius) { _Init(); }

	MATRICE_GLOBAL_INL value_t eval(const pointer _oth) const;

private:
	MATRICE_GLOBAL_INL void _Init();
	pointer _Data;
	size_t  _Radius, _Size;
	value_t _Option[2];
};

template<typename T, size_t _M, size_t _N> struct SMBase
{
	using value_t = T;
	SMBase(size_t m, size_t n) : m_data(m, n) {}
	const value_t avg() const { return m_data.sum()/m_data.size(); }
private:
	types::Matrix_<value_t, _M, _N> m_data;
};
template<metric_fn _Mety, typename T, size_t N> 
class Similarity
{
	using Op = Metric_<_Mety, T>;
public:
	Similarity() {}
	Similarity(size_t m, size_t n) {}

private:
	Op _Impl;
};
MATRICE_ALGS_END
#include "../private/_similarity.inl"