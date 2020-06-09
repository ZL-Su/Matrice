/***********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
***********************************************************************/
#pragma once
#include <core/matrix.h>
#include "../similarity.h"

MATRICE_ALGS_BEGIN
template<typename T>
concept real_type = requires {is_floating_point_v<T>; };

/// <summary>
/// ZNCC fitness function
/// </summary>
/// <typeparam name="T"></typeparam>
template<typename T> struct zncc_fitness {
	size_t m_size = 31;
	T* m_ref = nullptr;

	/// <summary>
	/// Evaluate the fitness value.
	/// </summary>
	/// <param name="target">Subset to be matched against</param>
	/// <returns>Fitness value</returns>
	auto operator() (T* const target) const noexcept {
		Metric_<metric_fn::ZNCC, T> metric(m_ref, m_size);
		return metric.eval(target);
	}
};

_DETAIL_BEGIN

/// <summary>
/// Particle class for PSO algorithm
/// </summary>
/// <typeparam name="_Ty"> float or double </typeparam>
template<typename _Ty> 
requires real_type<_Ty> class _Particle 
{
public:
	using value_type = _Ty;
	using pos_type = Vec2_<value_type>;
	struct options {
		value_type max_x = 0, min_x = 0;
		value_type max_y = 0, min_y = 0;
		value_type max_v = 0, min_v = 0;
		value_type c_1 = 0, c_2 = 0;
	};

	explicit _Particle(const pos_type& pos) noexcept
		:_Mypos(pos) {
	}


private:
	options _Myopt;
	value_type _Myfs;
	pos_type _Mypos, _Mybest;
	pos_type& _Myglobal;

};

template<class _Type>
MATRICE_GLOBAL_INL auto _Optimize(const _Type& t) noexcept {

}
_DETAIL_END
MATRICE_ALGS_END