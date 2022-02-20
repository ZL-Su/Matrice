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

MATRICE_ALG_BEGIN(optim)
_DETAIL_BEGIN

struct _Constant_tag {};
struct _Linear_decrease_tag {};

template<typename _Ty, class _Tag> class _Inertia_weight {
	static_assert(true, "Unsupported _Tag in _Inertia_weight<_Tag>");
};

/// <summary>
/// Particle class for PSO algorithm
/// </summary>
/// <typeparam name="_Ty"> float or double </typeparam>
template<typename _Ty, size_t _Dim> class _Particle{};
template<typename _Ty>
class _Particle<_Ty, 2> 
{
public:
	using value_type = _Ty;
	using pos_type = Vec2_<value_type>;
	struct options {
		value_type max_x = 0, min_x = 0; //max and min x-position
		value_type max_y = 0, min_y = 0; //max and min y-position
		value_type max_v = 0, min_v = 0; //max and min velocity
		value_type c_1 = 0, c_2 = 0;
	};

	explicit _Particle(const pos_type& pos, const options& opt) noexcept
		:_Mypos(pos), _Myopt(opt) {
	}

	template<class _Fty>
	auto eval_fval(_Fty&& fn) const {
		return _Myfval = fn(_Mypos);
	}


private:
	const options& _Myopt;
	pos_type _Mypos, _Mybest;
	mutable value_type _Myfval;
	mutable pos_type& _Mygbest;

};

template<class _Type> class _Optimizer {

};

template<typename _Ty> 
class _Inertia_weight<_Ty, _Constant_tag> {
public:
	using value_type = _Ty;
	MATRICE_HOST_INL _Inertia_weight() noexcept = default;
	MATRICE_HOST_INL _Inertia_weight(const value_type w) noexcept 
		: _Myval(w) {
	}
	
	MATRICE_HOST_FINL auto operator()(size_t, size_t) const noexcept {
		return _Myval;
	}

private:
	value_type _Myval{ 1 };
};
template<typename _Ty> 
class _Inertia_weight<_Ty, _Linear_decrease_tag> {
public:
	using value_type = _Ty;
	MATRICE_HOST_INL _Inertia_weight() noexcept = default;
	MATRICE_HOST_INL _Inertia_weight(const value_type min, const value_type max)noexcept 
		: _Mymin(min), _Mymax(max) {
	}
	
	MATRICE_HOST_FINL auto operator()(size_t it, size_t max_it) const noexcept {
		return _Mymin + (_Mymax - _Mymin) * it / value_type(max_it);
	}

private:
	value_type _Mymin{ 0.4 };
	value_type _Mymax{ 0.9 };
};

_DETAIL_END
MATRICE_ALG_END(pso)