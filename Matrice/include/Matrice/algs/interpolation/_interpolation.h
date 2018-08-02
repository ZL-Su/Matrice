/*********************************************************************
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
***********************************************************************/
#pragma once
#include <memory>
#include "_bicubic_spline_interp.h"

MATRICE_ALGS_BEGIN

template<typename _Ty, size_t _Options>
class Interpolation MATRICE_NONHERITABLE
{
	using Op_t = interpolation_traits_t<_Ty, _Options>;
public:
	template<typename... _Args> 
	MATRICE_GLOBAL_FINL Interpolation(const _Args&... args) 
		:m_op(std::make_unique<Op_t>(args...)) {}

private:
	std::unique_ptr<Op_t> m_op;
};

template<typename _Ty>
class Interpolation<_Ty, INTERP|BICUBIC|BSPLINE> MATRICE_NONHERITABLE
{
	using Op_t = interpolation_traits_t<_Ty,INTERP|BICUBIC|BSPLINE>;
public:
	template<typename... _Args>
	MATRICE_GLOBAL_FINL Interpolation(const _Args&... args)
		:m_op(std::make_unique<Op_t>(args...)) {}

private:
	std::unique_ptr<Op_t> m_op;
};

MATRICE_ALGS_END