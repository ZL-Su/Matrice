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
#include "_splineinterp.h"

MATRICE_ALGS_BEGIN

template<typename _Ty, std::size_t _Opt> class Interpolation MATRICE_NONHERITABLE {};

template<typename _Ty>
class Interpolation<_Ty, INTERP|BICUBIC|BSPLINE> MATRICE_NONHERITABLE
{
public:
	using kernel_type = auto_interp_dispatcher_t<_Ty, INTERP | BICUBIC | BSPLINE>;
	using value_type = typename kernel_type::value_type;
	static constexpr auto options = kernel_type::options;

	template<typename... _Args>
	MATRICE_GLOBAL_FINL Interpolation(const _Args&... args)
		:m_op(std::make_shared<kernel_type>(args...)) {}

	MATRICE_GLOBAL_FINL auto& operator()() {
		return (*m_op);
	}
	MATRICE_GLOBAL_FINL const auto& operator()() const {
		return (*m_op);
	}

private:
	mutable std::shared_ptr<kernel_type> m_op;
};

MATRICE_ALGS_END