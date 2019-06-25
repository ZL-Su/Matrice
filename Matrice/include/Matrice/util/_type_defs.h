/***************************************************************************
This file is part of Matrice, an effcient and elegant C++ library for SC.
Copyright(C) 2018-2019, Zhilong (Dgelom) Su, all rights reserved.

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
#include <functional>

namespace dgelom 
{
	using default_type = double;
	using index_t = long long;
	enum solver_type {
		AUTO = 0x0000, LUF = 0x0001, CHD = 0x0002, QRD = 0x0003,
		SVD = 0x0004, ITER = 0x0005, SPAR = 0x0006, GLS = 0x0007
	};
	enum {
		expr = 0x0099,
		rmaj = 101,
		cmaj = 102,
		gene = 103,
		symm = 104,
		diag = 105,
		band = 106,
		utri = 107,
		ltri = 108,
		spar = 109,
	};
	enum { SOBEL = 0, BSPL3 = 3, BSPL5 = 5, BSPL7 = 7 };
	enum class axis {all = 0x00, x = 0x01, y = 0x02, z = 0x02 };

	/**
	 * \enum type for matrix transpose tag
	 * ttag::A for automatic checking; 
	 * ttag::Y for transpose;
	 * ttag::N for no transpose
	 */
	enum class ttag {A = 000, Y = 111, N = 112};
	using transp = ttag;

	template<typename _Ty>
	using deleter_type = std::function<void(std::add_pointer_t<_Ty>)>;

	template<typename _Ty, typename _Tag> struct index {};
}