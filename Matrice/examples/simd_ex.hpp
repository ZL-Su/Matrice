/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2023, Zhilong(Dgelom) Su, all rights reserved.

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
**********************************************************************/
#pragma once
#include <arch/simd.h>
DGE_MATRICE_BEGIN
namespace example {
	// SIMD vector type deduction, pass
	using simd_vector_f32 = simd::packet_traits<uint8_t>::type;
	using simd_vector_f64 = simd::packet_traits<double>::type;
	auto a = simd_vector_f32{};
	auto b = simd_vector_f64{};
	auto c = decltype(b){};
	
	// Alias of SIMD type.
	using packet_t = simd::Packet_<double>;

	constexpr index_t packed_size = packet_t::size;

	packet_t add(const packet_t& _Left, const packet_t _Right) {
		return _Left + _Right;
	}
}
DGE_MATRICE_END