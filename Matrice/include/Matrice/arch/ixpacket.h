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
#include "_ixbase.h"
#include "./inl/_ixpackedops.inl.hpp"
#ifdef __AVX__
MATRICE_ARCH_BEGIN

template<typename T, int _Elems>
class Packet_ MATRICE_NONHERITABLE : public simd::simd_base_<T, _Elems>
{
	using myt = Packet_;
	using xbase_t = simd::simd_base_<T, _Elems>;
	using op_t = simd_op<T, _Elems>;
	using internal_t = typename xbase_t::internal_t;
	using initlist_t = typename xbase_t::initlist_t;
	using xbase_t::m_data;
public:
	using typename xbase_t::value_t;
	using pointer = value_t*;
	MATRICE_HOST_INL Packet_(const value_t _value) noexcept 
		: xbase_t(_value) {}
	MATRICE_HOST_INL Packet_(const pointer _data) noexcept
		: xbase_t(_data) {}
	MATRICE_HOST_INL Packet_(const initlist_t _list) noexcept
		: xbase_t(pointer(_list.begin())) {}

	MATRICE_HOST_INL auto& operator()() { return (m_data); }
	MATRICE_HOST_INL const auto& operator()() const { return (m_data); }
	MATRICE_HOST_INL myt& operator+ (const Packet_& _other);
	MATRICE_HOST_INL myt& operator- (const Packet_& _other);
	MATRICE_HOST_INL myt& operator* (const Packet_& _other);
	MATRICE_HOST_INL myt& operator/ (const Packet_& _other);
};
MATRICE_ARCH_END
#include "./inl/_ixpacket.inl"
#endif
