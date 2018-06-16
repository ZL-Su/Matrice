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

///<brief>
// \template param: T is a scalar type, _Elems is the number of packed elements
///</brief>
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
	enum {size = _Elems};
	using typename xbase_t::value_t;
	using pointer = value_t*;
	MATRICE_HOST_INL Packet_() noexcept : xbase_t() {}
	MATRICE_HOST_INL Packet_(const value_t _arg) noexcept : xbase_t(_arg) {}
	MATRICE_HOST_INL Packet_(const pointer _arg) noexcept : xbase_t(_arg) {}
	MATRICE_HOST_INL Packet_(const initlist_t _arg) noexcept : xbase_t(pointer(_arg.begin())) {}
	MATRICE_HOST_INL Packet_(const internal_t _arg) noexcept : xbase_t(_arg) {}

	MATRICE_HOST_INL auto& operator()() { return (m_data); }
	MATRICE_HOST_INL const auto& operator()() const { return (m_data); }
	MATRICE_HOST_INL myt& operator+ (const Packet_& _other);
	MATRICE_HOST_INL myt& operator- (const Packet_& _other);
	MATRICE_HOST_INL myt& operator* (const Packet_& _other);
	MATRICE_HOST_INL myt& operator/ (const Packet_& _other);
	MATRICE_HOST_INL myt  abs() const;
};
MATRICE_ARCH_END
#pragma region <!-- operators -->
template<typename T, int _N, typename Packet = dgelom::simd::Packet_<T, _N>>
MATRICE_HOST_FINL auto operator+ (const Packet& _Left, const Packet& _Right);
template<typename T, int _N, typename Packet = dgelom::simd::Packet_<T, _N>>
MATRICE_HOST_FINL auto operator- (const Packet& _Left, const Packet& _Right);
template<typename T, int _N, typename Packet = dgelom::simd::Packet_<T, _N>>
MATRICE_HOST_FINL auto operator* (const Packet& _Left, const Packet& _Right);
template<typename T, int _N, typename Packet = dgelom::simd::Packet_<T, _N>>
MATRICE_HOST_FINL auto operator/ (const Packet& _Left, const Packet& _Right);

MATRICE_NAMESPACE_BEGIN_
template<typename T, int _N, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
MATRICE_HOST_FINL auto abs (const dgelom::simd::Packet_<T, _N>& _Right);
template<typename T, int _N, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
MATRICE_HOST_FINL auto reduce(const dgelom::simd::Packet_<T, _N>& _Right);
_MATRICE_NAMESPACE_END
#pragma endregion
#include "./inl/_ixpacket.inl"
#endif