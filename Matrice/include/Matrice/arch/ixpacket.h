/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

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

#ifdef MATRICE_SIMD_ARCH
#include "_ixbase.h"
#include "./inl/_ixop_impls.hpp"

MATRICE_ARCH_BEGIN

/*<brief>
 *\class template for SIMD vectorization 
 *\param <T> must be a scalar type
 *\param <_Elems> the number of packet elements, which default value is auto-deduced
</brief>*/
template<typename T, int _Elems = packet_size_v<T>>
class Packet_ MATRICE_NONHERITABLE : public simd::simd_base_<T, _Elems>
{
	using Myt = Packet_;
	using _Mybase = simd::simd_base_<T, _Elems>;
	using typename _Mybase::initlist_t;
	using typename _Mybase::internal_t;
	using typename _Mybase::op_t;
	using _Mybase::m_data;
public:
	enum {size = _Mybase::size};
	using typename _Mybase::value_t;
	using pointer = std::add_pointer_t<value_t>;

	MATRICE_HOST_FINL Packet_() noexcept 
		: _Mybase() {}
	MATRICE_HOST_FINL Packet_(const value_t _arg) noexcept 
		: _Mybase(_arg) {}
	MATRICE_HOST_FINL Packet_(const pointer _arg) noexcept 
		: _Mybase(_arg) {}
	MATRICE_HOST_FINL Packet_(const initlist_t _arg) noexcept 
		: _Mybase(pointer(_arg.begin())) {}
	MATRICE_HOST_FINL Packet_(const internal_t _arg) noexcept 
		: _Mybase(_arg) {}
	template<typename _Fwdty, MATRICE_ENABLE_IF(has_data_v<_Fwdty>)>
	MATRICE_HOST_FINL Packet_(const _Fwdty& _arg) noexcept 
		: _Mybase(_arg.data()) {}

	MATRICE_HOST_FINL auto& operator()() { return (m_data); }
	MATRICE_HOST_FINL const auto& operator()() const { return (m_data); }
	MATRICE_HOST_FINL Myt& operator+ (const Packet_& _other);
	MATRICE_HOST_FINL Myt& operator- (const Packet_& _other);
	MATRICE_HOST_FINL Myt& operator* (const Packet_& _other);
	MATRICE_HOST_FINL Myt& operator/ (const Packet_& _other);
	MATRICE_HOST_FINL Myt  abs() const;

	template<typename T, int _Elems, typename Packet> friend
	MATRICE_HOST_FINL Packet abs(const Packet_<T, _Elems>& _Right);
	template<typename T, int _Elems, typename> friend
	MATRICE_HOST_FINL T reduce(const Packet_<T, _Elems>& _Right);
};

#pragma region <!-- operators -->
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator+ (const Packet_<T, _Elems>& _Left, const Packet_<T, _Elems>& _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator- (const Packet_<T, _Elems>& _Left, const Packet_<T, _Elems>& _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator* (const Packet_<T, _Elems>& _Left, const Packet_<T, _Elems>& _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator/ (const Packet_<T, _Elems>& _Left, const Packet_<T, _Elems>& _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator+ (const Packet_<T, _Elems>& _Left, T _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator- (const Packet_<T, _Elems>& _Left, T _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator* (const Packet_<T, _Elems>& _Left, T _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator/ (const Packet_<T, _Elems>& _Left, T _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator+ (T _Left, const Packet_<T, _Elems>& _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator- (T _Left, const Packet_<T, _Elems>& _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator* (T _Left, const Packet_<T, _Elems>& _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator/ (T _Left, const Packet_<T, _Elems>& _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator+ (const Packet_<T, _Elems>& _Left, T* _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator- (const Packet_<T, _Elems>& _Left, T* _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator* (const Packet_<T, _Elems>& _Left, T* _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator/ (const Packet_<T, _Elems>& _Left, T* _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator+ (T* _Left, const Packet_<T, _Elems>& _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator- (T* _Left, const Packet_<T, _Elems>& _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator* (T* _Left, const Packet_<T, _Elems>& _Right);
template<typename T, int _Elems, typename Packet = Packet_<T, _Elems>>
MATRICE_HOST_FINL Packet operator/ (T* _Left, const Packet_<T, _Elems>& _Right);
#pragma endregion

MATRICE_ARCH_END
#include "./inl/_ixpacket.inl"
#endif