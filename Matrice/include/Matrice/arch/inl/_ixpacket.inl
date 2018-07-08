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
#ifdef __AVX__
#include "../ixpacket.h"
#define _SIMD

template<typename T, int _Elems> MATRICE_HOST_FINL
Packet_<T, _Elems>& Packet_<T, _Elems>::operator+ (const Packet_<T, _Elems>& _other)
{
	m_data = transform(details::Op_<value_t>::plus<size>(), *this, _other);
	return (*this);
}
template<typename T, int _Elems> MATRICE_HOST_FINL 
Packet_<T, _Elems>& Packet_<T, _Elems>::operator- (const Packet_<T, _Elems>& _other)
{
	m_data = transform(details::Op_<value_t>::minus<size>(), *this, _other);
	return (*this);
}
template<typename T, int _Elems> MATRICE_HOST_FINL
Packet_<T, _Elems>& Packet_<T, _Elems>::operator* (const Packet_<T, _Elems>& _other)
{
	m_data = transform(details::Op_<value_t>::multiplies<size>(), *this, _other);
	return (*this);
}
template<typename T, int _Elems> MATRICE_HOST_FINL
Packet_<T, _Elems>& Packet_<T, _Elems>::operator/ (const Packet_<T, _Elems>& _other)
{
	m_data = transform(details::Op_<value_t>::divides<size>(), *this, _other);
	return (*this);
}
template<typename T, int _Elems> MATRICE_HOST_FINL
Packet_<T, _Elems> Packet_<T, _Elems>::abs() const
{
	return (transform(details::Op_<value_t>::abs<size>(), *this));
}

template<typename T, int _N, typename Packet = Packet_<T, _N>>
MATRICE_HOST_FINL Packet operator+ (const Packet_<T, _N>& _Left, const Packet_<T, _N>& _Right)
{
	return (transform(details::Op_<T>::plus<_N>(), _Left, _Right));
}
template<typename T, int _N, typename Packet = Packet_<T, _N>>
MATRICE_HOST_FINL Packet operator- (const Packet_<T, _N>& _Left, const Packet_<T, _N>& _Right)
{
	return (transform(details::Op_<T>::minus<_N>(), _Left, _Right));
}
template<typename T, int _N, typename Packet = Packet_<T, _N>>
MATRICE_HOST_FINL Packet operator* (const Packet_<T, _N>& _Left, const Packet_<T, _N>& _Right)
{
	return (transform(details::Op_<T>::multiplies<_N>(), _Left, _Right));
}
template<typename T, int _N, typename Packet = Packet_<T, _N>>
MATRICE_HOST_FINL Packet operator/ (const Packet_<T, _N>& _Left, const Packet_<T, _N>& _Right)
{
	return (transform(details::Op_<T>::divides<_N>(), _Left, _Right));
}
template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator+(const Packet_<T, _Elems>& _Left, T _Right)
{
	return (_Left + Packet(_Right));
}
template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator-(const Packet_<T, _Elems>& _Left, T _Right)
{
	return (_Left - Packet(_Right));
}
template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator*(const Packet_<T, _Elems>& _Left, T _Right)
{
	return (_Left * Packet(_Right));
}
template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator/(const Packet_<T, _Elems>& _Left, T _Right)
{
	return (_Left / Packet(_Right));
}
template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator+(T _Left, const Packet_<T, _Elems>& _Right)
{
	return (Packet(_Left) + _Right);
}
template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator-(T _Left, const Packet_<T, _Elems>& _Right)
{
	return (Packet(_Left) - _Right);
}
template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator*(T _Left, const Packet_<T, _Elems>& _Right)
{
	return (Packet(_Left) * _Right);
}
template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator/(T _Left, const Packet_<T, _Elems>& _Right)
{
	return (Packet(_Left) / _Right);
}

template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator+(const Packet_<T, _Elems>& _Left, T* _Right)
{
	return (_Left + Packet(_Right));
}
template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator-(const Packet_<T, _Elems>& _Left, T* _Right)
{
	return (_Left - Packet(_Right));
}
template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator*(const Packet_<T, _Elems>& _Left, T* _Right)
{
	return (_Left * Packet(_Right));
}
template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator/(const Packet_<T, _Elems>& _Left, T* _Right)
{
	return (_Left / Packet(_Right));
}
template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator+(T* _Left, const Packet_<T, _Elems>& _Right)
{
	return (Packet(_Left) + _Right);
}
template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator-(T* _Left, const Packet_<T, _Elems>& _Right)
{
	return (Packet(_Left) - _Right);
}
template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator*(T* _Left, const Packet_<T, _Elems>& _Right)
{
	return (Packet(_Left) * _Right);
}
template<typename T, int _Elems, typename Packet>
MATRICE_HOST_FINL Packet operator/(T* _Left, const Packet_<T, _Elems>& _Right)
{
	return (Packet(_Left) / _Right);
}

template<typename T, int _N, typename Packet = Packet_<T, _N>>
MATRICE_HOST_FINL Packet abs(const Packet_<T, _N>& _Right)
{
	return (transform(details::Op_<T>::abs<_N>(), _Right));
}
template<typename T, int _N, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
MATRICE_HOST_FINL T reduce(const Packet_<T, _N>& _Right)
{
	return (_Right.reduce());
}
template<typename _InIt, typename _OutIt, typename _Op, typename = std::enable_if_t<std::is_pointer_v<_InIt>&&std::is_pointer_v<_OutIt>>>
void transform(const _InIt _Src, _OutIt _Dst, size_t _N, _Op _Fn) {
	enum { _Elems = 4 };
	using value_t = std::remove_reference_t<decltype(*_First)>;
	size_t _Packed_size = _N / _Elems;
	for (size_t i = 0; i < _Packed_size; ++i) {
		_Fn(Packet_<value_t, _Elems>(_Src + i * _Elems)).unpack(_Dest + i * _Elems);
	}
	for (size_t i = _Packed_size * _Elems; i < _N; ++i) {
		*(_Dst + i) = _Fn(*(_Src + i));
	}
}

#endif