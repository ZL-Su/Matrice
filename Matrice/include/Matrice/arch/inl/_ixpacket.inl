/**************************************************************************
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
**************************************************************************/
#pragma once
#include "util/genalgs.h"
#ifdef MATRICE_SIMD_ARCH
#include "../simd.h"
MATRICE_ARCH_BEGIN
template<typename T> MATRICE_HOST_FINL
Packet_<T>& Packet_<T>::operator+(const Packet_<T>& _other)
{
	m_data = transform<op_t::plus<size>>(*this, _other);
	return (*this);
}
template<typename T> MATRICE_HOST_FINL 
Packet_<T>& Packet_<T>::operator-(const Packet_<T>& _other)
{
	m_data = transform<op_t::minus<size>>(*this, _other);
	return (*this);
}
template<typename T> MATRICE_HOST_FINL
Packet_<T>& Packet_<T>::operator*(const Packet_<T>& _other)
{
	m_data = transform<op_t::multiplies<size>>(*this, _other);
	return (*this);
}
template<typename T> MATRICE_HOST_FINL
Packet_<T>& Packet_<T>::operator/(const Packet_<T>& _other)
{
	m_data = transform<op_t::divides<size>>(*this, _other);
	return (*this);
}
template<typename T> MATRICE_HOST_FINL
Packet_<T> Packet_<T>::abs() const
{
	return (transform<op_t::abs<size>>(*this));
}

template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator+(const Packet_<T>& _Left, const Packet_<T>& _Right)
{
	return (transform<Packet::op_t::plus<Packet::size>>(_Left, _Right));
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator-(const Packet_<T>& _Left, const Packet_<T>& _Right)
{
	return (transform<Packet::op_t::minus<Packet::size>>(_Left, _Right));
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator*(const Packet_<T>& _Left, const Packet_<T>& _Right)
{
	return (transform<Packet::op_t::multiplies<Packet::size>>(_Left, _Right));
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator/(const Packet_<T>& _Left, const Packet_<T>& _Right)
{
	return (transform<Packet::op_t::divides<Packet::size>>(_Left, _Right));
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator+(const Packet_<T>& _Left, T _Right)
{
	return (_Left + Packet(_Right));
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator-(const Packet_<T>& _Left, T _Right)
{
	return (_Left - Packet(_Right));
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator*(const Packet_<T>& _Left, T _Right)
{
	return (_Left * Packet(_Right));
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator/(const Packet_<T>& _Left, T _Right)
{
	return (_Left / Packet(_Right));
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator+(T _Left, const Packet_<T>& _Right)
{
	return (Packet(_Left) + _Right);
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator-(T _Left, const Packet_<T>& _Right)
{
	return (Packet(_Left) - _Right);
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator*(T _Left, const Packet_<T>& _Right)
{
	return (Packet(_Left) * _Right);
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator/(T _Left, const Packet_<T>& _Right)
{
	return (Packet(_Left) / _Right);
}

template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator+(const Packet_<T>& _Left, T* _Right)
{
	return (_Left + Packet(_Right));
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator-(const Packet_<T>& _Left, T* _Right)
{
	return (_Left - Packet(_Right));
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator*(const Packet_<T>& _Left, T* _Right)
{
	return (_Left * Packet(_Right));
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator/(const Packet_<T>& _Left, T* _Right)
{
	return (_Left / Packet(_Right));
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator+(T* _Left, const Packet_<T>& _Right)
{
	return (Packet(_Left) + _Right);
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator-(T* _Left, const Packet_<T>& _Right)
{
	return (Packet(_Left) - _Right);
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator*(T* _Left, const Packet_<T>& _Right)
{
	return (Packet(_Left) * _Right);
}
template<typename T, typename Packet>
MATRICE_HOST_FINL Packet operator/(T* _Left, const Packet_<T>& _Right)
{
	return (Packet(_Left) / _Right);
}

template<typename T, typename Packet = Packet_<T>>
MATRICE_HOST_FINL Packet abs(const Packet_<T>& _Right)
{
	return (transform<Packet::op_t::abs<Packet::size>>(_Right));
}
template<typename T, typename = enable_if_t<is_arithmetic_v<T>>>
MATRICE_HOST_FINL T reduce(const Packet_<T>& _Right)
{
	return (dgelom::reduce_n_t<Packet_<T>::size>::value(_Right.begin()));
}
template<typename _InIt, typename _OutIt, typename _Op, 
	typename = enable_if_t<is_pointer_v<_InIt>&&is_pointer_v<_OutIt>>>
void transform(const _InIt _Begin, const _InIt _End, _OutIt _Dst, _Op _Fn) {
	enum { _Elems = 4 };
	using value_t = remove_reference_t<decltype(*_Begin)>;

	const auto _N = std::distance(_Begin, _End);
	const auto _Packed_size = vsize<_Elems>(_N);
	for (size_t i = 0; i < _Packed_size; i += _Elems) {
		_Fn(Packet_<value_t>(_Begin + i)).unpack(_Dst + i);
	}
	for (size_t i = _Packed_size; i < _N; ++i) {
		*(_Dst + i) = _Fn(*(_Begin + i));
	}
}
MATRICE_ARCH_END
#endif