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
#define _SIMD dgelom::simd::

template<typename T, int _Elems> MATRICE_HOST_INL
_SIMD Packet_<T, _Elems>& _SIMD Packet_<T, _Elems>::operator+ (const _SIMD Packet_<T, _Elems>& _other)
{
	m_data = op_t::sum(m_data, _other.m_data); return (*this);
}
template<typename T, int _Elems> MATRICE_HOST_INL 
_SIMD Packet_<T, _Elems>& _SIMD Packet_<T, _Elems>::operator- (const _SIMD Packet_<T, _Elems>& _other)
{
	m_data = op_t::sub(m_data, _other.m_data); return (*this);
}
template<typename T, int _Elems> MATRICE_HOST_INL
_SIMD Packet_<T, _Elems>& _SIMD Packet_<T, _Elems>::operator* (const _SIMD Packet_<T, _Elems>& _other)
{
	m_data = op_t::mul(m_data, _other.m_data); return (*this);
}
template<typename T, int _Elems> MATRICE_HOST_INL
_SIMD Packet_<T, _Elems>& _SIMD Packet_<T, _Elems>::operator/ (const _SIMD Packet_<T, _Elems>& _other)
{
	m_data = op_t::div(m_data, _other.m_data); return (*this);
}
template<typename T, int _Elems> MATRICE_HOST_INL
_SIMD Packet_<T, _Elems> _SIMD Packet_<T, _Elems>::abs() const
{
	return (_SIMD Packet_<T, _Elems>(op_t::abs(m_data)));
}

template<typename T, int _N, typename Packet>
MATRICE_HOST_FINL auto operator+ (const Packet& _Left, const Packet& _Right)
{
	return Packet(simd_op<T, _N>::sum(_Left(), _Right()));
}
template<typename T, int _N, typename Packet>
MATRICE_HOST_FINL auto operator- (const Packet& _Left, const Packet& _Right)
{
	return Packet(simd_op<T, _N>::sub(_Left(), _Right()));
}
template<typename T, int _N, typename Packet>
MATRICE_HOST_FINL auto operator* (const Packet& _Left, const Packet& _Right)
{
	return Packet(simd_op<T, _N>::mul(_Left(), _Right()));
}
template<typename T, int _N, typename Packet>
MATRICE_HOST_FINL auto operator/ (const Packet& _Left, const Packet& _Right)
{
	return Packet(simd_op<T, _N>::div(_Left(), _Right()));
}
template<typename T, int _N, typename>
MATRICE_HOST_FINL auto dgelom::abs(const _SIMD Packet_<T, _N>& _Right)
{
	return (_SIMD Packet_<T, _N>(_SIMD simd_op<T, _N>::abs(_Right())));
}
template<typename T, int _N, typename>
MATRICE_HOST_FINL auto dgelom::reduce(const _SIMD Packet_<T, _N>& _Right)
{
	return (_Right.reduce());
}
#endif