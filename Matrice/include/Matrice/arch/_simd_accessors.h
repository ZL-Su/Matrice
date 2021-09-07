/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

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
*********************************************************************/
#pragma once

#ifdef MATRICE_SIMD_ARCH
#include "_simd_traits.hpp"

MATRICE_ARCH_BEGIN
struct load {};
struct store {};
_DETAIL_BEGIN

// For internal use...
namespace internal {
template<typename... _Args>
MATRICE_HOST_FINL auto _Load(_Args... _) {};
template<typename... _Args>
MATRICE_HOST_FINL auto _Store(_Args... _) {};
template<typename... _Args>
MATRICE_HOST_FINL auto _Set_zero() {};
}

// Register storer
template<typename _Vec> struct _Storer {
	const _Vec& _Data;
	_Storer(const _Vec& _Src) : _Data(_Src) {}
	template<typename _Ty>
	MATRICE_HOST_FINL auto operator()(_Ty* _Dst) noexcept {
		return internal::_Store(_Data, _Dst);
	}
};

// Register loader
template<typename _Vec> struct _Loader {
	template<typename... _Args>
	MATRICE_HOST_FINL auto operator()(_Args... _Src) const noexcept {
		return internal::_Load(_Vec{}, _Src...);
	}
};

/// <summary>
/// \brief CLASS TEMPLATE, uniform interface for accessing the register.
/// </summary>
template<packetable_scalar _Ty, 
	uint8_t _N = packet_traits<_Ty>::size> 
class _Accessor : public packed_vector<_Ty, _N> {
	using _Mybase = packed_vector<_Ty, _N>;
	using _Native_type = typename _Mybase::type;
public:
	MATRICE_HOST_FINL static _Loader<_Mybase> fwd() noexcept {
		return _Loader<_Mybase>();
	}
	MATRICE_HOST_FINL static _Storer<_Native_type> bwd(const _Native_type& _Reg) noexcept {
		return _Storer<_Native_type>{_Reg};
	}
};

_DETAIL_END

template<typename _Ty, uint8_t _Size= packet_traits<_Ty>::size>
using accessor = detail::_Accessor<_Ty, _Size>;

MATRICE_ARCH_END
#include "inl/_simd_accessor.hpp"
#endif // MATRICE_SIMD_ARCH
