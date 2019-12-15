/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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

#include "_type_traits.h"

DGE_MATRICE_BEGIN

template<typename _Ty>
class Scalar {
	static_assert(is_scalar_v<remove_all_t<_Ty>>, 
		"_Ty in Scalar must be a primitive scalar type.");
	using _Myt = Scalar;
public:
	enum {
		Size = 1,
		rows_at_compiletime = 1,
		cols_at_compiletime = 1,
	};
	using value_type = remove_all_t<_Ty>;
	using reference = value_type&;
	using pointer = value_type*;

	MATRICE_GLOBAL_FINL Scalar() noexcept
		: m_value(0) {
	}
	template<typename _Uy>
	MATRICE_GLOBAL_FINL inline Scalar(const _Uy s) noexcept 
		: m_value(s) {
	}

	MATRICE_GLOBAL_FINL operator reference() const noexcept {
		return (m_value);
	}
	MATRICE_GLOBAL_FINL operator reference() noexcept {
		return (m_value);
	}

	MATRICE_GLOBAL_FINL operator pointer() const noexcept {
		return (&m_value);
	}
	MATRICE_GLOBAL_FINL operator pointer() noexcept {
		return (&m_value);
	}

private:
	value_type m_value;
};
template<typename _Ty>
struct is_scalar<Scalar<_Ty>> {
	static constexpr auto value = true;
};
template<typename _Ty> struct remove_all<Scalar<_Ty>> {
	using type = typename Scalar<_Ty>::value_type;
};
DGE_MATRICE_END