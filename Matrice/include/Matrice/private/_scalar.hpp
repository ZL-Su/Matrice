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
	static_assert(std::is_scalar_v<_Ty>, 
		"_Ty in Scalar must be a native scalar type.");
	using _Myt = Scalar;
public:
	using value_type = dgelom::remove_all_t<_Ty>;
	using refernce = std::add_lvalue_reference_t<value_type>;
	using pointer = std::add_pointer_t<value_type>;

	MATRICE_GLOBAL_FINL Scalar() noexcept {
	}
	template<typename _Uy>
	MATRICE_GLOBAL_FINL inline Scalar(const _Uy s) noexcept 
		: m_value(s) {
	}

	MATRICE_GLOBAL_FINL operator refernce() noexcept {
		return (m_value);
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
DGE_MATRICE_END