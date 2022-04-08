/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2020-2022, Zhilong(Dgelom) Su, all rights reserved.

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
#include "util/_macros.h"
#include "util/_type_defs.h"
#include "private/_tag_defs.h"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Derived> class _Index_base
{
	using _Myt = _Index_base;
public:
	MATRICE_GLOBAL_INL _Derived& operator++() noexcept {
		static_cast<_Derived*>(this)->_Myval += 1;
		return *static_cast<_Derived*>(this);
	}
	MATRICE_GLOBAL_INL _Derived operator++(int) noexcept {
		auto _Tmp = *static_cast<_Derived*>(this);
		static_cast<_Derived*>(this)->_Myval += 1;
		return (_Tmp);
	}
	MATRICE_GLOBAL_INL _Derived& operator--() noexcept {
		static_cast<_Derived*>(this)->_Myval -= 1;
		return *static_cast<_Derived*>(this);
	}
	MATRICE_GLOBAL_INL _Derived operator--(int) noexcept {
		auto _Tmp = *static_cast<_Derived*>(this);
		static_cast<_Derived*>(this)->_Myval -= 1;
		return (_Tmp);
	}
	MATRICE_GLOBAL_INL _Derived& operator+=(diff_t _Offset) noexcept {
		static_cast<_Derived*>(this)->_Myval += _Offset;
		return *static_cast<_Derived*>(this);
	}
	MATRICE_GLOBAL_INL _Derived& operator-=(diff_t _Offset) noexcept {
		static_cast<_Derived*>(this)->_Myval -= _Offset;
		return *static_cast<_Derived*>(this);
	}
};
_DETAIL_END

/**
 * \brief Specify index struct to signed 64-bit integer type (long long).
 */
template<>
struct Index<int64_t>: detail::_Index_base<Index<int64_t>> {
	using value_type = int64_t;
	using value_t = value_type;
	using category = tag::scalar;

	Index(value_type val = 0) noexcept :_Myval(val) {}
	template<typename _Uy>
	Index& operator=(const _Uy val)noexcept {
		_Myval = value_type(val);
	}

	MATRICE_GLOBAL_INL operator value_type() const noexcept {
		return _Myval;
	}

	value_type _Myval{ 0 };
};
using Index_t = Index<int64_t>;
using Index2_t = Index<int64_t, int64_t>;
using Index3_t = Index<int64_t, int64_t, int64_t>;
DGE_MATRICE_END