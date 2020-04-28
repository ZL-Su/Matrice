/**********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/
#pragma once
#include "util/utils.h"
#include "_type_traits.h"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
class _Bounded_pointer_base {
	using _Myt = _Bounded_pointer_base;
public:
	template<typename _Ty>
	constexpr _Bounded_pointer_base(_Ty* _Data, diff_t _Size) noexcept
		:_Mydata(_Data), _Mysize(_Size) {
	}

	template<typename _Ty>
	constexpr operator _Ty* () noexcept {
		return static_cast<_Ty*>(_Mydata);
	}

	template<typename _Ty>
	constexpr _Ty& operator*() noexcept {
		return *static_cast<_Ty*>(_Mydata);
	}

protected:
	void* _Mydata = nullptr;
	diff_t _Mysize;
};
_DETAIL_END
DGE_MATRICE_END