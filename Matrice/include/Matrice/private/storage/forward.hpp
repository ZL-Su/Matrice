/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
**************************************************************************/
#pragma once
#include "util/_macros.h"

namespace {
	static constexpr int dynamic = 0;
#ifdef MATRICE_ENABLE_CUDA
	static constexpr int device = -1;
	static constexpr int global = -2;
#endif
	static constexpr int extent_x = 0;
	static constexpr int extent_y = 1;
}

DGE_MATRICE_BEGIN
struct stack_alloc_tag {};
struct heap_alloc_tag {};
#ifdef MATRICE_ENABLE_CUDA
struct device_alloc_tag {};
struct global_alloc_tag {};
#endif

struct plain_layout {
	struct row_major { static constexpr auto value = 101; };
	struct col_major { static constexpr auto value = 102; };
};
_DETAIL_BEGIN
template<typename _Ty, diff_t _M, diff_t _N, size_t _Opt, typename _Layout> class _Allocator;
_DETAIL_END
template<typename _Ty, diff_t _RowsAtCT, diff_t _ColsAtCT, size_t _Options, class _Ly>
using Allocator = detail::_Allocator<_Ty, _RowsAtCT, _ColsAtCT, _Options, _Ly>;
DGE_MATRICE_END