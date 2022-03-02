/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2022, Zhilong(Dgelom) Su, all rights reserved.

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

MATRICE_NAMESPACE_BEGIN(placeholders)

struct all_t { all_t() {}; };
/// <summary>
/// \brief VAR all, can be used to Matrix, Array, or Vector slicing and indexing.
/// </summary>
inline static const auto all = all_t{};

MATRICE_NAMESPACE_END(placeholders)