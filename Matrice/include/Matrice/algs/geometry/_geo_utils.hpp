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
***********************************************************************/
#pragma once
#include "core/matrix.h"
#include "core/vector.h"

#ifdef _MSC_VER
#pragma warning(push)
#endif

DGE_MATRICE_BEGIN
/// <summary>
/// \brief datatype for defining dimensionality.
/// </summary>
enum dim_type : uint8_t {
	d1 = 1,
	d2 = 2,
	d3 = 3
};
using dim_t = dim_type;

namespace geo 
{

}
DGE_MATRICE_END

/// <summary>
/// \brief Lateral operator for specifying geometric dimensionality.
/// </summary>
/// <param name="_Dim">Dim with values of 1, 2, or 3.</param>
/// <returns>dgelom::dim_type</returns>
MATRICE_GLOBAL_FINL constexpr auto operator""D(size_t _Dim) noexcept {
	return dgelom::dim_type(static_cast<uint8_t>(_Dim));
}

#ifdef _MSC_VER
#pragma warning(pop)
#endif