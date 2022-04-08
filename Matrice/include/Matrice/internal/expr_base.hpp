/**********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2022, Zhilong(Dgelom) Su, all rights reserved.

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
**********************************************************************/
#pragma once

#include "util/_macros.h"

DGE_MATRICE_BEGIN
namespace xpr {
/// <summary>
/// \brief Provide a unified identity for expression.
/// </summary>
struct __xpr__{};

}

/// <summary>
/// \brief Expression concept. 
/// _Xpr is a concept iff it's derived from class xpr::__xpr__.
/// </summary>
template<typename _Xpr>
concept Expr = std::is_base_of_v<xpr::__xpr__, _Xpr>;
DGE_MATRICE_END