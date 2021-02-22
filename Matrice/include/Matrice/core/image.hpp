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
#include "matrix.h"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
/// <summary>
/// \brief CLASS TEMPLATE, Smooth image type
/// </summary>
/// <typeparam name="_Ty">Floating point type</typeparam>
/// <typeparam name="_Ip">Interpolation operator type</typeparam>
template<typename _Ty, class _Ip> requires is_floating_point_v<_Ty>
class _Image : public Matrix_<_Ty, ::dynamic, ::dynamic> {

};

/// <summary>
/// \brief CLASS TEMPLATE, General image type
/// </summary>
/// <typeparam name="_Ty">Scalar value type</typeparam>
template<typename _Ty>
class _Image : public Matrix_<_Ty, ::dynamic, ::dynamic> {

};
_DETAIL_END
DGE_MATRICE_END