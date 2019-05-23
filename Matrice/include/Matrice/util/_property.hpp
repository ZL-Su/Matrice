/**************************************************************************
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
**************************************************************************/
#pragma once
#include "_macros.h"

DGE_MATRICE_BEGIN

#define MATRICE_GET(T, NAME, EXPR) \
__declspec(property(get=_##NAME##_getter))T NAME;\
MATRICE_HOST_FINL T _##NAME##_getter() const noexcept { return EXPR; }

#define MATRICE_SET(T, NAME, EXPR) \
__declspec(property(put=_##NAME##_setter))T NAME;\
MATRICE_HOST_FINL void _##NAME##_setter(const T NAME) noexcept { EXPR = NAME; }

DGE_MATRICE_END