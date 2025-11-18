/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2025, Zhilong(Dgelom) Su, all rights reserved.

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
*********************************************************************/
#ifdef MATRICE_SIMD_ARCH
#include "arch/internal/_regix.hpp"

MATRICE_ARCH_BEGIN
_DETAIL_BEGIN

#define MATRICE_MAKE_REGIX(INTR, TYPE) \
template<> struct _Regix<INTR> {       \
using value_t = TYPE;                  \
using pointer = value_t*;              \
INTR _Myreg;
#define MATRICE_END_REGIX };

#undef MATRICE_MAKE_REGIX
#undef MATRICE_END_REGIX

_DETAIL_END
MATRICE_ARCH_END
#endif
