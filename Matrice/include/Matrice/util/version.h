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
*********************************************************************/
#pragma once
#include "_macros.h"

#define MATRICE_VERSION_MAJOR 2021
#define MATRICE_VERSION_MINOR 3
#define MATRICE_VERSION_REVISION 10

// Classic CPP stringifcation; the extra level of indirection allows the
// preprocessor to expand the macro before being converted to a string.
#define MATRICE_VER_STRING(x) MATRICE_STRINGFY(x)

// The Matrice version as a string; for example "2021.3.10".
#define MATRICE_VERSION_STRING \
MATRICE_VER_STRING(MATRICE_VERSION_MAJOR) "." \
MATRICE_VER_STRING(MATRICE_VERSION_MINOR) "." \
MATRICE_VER_STRING(MATRICE_VERSION_REVISION)