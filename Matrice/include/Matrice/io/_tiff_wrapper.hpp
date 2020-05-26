/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

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
#include "../util/_macros.h"
#include "../private/_type_traits.h"
#include "../forward.hpp"

// TIFF forward declaration
struct tiff;

DGE_MATRICE_BEGIN
using tiff_type = tiff;
using tiff_pointer = std::add_pointer_t<tiff_type>;
struct image_instance;

/**
 *\brief read tiff file to image_instance
 *\param [fpath] file path input.
 */
MATRICE_HOST_ONLY image_instance read_tiff_file(const char* fpath);
DGE_MATRICE_END