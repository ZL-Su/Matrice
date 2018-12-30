/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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
#include<experimental/filesystem>
#include "../util/_macros.h"

DGE_MATRICE_BEGIN namespace io {
/**
 * \for filesystem namespace in STD library
 */
namespace fs = std::experimental::filesystem;

namespace detail {

struct folder_tag {};

template<typename _Tag = folder_tag> class _Dir_impl {};

template<> class _Dir_impl<folder_tag> {
	using value_type = std::string_view;
	using path_type = fs::path;
public:
	_Dir_impl(value_type _Root, path_type&& _Workdir = fs::current_path()) {

	}

private:
	path_type _Mypath;
};
}


} DGE_MATRICE_END