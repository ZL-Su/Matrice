/*********************************************************************
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
***********************************************************************/
#pragma once
#include "../io.hpp"
#include "../_tiff_wrapper.hpp"
#include "../_png_wrapper.hpp"

DGE_MATRICE_BEGIN namespace io 
{ _DETAIL_BEGIN

template<typename T> MATRICE_HOST_INL
const std::string _Get_image_fmt(const T fname) noexcept {
	const std::string str = fname;
	return str.substr(str.size() - 3);
}

MATRICE_HOST_INL decltype(auto) _Imread(const char* path, uint8_t) {
	const auto fmt = _Get_image_fmt(path);

#ifdef MATRICE_DEBUG
	DGELOM_CHECK(fmt.size() == 3, "Unsupported image format.");
#endif // MATRICE_DEBUG

	if (fmt == "tif") {
		return read_tiff_file<uint8_t>(path);
	}
}

MATRICE_HOST_INL decltype(auto) _Imread(const char* path, float) {
	const auto fmt = _Get_image_fmt(path);

#ifdef MATRICE_DEBUG
	DGELOM_CHECK(fmt.size() == 3, "Unsupported image format.");
#endif // MATRICE_DEBUG

	if (fmt == "tif") {
		return read_tiff_file<float>(path);
	}
}

MATRICE_HOST_INL decltype(auto) _Imread(const char* path, double) {
	const auto fmt = _Get_image_fmt(path);

#ifdef MATRICE_DEBUG
	DGELOM_CHECK(fmt.size() == 3, "Unsupported image format.");
#endif // MATRICE_DEBUG

	if (fmt == "tif") {
		return read_tiff_file<double>(path);
	}
}

_DETAIL_END

template<typename _Ty, class _Pth>
MATRICE_HOST_INL decltype(auto) imread(const _Pth path) {
	//static_assert(is_any_of_v<std::decay_t<_Pth>, std::string, io::path_t, char*>,"Unknown path type in function imread().");

	if constexpr (is_same_v<remove_all_t<_Pth>, std::string>) {
		return (detail::_Imread(path.c_str(), _Ty()));
	}
	else if constexpr (is_same_v<remove_all_t<_Pth>, io::path_t>) {
		const auto str = path.generic_string();
		return (detail::_Imread(str.c_str(), _Ty()));
	}
	else return (detail::_Imread(path, _Ty()));
}
}
DGE_MATRICE_END