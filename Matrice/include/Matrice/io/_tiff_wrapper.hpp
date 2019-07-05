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
#include <libtiff/tiffio.h>
#include "../util/_macros.h"
#include "../util/_std_wrapper.h"

DGE_MATRICE_BEGIN
template<typename _Ty, int _M, int _N> class Matrix_;

using tiff_type = TIFF;
using tiff_pointer = std::add_pointer_t<tiff_type>;

template<typename _Ty> struct tiff_instance {
	using value_type = _Ty;

	std::add_pointer_t<value_type> m_data = nullptr;
	uint32_t m_rows = 0, m_cols = 0;
	std::size_t m_size = 0;
};

tiff_pointer open_tiff_file(const char* fname, const char* flag) {
	return TIFFOpen(fname, flag);
}

template<typename _Ty, typename _Ret = tiff_instance<_Ty>>
MATRICE_HOST_INL _Ret read_tiff_file(const char* fpath) {
	_Ret inst;
	if (const auto ptif = open_tiff_file(fpath, "r"); ptif) {
		uint32_t length;
		TIFFGetField(ptif, TIFFTAG_IMAGELENGTH, &length);

		const auto scanline = TIFFScanlineSize(ptif);
		auto buf = _TIFFmalloc(scanline);
		for (inst.m_rows = 0; inst.m_rows < length; ++inst.m_rows) {
			TIFFReadScanline(ptif, buf, inst.m_rows);
			for (inst.m_cols = 0; inst.m_cols < scanline; ++inst.m_cols) {
				inst.m_data[inst.m_size++] = static_cast<_Ty>(buf[inst.m_cols]);
			}
		}
		_TIFFfree(buf);
		TIFFClose(ptif);
	}
}

DGE_MATRICE_END