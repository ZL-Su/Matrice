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
#include <../addin/libtiff/tiffio.h>
#include "io/io.hpp"
#include "io/_tiff_wrapper.hpp"

DGE_MATRICE_BEGIN
inline tiff_pointer open_tiff_file(const char* fname, const char* flag) {
	return TIFFOpen(fname, flag);
}

MATRICE_HOST_ONLY image_instance read_tiff_file(const char* fpath) {
	image_instance inst;
	if (const auto ptif = open_tiff_file(fpath, "r"); ptif) {
		TIFFGetField(ptif, TIFFTAG_IMAGELENGTH, &inst.m_rows);
		TIFFGetField(ptif, TIFFTAG_IMAGEWIDTH, &inst.m_cols);
		const auto scanline = TIFFScanlineSize(ptif);
		inst.create(scanline / inst.m_cols);

		auto buf = (uint8_t*)_TIFFmalloc(scanline * sizeof(uint8_t));
		for (auto row = 0; row < inst.m_rows; ++row) {
			TIFFReadScanline(ptif, buf, row);
			auto pd = inst.m_data.data() + (row * inst.m_width);
			for (auto col = 0; col < inst.m_cols; ++col) {
				auto pb = buf + col * inst.m_nchs;
				for (auto chn = 0; chn < inst.m_nchs; ++chn) {
					pd[col + chn * inst.m_cols] = pb[chn];
				}
			}
		}

		_TIFFfree(buf);
		TIFFClose(ptif);
	}
	return forward<decltype(inst)>(inst);
}
DGE_MATRICE_END