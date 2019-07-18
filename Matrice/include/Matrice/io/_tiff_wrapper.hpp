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
#include "../../../addin/libtiff/tiffio.h"
#include "../util/_macros.h"
#include "../util/_std_wrapper.h"

DGE_MATRICE_BEGIN
using tiff_type = tiff;
using tiff_pointer = std::add_pointer_t<tiff_type>;
template<typename _Ty> using tiff_matrix_t = Matrix_<_Ty, 0, 0>;
template<typename _Ty> struct image_instance;

template<typename _Ty = uint8_t> 
struct tiff_instance : image_instance<_Ty> {

	template<typename _Uy = _Ty>
	MATRICE_HOST_INL tiff_matrix_t<_Uy> grayscale() const noexcept {
		if (this->m_nchs == 1) if constexpr (is_same_v<_Ty, _Uy>)
			return (this->m_data);

		tiff_matrix_t<_Uy> gc(this->m_rows, this->m_cols);
		if (this->m_nchs == 1) {
			gc.from(this->m_data.data());
		}
		else {
			for (auto r = 0; r < this->m_rows; ++r) {
				const auto pd = this->m_data[r];
				auto pg = gc[r];
				for (auto c = 0; c < this->m_cols; ++c) {
					pg[c] = pd[c] * 0.07 + pd[c + this->m_cols] * 0.72 + pd[c + 2 * this->m_cols] * 0.21;
				}
			}
		}
		if constexpr (is_same_v<_Ty, uint8_t>)
			if constexpr (is_floating_point_v<_Uy>)
				gc = gc / 255;

		return forward<decltype(gc)>(gc);
	}
};

tiff_pointer open_tiff_file(const char* fname, const char* flag) {
	return TIFFOpen(fname, flag);
}

template<typename _Ty=uint8_t, typename _Ret = tiff_instance<_Ty>>
MATRICE_HOST_INL _Ret read_tiff_file(const char* fpath) {
	_Ret inst;
	if (const auto ptif = open_tiff_file(fpath, "r"); ptif) {
		TIFFGetField(ptif, TIFFTAG_IMAGELENGTH, &inst.m_rows);
		TIFFGetField(ptif, TIFFTAG_IMAGEWIDTH, &inst.m_cols);
		const auto scanline = TIFFScanlineSize(ptif);
		inst.create(scanline/inst.m_cols);

		auto buf = (uint8_t*)_TIFFmalloc(scanline*sizeof(uint8_t));
		for (auto row = 0; row < inst.m_rows; ++row) {
			TIFFReadScanline(ptif, buf, row);
			auto pd = inst.m_data[row];
			for (auto col = 0; col < inst.m_cols; ++col) {
				auto pb = buf + col * inst.m_nchs;
				for (auto chn = 0; chn < inst.m_nchs; ++chn) {
					auto& val = pd[col + chn * inst.m_cols] = pb[chn];
					if constexpr (!is_same_v<_Ty, uint8_t>)
						val /= (_Ty)(255);
				}
			}
		}

		_TIFFfree(buf);
		TIFFClose(ptif);
	}
	return forward<_Ret>(inst);
}

DGE_MATRICE_END