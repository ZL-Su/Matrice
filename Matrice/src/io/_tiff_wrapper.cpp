#include <../addin/libtiff/tiffio.h>
#include <io/_tiff_wrapper.hpp>
#include <core/matrix.h>

DGE_MATRICE_BEGIN
tiff_pointer open_tiff_file(const char* fname, const char* flag) {
	return TIFFOpen(fname, flag);
}

template<typename _Ty>
MATRICE_HOST_ONLY tiff_instance<_Ty> read_tiff_file(const char* fpath) {
	tiff_instance<_Ty> inst;
	if (const auto ptif = open_tiff_file(fpath, "r"); ptif) {
		TIFFGetField(ptif, TIFFTAG_IMAGELENGTH, &inst.m_rows);
		TIFFGetField(ptif, TIFFTAG_IMAGEWIDTH, &inst.m_cols);
		const auto scanline = TIFFScanlineSize(ptif);
		inst.create(scanline / inst.m_cols);

		auto buf = (uint8_t*)_TIFFmalloc(scanline * sizeof(uint8_t));
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
	return forward<decltype(inst)>(inst);
}

template<> tiff_instance<uint8_t> read_tiff_file(const char*);
template<> tiff_instance<uint16_t> read_tiff_file(const char*);
template<> tiff_instance<uint32_t> read_tiff_file(const char*);
template<> tiff_instance<uint64_t> read_tiff_file(const char*);
template<> tiff_instance<float_t> read_tiff_file(const char*);
template<> tiff_instance<double_t> read_tiff_file(const char*);
DGE_MATRICE_END