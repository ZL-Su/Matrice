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
#include "../io.hpp"
#include "../_tiff_wrapper.hpp"
#include "../_png_wrapper.hpp"

DGE_MATRICE_BEGIN
struct image_constant {
	template<typename _Ty> static constexpr auto max_value = _Ty(255);
};

struct image_instance {
	using value_type = uint8_t;
	using data_type = std::vector<value_type>;

	template<typename _Ty>
	MATRICE_HOST_INL Matrix_<_Ty, ::dynamic> matrix()const noexcept {
		Matrix_<_Ty, ::dynamic> _Ret(m_rows, m_width);
		_Ret.from(m_data.data());
		if constexpr (is_floating_point_v<_Ty>)
			_Ret = _Ret / image_constant::max_value<_Ty>;
		return (_Ret);
	}
	MATRICE_HOST_INL image_instance& create(uint32_t nchs) {
		m_nchs = nchs;
		m_width = m_cols * m_nchs;
		m_data.resize(m_rows * m_width);
		return (*this);
	}
	template<typename _Op>
	MATRICE_HOST_INL auto operator()(_Op&& op) noexcept {
		return (op(m_data));
	}

	template<typename _Ty>
	MATRICE_HOST_INL Matrix_<_Ty, ::dynamic> grayscale() const noexcept {
		Matrix_<_Ty, ::dynamic> gc(this->m_rows, this->m_cols);
		if (this->m_nchs == 1) {
			gc.from(this->m_data.data());
		}
		else {
			for (auto r = 0; r < this->m_rows; ++r) {
				const auto pd = this->m_data.data() + r * m_width;
				auto pg = gc[r];
				for (auto c = 0; c < this->m_cols; ++c) {
					pg[c] = pd[c] * 0.07 + pd[c + m_cols] * 0.72 + pd[c + 2 * m_cols] * 0.21;
				}
			}
		}
		if constexpr (is_floating_point_v<_Ty>)
			gc = gc / image_constant::max_value<_Ty>;

		return forward<decltype(gc)>(gc);
	}

	MATRICE_HOST_INL size_t size() const noexcept {
		return m_data.size();
	}

	uint32_t m_rows = 0, m_cols = 0;
	uint32_t m_nchs = 1, m_width = 0;
	data_type m_data;
};
namespace io 
{ _DETAIL_BEGIN

template<typename T> MATRICE_HOST_INL
const std::string _Get_image_fmt(const T fname) noexcept {
	const std::string str = fname;
	return str.substr(str.size() - 3);
}

MATRICE_HOST_INL auto _Imread(const char* path, uint8_t) {
	const auto fmt = _Get_image_fmt(path);

#ifdef MATRICE_DEBUG
	DGELOM_CHECK(fmt.size() == 3, "Unsupported image format.");
#endif // MATRICE_DEBUG

	if (fmt == "tif") {
		return read_tiff_file(path);
	}
}

MATRICE_HOST_INL auto _Imread(const char* path, float) {
	const auto fmt = _Get_image_fmt(path);

#ifdef MATRICE_DEBUG
	DGELOM_CHECK(fmt.size() == 3, "Unsupported image format.");
#endif // MATRICE_DEBUG

	if (fmt == "tif") {
		return read_tiff_file(path);
	}
}

MATRICE_HOST_INL auto _Imread(const char* path, double) {
	const auto fmt = _Get_image_fmt(path);

#ifdef MATRICE_DEBUG
	DGELOM_CHECK(fmt.size() == 3, "Unsupported image format.");
#endif // MATRICE_DEBUG

	if (fmt == "tif") {
		return read_tiff_file(path);
	}
}

_DETAIL_END

template<typename _Ty, class _Pth>
MATRICE_HOST_INL auto imread(const _Pth path) {
	static_assert(is_any_of_v<remove_all_t<_Pth>, std::string, fs::path, const char*>,"Unknown path type in function imread().");

	if constexpr (is_same_v<remove_all_t<_Pth>, std::string>) {
		return (detail::_Imread(path.c_str(), _Ty()));
	}
	else if constexpr (is_same_v<remove_all_t<_Pth>, fs::path>) {
		const auto str = path.generic_string();
		return (detail::_Imread(str.c_str(), _Ty()));
	}
	else return (detail::_Imread(path, _Ty()));
}
}
DGE_MATRICE_END