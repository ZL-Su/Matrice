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
**********************************************************************/
#pragma once
#include <optional>
#include <io/io.hpp>
#include <core/matrix.h>
#include <algs/interpolation.h>

DGE_MATRICE_BEGIN
template<typename _Ty = float>
class speckle_generator {
	using _Myt = speckle_generator;
public:
	using value_type = _Ty;
	using matrix_type = Matrix<value_type>;
	using interp_type = interpolation<value_type>;
	using rect_type = rect<size_t>;

	speckle_generator(size_t rows, size_t cols) noexcept
		:m_rows(rows), m_cols(cols), m_srcinst(m_original) {
	}

	template<typename _Fn>
	inline void load_src(_Fn&& op) noexcept {
		const auto src = op();
		if constexpr (is_same_v<decltype(src), matrix_type>) {
			m_original = move(src);
		}
		else {
			m_original.create(src.shape());
			for (const auto& idx : range(0, src.size()))
				m_original(idx) = src(idx);
		}
		const auto x = (m_original.cols() - m_cols) / 2;
		const auto y = (m_original.rows() - m_rows) / 2;
		m_rect.set({ x, y }, m_cols, m_rows);
	}

	inline decltype(auto) original() noexcept {
		return (m_original);
	}

	inline decltype(auto) get_block_range() const noexcept {
		return std::make_tuple(m_rect.begin().x, m_rect.end().x, m_rect.begin().y, m_rect.end().y);
	}

	template<uint16_t _Order>
	inline std::optional<matrix_type> random_warp() noexcept {
		mt19937 eng;
		normal_distribution<value_type> ndist;
		interpolation<value_type, bilerp_tag> itp(m_original);

		matrix_type dst(m_rows, m_cols);
		if constexpr (_Order == static_cast<uint16_t>(0)) {
			for (const auto& idx : range(0, dst.size())) {
				auto [r, c] = unroll_linear_index<int>(idx, m_cols);
				const auto u = ndist(eng), v = ndist(eng);
				const auto x = c + m_rect.begin().x - u;
				const auto y = r + m_rect.begin().y - v;
				dst[r][c] = itp(x, y);
			}
		}
		else if constexpr(_Order == static_cast<uint16_t>(1)) {

		}
		else if constexpr (_Order == static_cast<uint16_t>(2)) {

		}
		else {
			static_assert("Unsupported deformation order.");
		}

		return std::make_optional(dst);
	}

private:
	size_t m_rows, m_cols;
	rect_type m_rect;
	matrix_type m_original;
	interp_type m_srcinst;
};
DGE_MATRICE_END