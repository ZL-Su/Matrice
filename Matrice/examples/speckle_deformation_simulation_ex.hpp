
#pragma once
#include "../include/io"
#include "../include/core"
#include "../include/Matrice/algs/interpolation.h"

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
	void load_src(_Fn&& op) noexcept {
		m_original = std::move(op);
		const auto x = (m_original.cols() - m_cols) / 2;
		const auto y = (m_original.rows() - m_rows) / 2;
		m_rect.set({ x, y }, m_cols, m_rows);
	}

	template<uint16_t _Order>
	matrix_type random_warp() noexcept {
		
		matrix_type dst(m_rows, m_cols);
		if constexpr (_Order == static_cast<uint16_t>(0)) {
			for (const auto& idx : range(0, dst.size())) {
				const auto r = idx / m_cols;
				const auto c = idx - r * m_cols;

			}
		}
		else if constexpr(_Order == static_cast<uint16_t>(1)) {

		}
		else if constexpr (_Order == static_cast<uint16_t>(2)) {

		}
		else {
			static_assert("Unsupported deformation order.")
		}

		return dst;
	}

private:
	size_t m_rows, m_cols;
	rect_type m_rect;
	matrix_type m_original;
	interp_type m_srcinst;
};
DGE_MATRICE_END