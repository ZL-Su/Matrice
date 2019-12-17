/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for
more detail.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once
#include "../_matrix_ops.hpp"
#include "../math/_linear_kernel.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<class _Mty>
class _Matrix_fact<_Mty, tag::_Linear_spd_tag> {
	using matrix_type = _Mty;
public:
	_Matrix_fact(_Mty& a) : m_ret(a) {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(a.rows() == a.cols(), "a must be a square matrix in _Matrix_linear_fact<>.");
#endif
		_Linear_spd_kernel(a.data(), a.rows());
	}

	MATRICE_GLOBAL_INL matrix_type& operator()() noexcept {
		return (m_ret);
	}
	MATRICE_GLOBAL_INL const matrix_type& operator()() const noexcept {
		return (m_ret);
	}

	MATRICE_GLOBAL_INL matrix_type inv() noexcept {
		matrix_type inv(m_ret.shape());
		_Linear_ispd_kernel(m_ret.data(), inv.data(), inv.rows());
		return inv;
	}

private:
	matrix_type& m_ret;
};
_DETAIL_END
DGE_MATRICE_END
