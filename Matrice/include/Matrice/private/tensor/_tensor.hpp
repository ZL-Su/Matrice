/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once

#include <valarray>
#include "../_matrix_base.hpp"
#include "../_matrix_exp.hpp"

DGE_MATRICE_BEGIN namespace detail {

template<typename _Ty, int _M = 0, int _N = 0, std::size_t _K = 0,
	typename matrix_type = types::Matrix_<_Ty, _M, _N>>
class _Tensor_impl MATRICE_NONHERITABLE : public std::valarray<matrix_type> {
	using base_t = std::valarray<matrix_type>;
	using _Myt_traits = matrix_traits<matrix_type>;
public:
	using type = matrix_type;
	using value_type = typename _Myt_traits::type;
	using value_t = value_type;

	_Tensor_impl(std::size_t _Rows)
		: base_t(m_size = _Rows), m_rows(_Rows), m_cols(1) {}
	_Tensor_impl(std::size_t _Rows, const matrix_type& _Mat)
		: base_t(_Mat, m_size = _Rows), m_rows(_Rows), m_cols(1) {}
	_Tensor_impl(std::size_t _Rows, std::size_t _Cols) 
		: base_t(m_size = _Rows*_Cols), m_rows(_Rows), m_cols(_Cols) {}
	_Tensor_impl(std::size_t _Rows, std::size_t _Cols, const matrix_type& _Mat) : base_t(_Mat, m_size = _Rows * _Cols), m_rows(_Rows), m_cols(_Cols) {}

private:
	std::size_t m_rows, m_cols, m_size;
};


} DGE_MATRICE_END