/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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
**************************************************************************/
#pragma once
#include <tuple>
#include "../_plain_base.hpp"
#include "../nonfree/_lnalge.h"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ty, 
	typename value_t = typename _Ty::value_t> 
class _Svd_impl
{
	using matrix_type = _Ty;
	using traits_type = matrix_traits<matrix_type>;
	using size_traits = typename traits_type::size;
	enum {
		CompileTimeRows = size_traits::rows::value,
		CompileTimeCols = size_traits::cols::value
	};

public:
	MATRICE_GLOBAL_INL _Svd_impl(const matrix_type& _A) : _A_view(_A) {
		auto [_A_view, _S, _Vt] = op(_A_view);
	}
	MATRICE_GLOBAL_INL auto operator()() {
		return std::make_tuple(_A_view, std::ref(_S), std::ref(_Vt));
	}

	/**
	 *\Static svd operator, which is thread-safe totally.
	 */
	MATRICE_GLOBAL_INL static auto op(const matrix_type& _A) {
		types::Matrix_<value_t, CompileTimeCols, 1> S;
		types::Matrix_<value_t, CompileTimeCols, CompileTimeCols> Vt;
		detail::_Lapack_kernel_impl<value_t>::svd(_A.data(), S.data(), Vt.data(), _A.shape());

		return std::make_tuple(std::ref(_A), std::ref(S), std::ref(Vt));
	}
private:
	const matrix_type& _A_view;
	types::Matrix_<value_t, CompileTimeCols, 1> _S;
	types::Matrix_<value_t, CompileTimeCols, CompileTimeCols> _Vt;

};
_DETAIL_END

DGE_MATRICE_END