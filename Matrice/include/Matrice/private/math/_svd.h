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
#include "../_matrix_base.hpp"
#ifdef __use_mkl__
#include <mkl.h>
#else
#include <fkl.h>
#endif

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
	MATRICE_GLOBAL_INL _Svd_impl(const matrix_type& _A) : _A_view(_A) {}
	MATRICE_GLOBAL_INL auto operator()() {
		return std::make_tuple(_A_view, std::ref(_S), std::ref(_V));
	}

private:
	const matrix_type& _A_view;
	types::Matrix_<value_t, CompileTimeCols, 1> _S;
	types::Matrix_<value_t, CompileTimeCols, CompileTimeCols> _V;

};
template<typename _Ty>
class _Svd_impl<_Ty, float>
{
	static_assert(std::is_same_v<typename _Ty::value_t, float >, "");
	using value_t = float;
	using pointer = std::add_pointer_t<value_t>;
	using matrix_type = _Ty;
	using traits_type = matrix_traits<matrix_type>;
	using size_traits = typename traits_type::size;
	enum {
		CompileTimeRows = size_traits::rows::value,
		CompileTimeCols = size_traits::cols::value
	};

public:
	MATRICE_GLOBAL_INL _Svd_impl(const matrix_type& _A) : _A_view(_A) {
		_S.create(_A_view.cols(), 1);
		_V.create(_A_view.cols(), _A_view.cols());
	}
	MATRICE_GLOBAL_INL auto operator()() {
		flapk::_sgesvd(_A_view.data(), _S.data(), _V.data(), _A_view.rows(), _A_view.cols());
		return std::make_tuple(_A_view, std::ref(_S), std::ref(_V));
	}
	/**
	 *\Static svd operator, which is thread-safe totally.
	 */
	MATRICE_GLOBAL_INL static auto op(const matrix_type& _A) {
		types::Matrix_<value_t, CompileTimeCols, 1> S(_A.rows(), _A.cols());
		types::Matrix_<value_t, CompileTimeCols, CompileTimeCols> V(_A.rows(), _A.cols());
		flapk::_sgesvd(_A.data(), S.data(), V.data(), _A.rows(), _A.cols());
		return std::make_tuple(_A, std::ref(S), std::ref(V));
	}

private:
	const matrix_type& _A_view;
	types::Matrix_<value_t, CompileTimeCols, 1> _S;
	types::Matrix_<value_t, CompileTimeCols, CompileTimeCols> _V;

};
template<typename _Ty>
class _Svd_impl<_Ty, double>
{
	static_assert(std::is_same_v<typename _Ty::value_t, double>, "");
	using value_t = double;
	using pointer = std::add_pointer_t<value_t>;
	using matrix_type = _Ty;
	using traits_type = matrix_traits<matrix_type>;
	using size_traits = typename traits_type::size;
	enum {
		CompileTimeRows = size_traits::rows::value,
		CompileTimeCols = size_traits::cols::value
	};

public:
	MATRICE_GLOBAL_INL _Svd_impl(const matrix_type& _A) : _A_view(_A) {
		_S.create(_A_view.cols(), 1);
		_V.create(_A_view.cols(), _A_view.cols());
	}
	MATRICE_GLOBAL_INL auto operator()() {
		flapk::_dgesvd(_A_view.data(), _S.data(), _V.data(), _A_view.rows(), _A_view.cols());
		return std::make_tuple(_A_view, std::ref(_S), std::ref(_V));
	}
	/**
	 *\Static svd operator, which is thread-safe totally.
	 */
	MATRICE_GLOBAL_INL static auto op(const matrix_type& _A) {
		types::Matrix_<value_t, CompileTimeCols, 1> S(_A.rows(), _A.cols());
		types::Matrix_<value_t, CompileTimeCols, CompileTimeCols> V(_A.rows(), _A.cols());
		flapk::_dgesvd(_A.data(), S.data(), V.data(), _A.rows(), _A.cols());
		return std::make_tuple(_A, std::ref(S), std::ref(V));
	}

private:
	const matrix_type& _A_view;
	types::Matrix_<value_t, CompileTimeCols, 1> _S;
	types::Matrix_<value_t, CompileTimeCols, CompileTimeCols> _V;

};
_DETAIL_END

DGE_MATRICE_END