/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU
General Public License for more detail.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once
#include "forward.hpp"
#include "../../math/_primitive_funcs.hpp"
#include "../../math/kernel_wrapper.hpp"

#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL

DGE_MATRICE_BEGIN
_INTERNAL_BEGIN
template<typename _Ptr>
lapack_int _lapack_potrf(int lyt, char ul, lapack_int n, _Ptr a, lapack_int lda);
template<typename _Ptr>
lapack_int _lapack_gesvd(int lyt, char jobu, char jobvt, lapack_int m, lapack_int n, _Ptr a, lapack_int lda, _Ptr s, _Ptr u, lapack_int ldu, _Ptr vt, lapack_int ldvt, _Ptr superb);
template<typename _Ptr>
lapack_int _lapack_syev(int lyt, char job, char ul, lapack_int n, _Ptr a, lapack_int lda, _Ptr w);
_INTERNAL_END

_DETAIL_BEGIN
struct _Lapack_kernel_wrapper {
	/**
	 *\brief computes the Cholesky factorization of a symmetric positive-definite matrix.
	 *\param [a] the matrix to be factorized;
	         [lu] optional par to indicate which part of a is stored, must be 'L' for the lower part or 'U' for the upper part. 
	 */
	template<class _Mty>
	static MATRICE_HOST_INL _Mty& spotrf(_Mty& a, char lu = 'L') {
		const auto lyt = a.allocator().fmt();
		auto sts = internal::_lapack_potrf(lyt, lu, a.rows(), a.data(), a.cols());
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(sts == 0, "Fail to decompose a with error code: " + cast_to_string(sts));
#endif
		return a;
	}

	/**
	 *\brief: computes SVD decomposition of a general matrix
	 *\param: [a] the matrix to be decomposed.
	          [eval_vt] specifies if Vt is evaluated or not. If true, the right singular vectors will be row-wisely stored in [a].
	 *\return: [s, sts] singular values and status flag.
	 */
	template<class _Mty>
	static MATRICE_HOST_INL auto gesvd(_Mty& a, bool eval_vt = false) {
		using value_type = typename remove_all_t<_Mty>::value_type;
		const auto lyt = a.allocator().fmt();
		const char jobu = 'N', jobvt = eval_vt? 'O' : 'N';
		constexpr auto null = (value_type*)(nullptr);
		Matrix_<value_type, _Mty::rows_at_compiletime, 1> s(diff_t(a.rows()));
		Matrix_<value_type, min(_Mty::rows_at_compiletime, _Mty::cols_at_compiletime), 1> superb(min<diff_t>(a.rows(), a.cols()));
		const auto sts = internal::_lapack_gesvd(lyt, jobu, jobvt, a.rows(), a.cols(), a.data(), a.cols(), s.data(), null, 1, null, 1, superb.data());
		return std::make_tuple(s, sts);
	}
	template<class _Mty>
	static MATRICE_HOST_INL auto gesvd(_Mty& a, bool eval_u, bool eval_vt) {
		using value_type = typename remove_all_t<_Mty>::value_type;
		const auto lyt = a.allocator().fmt();
		const char jobu = eval_u ? 'S' : 'N', jobvt = eval_vt ? 'O' : 'N';
		constexpr auto null = (value_type*)(nullptr);
		Matrix_<value_type, _Mty::rows_at_compiletime, 1> s(diff_t(a.rows()));
		Matrix_<value_type, min(_Mty::rows_at_compiletime, _Mty::cols_at_compiletime), 1> superb(min<diff_t>(a.rows(), a.cols()));
		Matrix_<value_type, _Mty::rows_at_compiletime> u(a.rows(), a.rows());
		const auto sts = internal::_lapack_gesvd(lyt, jobu, jobvt, a.rows(), a.cols(), a.data(), a.cols(), s.data(), u.data(), u.cols(), null, 1, superb.data());
		return std::make_tuple(s, u, sts);
	}
	template<class _Mty>
	static MATRICE_HOST_INL auto syev(_Mty& a, bool eval_vecs = false) {
		using value_type = typename remove_all_t<_Mty>::value_type;
		const auto jobz = eval_vecs ? 'V' : 'N';
		const auto lyt = a.allocator().fmt();
		Matrix_<value_type, _Mty::rows_at_compiletime, 1> w(diff_t(a.rows()));
		const auto sts = internal::_lapack_syev(lyt, jobz, 'L', a.rows(), a.data(), max(1, a.rows()), w.data());
		return forward<decltype(w)>(w);
	}
};
_DETAIL_END

_INTERNAL_BEGIN

_INTERNAL_END
DGE_MATRICE_END

#endif