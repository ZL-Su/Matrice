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
#include "forward.hpp"
#include "../../math/kernel_wrapper.hpp"

#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL

DGE_MATRICE_BEGIN
_INTERNAL_BEGIN
template<typename _Ptr>
lapack_int _lapack_syev(int lyt, char job, char ul, lapack_int n, _Ptr a, lapack_int lda, _Ptr w);
_INTERNAL_END

_DETAIL_BEGIN
struct _Lapack_kernel_wrapper {
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