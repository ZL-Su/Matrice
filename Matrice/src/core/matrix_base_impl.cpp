/**********************************************************************
    This file is part of Matrice, an effcient and elegant C++ library.
        Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

This program is free software :  you can redistribute it and/or modify
it  under the terms of the  GNU General Public License as published by
the  Free Software Foundation,  either  version  3  of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,    but
WITHOUT ANY WARRANTY;   without even the implied warranty of MERCHANT-
ABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public
License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/
#include "../../include/Matrice/core/matrix.h"
#ifdef __use_mkl__
#include <mkl.h>
#else
#include <fkl.h>
#endif

MATRICE_NAMESPACE_BEGIN_

template<typename _Rhs, typename value_t> value_t det_impl(const _Rhs & a)
{
	const int N = a.rows();
	const auto p = a.data();
	if (N == 2) return (p[0]*p[3] - p[1]*p[2]);
	if (N == 3) return (p[0]*(p[4]*p[8] - p[5]*p[7]) - p[1]*(p[3]*p[8] - p[5]*p[6]) + p[2]*(p[3]*p[7] - p[4]*p[6]));

	if constexpr (type_bytes<value_t>::value == 4)
		return fblas::_sdetm(fkl::sptr(p), a.rows());
	if constexpr (type_bytes<value_t>::value == 8)
		return fblas::_ddetm(fkl::dptr(p), a.rows());
	
	return std::numeric_limits<value_t>::infinity();
}
template float det_impl(const types::Matrix_<float, 2, 2>&);
template float det_impl(const types::Matrix_<float, 3, 3>&);
template float det_impl(const types::Matrix_<float, 0, 0>&);
template double det_impl(const types::Matrix_<double, 2, 2>&);
template double det_impl(const types::Matrix_<double, 3, 3>&);
template double det_impl(const types::Matrix_<double, 0, 0>&);

_MATRICE_NAMESPACE_END