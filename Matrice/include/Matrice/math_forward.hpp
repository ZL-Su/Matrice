/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

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

//\note: forward declarations of math operations and kernels

namespace dgelom {

namespace detail{template<class _Mty,typename _Op>class _Linear_solver;}

template<class _Mty, typename _Op>
using linear_solver_t = detail::_Linear_solver<_Mty, _Op>;

#define MATRICE_MAKE_LINEARSV_TYPE(OP) \
struct OP; \
template<class _Mty> \
using OP##_linear_solver_t = linear_solver_t<_Mty, OP>;

// \brief: Linear solver with SVD kernel
MATRICE_MAKE_LINEARSV_TYPE(svd);
// \brief: Linear solver with SPT kernel
MATRICE_MAKE_LINEARSV_TYPE(spt);

#undef MATRICE_MAKE_LINEARSV_TYPE
}
