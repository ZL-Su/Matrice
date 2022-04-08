/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2022, Zhilong(Dgelom) Su, all rights reserved.

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
********************************************************************/
#pragma once

#include "core/matrix.h"

DGE_MATRICE_BEGIN

namespace tag {
	struct lu{};
	struct qr{};
	struct spd{};
	struct svd{};
	struct eig{};
}

namespace detail {
	template<class _Exp> struct _Linear_decomp_base {};
	template<typename _Ty, size_t _M, size_t _N>
	struct _Linear_decomp_base<Matrix_<_Ty, _M, _N>> {
		using value_t = remove_all_t<_Ty>;
		constexpr static auto rows_at_compiletime = _M;
		constexpr static auto cols_at_compiletime = _N;
	};
}

/// <summary>
/// \brief CLASS TEMPLATE linear_decomposer, performs linear 
/// decomposition on a given expression with a specified
/// linear operator.
/// </summary>
/// <typeparam name="_Exp"></typeparam>
/// <typeparam name="_Op"></typeparam>
template<class _Exp, class _Op> class linear_decomposer {};

template<class _Exp>
class linear_decomposer<_Exp, tag::spd> 
	: public detail::_Linear_decomp_base<_Exp>
{
	using _Mybase = detail::_Linear_decomp_base<_Exp>;
public:
	using _Mybase::value_t;
	linear_decomposer(const _Exp& A) noexcept
		:_Myc(A) {
	};

	MATRICE_HOST_FINL ~linear_decomposer()noexcept {};



private:
	Matrix_<value_t, 
		_Mybase::rows_at_compiletime, 
		_Mybase::rows_at_compiletime> _Myc;
};
DGE_MATRICE_END