/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

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
*********************************************************************/
#pragma once
#include "matrix.h"
#include "algs/forward.hpp"

MATRICE_ALGS_BEGIN

MATRICE_ALGS_END

DGE_MATRICE_BEGIN
/// <summary>
/// \brief ENUM, imgae_func_type, for specify the image function.
/// </summary>
enum class image_func_type {
	identity,
	bilinear,
	bplineb3,
	bplineb5,
	bplineb7,
};

_DETAIL_BEGIN

/// <summary>
/// \brief CLASS TEMPLATE, Smooth image type
/// </summary>
/// <typeparam name="_Ty">Floating point type</typeparam>
/// <typeparam name="_Ip">Interpolation operator type</typeparam>
template<typename _Pre = bicerp_tag, typename _Ty = float>
class _Image {
	using _Myt = _Image;
	using _Myctx = typename algs::auto_interp_dispatcher<_Ty, _Pre>::type;
public:
	using matrix_type = typename _Myctx::matrix_type;
	using value_type = typename matrix_type::value_type;

	_Image(const matrix_type& _Src)
		:_Context(std::make_shared<_Myctx>(_Src)) {
	}

	/**
	 * \brief OP, eval image value at coordinates (x, y). 
	 */
	template<typename _Uy>
	decltype(auto) operator()(_Uy x, _Uy y) const {
		return (*_Context)({ value_type(x), value_type(y) });
	}
	template<typename _Uy>
	decltype(auto) operator()(_Uy x, _Uy y) {
		return (*_Context)({ value_type(x), value_type(y) });
	}

private:
	shared_ptr<_Myctx> _Context;
};
_DETAIL_END

/// <summary>
/// \brief ALIAS TEMPLATE, smooth image type
/// </summary>
/// <typeparam name="_Ty"></typeparam>
/// <typeparam name="_Func"></typeparam>
template<typename _Ty, typename _Func = bicerp_tag>
using image_t = detail::_Image<_Func, _Ty>;

DGE_MATRICE_END