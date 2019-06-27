/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
***********************************************************************/
#pragma once
#include "_layer_base.hpp"

DGE_MATRICE_BEGIN
namespace dnn {

	template<typename _Ity> class _Input_layer {
		static_assert(sizeof(_Ity) != sizeof(_Ity),
			"Unsupported type given to _Input_layer<>, which only supports dgelom::Matrix_, dgelom::tensor or array of dgelom::Matrix_.");
	};

	template<typename _Ty, int _M, int _N>
	class _Input_layer<Matrix_<_Ty, _M, _N>> {
		using _Myt = _Input_layer;
	public:
		using input_type = Matrix_<_Ty, _M, _N>;
		using category = typename _Layer_tag::input;

		_Input_layer() noexcept {}
		_Input_layer(const _Myt&) noexcept {}

		template<typename _Inty>
		_Input_layer(const _Inty&) noexcept {}

	};

	template<typename _Ty, int _M, int _N, size_t _D>
	class _Input_layer<std::array<Matrix_<_Ty, _M, _N>, _D>> {
		using _Myt = _Input_layer;
	public:
		using input_type = std::array<Matrix_<_Ty, _M, _N>, _D>;
		using category = typename _Layer_tag::input;

		_Input_layer() noexcept {}
		_Input_layer(const _Myt&) noexcept {}

		template<typename _Inty>
		_Input_layer(const _Inty&) noexcept {}

	};
}
DGE_MATRICE_END