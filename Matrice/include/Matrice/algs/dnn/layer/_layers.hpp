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
	struct convaxis {
		struct x { static constexpr auto value = 1; };
		struct y { static constexpr auto value = 0; };
	};

	struct conv_out_depth {
		MATRICE_HOST_INL conv_out_depth(uint32_t d) noexcept : value(d) {}
		uint32_t value = 1;
	};

	template<uint32_t _D,
		uint32_t _FH, uint32_t _FW = _FH,
		uint32_t _Sy = 1, uint32_t _Sx = _Sy,
		int _Py = conditional_size_v<_Sy != 1, 0, _FH / 2>,
		int _Px = conditional_size_v<_Sx != 1, 0, _FW / 2>>
	class _Conv2d_layer 
		: public _Layer<_Conv2d_layer<_D,_FH,_FW,_Sy,_Sx,_Py,_Px>>
	{
	public:
		static constexpr auto depth = _D;
		static constexpr auto k_rows = _FH, k_cols = _FW;
		static constexpr auto stride_y = _Sy, stride_x = _Sx;
		static constexpr auto padding_y = _Py, padding_x = _Px;

		using category = typename _Layer_tag::conv2d;
		using value_type = dnn_default_value_type;

		_Conv2d_layer(conv_out_depth cod) noexcept
			: m_nkernels(cod.value) {
			DGELOM_CHECK(m_nkernels > 0,
				"the number of kernels must be greater than 0.");
		}
		_Conv2d_layer() noexcept
			: _Conv2d_layer(conv_out_depth(depth)) {
		}



	private:
		uint32_t m_nkernels = depth;
		value_type m_lr_coef = 1;
		value_type m_decay_coef = 1;
		
		Tensor<value_type, depth> weights, biases;
	};
}
DGE_MATRICE_END