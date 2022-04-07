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
***********************************************************************/
#pragma once

#include "core/vector.h"
#include "algs/dnn/functions.h"
#include "algs/dnn/modules.h"

MATRICE_NAMESPACE_BEGIN(dnn)
namespace detail {

/// <summary>
/// \brief CLASS template, _PlainSelfAttention interface
/// 	-- a trivial implementation for self-attention.
/// </summary>
/// <typeparam name="_Input"></typeparam>
template<class _Input, size_t _Embed, size_t _Depth> 
class _PlainSelfAttention{
	static_assert(std::true_type::value, "Unsupported type '_Input'\
		in _PlainSelfAttention<_Input, _Embed, _Depth>.");
};

/// <summary>
/// \brief CLASS, specialization for a single vector input.
/// </summary>
/// <typeparam name="_Input">Vector or Matrix with cols or rows is one.</typeparam>
template<class _Input, size_t _Embed>
class _PlainSelfAttention<_Input, _Embed, 1>: public xpr::__xpr__ {
public:
	enum{depth = 1};
	using input_t = input_layer<_Input>;
	using linear_t = linear_layer<input_t, _Embed>;
	using value_t = _Input::value_t;
	using value_type = value_t;

	_PlainSelfAttention() noexcept 
		:_Myinput{} {
	}
	_PlainSelfAttention(const input_t& x) noexcept
		:_Myinput{ x } {
		_MyWQ = linear_t()(_Myinput);
		_MyWK = linear_t()(_Myinput);
		_MyWV = linear_t()(_Myinput);
	}

	MATRICE_GLOBAL_INL auto forward() const {
		const auto _Score = _MyWQ.dot(_MyWK);
		Vec_<value_t, depth> _Weights(_Score);
		using SoltMax_t = functional::softmax<decltype(_Weights)>;
		_Weights = SoltMax_t::forward(_Weights);
		return (_Weights(0) * _MyWV).eval();
	}

protected:
	input_t _Myinput;
	Matrix_<value_type, 1, _Embed> _MyWQ;
	Matrix_<value_type, 1, _Embed> _MyWK;
	Matrix_<value_type, 1, _Embed> _MyWV;
};
}

template<class _Input, size_t _Embed>
using plain_self_attention_t = detail::_PlainSelfAttention<_Input, _Embed,
	conditional_size_v<_Input::rows_at_compiletime==1, 1,
	conditional_size_v<_Input::cols_at_compiletime==1, 1,
	_Input::rows_at_compiletime>>>;
MATRICE_NAMESPACE_END(dnn)