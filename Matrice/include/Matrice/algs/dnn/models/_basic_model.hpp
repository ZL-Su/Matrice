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
#include "../../../core"

DGE_MATRICE_BEGIN
namespace dnn {

	using dnn_default_value_type = float;

_DETAIL_BEGIN
template<typename _Ty, uint32_t _D, uint32_t _E, bool _Req_biases> 
class _Model {
	static_assert("Undefined dnn::detail::_Model<>.");
};

template<typename _Ty, uint32_t _D, uint32_t _E>
class _Model<_Ty, _D, _E, std::false_type::value> {
public:
	using value_type = _Ty;
	using tensor_type = Tensor<value_type, _D, _E>;

protected:
	tensor_type m_weights;
};

template<typename _Ty, uint32_t _D, uint32_t _E>
class _Model<_Ty, _D, _E, std::true_type::value> {
public:
	using value_type = _Ty;
	using tensor_type = Tensor<value_type, _D, _E>;

protected:
	tensor_type m_weights, m_biases;
};
_DETAIL_END
}
DGE_MATRICE_END