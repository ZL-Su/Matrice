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
#include "../models/_basic_model.hpp"


DGE_MATRICE_BEGIN
namespace dnn {
/**
 *\brief network layer traits
 *\param <_Lyr> layer type
 */
template<typename _Lyr> struct _Layer_traits : traits<_Lyr> {};

/**
 *\brief define layer type
 */
struct _Layer_tag {
	struct linear {};  struct conv2d {};
	struct maxpool {}; struct avgpool {};
	struct input {};   struct output {};
};

/**
 *\brief network layer base class
 *\param <_Derived> derived layer type
 */
template<typename _Derived, typename _Traits = _Layer_traits<_Derived>>
class _Layer : public detail::_Model<typename _Traits::value_type, _Traits::depth, _Traits::extent, _Traits::has_bias> {
	using _Myt = _Layer;
	using _Mydt = _Derived;
	using _Mytraits = _Traits;
	using _Mybase = detail::_Model<typename _Traits::value_type, 
		_Traits::depth, _Traits::extent, _Traits::has_bias>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::tensor_type;

	MATRICE_HOST_INL auto forward() noexcept {
		//return static_cast<_Mydt*>(this)->_Forward_impl();
	}

	MATRICE_HOST_INL auto backward() noexcept {
		//return static_cast<_Mydt*>(this)->_Backward_impl();
	}

protected:
	tensor_type m_params;
};
}
DGE_MATRICE_END