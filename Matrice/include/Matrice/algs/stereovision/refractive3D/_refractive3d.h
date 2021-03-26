/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library for 
3D vision and photo-mechanics.
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
#include "../_camera.hpp"

MATRICE_ALG_BEGIN(vision)

inline _DETAIL_BEGIN
template<typename _Ty>
class _Refractive_reconstruction {
	using _Myt = _Refractive_reconstruction;
public:
	using value_type = _Ty;
	template<size_t _Dim>
	using vector_t = auto_vector_t<value_type, _Dim>;
	struct interface_type
	{
		vector_t<3> _Mynormal;
		value_type  _Mydist{ 0 };
		value_type  _Myshift{ 0 };
	};

	explicit _Refractive_reconstruction() noexcept {
	}

private:
	interface_type _Myinterface;
};
_DETAIL_END

template<typename _Ty>
using refractive_reconstruction = detail::_Refractive_reconstruction;

MATRICE_ALG_END(vision)