/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2021, Zhilong(Dgelom) Su, all rights reserved.

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
**********************************************************************/
#pragma once

#include <core/vector.h>
#include <algs/geometry.hpp>

DGE_MATRICE_BEGIN
namespace example {
template<typename _Ty, size_t _Dim>
class rotation {
public:
	using value_t = _Ty;
	using vector_t = auto_vector_t<value_t, _Dim>;

	explicit rotation(const vector_t& data) noexcept {
		_Data = data;
	}
	explicit rotation(const vector_t& data, value_t scale) noexcept {
		_Data = data*scale;
	}

	inline auto matrix() const noexcept {
		return rodrigues(_Data);
	}

private:
	vector_t _Data;
};
}
DGE_MATRICE_END