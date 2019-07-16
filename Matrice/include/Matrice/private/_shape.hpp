/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
#include <array>
#include "../util/_macros.h"
#include "../util/_std_wrapper.h"
#include "../util/_exception.h"

DGE_MATRICE_BEGIN
template<typename _Ity = size_t>
class shape_ : public std::array<_Ity, 4>
{
	using _Mybase = std::array<_Ity, 4>;
public:
	using index_type = typename _Mybase::value_type;

	shape_(initlist<index_type> shape) noexcept {
	}

	const index_type& height() const noexcept {
		return (*this)[0];
	}
	index_type& height() noexcept {
		return (*this)[0];
	}
	const index_type& width() const noexcept {
		return (*this)[1];
	}
	index_type& width() noexcept {
		return (*this)[1];
	}
	const index_type& depth() const noexcept {
		return (*this)[2];
	}
	index_type& depth() noexcept {
		return (*this)[2];
	}
};
DGE_MATRICE_END
