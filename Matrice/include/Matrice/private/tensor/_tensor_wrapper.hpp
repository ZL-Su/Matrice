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
#include "_tensor.hpp"

DGE_MATRICE_BEGIN

/* tensor instance */
template<typename _Ty, typename dev_tag> 
struct tensor_inst {};

/* tensor descriptor */
struct tensor_desc {
	uint32_t depth;
	uint32_t extent;
	size_t height, width;
};

template<typename _Ty>
class tensor_view {
public:
	using value_type = _Ty;
	using pointer = std::add_pointer_t<value_type>;

	MATRICE_GLOBAL_INL const pointer data() const noexcept {
		return m_data;
	}
	MATRICE_GLOBAL_INL pointer data() noexcept {
		return m_data;
	}
	MATRICE_GLOBAL_INL const tensor_desc& desc() const noexcept {
		return m_desc;
	}
	MATRICE_GLOBAL_INL tensor_desc& desc() noexcept {
		return m_desc;
	}
private:
	tensor_desc m_desc;
	pointer m_data;
};
DGE_MATRICE_END