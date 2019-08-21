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
#include "algs/correlation/_optim.h"
#include "_projection.hpp"

MATRICE_ALGS_BEGIN
_DETAIL_BEGIN
template<typename _Ty>
class _Ewise_depth_estimation {
	using _Myt = _Ewise_depth_estimation;
	using _Projection = _Aligned_projection<_Ty, cs_alignment_tag::left>;
public:
	using point_type = typename _Projection::point_type;
	_Ewise_depth_estimation(const point_type& r, const point_type& t)
		:m_projection(r,t){
	}

private:
	_Projection m_projection;
};
_DETAIL_END
MATRICE_ALGS_END
