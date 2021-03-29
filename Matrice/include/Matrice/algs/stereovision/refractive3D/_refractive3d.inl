/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library for
3D Vision and Photo-Mechanics.
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
#include "_refractive3d.h"

MATRICE_ALG_BEGIN(vision)

template<typename _Ty> MATRICE_HOST_FINL
auto detail::_Refractive_reconstruction<_Ty>::
_Eval_incident_directions(const vector_t<3>& _P) const noexcept 
{
	vector_t<3> _Diff_po = _P - _Myorigin;
	vector_t<3> _Diff_pt = _P - _Mytrans;

	_Diff_po = _Diff_po.normalize(_Diff_po.norm());
	_Diff_pt = _Diff_pt.normalize(_Diff_pt.norm());

	return tuple{ _Diff_po, _Diff_pt };
}

template<typename _Ty> template<typename _Tag, typename> MATRICE_HOST_FINL
auto detail::_Refractive_reconstruction<_Ty>::_Eval_intersections(_Tag, const ray_type& l) const noexcept
{

}

MATRICE_ALG_END(vision)