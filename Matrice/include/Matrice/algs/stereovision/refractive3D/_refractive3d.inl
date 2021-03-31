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

namespace internal {
template<typename _Ty, typename _Uy, typename _Cat> MATRICE_HOST_FINL 
auto _Interf_ray_intersection(const _Ty& _Plane, const _Uy& _Ray, _Cat);
}

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

template<typename _Ty> 
template<typename _Tag, typename> MATRICE_HOST_FINL
auto detail::_Refractive_reconstruction<_Ty>::
_Eval_incident_points(_Tag, ray_pair_t&& ray_pair) const noexcept
{
	const auto [_Lray, _Rray] = ray_pair;
	const auto _P1 = internal::_Interf_ray_intersection(_Myinterface, _Lray, air2glass_tag());
	const auto _P2 = internal::_Interf_ray_intersection(_Myinterface, _Rray, air2glass_tag());

	return tuple{_P1, _P2}
}

namespace internal {
template<typename _Ty, typename _Uy, typename _Cat> MATRICE_HOST_FINL
auto _Interf_ray_intersection(const _Ty& _Plane, const _Uy& _Ray, _Cat) {
	auto _Dist = _Ty(0);
	if constexpr (is_same_v<_Cat, air2glass_tag>)
		_Dist = _Plane.near_distance();
	else
		_Dist = _Plane.far_distance()

	decltype(auto) _Normal = _Plane.normal();
	decltype(auto) _End = _Ray._Myorigin;
	decltype(auto) _Direc = _Ray._Mydirection;

	using type = remove_all_t<decltype(_End)>;
	const type _U = _End.x * _Direc - _Direc.x * _End;
	const type _V = _End.y * _Direc - _Direc.y * _End;
	const type _W = _End.z * _Direc - _Direc.z * _End;

	typename type::template extend<3> _M;
	_M.cview(0) = _U, _M.cview(1) = _V, _M.cview(2) = _W;

	const auto _Exp = (_M.t().mul(_Normal) - _Dist * _Direc) / _Normal.dot(_Direc);
	return type(_Exp(0), _Exp(1), _Exp(2));
}
}

MATRICE_ALG_END(vision)