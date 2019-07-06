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
struct cs_alignment_tag { 
	struct left {}; struct middle {}; struct right{};
};
_DETAIL_BEGIN
template<typename _Ty> class _Projection_base {
	using _Myt = _Projection_base<_Ty>;
public:
	static constexpr auto dim = 3;
	using value_type = _Ty;
	using point_type = Vec3_<value_type>;
	using matrix_type = Matrix_<value_type, dim, dim>;
	template<typename _An> struct depth_paramterized_projection {
		point_type point;
		MATRICE_HOST_INL point_type grad() noexcept {
			return point_type{};
		}
	};
	/**
	 *\brief initialize the rotation and translation
	 *\param [_Ext] external geometry parameters: rx, ry, rz, tx, ty, tz
	 */
	_Projection_base(const array_n<value_type, dim<<1>& _Ext) noexcept
		:_MyT(_Ext(3), _Ext(4), _Ext(5)) {
		point_type r{ _Ext(0), _Ext(1), _Ext(2) };
		const auto theta = sqrt(r.dot(r));
		const auto sval = sin(theta);
		auto cval = cos(theta);

		_MyR.identity(); cval = 1 - cval; r = r / theta;

		_MyR(0) += cval * r[0] * r[0];
		_MyR(1) += cval * r[0] * r[1];
		_MyR(2) += cval * r[0] * r[2];
		_MyR(3) += cval * r[1] * r[0];
		_MyR(4) += cval * r[1] * r[1];
		_MyR(5) += cval * r[1] * r[2];
		_MyR(6) += cval * r[2] * r[0];
		_MyR(7) += cval * r[2] * r[1];
		_MyR(8) += cval * r[2] * r[2];

		_MyR(1) -= sval * r[2], _MyR(2) += sval * r[1];
		_MyR(3) += sval * r[2], _MyR(5) -= sval * r[0];
		_MyR(6) -= sval * r[1], _MyR(7) += sval * r[0];
	}

	/**
	 *\brief rotation X with [x, y, z]^T = RX
	 *\param [_X] input 3d point
	 */
	MATRICE_HOST_INL point_type rotation(const point_type& _X)noexcept{
		return (_MyR.mul(_X) + _MyT);
	}
	/**
	 *\brief transform and normalize with [x, y, 1]^T = <RX + T>
	 *\param [_X] input 3d point
	 */
	MATRICE_HOST_INL point_type forward(const point_type& _X)noexcept{
		point_type p = this->rotation(_X) + _MyT;
		return (p.normalize(p.z));
	}

protected:
	matrix_type _MyR;
	point_type _MyT;
};

template<typename _Ty, typename _An> class _Aligned_projection {
	static_assert("Unsupported coordinate system alignment type.");
};

template<typename _Ty>
class _Aligned_projection<_Ty, cs_alignment_tag::left>
	:public _Projection_base<_Ty> {
	using _Myt = _Aligned_projection;
	using _Mybase = _Projection_base<_Ty>;
public:
	using _Mybase::_Projection_base;
	using typename _Mybase::point_type;
	using typename _Mybase::value_type;

	/**
	 *\brief compute the gradient of projected point w.r.t. depth
	 *\param [pd] the first two elements are coordinates of a normalized, distortion-rectified image point, the last one is depth. 
	 */
	MATRICE_HOST_INL point_type backward(const point_type& pd) noexcept {
		const auto d = pd.z; //depth
		const auto X = _Mybase::rotation({ pd.x, pd.y, 1 });
		const auto s = 1 / sqr(X.z*d + _MyT.z);
		const auto gx = (X.x*_MyT.z - X.z*_MyT.x)*s;
		const auto gy = (X.y*_MyT.z - X.z*_MyT.y)*s;
		return point_type{ gx, gy, 0 };
	}
	MATRICE_HOST_INL point_type backward(value_type x, value_type y, value_type depth) noexcept {
		const auto X = _Mybase::rotation({ x, y, 1 });
		const auto s = 1 / sqr(X.z*depth + _MyT.z);
		const auto gx = (X.x*_MyT.z - X.z*_MyT.x)*s;
		const auto gy = (X.y*_MyT.z - X.z*_MyT.y)*s;
		return point_type{ gx, gy, 0 };
	}
private:
	using _Mybase::_MyT;
};
_DETAIL_END

DGE_MATRICE_END