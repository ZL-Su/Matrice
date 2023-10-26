#pragma once

#include "core/matrix.h"
#include "core/vector.h"

MATRICE_ALGS_BEGIN
template<typename _Ty>
class PlaneFitting {
public:
	using value_t = _Ty;
	using matrix_t = Matrix<value_t>;
	using point_t = Vector3<value_t>;

	PlaneFitting(const matrix_t& data) noexcept 
		: _Mydata(data) {
	}

private:
	auto find_plane(const point_t& x, const point_t& y, const point_t& z);

	const matrix_t& _Mydata;
};

template<typename _Ty> MATRICE_HOST_FINL
auto PlaneFitting<_Ty>::find_plane(const point_t& x, const point_t& y, const point_t& z)
{
	const point_t yx = y - x, yz = y - z;
	auto _Norm = yx.cross(yz);
	_Norm = _Norm.normalize(_Norm.norm());
}

MATRICE_ALGS_END
