
#pragma once

#include "../../core/matrix.h"
#include "../../core/vector.h"

DGE_MATRICE_BEGIN _DETAIL_BEGIN
template<typename _Ty, typename = std::enable_if_t<std::is_arithmetic_v<_Ty>>>
MATRICE_HOST_INL Vec3_<_Ty> _Rodrigues_impl(const Matrix_<_Ty, 3, 3>& _R) {
	using value_t = _Ty;
	constexpr const auto _Unit = value_t(1.0);
	constexpr const auto _Half = value_t(0.5);
	constexpr const auto _Zero = zero<value_t>::value;

	Vec3_<value_t> _Ret{_R[2][1]-_R[1][2], _R[0][2]-_R[2][0], _R[1][0]-_R[0][1]};

	auto s = _Ret.norm<2>() * _Half;
	auto c = (_R.trace() - _Unit)*_Half;
	c = c > _Unit ? _Unit : c < -_Unit ? -_Unit : c;
	auto a = std::acos(c);

	if (s < 1e-5) {
		if (c > _Zero) return (_Ret = { _Zero });

		auto t = (_R[0][0] + _Unit)*_Half;
		_Ret.x = sqrt(max(t, _Zero));
		t = (_R[1][1] + _Unit)*_Half;
		_Ret.y = sqrt(max(t, _Zero))*(_R[0][1] < _Zero ? -_Unit : _Unit);
		t = (_R[2][2] + _Unit)*_Half;
		_Ret.z = sqrt(max(t, _Zero))*(_R[0][2] < _Zero ? -_Unit : _Unit);

		if (abs(_Ret.x)<abs(_Ret.y)&&abs(_Ret.x)<abs(_Ret.z)) {
			if ((_R[1][2] > _Zero) != (_Ret.y*_Ret.z > _Zero)) {
				_Ret.z = -_Ret.z;
			}
		}
		a /= _Ret.norm<2>();
		return (_Ret = { a*_Ret.x, a*_Ret.y, a*_Ret.z });
	}
	
	a = _Half * a / s;
	return (_Ret = { a*_Ret.x, a*_Ret.y, a*_Ret.z });
}
_DETAIL_END

template<typename _Ty, size_t _Cols, typename _Outty>
MATRICE_HOST_INL auto rodrigues(const Matrix_<_Ty, 3, _Cols>& _Left, _Outty _Right) {
	auto _Ret = detail::_Rodrigues_impl(_Left);
	if constexpr (std::is_pointer_v<_Outty>)
		dgelom::transform(_Ret.begin(), _Ret.end(), _Right);
	else _Right = _Ret;
}
DGE_MATRICE_END
