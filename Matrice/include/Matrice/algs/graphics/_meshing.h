/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include <memory>
#include "../../core/matrix.h"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Pty = Vec3_<float>>
class _Meshgrid_impl 
	: std::enable_shared_from_this<_Meshgrid_impl<_Pty>> {
public:
	using point_type = _Pty;
	using value_type = typename point_type::value_t;
	using ptlist_type = std::vector<point_type>;

	_Meshgrid_impl(const ptlist_type& _Pts)
		:_Mypts(std::make_shared(_Pts)) {}

	const decltype(auto) mesh_size() const {
		return (_Mymesz);
	}
	decltype(auto) mesh_size() {
		return (_Mymesz);
	}

private:
	MATRICE_HOST_INL void _Find_bound() {
		const auto _Minx = std::min(_Mypts.begin(), _Mypts, end(), 
			[](const auto& _Left, const auto& _Right) {
			return (_Left.x < _Right.x);
		}).x;
		const auto _Miny = std::min(_Mypts.begin(), _Mypts, end(),
			[](const auto& _Left, const auto& _Right) {
			return (_Left.y < _Right.y);
		}).y;
		const auto _Maxx = std::max(_Mypts.begin(), _Mypts, end(),
			[](const auto& _Left, const auto& _Right) {
			return (_Left.x > _Right.x);
		}).x;
		const auto _Maxy = std::max(_Mypts.begin(), _Mypts, end(),
			[](const auto& _Left, const auto& _Right) {
			return (_Left.y > _Right.y);
		}).y;
		_Mybound = { _Minx, _Miny, _Maxx, _Maxy };
	}
	std::shared_ptr<ptlist_type> _Mypts;
	Vec2_<int32_t> _Mymesz;
	Matrix_<value_type, 2, 2> _Mybound;
};
_DETAIL_END


DGE_MATRICE_END