/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
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
#include <memory>
#include "core/matrix.h"
#include "core/vector.h"

#ifdef _MSC_VER
#pragma warning(push)
#endif

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Pty = Vec3_<float>>
class _Meshgrid_based
	: std::enable_shared_from_this<_Meshgrid_based<_Pty>> {
	using _Myt = _Meshgrid_based<_Pty>;
public:
	using point_type = _Pty;
	using value_type = typename point_type::value_t;
	using ptlist_type = std::vector<point_type>;

	_Meshgrid_based(const ptlist_type& _Pts)
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

namespace graph {
/**
 *\brief <class> linspace<_Ty> </class>
 *\typen <_Ty> a scalar type
*/
template<typename _Ty = float> class linspace {
	static_assert(is_scalar_v<_Ty>, "_Ty in linspace<> must be a scalar.");
public:
	using value_type = _Ty;

	template<typename _Vt = size_t> MATRICE_GLOBAL_INL
		linspace(value_type begin, value_type end, _Vt span) noexcept
		: m_begin(begin), m_end(end), m_span(static_cast<value_type>(span)),
		m_step(safe_div(sub(m_end, m_begin), m_span)) {
	}

	MATRICE_GLOBAL_INL const value_type& begin() const noexcept {
		return (m_begin);
	}
	MATRICE_GLOBAL_INL value_type& begin() noexcept {
		return (m_begin);
	}
	MATRICE_GLOBAL_INL const value_type& end() const noexcept {
		return (m_end);
	}
	MATRICE_GLOBAL_INL value_type& end() noexcept {
		return (m_end);
	}

	MATRICE_GLOBAL_INL value_type& operator++() {
		return (m_current += m_step);
	}
	MATRICE_GLOBAL_INL value_type& operator++(int) {
		return (m_current += m_step);
	}
	MATRICE_GLOBAL_INL operator value_type() const noexcept {
		return m_current;
	}
	MATRICE_GLOBAL_INL operator bool() const noexcept {
		return (m_current <= m_end);
	}
private:
	value_type m_begin = 0, m_end = 1;
	value_type m_span, m_step = 1, m_current = m_begin;
};
}

DGE_MATRICE_END

#ifdef _MSC_VER
#pragma warning(pop)
#endif