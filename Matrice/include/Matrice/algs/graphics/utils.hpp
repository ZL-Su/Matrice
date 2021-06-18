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
template<typename _Ty> class linspace;

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

	decltype(auto) mesh_size() const {
		return (_Mymesz);
	}
	decltype(auto) mesh_size() {
		return (_Mymesz);
	}

private:
	MATRICE_HOST_INL void _Find_bound() {
		const auto _Minx = std::min(_Mypts.begin(), _Mypts.end(), 
			[](const auto& _Left, const auto& _Right) {
			return (_Left.x < _Right.x);
		}).x;
		const auto _Miny = std::min(_Mypts.begin(), _Mypts.end(),
			[](const auto& _Left, const auto& _Right) {
			return (_Left.y < _Right.y);
		}).y;
		const auto _Maxx = std::max(_Mypts.begin(), _Mypts.end(),
			[](const auto& _Left, const auto& _Right) {
			return (_Left.x > _Right.x);
		}).x;
		const auto _Maxy = std::max(_Mypts.begin(), _Mypts.end(),
			[](const auto& _Left, const auto& _Right) {
			return (_Left.y > _Right.y);
		}).y;
		_Mybound = { _Minx, _Miny, _Maxx, _Maxy };
	}

	std::shared_ptr<ptlist_type> _Mypts;
	Matrix_<value_type, 2, 2> _Mybound;
	Vec2_<int32_t> _Mymesz;
};
_DETAIL_END

/**
 *\brief CLASS TEMPLATE linspace<_Ty> </class>
 *\typename <_Ty> a scalar type
*/
template<typename _Ty = float> 
class linspace {
	static_assert(is_scalar_v<_Ty>, 
		"_Ty in linspace<> must be a scalar.");
	using _Myt = linspace<_Ty>;
public:
	using value_type = _Ty;
	using row_vector_type = Matrix_<value_type, 1, ::dynamic>;

	/// <summary>
	/// \brief Ctor of class linspace.
	/// </summary>
	/// <param name="'begin'"> The staring value of the sequence </param>
	/// <param name="'end'"> The ending value of the sequence </param>
	/// <param name="'num'"> Number of samples to generate. Default is 51.</param>
	linspace(value_type begin, value_type end, size_t num = 51) noexcept
		: m_begin(begin), m_end(end), m_num(num) {
		m_step = num == 1 ? sub(m_end, m_begin) : safe_div(sub(m_end, m_begin), num - 1);
	}

	MATRICE_GLOBAL_INL decltype(auto) begin() const noexcept {
		return (m_begin);
	}
	MATRICE_GLOBAL_INL decltype(auto) begin() noexcept {
		return (m_begin);
	}
	MATRICE_GLOBAL_INL decltype(auto) end() const noexcept {
		return (m_end);
	}
	MATRICE_GLOBAL_INL decltype(auto) end() noexcept {
		return (m_end);
	}

	MATRICE_GLOBAL_INL decltype(auto) operator++() {
		return (m_current += m_step);
	}
	MATRICE_GLOBAL_INL decltype(auto) operator++(int) {
		return (m_current += m_step);
	}
	MATRICE_GLOBAL_INL operator value_type() const noexcept {
		return m_current;
	}
	MATRICE_GLOBAL_INL operator bool() const noexcept {
		return (m_current <= m_end);
	}

	/// <summary>
	/// \brief Return evenly spaced numbers over the specified interval [begin, end].
	/// If paramerter 'endpoint' is false, then the end point of the interval is excluded.
	/// </summary>
	MATRICE_GLOBAL_INL row_vector_type operator()(bool endpoint = true) const {
		row_vector_type _Ret(endpoint ? m_num : m_num - 1);
		for (auto i = 0; i < _Ret.size(); ++i)
			_Ret(i) = i * m_step + m_begin;
		return _Ret;
	}

	/// <summary>
	/// \brief Return evenly spaced numbers over the specified interval [start, stop].
	/// </summary>
	/// <param name="'start'"> The staring value of the sequence. </param>
	/// <param name="'stop'"> The ending value of the sequence. </param>
	/// <param name="'num'"> Number of samples to generate. Default is 51. </param>
	/// <param name="'endpoint'"> If True, stop is the last sample. Otherwise, it is not included. Default is True. </param>
	/// <returns> A row vector holds the generated number sequence. </returns>
	static row_vector_type _(_Ty start, _Ty stop, size_t num = 51, bool endpoint = true) {
		_Myt linspace(start, stop, num);
		return linspace(endpoint);
	}

private:
	value_type m_begin = 0, m_end = 1;
	value_type m_step = 1, m_current = m_begin;
	size_t m_num;
};

/// <summary>
/// \brief Factory function to create evenly spaced numbers over the specified interval [begin, end].
/// </summary>
/// <typeparam name="_Ty"> Any scalar type. The type of generated numbers is the same as '_Ty'. </typeparam>
/// <param name="'start'"> The staring value of the sequence. </param>
/// <param name="'stop'"> The ending value of the sequence. </param>
/// <param name="'num'"> Number of samples to generate. Default is 51.</param>
/// <param name="'endpoint'">: Conditional argument. If endpoint is false, then the end point of the interval is excluded.</param>
/// <returns>dgelom::Matrix_ with size of 1-by-::dynamic</returns>
template<typename _Ty> requires is_scalar_v<_Ty>
MATRICE_GLOBAL_INL auto make_linspace(_Ty start, _Ty stop, size_t num = 51, bool endpoint = true) noexcept {
	return linspace<conditional_t<std::is_unsigned_v<_Ty>, int, _Ty>>::
		_(start, stop, num, endpoint);
}

DGE_MATRICE_END

#ifdef _MSC_VER
#pragma warning(pop)
#endif