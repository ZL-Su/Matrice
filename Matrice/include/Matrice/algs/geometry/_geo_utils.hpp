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

#ifdef _MSC_VER
#pragma warning(push)
#endif

DGE_MATRICE_BEGIN namespace geo 
{

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
		m_step(safe_div(sub(m_end,m_begin), m_span)) {
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