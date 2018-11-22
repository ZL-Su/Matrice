/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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
#include "_splineinterp.h"

MATRICE_ALGS_BEGIN

template<typename _Ty, std::size_t _Opt>
class _Interpolation_wrapper MATRICE_NONHERITABLE
{
public:
	using type = auto_interp_dispatcher_t<_Ty, _Opt>;
	using value_type = typename type::value_type;
	static constexpr auto option = type::option;

	template<typename... _Args>
	MATRICE_GLOBAL_FINL _Interpolation_wrapper(const _Args&... args)
		:m_op(std::make_shared<type>(args...)) {}

	/**
	 * \Get reference of interpolation kernel operator.
	 * \Example:
	 *		_Interpolation_wrapper<...> _Myname(...);
	 *		const auto& _Itp = _Myname(); \\return interpolation kernel.
	 *    const auto& _Coeff = _Itp(); \\return interpolation coeff.
	 *    auto _Value = _Itp(_Pos); \\return interpolated value at _Pos.
	 *    auto[_Gx, _Gy] = _Itp.grad(_Pos); \\return iterpolated grad. at _Pos.
	 *		_Gx = _Itp.grad<axis::x>(_Pos); \\return interpolated dI/dx.
	 *		_Gy = _Itp.grad<axis::y>(_Pos); \\return interpolated dI/dy.
	 */
	MATRICE_GLOBAL_FINL auto& operator()() {
		return (m_op);
	}
	MATRICE_GLOBAL_FINL const auto& operator()() const {
		return (m_op);
	}
	MATRICE_GLOBAL_INL operator type() const {
		return (*m_op);
	}

private:
	mutable std::shared_ptr<type> m_op;
};

MATRICE_ALGS_END