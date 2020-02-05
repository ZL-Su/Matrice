/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once
namespace dgelom {
namespace detail {
	template<typename _Ty> class _SO2;
	template<typename _Ty> class _SO3;
	template<typename _Ty> class _SE2;
	template<typename _Ty> class _SE3;
	template<typename _Ty> class _se2;
	template<typename _Ty> class _se3;
}

namespace internal {
template<class _Group> struct _Lie_group_prop {};
template<class _Algeb> struct _Lie_algebra_prop {};

template<typename _Ty>
struct _Lie_group_prop<detail::_SO2<_Ty>> {
	static constexpr auto dim = 2;
	static constexpr auto dof = 2;
	using type = detail::_SO2<_Ty>;
};

template<typename _Ty>
struct _Lie_group_prop<detail::_SO3<_Ty>>{
	static constexpr auto dim = 3;
	static constexpr auto dof = 3;
	using type = detail::_SO3<_Ty>;
};

template<typename _Ty>
struct _Lie_group_prop<detail::_SE2<_Ty>> {
	static constexpr auto dim = 2;
	static constexpr auto dof = 2;
	using type = detail::_SE2<_Ty>;
};

template<typename _Ty>
struct _Lie_group_prop<detail::_SE3<_Ty>> {
	static constexpr auto dim = 3;
	static constexpr auto dof = 3;
	using type = detail::_SE2<_Ty>;
};

}
}