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
	template<size_t _Dim, typename _Ty> class _SO;
	template<size_t _Dim, typename _Ty> class _SE;
}

namespace internal {
template<class _Group> struct _Lie_group_prop {};

template<typename _Ty>
struct _Lie_group_prop<detail::_SO<2, _Ty>> {
	static constexpr auto dim = 2;
	static constexpr auto dof = 2;
};

template<typename _Ty>
struct _Lie_group_prop<detail::_SO<3, _Ty>>{
	static constexpr auto dim = 3;
	static constexpr auto dof = 3;
};

template<typename _Ty>
struct _Lie_group_prop<detail::_SE<2, _Ty>> {
	static constexpr auto dim = 2;
	static constexpr auto dof = 2;
};

template<typename _Ty>
struct _Lie_group_prop<detail::_SE<3, _Ty>> {
	static constexpr auto dim = 3;
	static constexpr auto dof = 3;
};

}
}