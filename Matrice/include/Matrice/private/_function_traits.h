/*  *************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
*	*************************************************************************/
#pragma once
#include <functional>
#include "_type_traits.h"
#include "_size_traits.h"
#include "_typelist.h"

DGE_MATRICE_BEGIN
/**
 * is_function_type<T> is true_type iff T is a plain function type (i.e. "Result(Args...)")
 */
template<typename _Ty>
struct is_function_type : std::false_type {};
template<typename _Ret, typename... _Args>
struct is_function_type<_Ret(_Args...)> : std::true_type {};
template<typename _Ty> using is_function_type_t = typename is_function_type<_Ty>::type;

/**
 * Access information about result type or arguments from a function type.
 * Example:
 * using A = function_traits<int (float, double)>::return_type // A == int
 * using A = function_traits<int (float, double)>::parameter_types::tuple_type // A == tuple<float, double>
 */
template<typename _Fty, MATRICE_ENABLE_IF(is_function_v<_Fty>)> struct function_traits {};
template<typename _Ret, typename... _Args> 
struct function_traits<_Ret(_Args...)> {
	using plain_type = _Ret(_Args...);
	using function_type = std::function<_Ret(_Args...)>;
	using return_type = _Ret;
	using argument_type_list = typelist<_Args...>;
	static constexpr auto nargts = argument_type_list::ntypes;
};
DGE_MATRICE_END
