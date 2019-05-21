/*  *************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-, Zhilong(Dgelom) Su, all rights reserved.

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

#include "_type_traits.h"

DGE_MATRICE_BEGIN namespace tl {
namespace detail { 
	template<typename... _Args> struct false_t : std::false_type {};
	template<typename... _Args> MATRICE_GLOBAL_INL constexpr auto false_t_v = false_t<_Args...>::value;
}
/**
 * Forward declarations
 */
template<typename _TyList> struct to_tuple;
template<typename _TyList> struct size;
template<size_t _Idx, typename _TyList> struct element;
/**
 * Type holding a list of types for compile time type computations
 */
template<typename... _Args> 
struct typelist MATRICE_NONHERITABLE {

	static constexpr size_t ntypes = size<typelist>::value;
	using tuple_type = typename to_tuple<typelist>::type;
	template<size_t _Idx>
	using type = typename element<_Idx, typelist>::type;

private:
	typelist() = delete; // not for instantiation
};

/**
 * Returns the number of types in a typelist
 * Example:
 *   3  ==  size<typelist<int, int, double>>::value
 */
template<typename _TyList> 
struct size MATRICE_NONHERITABLE {
	static_assert(detail::false_t_v<_TyList>, "In tl::size<T>, T must be typelist<...>.");
};
template<typename... _Tys> 
struct size<typelist<_Tys...>> MATRICE_NONHERITABLE {
	static constexpr size_t value = sizeof...(_Tys);
};
template<typename... _Tys> 
MATRICE_GLOBAL_INL constexpr auto size_v = size<_Tys>::value;

/**
 * Transforms a list of types into a tuple holding these types.
 * Example:
 *   std::tuple<int, string>  ==  to_tuple_t<typelist<int, string>>
 */
template<typename _TyList> 
struct to_tuple MATRICE_NONHERITABLE {
	static_assert(detail::false_t_v<_TyList>, "In tl::to_tuple<T>, T must be typelist<...>.");
};
template<typename... _Tys> 
struct to_tuple<typelist<_Tys...>> MATRICE_NONHERITABLE {
	using type = tuple<_Tys...>;
};
template<typename _TyList> using to_tuple_t = typename to_tuple<_TyList>::type;


/**
 * Creates a typelist containing the types of a given tuple.
 * Example:
 *   typelist<int, string>  ==  from_tuple_t<std::tuple<int, string>>
 */
template<typename _Tuple> 
struct from_tuple MATRICE_NONHERITABLE {
	static_assert(detail::false_t_v<_Tuple>, "In tl::from_tuple<T>, T must be std::tuple<...>.");
};
template<typename... _Tys> 
struct from_tuple<std::tuple<_Tys...>> MATRICE_NONHERITABLE {
	using type = typelist<_Tys...>;
};
template<class Tuple> using from_tuple_t = typename from_tuple<Tuple>::type;


/**
 * Concatenates multiple type lists.
 * Example:
 *   typelist<int, string, int>  ==  concat_t<typelist<int, string>, typelist<int>>
 */
template<typename... _TyLists> 
struct concat MATRICE_NONHERITABLE {
	static_assert(detail::false_t_v<_TyLists...>, "In tl::concat<T1, ...>, the T arguments each must be typelist<...>.");
};
template<typename... _Head1Tys, typename... _Head2Tys, typename... _TailLists>
struct concat<typelist<_Head1Tys...>, typelist<_Head2Tys...>, _TailLists...> MATRICE_NONHERITABLE {
	using type = typename concat<typelist<_Head1Tys..., _Head2Tys...>, _TailLists...>::type;
};
template<typename... _HeadTys>
struct concat<typelist<_HeadTys...>> MATRICE_NONHERITABLE { using type = typelist<_HeadTys...>; };
template<> struct concat<> MATRICE_NONHERITABLE { using type = typelist<>; };
template<typename... _TyLists> using concat_t = typename concat<_TyLists...>::type;

/**
 * Filters the types in a type list by a type trait.
 * Examples:
 *   typelist<int&, const string&&>  ==  filter_t<std::is_reference, typelist<void, string, int&, bool, const string&&, int>>
 */
template<template <typename> class _Cond, typename _TyList> 
struct filter MATRICE_NONHERITABLE {
	static_assert(detail::false_t_v<_TyList>, "In tl::filter<_Cond, _TyList>, the _TyList argument must be typelist<...>.");
};
template<template <typename> class _Cond, class _Head, class... _Tails>
struct filter<_Cond, typelist<_Head, _Tails...>> MATRICE_NONHERITABLE {
	static_assert(is_type_condition_v<_Cond>, "In tl::filter<_Cond, _TyList>, the _Cond argument must be a condition type trait, i.e. have a static constexpr bool ::value member.");
	using type = conditional_t<_Cond<_Head>::value,
		concat_t<typelist<_Head>, typename filter<_Cond, typelist<_Tails...>>::type>,
		typename filter<_Cond, typelist<_Tails...>>::type>;
};
template<template <typename> class _Cond>
struct filter<_Cond, typelist<>> MATRICE_NONHERITABLE {
	static_assert(is_type_condition_v<_Cond>, "In tl::filter<_Cond, _TyList>, the _Cond argument must be a condition type trait, i.e. have a static constexpr bool ::value member.");
	using type = typelist<>;
};
template<template <typename> class _Cond, typename _TyList>
using filter_t = typename filter<_Cond, _TyList>::type;

/**
 * Counts how many types in the list fulfill a type trait
 * Examples:
 *   2  ==  count_if<std::is_reference, typelist<void, string, int&, bool, const string&&, int>>
 */
template<template <typename> class _Cond, typename _TyList>
struct count_if final {
	static_assert(is_type_condition_v<_Cond>, "In tl::count_if<_Cond, _TyList>, the _Cond argument must be a condition type trait, i.e. have a static constexpr bool ::value member.");
	static_assert(is_instantiation_of<typelist, _TyList>::value, "In tl::count_if<_Cond, _TyList>, the _TyList argument must be typelist<...>.");
	// TODO Direct implementation might be faster
	static constexpr size_t value = size_v<filter_t<_Cond, _TyList>>;
};
template<template <typename> class _Cond, typename _TyList> MATRICE_GLOBAL_INL constexpr size_t count_if_v = count_if<_Cond, _TyList>::value;

/**
 * Returns true iff the type trait is true for all types in the type list
 * Examples:
 *   true   ==  true_for_each_type<std::is_reference, typelist<int&, const float&&, const MyClass&>>::value
 *   false  ==  true_for_each_type<std::is_reference, typelist<int&, const float&&, MyClass>>::value
 */
template<template <typename> class _Cond, typename _TyList> 
struct true_for_each_type MATRICE_NONHERITABLE {
	static_assert(detail::false_t_v<_TyList>, "In tl::true_for_each_type<_Cond, _TyList>, the _TyList argument must be typelist<...>.");
};
template<template <typename> class _Cond, class... _Tys>
struct true_for_each_type<_Cond, typelist<_Tys...>> MATRICE_NONHERITABLE 
	: std::conjunction<_Cond<_Tys>...> {
	static_assert(is_type_condition_v<_Cond>, "In tl::true_for_each_type<_Cond, _TyList>, the _Cond argument must be a condition type trait, i.e. have a static constexpr bool ::value member.");
};

/**
 * Maps types of a type list using a type trait
 * Example:
 *  typelist<int&, double&, string&>  ==  map_t<std::add_lvalue_reference_t, typelist<int, double, string>>
 */
template<template <typename> class _Mapper, typename _TyList> struct map final {
	static_assert(detail::false_t_v<_TyList>, "In tl::map<_Mapper, _TyList>, the _TyList argument must be typelist<...>.");
};
template<template <typename> class _Mapper, typename... _Tys>
struct map<_Mapper, typelist<_Tys...>> final { using type = typelist<_Mapper<_Tys>...>; };
template<template <class> class Mapper, typename _TyList>
using map_t = typename map<Mapper, _TyList>::type;

/**
 * Returns the first element of a type list.
 * Example:
 *   int  ==  head_t<typelist<int, string>>
 */
template<typename _TyList> struct head final {
	static_assert(detail::false_t_v<_TyList>, "In tl::head<T>, the T argument must be typelist<...>.");
};
template<typename _Head, typename... _Tails> struct head<typelist<_Head, _Tails...>> final {
	using type = _Head;
};
template<typename _TyList> using head_t = typename head<_TyList>::type;

/**
 * Returns the N-th element of a type list.
 * Example:
 * int == element_t<1, typelist<float, int, char>>
 */

 /// Base template.
template<size_t _Idx, typename _TyList> 
struct element MATRICE_NONHERITABLE {
	static_assert(detail::false_t_v<_TyList>, "In tl::element<T>, the T argument must be typelist<...>.");
};

/// Successful case, we have reached the zero index and can "return" the head type.
template<typename _Head, typename... _Tails> struct element<0, typelist<_Head, _Tails...>> { using type = _Head; };

/// Error case, we have an index but ran out of types! It will only be selected
/// if `Ts...` is actually empty!
template <size_t _Idx, class... _Tys>
struct element<_Idx, typelist<_Tys...>> {
	static_assert(_Idx < sizeof...(_Tys), "Index is out of bounds in tl::element");
};

/// Shave off types until we hit the <0, Head, Tail...> or <Index> case.
template<size_t _Idx, typename _Head, typename... _Tails> 
struct element<_Idx, typelist<_Head, _Tails...>> : element<_Idx - 1, typelist<_Tails...>> { };

/// Convenience alias.
template<size_t _Idx, typename _TyList>
using element_t = typename element<_Idx, _TyList>::type;

/**
 * Returns the last element of a type list.
 * Example:
 *   int  ==  last_t<typelist<int, string>>
 */
template <typename _TyList>
struct last final {
	static_assert(detail::false_t_v<_TyList>, "In tl::last<T>, the T argument must be typelist<...>.");
};
template <typename _Head, typename... _Tails>
struct last<typelist<_Head, _Tails...>> final { using type = typename last<typelist<_Tails...>>::type; };
template <typename _Head>
struct last<typelist<_Head>> final { using type = _Head; };
template <typename _TyList>
using last_t = typename last<_TyList>::type;
static_assert(std::is_same_v<int, last_t<typelist<double, float, int>>>, "");

/**
 * Reverses a typelist.
 * Example:
 *   typelist<int, string>  == reverse_t<typelist<string, int>>
 */
template<typename _TyList> 
struct reverse MATRICE_NONHERITABLE {
	static_assert(detail::false_t_v<_TyList>, "In dgelom::reverse<T>, the T argument must be typelist<...>.");
};
template<typename _Head, typename... _Tails> 
struct reverse<typelist<_Head, _Tails...>> MATRICE_NONHERITABLE {
	using type = concat_t<typename reverse<typelist<_Tails...>>::type, typelist<_Head>>;
};
template<> struct reverse<typelist<>> final { using type = typelist<>; };
template<typename _TyList> using reverse_t = typename reverse<_TyList>::type;

/**
 * Maps a list of types into a list of values.
 * Examples:
 *   // C++14 example
 *   auto sizes =
 *     map_types_to_values<typelist<int64_t, bool, uint32_t>>(
 *       [] (auto t) { return sizeof(decltype(t)::type); }
 *     );
 *   //  sizes  ==  std::tuple<size_t, size_t, size_t>{8, 1, 4}
 *
 *   // C++14 example
 *   auto shared_ptrs =
 *     map_types_to_values<typelist<int, double>>(
 *       [] (auto t) { return make_shared<typename decltype(t)::type>(); }
 *     );
 *   // shared_ptrs == std::tuple<shared_ptr<int>, shared_ptr<double>>()
 *
 *   // C++11 example
 *   struct map_to_size {
 *     template<class T> constexpr size_t operator()(T) {
 *       return sizeof(typename T::type);
 *     }
 *   };
 *   auto sizes =
 *     map_types_to_values<typelist<int64_t, bool, uint32_t>>(
 *       map_to_size()
 *     );
 *   //  sizes  ==  std::tuple<size_t, size_t, size_t>{8, 1, 4}
 */
namespace detail {
	template<class T> struct type_ final { using type = T; };

	template<typename _TyList> struct map_types_to_values final {
		static_assert(detail::false_t_v<_TyList>, "In tl::map_types_to_values<T>, the T argument must be typelist<...>.");
	};
	template<typename... _Tys> struct map_types_to_values<typelist<_Tys...>> final {
		template<typename _Fty>
		static tuple<std::invoke_result_t<_Fty(type_<_Tys>)>...> call(_Fty&& func) {
			return tuple<std::invoke_result_t<_Fty(type_<_Tys>)>...> { forward<_Fty>(func)(type_<_Tys>())... };
		}
	};
}
}

template<typename... _Args> using typelist = tl::typelist<_Args...>;

template<typename _TyList, typename _Fty> 
auto map_types_to_values(_Fty&& func)->
decltype(tl::detail::map_types_to_values<_TyList>::call(forward<_Fty>(func))) {
	return tl::detail::map_types_to_values<_TyList>::call(forward<_Fty>(func));
} 

DGE_MATRICE_END