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
#include <type_traits>
#include "../util/_macros.h"
#include "../util/_type_defs.h"
#include "../util/_std_wrapper.h"

DGE_MATRICE_BEGIN

template<typename T> struct remove_reference { using type = typename std::remove_reference<T>::type; };
template<typename T> using remove_reference_t = typename remove_reference<T>::type;

/**
 *\brief retrieve plain type of type T.
 */
template<typename T> struct remove_all { using type = std::remove_cv_t<remove_reference_t<T>>; };
template<typename T> using remove_all_t = typename remove_all<T>::type;

/**
 *\brief get bytes of type T.
 */
template<typename T> struct type_bytes { enum { value = sizeof(T) }; };
template<typename T> MATRICE_GLOBAL_INL constexpr int type_bytes_v = type_bytes<T>::value;

/**
 *\brief conditional_t<_Cond, T, U> is T if _Cond is true, else it is U.  
 */
template<bool _Cond, typename T, typename U> struct conditional {};
template<typename T, typename U> struct conditional<true, T, U> { using type = T; };
template<typename T, typename U> struct conditional<false, T, U> { using type = U; };
template<bool _Cond, typename T, typename U> using conditional_t = typename conditional<_Cond, T, U>::type;

/**
 *\brief has_value_t<T> is true iff T has value_t.
 */
template<typename T, typename Enable = void> struct has_value_t : std::false_type {};
template<typename T> MATRICE_GLOBAL_INL constexpr auto has_value_t_v = has_value_t<T>::value;

template<typename T, typename = enable_if_t<has_value_t_v<T>>> struct value_type { using type = conditional_t<is_arithmetic_v<T>, remove_all_t<T>, typename T::value_t>; };
template<typename T> using value_type_t = typename value_type<T>::type;

template<typename T, typename U> struct common_value_type { using type = common_type_t<value_type_t<T>, value_type_t<U>>; };
template<typename T, typename U> using common_value_t = typename common_value_type<T, U>::type;

/**
 *\brief is_zero_v<_Val> is true iff _Val == 0.
 */
template<int _Val> struct is_zero { enum { value = _Val == 0 }; };
template<int _Val> MATRICE_GLOBAL_INL constexpr auto is_zero_v = is_zero<_Val>::value;

template<typename T> struct is_scalar {
	constexpr static auto value = std::is_scalar_v<T>;
};
template<typename T> inline constexpr auto is_scalar_v = is_scalar<T>::value;

/**
 *\brief is_static_v<_R, _C> is true iff both _R and _C is greater than 0.
 */
template<int _R, int _C> struct is_static {enum {value = _R > 0 && _C >0 }; };
template<int _R, int _C> MATRICE_GLOBAL_INL constexpr auto is_static_v = is_static<_R, _C>::value;

template<typename T> struct is_common_int64 { enum { value = is_integral_v<T> && sizeof(T) == 8 }; };
template<typename T> MATRICE_GLOBAL_INL constexpr auto is_common_int64_v = is_common_int64<T>::value;

template<typename T> struct is_int64 { enum { value = std::is_signed_v<T> && is_common_int64_v<T> }; };
template<typename T> MATRICE_GLOBAL_INL constexpr auto is_int64_v = is_int64<T>::value;

template<typename T> struct is_uint64 { enum { value = std::is_unsigned_v<T> && is_common_int64_v<T> }; };
template<typename T> MATRICE_GLOBAL_INL constexpr auto is_uint64_v = is_uint64<T>::value;

template<typename T> struct is_float32 { enum { value = std::is_floating_point<T>::value && (sizeof(T) == 4) }; };
template<typename T> MATRICE_GLOBAL_INL constexpr auto is_float32_v = is_float32<T>::value;

template<typename T> struct is_float64 { enum { value = std::is_floating_point<T>::value && (sizeof(T) == 8) }; };
template<typename T> MATRICE_GLOBAL_INL constexpr auto is_float64_v = is_float64<T>::value;

template<typename T> MATRICE_GLOBAL_INL constexpr auto is_floating_point_v = is_float64_v<T> || is_float32_v<T>;

template<typename T> struct add_const_reference {
	using type = std::add_lvalue_reference_t<std::add_const_t<T>>;
};
template<typename T> using add_const_reference_t = typename add_const_reference<T>::type;

template<typename T> struct add_const_pointer {
	using type = std::add_pointer_t<std::add_const_t<T>>;
};
template<typename T> using add_const_pointer_t = typename add_const_pointer<T>::type;

template<typename T> struct add_pointer_const {
	using type = std::add_const_t<std::add_pointer_t<T>>;
};
template<typename T> using add_pointer_const_t = typename add_pointer_const<T>::type;

/**
 *\brief is_iterator_v<T> is true iff T is an iterator or pointer.
 */
template<typename T> constexpr bool is_iterator_v = std::_Is_iterator_v<T>;
/**
 *\brief retrieve iterator traits, including categary, value, difference, pointer, and reference types, etc.
 */
template<typename T> using iterater_traits = std::iterator_traits<T>;

template<typename T> struct _View_trait { enum { value = 0x0008*sizeof(T) }; };
template<> struct _View_trait<unsigned char> { enum { value = 0x0008 }; };
template<> struct _View_trait<int> { enum { value = 0x0016 }; };
template<> struct _View_trait<size_t> { enum { value = 0x0016 }; };
template<> struct _View_trait<float> { enum { value = 0x0032 }; };
template<> struct _View_trait<double> { enum { value = 0x0064 }; };
template<typename T> MATRICE_GLOBAL_INL constexpr auto plane_view_v = _View_trait<T>::value;

template<typename T>
struct traits { 
	using type =
		conditional_t<is_class_v<T>, typename T::value_type, remove_all_t<T>>;
};

template<typename T> struct is_matrix : std::false_type {};
template<typename T> MATRICE_GLOBAL_INL constexpr 
bool is_matrix_v = is_matrix<T>::value;
template<typename Mty>
struct matrix_traits : traits<Mty> {};

template<typename Exp> struct expression_options { enum { value = Exp::flag | expr }; };
template<typename Exp> struct is_expression :std::false_type {};
template<typename Exp> MATRICE_GLOBAL_INL constexpr 
bool is_expression_v = is_expression<Exp>::value;
template<class Exp, typename = enable_if_t<is_expression_v<Exp>>>
struct expression_traits : traits<Exp> {};

template<typename T> struct is_mtxview : std::false_type {};
/**
 *\brief is_mtxview_v<T> is true iff T is dgelom::detail::_Matrix_view or its derived type.
 */
template<typename T> MATRICE_GLOBAL_INL constexpr bool is_mtxview_v = is_mtxview<T>::value;
template<typename View> struct mtxview_traits : traits<View> {};

template<typename T> struct is_fxdvector : std::false_type {};
/**
 *\brief is_fxdvector_v<T> is true iff T is dgelom::Vec_ or its derived type.
 */
template<typename T> MATRICE_GLOBAL_INL constexpr bool is_fxdvector_v = is_fxdvector<T>::value;

/**
 *\brief is_matrix_convertible_v<T> is true type iff T is dgelom::Matrix_<...> or related matrix expression or view type.  
 */
template<typename T> MATRICE_GLOBAL_INL constexpr bool is_matrix_convertible_v = is_matrix_v<T> || is_expression_v<T> || is_mtxview_v<T>;

template<int _M, int _N> struct allocator_traits;
template<int _M, int _N> MATRICE_GLOBAL_INL constexpr auto allocator_traits_v = allocator_traits<_M, _N>::value;

/**
 *\brief internal type for accessing allocator traits.
 */
template<typename Al> struct _Allocator_traits; /*{
	static_assert(false, "Unknown allocator Al in _Allocator_traits.");
};*/

template<class T, typename = enable_if_t<is_matrix_v<T> || is_expression_v<T>>> 
struct layout_traits : traits<T> {
	MATRICE_GLOBAL_FINL static auto layout_type(size_t _format) {
		return (_format & rmaj == rmaj) ? rmaj : cmaj;
	}
	MATRICE_GLOBAL_FINL static auto storage_type(size_t _format) {
		return (_format & symm == symm) ? symm : (_format & diag == diag) ? diag :
			    (_format & band == band) ? band : (_format & utri == utri) ? utri :
				 (_format & ltri == ltri) ? ltri : (_format & spar == spar) ? spar :
				 gene;
	}
	MATRICE_GLOBAL_FINL static bool is_rmajor(size_t _format) {
		return _format & rmaj == rmaj;
	}
	MATRICE_GLOBAL_FINL static bool is_cmajor(size_t _format) {
		return _format & cmaj == cmaj;
	}
	MATRICE_GLOBAL_FINL static bool is_symmetric(size_t _format) {
		return _format & symm == symm;
	}
	MATRICE_GLOBAL_FINL static bool is_diagnal(size_t _format) {
		return _format & diag == diag;
	}
	MATRICE_GLOBAL_FINL static bool is_banded(size_t _format) {
		return _format & band == band;
	}
	MATRICE_GLOBAL_FINL static bool is_sparse(size_t _format) {
		return _format & spar == spar;
	}
	MATRICE_GLOBAL_FINL static bool is_uppertr(size_t _format) {
		return _format & utri == utri;
	}
	MATRICE_GLOBAL_FINL static bool is_lowertr(size_t _format) {
		return _format & ltri == ltri;
	}
};

/**
 *\brief is_equality_comparable<T> is true_type iff the equality operator is defined for T.
 */
template<typename T, typename Enable = void> struct is_equality_comparable : std::false_type {};
template<typename T> struct is_equality_comparable<T, std::void_t<decltype(std::declval<T&>() == std::declval<T&>())>> : std::true_type {};
template<typename T> MATRICE_GLOBAL_INL constexpr
auto is_equality_comparable_v = is_equality_comparable<T>::value;

/**
 *\brief is_hashable<T> is true_type iff std::hash is defined for T
 */
template<typename T, typename Enable = void> struct is_hashable : std::false_type {};
template<typename T> struct is_hashable<T, std::void_t<decltype(std::hash<T>()(std::declval<T&>()))>> : std::true_type {};
template<typename T> MATRICE_GLOBAL_INL constexpr 
auto is_hashable_v = is_hashable<T>::value;

/**
 *\brief is_instantiation_of_v<T, I> is true_type iff I is a template instantiation of T (e.g. vector<int> is an instantiation of vector)
 *  Example:
 *    is_instantiation_of_v<vector, vector<int>> // true
 *    is_instantiation_of_v<pair, pair<int, string>> // true
 *    is_instantiation_of_v<vector, pair<int, string>> // false
 */
template<template<class...> class Template, class T, int... sizes>
struct is_instantiation_of : std::false_type {};
template<template<class...> class Template, class... Args>
struct is_instantiation_of<Template, Template<Args...>> : std::true_type {};
template<template<class...> class Template, class T> 
MATRICE_GLOBAL_INL constexpr auto is_instantiation_of_v = is_instantiation_of<Template, T>::value;

/**
 *\brief is_type_condition<C> is true_type iff C<...> is a type trait representing a condition (i.e. has a constexpr static bool ::value member)
 * Example:
 *   is_type_condition<std::is_reference>  // true
 */
template<template<typename> typename C, typename Enable = void>
struct is_type_condition : std::false_type {};
template<template<typename> typename C>
struct is_type_condition<C, std::enable_if_t<std::is_same<bool, std::remove_cv_t<decltype(C<int>::value)>>::value>> : std::true_type {};
template<template<typename> typename C> 
MATRICE_GLOBAL_INL constexpr auto is_type_condition_v = is_type_condition<C>::value;

/**
 *\brief retrieve solver traits
 */
template<typename _Ty> struct solver_traits {};

/**
 *\brief has_method_Fn<T, _Args> is true_type iff T has method T.Fn(_Args...)
 */
#define _HAS_METHOD(_Name) \
template<typename T, typename... _Args> struct has_method_##_Name { \
	static constexpr bool value = std::is_same_v<decltype(_Check<T>(0)), std::true_type::value>; \
private: \
	template<typename _C> static auto _Check(int)->decltype(std::declval<_C>()._Name(std::declval<_Args>()...), std::true_type()); \
	template<typename _C> static std::false_type _Check(...); \
};

 /**
  *\brief has_value_t<T> is true_type iff T has member value_t
  */
template<typename T>
struct has_value_t<T> { static constexpr auto value = is_expression_v<T> || is_matrix_v<T> || is_mtxview_v<T>; };

/**
 *\brief is_tensor<T> is true_type iff T is dgelom::tensor<...>
 */
template<typename T> struct is_tensor : std::false_type {};
template<typename T>
MATRICE_GLOBAL_INL constexpr auto is_tensor_v = is_tensor<T>::value;

/**
 *\brief tensor_traits<T> for accessing tensor members
 */
template<typename T> struct tensor_traits {};

/**
 *\brief category_type<T> for accessing category member of T
 */
template<typename T> struct category_type { using type = typename T::category; };
template<typename T> using category_type_t = typename category_type<T>::type;

/**
 *\brief is_same_v<T,U> is true iff. T and U are same type.
 */
template<typename T, typename U> 
MATRICE_GLOBAL_INL constexpr auto is_same_v = std::is_same<T, U>::value;

/**
 *\brief is_same_v<T,U> is true iff. T and U are not same type.
 */
template<typename T, typename U>
MATRICE_GLOBAL_INL constexpr auto is_not_same_v = !std::is_same<T, U>::value;

/**
 *\brief is_any_of_v<T, Ts...> is true iff. T is one of types encapsulated in Ts....
 */
template<typename T, typename... Ts>
MATRICE_GLOBAL_INL constexpr auto is_any_of_v = std::disjunction_v<std::is_same<T, Ts>...>;

DGE_MATRICE_END