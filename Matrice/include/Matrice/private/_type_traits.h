/*  *************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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
#include "../util/_type_defs.h"
#include "_memory.h"

MATRICE_NAMESPACE_BEGIN_

template<typename T> struct remove_reference { using type = typename std::remove_reference<T>::type; };
template<typename T> using remove_reference_t = typename remove_reference<T>::type;

template<typename T> struct type_bytes { enum { value = sizeof(T) }; };
template<typename T> MATRICE_GLOBAL_INL constexpr int type_bytes_v = type_bytes<T>::value;

template<bool _Test, typename T1, typename T2> struct conditonal {};
template<typename T1, typename T2> struct conditonal<true, T1, T2> { using type = T1; };
template<typename T1, typename T2> struct conditonal<false, T1, T2> { using type = T2; };
template<bool _Test, typename T1, typename T2> struct conditional {};
template<typename T1, typename T2> struct conditional<true, T1, T2> { using type = T1; };
template<typename T1, typename T2> struct conditional<false, T1, T2> { using type = T2; };
template<bool _Test, typename T1, typename T2> using conditional_t = typename conditional<_Test, T1, T2>::type;

template<int _Val> struct is_zero { enum { value = _Val == 0 }; };
template<int _R, int _C> struct is_static {enum {value = _R > 0 && _C >0 }; };

template<typename T> struct is_common_int64 { enum { value = std::is_integral_v<T> && sizeof(T) == 8 }; };
template<typename T> struct is_int64 { enum { value = std::is_signed_v<T> && std::is_integral_v<T> && sizeof(T) == 8 }; };
template<typename T> struct is_uint64 { enum { value = std::is_unsigned_v<T> && std::is_integral_v<T> && sizeof(T) == 8 }; };

template<typename T> struct is_float32 { enum { value = std::is_floating_point<T>::value && (sizeof(T) == 4) }; };
template<typename T> struct is_float64 { enum { value = std::is_floating_point<T>::value && (sizeof(T) == 8) }; };

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

template<typename T> struct _View_trait { enum { value = 0x0008*sizeof(T) }; };
template<> struct _View_trait<unsigned char> { enum { value = 0x0008 }; };
template<> struct _View_trait<int> { enum { value = 0x0016 }; };
template<> struct _View_trait<float> { enum { value = 0x0032 }; };
template<> struct _View_trait<double> { enum { value = 0x0064 }; };

template<typename T, typename = std::enable_if_t<std::is_class_v<T>>>
struct traits { using type = typename T::value_t; };

template<typename _Ty> struct is_matrix : std::false_type {};
template<typename _Ty> MATRICE_GLOBAL_INL constexpr bool is_matrix_v = is_matrix<_Ty>::value;
template<typename Mty, typename = std::enable_if_t<is_matrix_v<Mty>>>
struct matrix_traits : traits<Mty> {};

template<typename Exp> struct expression_options { enum { value = Exp::flag | expr }; };
template<typename Exp> struct is_expression :std::false_type {};
template<typename Exp> MATRICE_GLOBAL_INL constexpr bool is_expression_v = is_expression<Exp>::value;
template<class Exp, typename = std::enable_if_t<is_expression_v<Exp>>>
struct expression_traits : traits<Exp> {};

template<typename _Ty> struct is_iterator: std::false_type {};
template<typename _Ty> MATRICE_GLOBAL_INL constexpr bool is_iterator_v = is_iterator<_Ty>::value;
template<typename Itr> struct iterator_traits : traits<Itr> {};

template<typename _Ty> struct is_mtxview : std::false_type {};
template<typename _Ty> MATRICE_GLOBAL_INL constexpr bool is_mtxview_v = is_mtxview<_Ty>::value;
template<typename View> struct mtxview_traits : traits<View> {};

template<int _M, int _N> struct allocator_traits {
	enum {
		value = _M > 0 && _N > 0 ? LINEAR + COPY :  // stack allocator
		_M == 0 && _N == -1 ? LINEAR :  // linear device allocator
		_M == -1 && _N == -1 ? PITCHED :  // pitched device allocator
#ifdef __CXX11_SHARED__
		LINEAR + SHARED  // smart heap or global allocator
#else
		LINEAR + COPY    // deep heap or global allocator
#endif      
	};
};

template<class _Ty, typename = std::enable_if_t<is_matrix_v<_Ty> || is_expression_v<_Ty>>> 
struct layout_traits : traits<_Ty> {
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
 * is_equality_comparable<T> is true_type iff the equality operator is defined for T.
 */
template<typename T, typename Enable = void> struct is_equality_comparable : std::false_type {};
template<typename T> struct is_equality_comparable<T, std::void_t<decltype(std::declval<T&>() == std::declval<T&>())>> : std::true_type {};
template<typename T> using is_equality_comparable_t = typename is_equality_comparable<T>::type;

/**
 * is_hashable<T> is true_type iff std::hash is defined for T
 */
template<typename T, typename Enable = void> struct is_hashable : std::false_type {};
template<typename T> struct is_hashable<T, std::void_t<decltype(std::hash<T>()(std::declval<T&>()))>> : std::true_type {};
template<typename T> using is_hashable_t = typename is_hashable<T>::type;

/**
 * is_instantiation_of<T, I> is true_type iff I is a template instantiation of T (e.g. vector<int> is an instantiation of vector)
 *  Example:
 *    is_instantiation_of_t<vector, vector<int>> // true
 *    is_instantiation_of_t<pair, pair<int, string>> // true
 *    is_instantiation_of_t<vector, pair<int, string>> // false
 */
template <template<class...> class Template, class T>
struct is_instantiation_of : std::false_type {};
template <template <class...> class Template, class... Args>
struct is_instantiation_of<Template, Template<Args...>> : std::true_type {};
template<template<class...> class Template, class T> using is_instantiation_of_t = typename is_instantiation_of<Template, T>::type;

/**
 * is_type_condition<C> is true_type iff C<...> is a type trait representing a condition (i.e. has a constexpr static bool ::value member)
 * Example:
 *   is_type_condition<std::is_reference>  // true
 */
template<template<typename> typename C, typename Enable = void>
struct is_type_condition : std::false_type {};
template<template<typename> typename C>
struct is_type_condition<C, std::enable_if_t<std::is_same<bool, std::remove_cv_t<decltype(C<int>::value)>>::value>> : std::true_type {};

_MATRICE_NAMESPACE_END