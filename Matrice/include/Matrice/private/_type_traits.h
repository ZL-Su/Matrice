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
template<typename T> using remove_reference_t = remove_reference<T>;

template<typename T> struct type_bytes { enum { value = sizeof(T) }; };

template<bool _Test, typename T1, typename T2> struct conditonal {};
template<typename T1, typename T2> struct conditonal<true, T1, T2> { using type = T1; };
template<typename T1, typename T2> struct conditonal<false, T1, T2> { using type = T2; };
template<bool _Test, typename T1, typename T2> struct conditional {};
template<typename T1, typename T2> struct conditional<true, T1, T2> { using type = T1; };
template<typename T1, typename T2> struct conditional<false, T1, T2> { using type = T2; };
template<bool _Test, typename T1, typename T2> using conditional_t = typename conditional<_Test, T1, T2>::type;

template<typename Exp> struct expression_options { enum { value = Exp::flag | expr }; };

template<int _Opt> struct is_expression { enum { value = _Opt & expr == expr ? true : false }; };
template<typename Exp> constexpr bool is_expression_v = is_expression<Exp::options>::value;

template<int _Val> struct is_zero { enum { value = _Val == 0 ? true : false }; };

template<int _R, int _C> struct is_static {enum {value = _R > 0 && _C >0 ? true : false}; };

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

template<typename Mty, typename = std::enable_if_t<std::is_class_v<Mty>>>
struct matrix_traits : traits<Mty> {
	enum { M = Mty::CompileTimeRows };
	enum { N = Mty::CompileTimeRows };
};

template<typename Exp, typename = std::enable_if_t<std::is_class_v<Exp>>>
struct expression_traits : traits<Exp> {
	enum { options = Exp::options };
	static_assert(is_expression<options>::value, "Not expression type.");
};

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

template<typename Mty, typename = std::enable_if_t<std::is_class_v<Mty>>> 
struct layout_traits : traits<Mty> {
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

template<int _Rows = 0, int _Cols = 0> struct compile_time_size {
	enum { val_1 = 0x0001, val_2 = 0x0002, val_3 = 0x0003, val_4 = 0x0004 };
	enum {
		RunTimeDeduceInHost = 0,
		RunTimeDeduceInDevi = -1,
		CompileTimeRows = _Rows,
		CompileTimeCols = _Cols,
	};
	static const int RunTimeDeducedInHost = 0, RunTimeDeducedInDevice = -1;
};
_MATRICE_NAMESPACE_END