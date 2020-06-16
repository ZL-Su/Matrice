/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

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
#include "_ad_utils.h"

DGE_MATRICE_BEGIN
namespace ade {
template<typename _Lty> class auto_diff_variable;
template<typename _Lty, typename _Rty> class diff_add_exp;
template<typename _Lty, typename _Rty> class diff_sub_exp;
template<typename _Lty, typename _Rty> class diff_mul_exp;

template<class _Derived, class _Traits = autodiff_exp_traits<_Derived>>
class auto_diff_exp {
	using _Myt = auto_diff_exp;

protected:
	struct _Myval_tag {};
	struct _Mydiff_tag {};

public:
	using derived_type = _Derived;
	using value_type = typename _Traits::value_type;

	/// <summary>
	/// Evaluate the value of this expression.
	/// </summary>
	MATRICE_GLOBAL_FINL constexpr value_type operator()()const noexcept {
		return this->value();
	}
	MATRICE_GLOBAL_FINL constexpr value_type value()const noexcept {
		return _Dptr()->_Eval<_Myval_tag>();
	}
	
	/// <summary>
	/// Evaluate the derivative of this expression
	/// </summary>
	MATRICE_GLOBAL_FINL constexpr const value_type deriv()const noexcept {
		return _Dptr()->_Eval<_Mydiff_tag>();
	}
	MATRICE_GLOBAL_FINL constexpr value_type deriv() noexcept {
		return _Dptr()->_Eval<_Mydiff_tag>();
	}

private:
	MATRICE_GLOBAL_FINL decltype(auto) _Dptr()const noexcept {
		return static_cast<const derived_type*>(this);
	}
	MATRICE_GLOBAL_FINL decltype(auto) _Dptr() noexcept {
		return static_cast<derived_type*>(this);
	}
};
template<class T>
struct is_autodiff_exp<auto_diff_exp<T>> : std::true_type {};

template<typename _Ty>
class auto_diff_variable : public auto_diff_exp<auto_diff_variable<_Ty>> {
	using _Mybase = auto_diff_exp<auto_diff_variable<_Ty>>;
public:
	using typename _Mybase::value_type;

	MATRICE_GLOBAL_INL auto_diff_variable() noexcept = default;
	MATRICE_GLOBAL_INL auto_diff_variable(value_type x) noexcept
		:_Myval(x) {}

	template<class _Ety>
	MATRICE_GLOBAL_FINL constexpr auto _Eval() const noexcept {
		if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>)
			return _Myval;
		if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>)
			return value_type(1);
		static_assert(is_any_of_v<_Ety, typename _Mybase::_Myval_tag, typename _Mybase::_Mydiff_tag>, 
			"Invalid type parameter in ::_Eval<_Ety>.");
	}
private:
	value_type _Myval = 0;
};
template<typename _Ty>
struct is_autodiff_exp<auto_diff_variable<_Ty>> : std::true_type {};
template<typename _Ty>
struct autodiff_exp_traits<auto_diff_variable<_Ty>> { using value_type = remove_all_t<_Ty>; };

/// <summary>
/// Auto differential expression for binary addition.
/// </summary>
/// <typeparam name="_Lty">Type for left hand side operand</typeparam>
/// <typeparam name="_Rty">Type for right hand side operand</typeparam>
template<class _Lty, class _Rty>
class diff_add_exp : public auto_diff_exp<diff_add_exp<_Lty, _Rty>> {
	using _Mybase = auto_diff_exp<diff_add_exp<_Lty, _Rty>>;
public:
	using typename _Mybase::value_type;

	MATRICE_GLOBAL_INL diff_add_exp(const _Lty& x, const _Rty& y)noexcept
		:m_varx(x), m_vary(y) {
	}

	template<class _Ety>
	MATRICE_GLOBAL_FINL constexpr auto _Eval() const noexcept {
		if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>)
			return  m_varx.value() + m_vary.value();;
		if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>)
			return m_varx.deriv() + m_vary.deriv();
		static_assert(is_any_of_v<_Ety, typename _Mybase::_Myval_tag, typename _Mybase::_Mydiff_tag>, 
			"Invalid type parameter in ::_Eval<_Ety>.");
	}

private:
	const _Lty& m_varx;
	const _Rty& m_vary;
};
template<class _Lty>
class diff_add_exp<_Lty, typename _Lty::value_type> : public auto_diff_exp<diff_add_exp<_Lty, typename _Lty::value_type>> {
	using _Mybase = auto_diff_exp<diff_add_exp<_Lty, typename _Lty::value_type>>;
public:
	using typename _Mybase::value_type;
	MATRICE_GLOBAL_INL diff_add_exp(const _Lty& x, const value_type s)noexcept
		:m_varx(x), m_const(s) {
	}

	template<class _Ety>
	MATRICE_GLOBAL_FINL constexpr auto _Eval() const noexcept {
		if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>)
			return  m_varx.value() + m_const;
		if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>)
			return m_varx.deriv();
		static_assert(is_any_of_v<_Ety, typename _Mybase::_Myval_tag, typename _Mybase::_Mydiff_tag>,
			"Invalid type parameter in ::_Eval<_Ety>.");
	}
private:
	const _Lty& m_varx;
	value_type m_const;
};
template<class _Rty>
struct diff_add_exp<typename _Rty::value_type, _Rty>
	: public diff_add_exp<_Rty, typename _Rty::value_type> {
	using _Mybase = diff_add_exp<_Rty, typename _Rty::value_type>;
public:
	using typename _Mybase::value_type;

	MATRICE_GLOBAL_INL diff_add_exp(const value_type s, const _Rty& x)noexcept
		:_Mybase(x, s) {
	}
};
template<typename T, typename U>
struct is_autodiff_exp<diff_add_exp<T, U>> : std::true_type {};
template<typename T, typename U>
struct autodiff_exp_traits<diff_add_exp<T, U>> { 
	using value_type = conditional_t<is_scalar_v<T>, T, typename T::value_type>;
};

template<class _Lty, class _Rty>
class diff_sub_exp : public auto_diff_exp<diff_sub_exp<_Lty, _Rty>> {
	using _Mybase = auto_diff_exp<diff_sub_exp<_Lty, _Rty>>;
public:
	using typename _Mybase::value_type;
	
	MATRICE_GLOBAL_INL diff_sub_exp(const _Lty& x, const _Rty& y)noexcept
		:m_varx(x), m_vary(y) {
	}

	template<class _Ety>
	MATRICE_GLOBAL_FINL constexpr auto _Eval() const noexcept {
		if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>)
			return  m_varx.value() - m_vary.value();
		if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>)
			return m_varx.deriv() - m_vary.deriv();
		static_assert(is_any_of_v<_Ety, typename _Mybase::_Myval_tag, typename _Mybase::_Mydiff_tag>,
			"Invalid type parameter in ::_Eval<_Ety>.");
	}

private:
	const _Lty& m_varx;
	const _Rty& m_vary;
};
template<class _Lty>
class diff_sub_exp<_Lty, typename _Lty::value_type> : public auto_diff_exp<diff_sub_exp<_Lty, typename _Lty::value_type>> {
	using _Mybase = auto_diff_exp<diff_sub_exp<_Lty, typename _Lty::value_type>>;
public:
	using typename _Mybase::value_type;

	MATRICE_GLOBAL_INL diff_sub_exp(const _Lty& x, const value_type s)noexcept
		:m_varx(x), m_const(s) {
	}

	template<class _Ety>
	MATRICE_GLOBAL_FINL constexpr auto _Eval() const noexcept {
		if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>)
			return  m_varx.value() - m_const;
		if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>)
			return m_varx.deriv();
		static_assert(is_any_of_v<_Ety, typename _Mybase::_Myval_tag, typename _Mybase::_Mydiff_tag>,
			"Invalid type parameter in ::_Eval<_Ety>.");
	}
private:
	const _Lty& m_varx;
	value_type m_const;
};
template<class _Rty>
struct diff_sub_exp<typename _Rty::value_type, _Rty>
	: public diff_sub_exp<_Rty, typename _Rty::value_type> {
	using _Mybase = diff_sub_exp<_Rty, typename _Rty::value_type>;
public:
	using typename _Mybase::value_type;

	MATRICE_GLOBAL_INL diff_sub_exp(const value_type s, const _Rty& x)noexcept
		:_Mybase(x, s) {
	}
};
template<typename T, typename U>
struct is_autodiff_exp<diff_sub_exp<T, U>> : std::true_type {};
template<typename T, typename U>
struct autodiff_exp_traits<diff_sub_exp<T, U>> {
	using value_type = conditional_t<is_scalar_v<T>, T, typename T::value_type>;
};

template<class _Lty, class _Rty>
class diff_mul_exp : public auto_diff_exp<diff_mul_exp<_Lty, _Rty>> {
	using _Mybase = auto_diff_exp<diff_mul_exp<_Lty, _Rty>>;
public:
	using typename _Mybase::value_type;

	MATRICE_GLOBAL_INL diff_mul_exp(const _Lty& x, const _Rty& y)noexcept
		:m_varx(x), m_vary(y) {
	}

	template<class _Ety>
	MATRICE_GLOBAL_FINL constexpr auto _Eval() const noexcept {
		if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>)
			return  m_varx.value() * m_vary.value();
		if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>)
			return m_varx.deriv() * m_vary.value() + m_varx.value() * m_vary.deriv();
		static_assert(is_any_of_v<_Ety, typename _Mybase::_Myval_tag, typename _Mybase::_Mydiff_tag>,
			"Invalid type parameter in ::_Eval<_Ety>.");
	}

private:
	const _Lty& m_varx;
	const _Rty& m_vary;
};
template<class _Lty>
class diff_mul_exp<_Lty, typename _Lty::value_type> 
	: public auto_diff_exp<diff_mul_exp<_Lty, typename _Lty::value_type>> {
	using _Mybase = auto_diff_exp<diff_mul_exp<_Lty, typename _Lty::value_type>>;
public:
	using typename _Mybase::value_type;
	MATRICE_GLOBAL_FINL diff_mul_exp(const _Lty& x, const value_type s)noexcept
		:m_varx(x), m_const(s) {
	}

	template<class _Ety>
	MATRICE_GLOBAL_FINL constexpr auto _Eval() const noexcept {
		if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>)
			return  m_varx.value() * m_const;
		if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>)
			return m_varx.deriv() * m_const;
		static_assert(is_any_of_v<_Ety, typename _Mybase::_Myval_tag, typename _Mybase::_Mydiff_tag>,
			"Invalid type parameter in ::_Eval<_Ety>.");
	}

private:
	const _Lty& m_varx;
	value_type m_const;
};
template<class _Rty>
struct diff_mul_exp<typename _Rty::value_type, _Rty>
	: public diff_mul_exp<_Rty, typename _Rty::value_type> {
	using _Mybase = diff_mul_exp<_Rty, typename _Rty::value_type>;
public:
	using typename _Mybase::value_type;
	MATRICE_GLOBAL_INL diff_mul_exp(const value_type s, const _Rty& x)noexcept
		:_Mybase(x, s) {
	}
};
template<typename T, typename U>
struct is_autodiff_exp<diff_mul_exp<T, U>> : std::true_type {};
template<typename T, typename U>
struct autodiff_exp_traits<diff_mul_exp<T, U>> {
	using value_type = conditional_t<is_scalar_v<T>, T, typename T::value_type>;
};

#define DGELOM_MAKE_UNARY_AUTODIFF_EXP(NAME, VEXP, DEXP) \
	template<class _Lty> class diff_##NAME##_exp \
	: public auto_diff_exp<diff_##NAME##_exp<_Lty>> { \
		using _Mybase = auto_diff_exp<diff_##NAME##_exp<_Lty>>; \
	public: \
		using _Mybase::value_type; \
		MATRICE_GLOBAL_INL diff_##NAME##_exp(const _Lty& x) noexcept \
			:m_varx(x) {} \
		template<class _Ety> \
	    MATRICE_GLOBAL_FINL constexpr auto _Eval() const noexcept { \
	    if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>) \
		    VEXP; \
	    if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>) \
		    DEXP; \
	    static_assert(is_any_of_v<_Ety, typename _Mybase::_Myval_tag, \
			typename _Mybase::_Mydiff_tag>, \
		    "Invalid type parameter in ::_Eval<_Ety>."); \
		}\
	private: \
		const _Lty& m_varx; \
	}; \
	template<typename T> \
	struct is_autodiff_exp<diff_##NAME##_exp<T>> : std::true_type {}; \
	template<typename T> \
	struct autodiff_exp_traits<diff_##NAME##_exp<T>> { \
		using value_type = conditional_t<is_scalar_v<T>, T, typename T::value_type>; \
	}; \
	template<class _Lty, MATRICE_ENABLE_IF(is_autodiff_exp_v<_Lty>)> \
	MATRICE_GLOBAL_FINL auto NAME(const _Lty& x) noexcept { \
		return diff_##NAME##_exp{ x }; \
	}

DGELOM_MAKE_UNARY_AUTODIFF_EXP(sin,
	return std::sin(m_varx.value()),
	return this->value() * m_varx.deriv()
);
DGELOM_MAKE_UNARY_AUTODIFF_EXP(cos,
	return std::cos(m_varx.value()),
	return -std::sin(m_varx.value()) * m_varx.deriv()
);
DGELOM_MAKE_UNARY_AUTODIFF_EXP(exp,
	return std::exp(m_varx.value()),
	return this->value() * m_varx.deriv()
);
DGELOM_MAKE_UNARY_AUTODIFF_EXP(cube,
	return std::pow(m_varx.value(), 3),
	return 3 * m_varx.value() * m_varx.value() * m_varx.deriv()
);
DGELOM_MAKE_UNARY_AUTODIFF_EXP(sq,
	return m_varx.value() * m_varx.value(),
	return 2 * m_varx.value() * m_varx.deriv()
);

template<class _Lty, class _Rty>
MATRICE_GLOBAL_FINL auto operator+(const _Lty& l, const _Rty& r)noexcept {
	return diff_add_exp{ l, r };
}
template<class _Lty, class _Rty>
MATRICE_GLOBAL_FINL auto operator-(const _Lty& l, const _Rty& r)noexcept {
	return diff_sub_exp{ l, r };
}
template<class _Lty, class _Rty>
MATRICE_GLOBAL_FINL auto operator*(const _Lty& l, const _Rty& r)noexcept {
	return diff_mul_exp{ l, r };
}
#undef DGELOM_MAKE_UNARY_AUTODIFF_EXP
}
DGE_MATRICE_END