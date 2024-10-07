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
#include "core/vector.h"

DGE_MATRICE_BEGIN
namespace ade {
template<typename _Lty, int32_t _N> class auto_diff_variable;
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
	using variable_type = typename _Traits::variable_type;
	using scalar_type = typename _Traits::scalar_type;

	/// <summary>
	/// Evaluate the value of this expression.
	/// </summary>
	MATRICE_GLOBAL_FINL constexpr auto operator()()const noexcept {
		return this->value();
	}
	MATRICE_GLOBAL_FINL constexpr auto value()const noexcept {
		return _Dptr()->_Eval<_Myval_tag>();
	}
	
	/// <summary>
	/// Evaluate the gradient of this expression w.r.t. all variable elements.
	/// </summary>
	MATRICE_GLOBAL_FINL constexpr const auto grad()const noexcept {
		return _Dptr()->_Eval<_Mydiff_tag>();
	}
	MATRICE_GLOBAL_FINL constexpr auto grad() noexcept {
		return _Dptr()->_Eval<_Mydiff_tag>();
	}

	/**
	 * Evaluate the autodiff expression to a variable.
	 * @return An ```auto_diff_variable``` object.
	 */
	MATRICE_GLOBAL_FINL constexpr auto eval()const noexcept {
		variable_type var;
		var.value() = value();
		var.deriv() = grad();
		return var;
	}

private:
	MATRICE_GLOBAL_FINL decltype(auto)_Dptr()const noexcept {
		return static_cast<const derived_type*>(this);
	}
	MATRICE_GLOBAL_FINL decltype(auto)_Dptr() noexcept {
		return static_cast<derived_type*>(this);
	}
};
template<class T>
struct is_autodiff_exp<auto_diff_exp<T>> : std::true_type {};

/**
 * CLASS TEMPLATE, declare a vector variable for auto-diff.
 * @param <_Ty> Template type, any scalar type is allowed.
 * @param <_N> Template size, dimension of the vector variable.
 */
template<typename _Ty, int32_t _N=1>
class auto_diff_variable : public auto_diff_exp<auto_diff_variable<_Ty, _N>> 
{
	using _Myt = auto_diff_variable;
	using _Mybase = auto_diff_exp<auto_diff_variable<_Ty, _N>>;
public:
	/// scalar type
	using typename _Mybase::scalar_type;
	/// variable type
	using typename _Mybase::variable_type;
	/// vector type
	using value_type=detail::Vector<scalar_type, _N>;

	/// @brief Default ctor.
	MATRICE_GLOBAL_INL auto_diff_variable() noexcept = default;
	/// @brief Ctor from a vector ```x```.
	MATRICE_GLOBAL_INL auto_diff_variable(const value_type& x) noexcept
		:_Myval(x) {
	}
	/// @brief Ctor from a scalar ```x```.
	MATRICE_GLOBAL_INL auto_diff_variable(const scalar_type& x) noexcept
		:_Myval(x) {
	}
	/// @brief Ctor from a value vector `x` and grad vector `d`.
	MATRICE_GLOBAL_INL auto_diff_variable(const value_type& x, const value_type& d) noexcept
		:_Myval(x), _Mydev(d) {
	}
	/// @brief Ctor from a scalar `x`.
	MATRICE_GLOBAL_INL auto_diff_variable(const _Myt& other) noexcept
		:_Myval(other._Myval), _Mydev(other._Mydev) {
	}
	/// @brief Ctor from an autodiff expression `exp`.
	template<class Op>
	MATRICE_GLOBAL_INL auto_diff_variable(auto_diff_exp<Op> exp) noexcept
	{
		_Myval = exp.value();
		_Mydev = exp.grad();
	}
	
	/// @brief Copy from a scalar value `x`.
	MATRICE_GLOBAL_INL _Myt& operator=(const scalar_type& x) noexcept {
		_Myval = x;
		return (*this);
	}
	/// @brief Copy from another value `x`.
	MATRICE_GLOBAL_INL _Myt& operator=(const value_type& x) noexcept {
		_Myval = x;
		return (*this);
	}
	/// @brief Copy from another variable `other`.
	MATRICE_GLOBAL_INL _Myt& operator=(const _Myt& other) noexcept {
		_Myval = other._Myval;
		_Mydev = other._Mydev;
		return (*this);
	}

	template<class _Ety>
	MATRICE_GLOBAL_FINL constexpr auto& _Eval() const noexcept {
		if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>)
			return _Myval;
		if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>)
			return _Mydev;
		static_assert(
			is_any_of_v<_Ety, typename _Mybase::_Myval_tag, typename _Mybase::_Mydiff_tag>, 
			"Invalid type parameter in ::_Eval<_Ety>.");
	}
	template<class _Ety>
	MATRICE_GLOBAL_FINL constexpr auto& _Eval() noexcept {
		if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>)
			return _Myval;
		if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>)
			return _Mydev;
		static_assert(
			is_any_of_v<_Ety, typename _Mybase::_Myval_tag, typename _Mybase::_Mydiff_tag>, 
			"Invalid type parameter in ::_Eval<_Ety>.");
	}
private:
	value_type _Myval{0};
	value_type _Mydev{1};
};
template<class T> struct is_autodiff_var : std::false_type{};
template<typename _Ty, int32_t _N>
struct is_autodiff_exp<auto_diff_variable<_Ty, _N>> : std::true_type {};
template<typename _Ty, int32_t _N>
struct is_autodiff_var<auto_diff_variable<_Ty, _N>> : is_autodiff_exp<auto_diff_variable<_Ty, _N>>{};
template<typename _Ty, int32_t _N>
struct autodiff_exp_traits<auto_diff_variable<_Ty, _N>> {
	inline static constexpr auto dim_v = _N;
	using scalar_type = remove_all_t<_Ty>;
	using value_type = detail::Vector<scalar_type, dim_v>;
	using variable_type = ade::auto_diff_variable<scalar_type, dim_v>;
};

/**
 * Factory function for declaring an auto-differential N-dim variable.
 * @tparam N dimension of the variable, default is `1`.
 * @tparam T scalar type used in the variable, default is `float`.
 * @param scalar A scalar value for initializing the variable.
 */
template<int32_t N=1, typename T=float_t>
MATRICE_GLOBAL_INL auto make_variable(T scalar) noexcept {
	return auto_diff_variable<T, N>(scalar);
}

/** <summary>
 * \brief Auto differential expression for binary addition.
 * </summary>
 * <typeparam name="_Lty">Left-hand-side operand type, can be scalar, variable, and exp.</typeparam>
 * <typeparam name="_Rty">Type for right hand side operand</typeparam>
 */
template<class _Lty, class _Rty>
class diff_add_exp : public auto_diff_exp<diff_add_exp<_Lty, _Rty>> {
	using _Mybase = auto_diff_exp<diff_add_exp<_Lty, _Rty>>;
public:
	using typename _Mybase::scalar_type;
	using typename _Mybase::variable_type;

	MATRICE_GLOBAL_INL diff_add_exp(const _Lty& x, const _Rty& y)noexcept
		:m_varx(x), m_vary(y) {
	}

	template<class _Ety>
	MATRICE_GLOBAL_FINL constexpr auto _Eval() const noexcept {
		if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>)
			return (m_varx.value() + m_vary.value()).eval();
		if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>)
			return (m_varx.grad() + m_vary.grad()).eval();
		static_assert(is_any_of_v<_Ety, typename _Mybase::_Myval_tag, typename _Mybase::_Mydiff_tag>, 
			"Invalid type parameter in ::_Eval<_Ety>.");
	}

private:
	const _Lty& m_varx;
	const _Rty& m_vary;
};
template<class _Lty>
class diff_add_exp<_Lty, typename _Lty::scalar_type> 
    : public auto_diff_exp<diff_add_exp<_Lty, typename _Lty::scalar_type>> {
	using _Mybase = auto_diff_exp<diff_add_exp<_Lty, typename _Lty::scalar_type>>;
public:
	using typename _Mybase::scalar_type;
	using typename _Mybase::variable_type;

	MATRICE_GLOBAL_INL diff_add_exp(const _Lty& x, scalar_type s)noexcept
		:m_varx(x), m_const(s) {
	}

	template<class _Ety>
	MATRICE_GLOBAL_FINL constexpr auto _Eval() const noexcept {
		if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>)
			return (m_varx.value() + m_const).eval();
		if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>)
			return m_varx.grad();
		static_assert(is_any_of_v<_Ety, typename _Mybase::_Myval_tag, typename _Mybase::_Mydiff_tag>,
			"Invalid type parameter in ::_Eval<_Ety>.");
	}
private:
	const _Lty& m_varx;
	scalar_type m_const;
};
template<class _Rty>
struct diff_add_exp<typename _Rty::scalar_type, _Rty>
	: public diff_add_exp<_Rty, typename _Rty::scalar_type> {
	using _Mybase = diff_add_exp<_Rty, typename _Rty::scalar_type>;
public:
	using typename _Mybase::scalar_type;
	using typename _Mybase::variable_type;

	MATRICE_GLOBAL_INL diff_add_exp(scalar_type s, const _Rty& x)noexcept
		:_Mybase(x, s) {
	}
};
template<typename T, typename U>
struct is_autodiff_exp<diff_add_exp<T, U>> : std::true_type {};
template<typename T, typename U>
struct autodiff_exp_traits<diff_add_exp<T, U>> {
	using _Adexp_type = conditional_t<is_scalar_v<remove_all_t<T>>, remove_all_t<U>, remove_all_t<T>>;
	using variable_type = typename autodiff_exp_traits<_Adexp_type>::variable_type;
	using scalar_type = typename autodiff_exp_traits<variable_type>::scalar_type;
};

template<class _Lty, class _Rty>
class diff_sub_exp : public auto_diff_exp<diff_sub_exp<_Lty, _Rty>> {
	using _Mybase = auto_diff_exp<diff_sub_exp<_Lty, _Rty>>;
public:
	using typename _Mybase::variable_type;
	
	MATRICE_GLOBAL_INL diff_sub_exp(const _Lty& x, const _Rty& y)noexcept
		:m_varx(x), m_vary(y) {
	}

	template<class _Ety>
	MATRICE_GLOBAL_FINL constexpr auto _Eval() const noexcept {
		if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>)
			return  (m_varx.value() - m_vary.value()).eval();
		if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>)
			return (m_varx.grad() - m_vary.grad()).eval();
		static_assert(is_any_of_v<_Ety, typename _Mybase::_Myval_tag, typename _Mybase::_Mydiff_tag>,
			"Invalid type parameter in ::_Eval<_Ety>.");
	}

private:
	const _Lty& m_varx;
	const _Rty& m_vary;
};
template<class _Lty>
class diff_sub_exp<_Lty, typename _Lty::scalar_type> 
    : public auto_diff_exp<diff_sub_exp<_Lty, typename _Lty::scalar_type>> {
	using _Mybase = auto_diff_exp<diff_sub_exp<_Lty, typename _Lty::scalar_type>>;
public:
	using typename _Mybase::variable_type;
	using typename _Mybase::scalar_type;

	MATRICE_GLOBAL_INL diff_sub_exp(const _Lty& x, scalar_type s)noexcept
		:m_varx(x), m_const(s) {
	}

	template<class _Ety>
	MATRICE_GLOBAL_FINL constexpr auto _Eval() const noexcept {
		if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>)
			return (m_varx.value() - m_const).eval();
		if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>)
			return m_varx.grad();
		static_assert(is_any_of_v<_Ety, typename _Mybase::_Myval_tag, typename _Mybase::_Mydiff_tag>,
			"Invalid type parameter in ::_Eval<_Ety>.");
	}
private:
	const _Lty& m_varx;
	scalar_type m_const;
};
template<class _Rty>
struct diff_sub_exp<typename _Rty::scalar_type, _Rty>
	: public diff_sub_exp<_Rty, typename _Rty::scalar_type> {
	using _Mybase = diff_sub_exp<_Rty, typename _Rty::scalar_type>;
public:
	using typename _Mybase::scalar_type;

	MATRICE_GLOBAL_INL diff_sub_exp(scalar_type s, const _Rty& x)noexcept
		:_Mybase(x, s) {
	}
};
template<typename T, typename U>
struct is_autodiff_exp<diff_sub_exp<T, U>> : std::true_type {};
template<typename T, typename U>
struct autodiff_exp_traits<diff_sub_exp<T, U>> {
	using _Adexp_type = conditional_t<is_scalar_v<remove_all_t<T>>, remove_all_t<U>, remove_all_t<T>>;
	using variable_type = typename autodiff_exp_traits<_Adexp_type>::variable_type;
	using scalar_type = typename autodiff_exp_traits<variable_type>::scalar_type;
};

template<class _Lty, class _Rty>
class diff_mul_exp : public auto_diff_exp<diff_mul_exp<_Lty, _Rty>> {
	using _Mybase = auto_diff_exp<diff_mul_exp<_Lty, _Rty>>;
public:
	using typename _Mybase::variable_type;

	MATRICE_GLOBAL_INL diff_mul_exp(const _Lty& x, const _Rty& y)noexcept
		:m_varx(x), m_vary(y) {
	}

	template<class _Ety=typename _Mybase::_Myval_tag>
	MATRICE_GLOBAL_FINL constexpr auto _Eval() const noexcept {
		if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>)
			return (m_varx.value() * m_vary.value()).eval();
		if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>)
			return (m_varx.grad() * m_vary.value() + m_varx.value() * m_vary.grad()).eval();
		static_assert(is_any_of_v<_Ety, typename _Mybase::_Myval_tag, typename _Mybase::_Mydiff_tag>,
			"Invalid type parameter in ::_Eval<_Ety>.");
	}

private:
	const _Lty& m_varx;
	const _Rty& m_vary;
};
template<class _Lty>
class diff_mul_exp<_Lty, typename _Lty::scalar_type> 
	: public auto_diff_exp<diff_mul_exp<_Lty, typename _Lty::scalar_type>> {
	using _Mybase = auto_diff_exp<diff_mul_exp<_Lty, typename _Lty::scalar_type>>;
public:
	using typename _Mybase::scalar_type;
	MATRICE_GLOBAL_FINL diff_mul_exp(const _Lty& x, scalar_type s)noexcept
		:m_varx(x), m_const(s) {
	}

	template<class _Ety>
	MATRICE_GLOBAL_FINL constexpr auto _Eval() const noexcept {
		if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>)
			return (m_varx.value() * m_const).eval();
		if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>)
			return (m_varx.grad() * m_const).eval();
		static_assert(is_any_of_v<_Ety, typename _Mybase::_Myval_tag, typename _Mybase::_Mydiff_tag>,
			"Invalid type parameter in ::_Eval<_Ety>.");
	}

private:
	const _Lty& m_varx;
	scalar_type m_const;
};
template<class _Rty>
struct diff_mul_exp<typename _Rty::scalar_type, _Rty>
	: public diff_mul_exp<_Rty, typename _Rty::scalar_type> {
	using _Mybase = diff_mul_exp<_Rty, typename _Rty::scalar_type>;
public:
	using typename _Mybase::scalar_type;
	MATRICE_GLOBAL_INL diff_mul_exp(scalar_type s, const _Rty& x)noexcept
		:_Mybase(x, s) {
	}
};
template<typename T, typename U>
struct is_autodiff_exp<diff_mul_exp<T, U>> : std::true_type {};
template<typename T, typename U>
struct autodiff_exp_traits<diff_mul_exp<T, U>> {
	using _Adexp_type = conditional_t<is_scalar_v<remove_all_t<T>>, remove_all_t<U>, remove_all_t<T>>;
	using variable_type = typename autodiff_exp_traits<_Adexp_type>::variable_type;
	using scalar_type = typename autodiff_exp_traits<variable_type>::scalar_type;
};

#define DGELOM_MAKE_UNARY_AUTODIFF_EXP(NAME, VEXP, DEXP) \
template<class _Lty> class diff_##NAME##_exp \
: public auto_diff_exp<diff_##NAME##_exp<_Lty>> { \
	using _Mybase = auto_diff_exp<diff_##NAME##_exp<_Lty>>; \
public: \
	using typename _Mybase::variable_type; \
	using typename _Mybase::scalar_type; \
	MATRICE_GLOBAL_INL diff_##NAME##_exp(const _Lty& x) noexcept \
		:m_varx(x) {} \
	template<class _Ety> \
	MATRICE_GLOBAL_FINL constexpr auto _Eval() const noexcept { \
	if constexpr (is_same_v<_Ety, typename _Mybase::_Myval_tag>) \
	{	VEXP;  } \
	if constexpr (is_same_v<_Ety, typename _Mybase::_Mydiff_tag>) \
	{   DEXP;  }\
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
	using variable_type = typename autodiff_exp_traits<remove_all_t<T>>::variable_type;\
	using scalar_type = typename autodiff_exp_traits<variable_type>::scalar_type;\
}; \
template<class _Lty, MATRICE_ENABLE_IF(is_autodiff_exp_v<_Lty>)> \
MATRICE_GLOBAL_FINL auto NAME(const _Lty& x) noexcept { \
	return diff_##NAME##_exp{ x }; \
}

DGELOM_MAKE_UNARY_AUTODIFF_EXP(sin,
	MATRICE_USE_STD(sin)
	return sin(m_varx.value()),
	MATRICE_USE_STD(cos)
	return cos(m_varx.value()) * m_varx.grad()
);
DGELOM_MAKE_UNARY_AUTODIFF_EXP(cos,
	MATRICE_USE_STD(cos)
	return cos(m_varx.value()),
	MATRICE_USE_STD(sin)
	return -sin(m_varx.value()) * m_varx.grad()
);
DGELOM_MAKE_UNARY_AUTODIFF_EXP(exp,
	MATRICE_USE_STD(exp)
	return exp(m_varx.value()),
	return this->value() * m_varx.grad()
);
DGELOM_MAKE_UNARY_AUTODIFF_EXP(cube,
	MATRICE_USE_STD(pow)
	return pow(m_varx.value(), 3),
	return 3 * m_varx.value() * m_varx.value() * m_varx.grad()
);
DGELOM_MAKE_UNARY_AUTODIFF_EXP(sq,
	return m_varx.value() * m_varx.value(),
	return 2 * m_varx.value() * m_varx.grad()
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