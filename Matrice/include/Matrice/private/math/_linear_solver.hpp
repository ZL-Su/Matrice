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
#include "_type_traits.h"
#include "core/matrix.h"

DGE_MATRICE_BEGIN
struct svd { enum { singular_value_decomposition = 4 }; };
struct spt { enum { spdtr_cholesky_decomposition = 2 }; };
struct luf { enum { lower_upper_tr_decomposition = 1 }; };
struct qrd { enum { ortho_upper_tr_decomposition = 3 }; };
_INTERNAL_BEGIN
// \brief parse and invoke the linear algebra kernel to _Op.
template<typename _Op, typename... _Ts>
MATRICE_GLOBAL void _Lak_adapter(_Ts&... args);
_INTERNAL_END

_DETAIL_BEGIN
// \brief solver traits class
template<class _Sty> struct _Solver_traits {};

// \brief CLASS TEMPLATE for linear solver base
template<typename _Derived>
class _Linear_solver_base {
	using _Myderived = _Derived;
	using _Mytraits = _Solver_traits<_Myderived>;
public:
	using matrix_type = typename _Mytraits::matrix_type;
	using value_type = typename _Mytraits::value_type;

	_Linear_solver_base(matrix_type& coeff) noexcept
		:_Mycoeff(coeff) {
		this->_Forward();
	}

	template<class _Rty>
	MATRICE_GLOBAL_INL auto solve(const _Rty& rhs) const {
		return static_cast<_Myderived*>(this)->backward(rhs);
	}

private:
	/**
	 *\brief perf. coeff matrix decomposition
	 */
	MATRICE_GLOBAL_INL void _Forward() {
		auto _Dptr = static_cast<_Myderived*>(this);
		using kernel_t = typename _Mytraits::op_type;
		if constexpr (is_same_v<kernel_t, void>) {
			//_Mycoeff is a general square matrix, _Forward does nothing.
		}
		if constexpr (is_same_v<kernel_t, svd>) {
			internal::_Lak_adapter<svd>(_Mycoeff, _Dptr->s(), _Dptr->v());
		}
		if constexpr (is_same_v<kernel_t, spt>) {
			internal::_Lak_adapter<spt>(_Mycoeff);
		}
	}

protected:
	matrix_type& _Mycoeff;
};

// \brief CLASS TEMPLATE for linear solver
// \param <_Mty> Matrice compatible matrix type
// \param <_Op> lapack kernel operator, default is void
template<class _Mty, typename _Op = void> 
class _Linear_solver : public _Linear_solver_base<_Linear_solver<_Mty, _Op>> {
	using _Myt = _Linear_solver;
	using _Mybase = _Linear_solver_base <_Linear_solver<_Mty, _Op>>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::matrix_type;
	_Linear_solver(matrix_type& coeff) : _Mybase(coeff) {
	}

private:

};

// \brief CLASS TEMPLATE for SVD based linear solver
// \param <_Mty> Matrice compatible matrix type
template<class _Mty>
class _Linear_solver<_Mty, svd> : 
	public _Linear_solver_base<_Linear_solver<_Mty, svd>> {
	using _Myt = _Linear_solver;
	using _Mybase = _Linear_solver_base<_Linear_solver<_Mty, svd>>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::matrix_type;

	// \note: the input coeff will be overwritten by U
	_Linear_solver(matrix_type& coeff) 
		: _Mybase(coeff), 
		_Mys(coeff.cols()), 
		_Myv(coeff.cols(), coeff.cols()) {
	}

	MATRICE_GLOBAL_INL decltype(auto) u() const noexcept {
		return (_Mybase::_Mycoeff);
	}
	MATRICE_GLOBAL_INL decltype(auto) u() noexcept {
		return (_Mybase::_Mycoeff);
	}
	MATRICE_GLOBAL_INL decltype(auto) s() const noexcept {
		return (_Mys);
	}
	MATRICE_GLOBAL_INL decltype(auto) s() noexcept {
		return (_Mys);
	}
	MATRICE_GLOBAL_INL decltype(auto) v() const noexcept {
		return (_Myv);
	}
	MATRICE_GLOBAL_INL decltype(auto) v() noexcept {
		return (_Myv);
	}
private:
	Matrix_<value_type, _Mty::cols_at_compiletime, 1> _Mys;
	Matrix_<value_type, _Mty::cols_at_compiletime> _Myv;
};

template<class _Mty, typename _Op> 
struct _Solver_traits<_Linear_solver<_Mty, _Op>>{
	using matrix_type = _Mty;
	using value_type = typename matrix_type::value_type;
	using op_type = _Op;
};
_DETAIL_END

// \brief: linear solver factory function
template<class _Mty, typename _Op>
auto make_linear_solver(_Mty& A, _Op) noexcept {
	return detail::_Linear_solver<_Mty, _Op>(A);
}
DGE_MATRICE_END