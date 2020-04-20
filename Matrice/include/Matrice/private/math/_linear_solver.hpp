/**********************************************************************
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
***********************************************************************/
#pragma once
#include "../_type_traits.h"
#include "core/matrix.h"

DGE_MATRICE_BEGIN
struct svd { enum { singular_value_decomposition = 4 }; };
struct spt { enum { spdtr_cholesky_decomposition = 2 }; };
struct lud { enum { lower_upper_tr_decomposition = 1 }; };
struct qrd { enum { ortho_upper_tr_decomposition = 3 }; };
struct solver_status {
	int value = 1;
	MATRICE_GLOBAL_FINL bool success() noexcept {
		return value == 0;
	}
};
_INTERNAL_BEGIN
// \brief parse and invoke the linear algebra kernel to _Op.
template<typename _Op, typename... _Ts>
MATRICE_GLOBAL solver_status _Lak_adapter(_Ts... args);
template<typename _Op, typename... _Ts>
MATRICE_GLOBAL solver_status _Bwd_adapter(_Ts... args);
template<typename _Op, typename... _Ts>
MATRICE_GLOBAL solver_status _Inv_adapter(_Ts... args);
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
		:_Mycoeff(coeff), _Myeps(matrix_type::eps) {
	}

	/**
	 *\brief Solve AX = B. The solver kernel is dispatched according to the solver operator specified to _Op.
	 */
	template<class _Rty>
	MATRICE_GLOBAL_INL decltype(auto) solve(_Rty& b) {
		return ((_Myderived*)(this))->backward(b);
	}

	/**
	 *\brief Compute inverse of the coeff. matrix.
	 */
	MATRICE_GLOBAL_INL auto inv() {
		return this->_Inverse();
	}

protected:
	/**
	 *\brief Perform coeff matrix decomposition
	 */
	MATRICE_GLOBAL_INL void _Forward() {
		auto _Dptr = static_cast<_Myderived*>(this);
		using kernel_t = typename _Mytraits::op_type;
		if constexpr (is_same_v<kernel_t, void>) {
			//_Mycoeff is a general square matrix, _Forward does nothing.
		}
		if constexpr (is_same_v<kernel_t, svd>) {
			internal::_Lak_adapter<svd>(view(_Mycoeff), view(_Dptr->s()), view(_Dptr->vt()));
		}
		if constexpr (is_same_v<kernel_t, spt>) {
			internal::_Lak_adapter<spt>(view(_Mycoeff));
		}
		if constexpr (is_same_v<kernel_t, lud>) {
			internal::_Lak_adapter<spt>(view(_Mycoeff));
		}
	}

	MATRICE_GLOBAL_INL auto _Inverse() {
		auto _Dptr = (_Myderived*)(this);
		using kernel_t = typename _Mytraits::op_type;
		if constexpr (is_same_v<kernel_t, void>) {
			//general matrix inverse
			matrix_type _Ret(_Mycoeff.rows(), _Mycoeff.cols());

			return move(_Ret);
		}
		if constexpr (is_same_v<kernel_t, svd>) {
			Matrix_<value_type, matrix_type::cols_at_compiletime, matrix_type::rows_at_compiletime> _Ret(_Mycoeff.cols(), _Mycoeff.rows());
			_Ret.cview(0) = _Dptr->s();
			internal::_Inv_adapter<svd>(view(_Mycoeff), view(_Ret), view(_Dptr->vt()));
			return move(_Ret);
		}
		if constexpr (is_same_v<kernel_t, spt>) {
			matrix_type _Ret(_Mycoeff.rows(), _Mycoeff.cols());
			internal::_Inv_adapter<spt>(view(_Mycoeff), view(_Ret));
			return move(_Ret);
		}
	}

protected:
	matrix_type& _Mycoeff; //ref to the given coefficient matrix
	value_type _Myeps;     //default roundoff threshold
};

#define MATRICE_MAKE_LINEAR_SOLVER_SPEC_BEGIN(OP) \
template<class _Mty> \
class _Linear_solver<_Mty, OP> : \
	public _Linear_solver_base<_Linear_solver<_Mty, OP>> { \
	using _Myt = _Linear_solver; \
	using _Mybase = _Linear_solver_base<_Linear_solver<_Mty, OP>>; \
	using _Myop = OP; \
public: \
	using typename _Mybase::value_type; \
	using typename _Mybase::matrix_type;

#define MATRICE_MAKE_LINEAR_SOLVER_SPEC_END(OP) };

// \brief CLASS TEMPLATE for linear solver with Gauss-Jordan alg.
// \param <_Mty> Matrice compatible matrix type
// \param <_Op> lapack kernel operator, default is void
template<class _Mty, typename _Op = void> 
class _Linear_solver : public _Linear_solver_base<_Linear_solver<_Mty, _Op>> {
	using _Myt = _Linear_solver;
	using _Mybase = _Linear_solver_base <_Linear_solver<_Mty, _Op>>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::matrix_type;
	_Linear_solver(matrix_type& coeff) noexcept 
		: _Mybase(coeff) {
		_Mybase::_Forward();
	}

private:

};

// \brief CLASS TEMPLATE for Chelosky decomposition based linear solver
// \param <_Mty> Matrice compatible matrix type
// \note This solver is only for symmetric positive definite systems. 
MATRICE_MAKE_LINEAR_SOLVER_SPEC_BEGIN(spt)
public:
	// \note The input coeff will be overwritten by $L$
	_Linear_solver(matrix_type& coeff) noexcept
		: _Mybase(coeff) {
		_Mybase::_Forward();
	}

	// \brief Perform back substitution to solve ${AX = B}$, where B allowed to have multi-cols and will be overwritten by X.
	template<class _Rty>
	MATRICE_GLOBAL_INL _Rty& backward(_Rty& B) noexcept {
		decltype(auto) _L = _Mybase::_Mycoeff;
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_L.rows() == B.rows(),
			"The number of rows of the right-hand vector(s) is not identical to that of the coeff. matrix.");
#endif
		internal::_Bwd_adapter<_Myop>(view(_L), view(B));
		return (B);
	}

	// \brief Returns the determinant of the given coeff matrix.
	MATRICE_GLOBAL_INL auto (det)() const noexcept {
		decltype(auto) _L = _Mybase::_Mycoeff;
		auto _Ret = one<value_type>;
		for (auto _Idx = 0; _Idx < _L.rows(); ++_Idx)
			_Ret *= _L[_Idx][_Idx];
		return sqr(_Ret);
	}
MATRICE_MAKE_LINEAR_SOLVER_SPEC_END(spt)

// \brief CLASS TEMPLATE for SVD based linear solver
// \param <_Mty> Matrice compatible matrix type
MATRICE_MAKE_LINEAR_SOLVER_SPEC_BEGIN(svd)
public:
	// \note: the input coeff will be overwritten by U
	_Linear_solver(matrix_type& coeff) 
		: _Mybase(coeff), 
		_Mys(coeff.cols()), 
		_Myvt(coeff.cols(), coeff.cols()) {
		_Mybase::_Forward();
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
	MATRICE_GLOBAL_INL decltype(auto) vt() const noexcept {
		return (_Myvt);
	}
	MATRICE_GLOBAL_INL decltype(auto) vt() noexcept {
		return (_Myvt);
	}

	template<typename _Rty>
	MATRICE_GLOBAL_INL _Rty& backward(_Rty& b) noexcept {

	}

	/**
	 *\brief Solve $\mathbf{A}\mathbf{x} = \mathbf{0}$, but the method returns the view of the solution rather than a vector.
	 */
	MATRICE_GLOBAL_INL auto solve() noexcept {
		return _Myvt.rview(size_t(_Myvt.rows()) - 1);
	}
private:
	Matrix_<value_type, _Mty::cols_at_compiletime, 1> _Mys;
	Matrix_<value_type, _Mty::cols_at_compiletime> _Myvt;

MATRICE_MAKE_LINEAR_SOLVER_SPEC_END(svd);

template<class _Mty, typename _Op> 
struct _Solver_traits<_Linear_solver<_Mty, _Op>>{
	using matrix_type = _Mty;
	using value_type = typename matrix_type::value_type;
	using op_type = _Op;
};

#undef MATRICE_MAKE_LINEAR_SOLVER_SPEC_BEGIN
#undef MATRICE_MAKE_LINEAR_SOLVER_SPEC_END
_DETAIL_END

// \brief: linear solver factory function
template<typename _Op, class _Mty>
auto make_linear_solver(_Mty& A) noexcept {
	return detail::_Linear_solver<_Mty, _Op>(A);
}
DGE_MATRICE_END