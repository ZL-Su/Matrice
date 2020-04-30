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
struct svd { enum { singular_value_decomp = 4 }; };
struct spt { enum { spdtr_cholesky_decomp = 2 }; };
struct lud { enum { lower_upper_tr_decomp = 1 }; };
struct qrd { enum { ortho_upper_tr_decomp = 3 }; };
struct lls { enum { linear_least_squares = 0 }; };
struct solver_status {
	int value = 1;
	MATRICE_GLOBAL_FINL bool success() noexcept {
		return value == 0;
	}
};
_INTERNAL_BEGIN
// \brief Parse and invoke the linear algebra kernel according to _Op.
template<typename _Op, typename... _Ts>
MATRICE_GLOBAL solver_status _Lak_adapter(_Ts... args);
template<typename _Op, typename... _Ts>
MATRICE_GLOBAL solver_status _Bwd_adapter(_Ts... args);
template<typename _Op, typename... _Ts>
MATRICE_GLOBAL solver_status _Inv_adapter(_Ts... args);
template<typename... _Ts>
MATRICE_GLOBAL solver_status _Imp_adapter(_Ts... args);
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
	 *\brief Solve $\mathbf{AX} = \mathbf{B}$. The solver kernel is dispatched according to the solver operator specified to _Op.
	   The right-hand term $\mathbf{B}$ is allowed to hold more than one column when the linear kernel _Op is dgelom::spt.
	 *\param [args...] variadic arguments that may empty or more than one inputs. If args... is empty, the method solves a homogeneous linear system with form of $\mathbf{Ax} = \mathbf{0}$.
	 */
	template<class... _Args>
	MATRICE_GLOBAL_INL decltype(auto)solve(_Args&&... args) {
		return ((_Myderived*)(this))->_Xbackward(args...);
	}

	/**
	 *\brief Improve the solution of $\mathbf{Ax} = \mathbf{b}$ with an iterative error reduction. The method is recommended calling several times to refine the solution.
	 *\param [args] := $[\mathbf{A, b, x}]$ where $\mathbf{A}$ must be the original coeff matrix, $\mathbf{b}$ the original RHS vector and $\mathbf{x}$ the solved solution to be refined in the method.
	 */
	template<class... _Args>
	MATRICE_GLOBAL_INL decltype(auto)refine(_Args&... args) {
		return (this->_Improve(args...));
	}

	/**
	 *\brief Compute inverse of the coeff. matrix.
	 */
	MATRICE_GLOBAL_INL auto inv() {
		return (this->_Inverse());
	}

protected:
	/**
	 *\brief Perform coeff matrix decomposition
	 */
	MATRICE_GLOBAL_INL void _Forward() {
		using kernel_t = typename _Mytraits::op_type;

		auto _Mdp = static_cast<_Myderived*>(this);
		if constexpr (is_same_v<kernel_t, void>) {
			//_Mycoeff is a general square matrix, _Forward does nothing.
		}
		if constexpr (is_same_v<kernel_t, svd>) {
			internal::_Lak_adapter<svd>(view(_Mycoeff), view(_Mdp->s()),
				view(_Mdp->vt()));
		}
		if constexpr (is_same_v<kernel_t, spt>) {
			internal::_Lak_adapter<spt>(view(_Mycoeff));
		}
		if constexpr (is_same_v<kernel_t, lud>) {
			internal::_Lak_adapter<lud>(view(_Mycoeff), _Mdp->permut().data());
		}
	}

	MATRICE_GLOBAL_INL auto _Inverse() {
		using kernel_t = typename _Mytraits::op_type;

		auto _Mdp = (_Myderived*)(this);
		if constexpr (is_same_v<kernel_t, void>) {
			//general matrix inverse
			matrix_type _Ret(_Mycoeff.rows(), _Mycoeff.cols());

			return move(_Ret);
		}
		if constexpr (is_same_v<kernel_t, svd>) {
			Matrix_<value_type, matrix_type::cols_at_compiletime, matrix_type::rows_at_compiletime> _Ret(_Mycoeff.cols(), 
				_Mycoeff.rows());
			internal::_Inv_adapter<kernel_t>(view(_Mycoeff), view(_Mdp->s()),
				view(_Mdp->vt()), view(_Ret));
			return move(_Ret);
		}
		if constexpr (is_same_v<kernel_t, spt>) {
			matrix_type _Ret(_Mycoeff.rows(), _Mycoeff.cols());
			internal::_Inv_adapter<kernel_t>(view(_Mycoeff), view(_Ret));
			return move(_Ret);
		}
	}

private:
	template<typename _Bty, typename _Xty>
	decltype(auto) _Improve(matrix_type& _A, _Bty& b, _Xty& x)noexcept {
		internal::_Imp_adapter(view(_A), view(b), view(x));
		const auto _Res = this->solve(b);
		return (x = x - _Res);
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

};

// \brief CLASS TEMPLATE for linear least sqaure solver
// \param <_Mty> Matrice compatible matrix type
MATRICE_MAKE_LINEAR_SOLVER_SPEC_BEGIN(lls)
	_Linear_solver(matrix_type& coeff) noexcept
		:_Mybase(coeff) {
	}

	/**
	 *\brief Compute inverse of the coeff. matrix.
	 */
	MATRICE_GLOBAL_INL auto inv() const {
		auto _AtA = _Mybase::_Mycoeff.t().mul(_Mybase::_Mycoeff).eval();
		return (_AtA.inv().eval());
	}

	/**
	 *\brief Solve the linear system in form of
	 //tex: $\mathbf{AX} = \mathbf{B}$
	 */
	template<class _Bty>
	MATRICE_GLOBAL_INL auto solve(const _Bty& b) const {
		const auto _Atb = _Mybase::_Mycoeff.t().mul(b).eval();
		return (this->inv().mul(_Atb).eval());
	}
MATRICE_MAKE_LINEAR_SOLVER_SPEC_END(lls)

// \brief CLASS TEMPLATE for Chelosky decomposition based linear solver
// \param <_Mty> Matrice compatible matrix type
// \note This solver is only for symmetric positive definite systems. 
MATRICE_MAKE_LINEAR_SOLVER_SPEC_BEGIN(spt)
	/**
	 *\note The input coeff will be overwritten by $L$
	 */
	_Linear_solver(matrix_type& coeff) noexcept
		: _Mybase(coeff) {
		_Mybase::_Forward();
	}

	/**
	 *  \brief Returns the determinant of the given coeff matrix.
	 */
	MATRICE_GLOBAL_INL auto(det)()const noexcept {
		decltype(auto) _L = _Mybase::_Mycoeff;
		auto _Ret = one<value_type>;
		for (auto _Idx = 0; _Idx < _L.rows(); ++_Idx)
			_Ret *= _L[_Idx][_Idx];
		return sqr(_Ret);
	}

	/**
	 * \brief Perform back substitution to solve $\mathbf{AX} = \mathbf{B}$, where $\mathbf{B}$ allowed to have multi-cols and will be overwritten by $\mathbf{X}$.
	 * \note The method is for parsing by the base class of the linear solver rather than for external calling.
	 */
	template<class _Rty> [[non_external_callable]]
	MATRICE_GLOBAL_INL _Rty& _Xbackward(_Rty& B)const noexcept {
		decltype(auto) _L = _Mybase::_Mycoeff;
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_L.rows() == B.rows(),
			"The number of rows of the right-hand vector(s) is not identical to that of the coeff. matrix.");
#endif
		internal::_Bwd_adapter<_Myop>(view(_L), view(B));
		return (B);
	}
MATRICE_MAKE_LINEAR_SOLVER_SPEC_END(spt)

// \brief CLASS TEMPLATE for SVD based linear solver
// \param <_Mty> Matrice compatible matrix type
MATRICE_MAKE_LINEAR_SOLVER_SPEC_BEGIN(lud)
    /**
	 *\note The input coeff will be overwritten by $L\U$
	 */
	_Linear_solver(matrix_type& coeff) noexcept
	    : _Mybase(coeff), _Myidx(coeff.rows()+1){
	    _Mybase::_Forward();
    }
    
    /**
	 *\brief Get the row permutation produced by the partial pivoting.
	 */
    MATRICE_GLOBAL_INL decltype(auto) permut() const noexcept{
		return (_Myidx);
	}
	MATRICE_GLOBAL_INL decltype(auto) permut() noexcept {
		return (_Myidx);
	}

	/**
	 *\brief Get the parity produced by the partial pivoting.
	 */
	MATRICE_GLOBAL_INL decltype(auto) parity() const noexcept {
		return (_Myidx(_Myidx.rows()-1));
	}
	
	/**
	 *\brief Perform back substitution to solve
	 //tex:$\mathbf{Ax} = \mathbf{b}, \mathbf{b} \leftarrow \mathbf{x}$
	 *\note The method is auto-parsed by the solver engine.
	 */
    template<typename _Bty> [[non_external_callable]]
	MATRICE_GLOBAL_INL decltype(auto)_Xbackward(_Bty& b)const noexcept {

		return (b);
	}
private:
	array_n<int, _Mty::rows_at_compiletime + 1> _Myidx;
MATRICE_MAKE_LINEAR_SOLVER_SPEC_END(lud)

// \brief CLASS TEMPLATE for SVD based linear solver
// \param <_Mty> Matrice compatible matrix type
MATRICE_MAKE_LINEAR_SOLVER_SPEC_BEGIN(svd)
	/**
	 *  \note: the input coeff will be overwritten by $\mathbf{U}$
	 */
	_Linear_solver(matrix_type& coeff) 
		: _Mybase(coeff), 
		_Mys(coeff.cols()), _Myvt(coeff.cols(), coeff.cols()) {
		_Mybase::_Myeps *= value_type(0.5) * sqrt<value_type>(coeff.size());
		_Mybase::_Forward();
	}

	MATRICE_GLOBAL_INL decltype(auto)u()const noexcept {
		return (_Mybase::_Mycoeff);
	}
	MATRICE_GLOBAL_INL decltype(auto)u()noexcept {
		return (_Mybase::_Mycoeff);
	}
	MATRICE_GLOBAL_INL decltype(auto)s()const noexcept {
		return (_Mys);
	}
	MATRICE_GLOBAL_INL decltype(auto)s()noexcept {
		return (_Mys);
	}
	MATRICE_GLOBAL_INL decltype(auto)vt()const noexcept {
		return (_Myvt);
	}
	MATRICE_GLOBAL_INL decltype(auto)vt()noexcept {
		return (_Myvt);
	}

	/**
	 *\brief Return the rank of the given source matrix.
	 */
	MATRICE_GLOBAL_INL size_t rank()const noexcept {
		const auto _Thresh = _Mys(0) * _Mybase::_Myeps;
		auto _Ret = size_t(0);
		for (const auto& _Val : _Mys) {
			(_Val > _Thresh) ? ++_Ret : _Ret;
		}
		return _Ret;
	}
	/**
	 *\brief Return the nullity of the given source matrix.
	 */
	MATRICE_GLOBAL_INL size_t nullity()const noexcept {
		const auto _Thresh = _Mys(0) * _Mybase::_Myeps;
		auto _Ret = size_t(0);
		for (const auto& _Val : _Mys) {
			(_Val <= _Thresh) ? ++_Ret : _Ret;
		}
		return _Ret;
	}

	/**
     * \brief Solve the system $\mathbf{A}\mathbf{x} = \mathbf{b}$ with the pre-computed coeff. matrix, and $\mathbf{b}$ is overwritten by the solution $\mathbf{x}$.
	 * \note The method is for parsing by the base class of the linear solver rather than for external calling.
     */
	template<typename _Bty> [[non_external_callable]]
	MATRICE_GLOBAL_INL auto _Xbackward(_Bty& b) noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(b.size() == _Mybase::_Mycoeff.rows(), "Bad size in _Linear_solver<_Mty, svd>::backward(...).");
#endif
		decltype(_Mys) x;
		internal::_Bwd_adapter<_Myop>(view(_Mybase::_Mycoeff), view(_Mys), view(_Myvt), view(b), view(x));
		return forward<decltype(x)>(x);
	}
	/**
	 * \brief Solve $\mathbf{A}\mathbf{x} = \mathbf{0}$, but the method returns the view of the solution rather than a vector.
	 * \note The method is for parsing by the base class of the linear solver rather than for external calling.
	 */
	[[non_external_callable]]
	MATRICE_GLOBAL_INL auto _Xbackward()const noexcept {
		return _Myvt.rview(size_t(_Myvt.rows()) - 1);
	}
private:
	array_n<value_type, _Mty::cols_at_compiletime> _Mys;
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