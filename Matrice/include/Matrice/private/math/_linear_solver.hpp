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
#include "../_type_traits.h"
#include "../../forward.hpp"

DGE_MATRICE_BEGIN
struct svd { static constexpr auto value = solver_type::SVD; };
struct spt { static constexpr auto value = solver_type::CHD; };
struct luf { static constexpr auto value = solver_type::LUF; };
struct qrd { static constexpr auto value = solver_type::QRD; };

_INTERNAL_BEGIN
	// \brief parse and invoke the linear algebra kernel to _Op.
	template<typename _Op, typename... _Ts>
	MATRICE_GLOBAL void _Lak_adapter(_Ts&... args);
_INTERNAL_END

_DETAIL_BEGIN
// \brief 
template<class _Sty> struct _Solver_traits {};

// \brief CLASS TEMPLATE for linear solver base
template<typename _Derived>
class _Linear_solver_base {
	using _Myderived = _Derived;
	using _Mytraits = _Solver_traits<_Myderived>;
public:
	using matrix_type = typename _Mytraits::matrix_type;
	_Linear_solver_base(matrix_type& coeff) noexcept
		:_Mycoeff(coeff) {
		this->_Forward();
	}

	template<class _Rty>
	MATRICE_GLOBAL_INL auto solve(const _Rty& rhs) const {
		return static_cast<_Myderived*>(this)->backward(rhs);
	}

protected:
	/**
	 *\brief perf. coeff matrix decomposition
	 */
	MATRICE_GLOBAL_INL void _Forward() {
		auto derived = static_cast<_Myderived*>(this);
		using kernel_t = typename _Mytraits::op_type;
		if constexpr (is_same_v<kernel_t, void>) {
			//_Mycoeff is a general square matrix, _Forward does nothing.
		}
		if constexpr (is_same_v<kernel_t, svd>) {
			internal::_Lak_adapter<svd>(_Mycoeff, derived->s(), derived()->vt());
		}
		if constexpr (is_same_v<kernel_t, spt>) {
			internal::_Lak_adapter<spt>(_Mycoeff);
		}
	}

	matrix_type& _Mycoeff;
};

// \brief CLASS TEMPLATE for linear solver
// \param <_Mty> Matrice compatible matrix type
// \param <_Op> lapack kernel operator
template<class _Mty, typename _Op = void> 
class _Linear_solver : public _Linear_solver_base<_Linear_solver<_Mty, _Op>> {
	using _Myt = _Linear_solver;
	using _Mybase = _Linear_solver_base <_Linear_solver<_Mty, _Op>>;
public:
	_Linear_solver(_Mty& coeff) : _Mybase(coeff) {
	}

private:

};

template<class _Mty>
class _Linear_solver<_Mty, svd> : 
	public _Linear_solver_base<_Linear_solver<_Mty, svd>> {
	using _Myt = _Linear_solver;
	using _Mybase = _Linear_solver_base <_Linear_solver<_Mty, svd>>;
public:
	_Linear_solver(_Mty& coeff) : _Mybase(coeff) {
	}

	MATRICE_GLOBAL_INL decltype(auto) u() const noexcept {
		return _Mybase::_Mycoeff;
	}
	MATRICE_GLOBAL_INL decltype(auto) u() noexcept {
		return _Mybase::_Mycoeff;
	}
	MATRICE_GLOBAL_INL decltype(auto) s() const noexcept {

	}
	MATRICE_GLOBAL_INL decltype(auto) s() noexcept {

	}
	MATRICE_GLOBAL_INL decltype(auto) vt() const noexcept {

	}
	MATRICE_GLOBAL_INL decltype(auto) vt() noexcept {

	}
private:

};

template<class _Mty, typename _Op> 
struct _Solver_traits<_Linear_solver<_Mty, _Op>>{
	using matrix_type = _Mty;
	using op_type = _Op;
};
_DETAIL_END

DGE_MATRICE_END