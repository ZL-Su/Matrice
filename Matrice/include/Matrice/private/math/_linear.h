/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

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
********************************************************************/
#pragma once
#include "../_tag_defs.h"
#include "../nonfree/blas_lapack_kernel.h"

DGE_MATRICE_BEGIN _DETAIL_BEGIN

// Forward declaration
template<typename _Mty, typename _Tag> 
class _Matrix_decomposition {};

/// <summary>
/// \brief CLASS template for matrix decomposition.
/// </summary>
/// <typeparam name="_Mty">matrix type</typeparam>
template<typename _Mty> 
class _Matrix_decomposition<_Mty, tag::_Linear_spd_tag> {
	static_assert(is_matrix_v<_Mty>, 
		"_Mty in _Matrix_decomposition must be a Matrix_ type.");
public:
	using value_type = typename _Mty::value_type;
	using category = tag::_Linear_spd_tag;

	/**
	 * \brief Declare a matrix decomposition instance.
	 * \param "_A" The coeff matrix to be processed. 
	   We take a reference of the matrix such that the decomposition is built inplace.
	 */
	MATRICE_HOST_INL _Matrix_decomposition(const _Mty& _A) 
		: _Mycoef(_A) {}

	/**
	 * \brief Forward to perform matrix decomposition with a specified linear kernel.
	 * \return A reference to the overwritten coefficient matrix A.
	 */
	MATRICE_HOST_INL decltype(auto) forward() {
		_Lapack_kernel_impl<value_type>::spd(_Mycoef);
		return (_Mycoef);
	}

	/**
	 * \brief In-place solve linear system AX=B with the factorized A.
	 * \param "_X" A vector or matrix holds the right-hand side B of the linear system.
	 * \return The solution, which is stored in _X.
	 */
	template<typename _Rhs>
	MATRICE_HOST_INL decltype(auto) backward(const _Rhs& _X) const {
		return _Lapack_backward_impl<CHD>::eval(_Mycoef, _X);
	}

	/**
	 * \brief Accessor for the multiplier or damping factor $\lambda$.
	 */
	MATRICE_HOST_INL decltype(auto) lambda() const noexcept {
		return _Mylamb;
	}
	MATRICE_HOST_INL decltype(auto) lambda() noexcept {
		return _Mylamb;
	}

private:
	const _Mty& _Mycoef;
	typename _Mty::value_type _Mylamb{ 0 };
};

struct _Linear {

	template<typename _Maty> class _Base {
		static_assert(is_matrix_v<_Maty>, "");
	public:
		using matrix_type = _Maty;
		using value_type = typename matrix_type::value_type;

		MATRICE_HOST_INL _Base(matrix_type& _Data) noexcept
			: _Mydata(_Data) {}

		/**
		 * \Return precomputed matrix.
		 */
		MATRICE_HOST_INL decltype(auto)operator()()const noexcept { 
			return (_Mydata); 
		}
		MATRICE_HOST_INL decltype(auto)operator()()noexcept {
			return (_Mydata);
		}

	protected:
		matrix_type& _Mydata;
	};
	template<typename _Maty, solver_type _Tag> class _Decomp {};

	template<typename _Maty>
	class _Decomp<_Maty, solver_type::CHD> 
		: _Base<_Maty> {
		using _Mybase = _Base<_Maty>;
	public:
		using typename _Mybase::matrix_type;
		using typename _Mybase::value_type;
		static constexpr auto solver_tag = solver_type::CHD;

		MATRICE_HOST_INL _Decomp(matrix_type& _Data) 
			: _Mybase(_Data) {
			this->operator();
		}

	private:
		MATRICE_HOST_INL auto operator()() {
			_Lapack_kernel_impl<value_type>::spd(_Mybase::_Mydata.plvt());
		};
	};
};

_DETAIL_END 
template<typename... _Args>
using matrix_decomp = detail::_Matrix_decomposition<_Args...>;
DGE_MATRICE_END