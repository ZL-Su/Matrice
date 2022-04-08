/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
#include "_layer_base.hpp"
#include "_input_layer.hpp"

MATRICE_NAMESPACE_BEGIN(dnn)
namespace detail {
/// <summary>
/// \brief CLASS interface, _Linear_layer<_Input, _OutSize, _HasBias>.
/// </summary>
/// <typeparam name="_Input"></typeparam>
template<typename _Input, size_t _Out, bool _HasBias> class _Linear_layer{};
/// <summary>
/// \brief TRAITS, for _Linear_layer type.
/// </summary>
/// <typeparam name="_Input"></typeparam>
template<typename _Input, size_t _Out, bool _HasBias>
struct _Layer_traits<_Linear_layer<_Input, _Out, _HasBias>> {
	static constexpr auto depth = _Layer_traits<_Input>::depth;
	static constexpr auto has_bias = _HasBias;
	using value_type = typename _Input::value_type;
};

/// <summary>
/// \brief Specialize to the input with type of 1-row and fixed-column matrix.
/// </summary>
/// <typeparam name="_Ty"></typeparam>
template<typename _Ty, int _N, size_t _Out>
class _Linear_layer<_Input_layer<Matrix_<_Ty, 1, _N>>, _Out, false> 
	: public _Layer<_Linear_layer<_Input_layer<Matrix_<_Ty, 1, _N>>, _Out, false>> {
public:
	using category = _Layer_tag::linear;
	using input_t = _Input_layer<Matrix_<_Ty, 1, _N>>;
	using value_t = _Ty;

	_Linear_layer() noexcept {
		_Myweights = decltype(_Myweights)::randn();
	}
	
	MATRICE_GLOBAL_INL constexpr auto insize() const noexcept {
		return _N;
	}
	MATRICE_GLOBAL_INL constexpr auto outsize() const noexcept {
		return _Out;
	}
	MATRICE_GLOBAL_INL constexpr auto depth() const noexcept {
		return 1;
	}

	/**
	 * \brief Forward computation with lazy evaluation.
	 * \return Expression of $matmul(X,W)$, and requires
	           to call the method 'eval()' to collect the
			   output feature vector.
	 */
	MATRICE_GLOBAL_INL auto operator()(const input_t& _X) const noexcept {
		return _X.data().mul(_Myweights);
	}
	/**
	 * \brief Inplace forward evaluation.
	 * \return The resulted feature vector. 
	 */
	MATRICE_GLOBAL_INL auto forward(const input_t& _X) const {
		return this->operator()(_X).eval();
	}

protected:
	Matrix_<value_t, _N, _Out> _Myweights;
};

/// <summary>
/// \brief Specialize to the input with type of fixed vector.
/// </summary>
template<typename _Ty, int _N, size_t _Out>
class _Linear_layer<_Input_layer<Vec_<_Ty, _N>>, _Out, false> :
	public _Linear_layer<_Input_layer<Matrix_<_Ty, 1, _N>>, _Out, false> {
};

}
template<typename _Input, size_t _Out, bool _HasBias=false>
using linear_layer = detail::_Linear_layer<_Input, _Out, _HasBias>;
MATRICE_NAMESPACE_END(dnn)
