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
#include "../../../core"

DGE_MATRICE_BEGIN 
namespace dnn {

///<functional> implementations for DNNs </functional>
struct functional {

/* Sigmoid activation function */
struct sigmoid {
	/**
	 *\brief for executing nodal activation(s):
	  //tex:
	  //$$\sigma(x) = \dfrac{1}{1+e^{-x}}$$
	 *\param [x] input data with a type of scalar, ::Matrix or ::tensor.
	 */
	template<typename _Ty> 
	MATRICE_GLOBAL_INL static _Ty forward(const _Ty& x) {
		using value_t = conditional_t<is_scalar_v<_Ty>, _Ty, 
			typename _Ty::value_t>;
		return (one<value_t>/(one<value_t> + exp(-1*x)));
	}
	/**
	 *\brief for evaluating gradient(s) w.r.t. the input:
	  //tex:
	  //$$\sigma'(x) = \sigma(x)(1-\sigma(x))$$
	 *\param [y] data produced by sigmoid::forward().
	 */
	template<typename _Ty>
	MATRICE_GLOBAL_INL static _Ty grad(const _Ty& y) {
		using value_t = conditional_t<is_scalar_v<_Ty>, _Ty, 
			typename _Ty::value_t>;
		return (y * (one<value_t> - y));
	}
	/**
	 *\brief backward update gradient w.r.t. input params in-place.
	 *\param [g] input gradients, [y] tensor produced by forward().
	 */
	template<typename _Ty>
	MATRICE_GLOBAL_INL static _Ty& backward(_Ty& g, const _Ty& y) {
		const auto jac = sigmoid::grad(y);
		return (g);
	}
};

/* Tanh activation function */
struct tanh {
	/**
	 *\brief for executing nodal activation(s):
	  //tex:
	  //$$\sigma(x) = \dfrac{2}{1+e^{-2x}}-1$$
	 *\param [x] input data with a type of scalar, ::Matrix or ::tensor.
	 */
	template<typename _Ty>
	MATRICE_GLOBAL_INL static _Ty forward(const _Ty& x) {
		using value_t = conditional_t<is_scalar_v<_Ty>, 
			_Ty, typename _Ty::value_t>;
		return two<value_t>/(one<value_t> + exp(-2*x)) - one<value_t>;
	}
	/**
	 *\brief evaluating gradient(s) w.r.t. the input:
	  //tex:
	  //$$\sigma'(x) = 1-\sigma^2(x)$$
	 *\param [y] data produced by tanh::forward().
	 */
	template<typename _Ty>
	MATRICE_GLOBAL_INL static _Ty grad(const _Ty& y) {
		using value_t = conditional_t<is_scalar_v<_Ty>, 
			_Ty, typename _Ty::value_t>;
		return (one<value_t> - sqr(y));
	}
	/**
	 *\brief backward update gradient w.r.t. input params in-place.
	 *\param [g] input gradients, [y] tensor produced by forward().
	 */
	template<typename _Ty>
	MATRICE_GLOBAL_INL static _Ty& backward(_Ty& g, const _Ty& y) {
		const auto jac = tanh::grad(y);
		return (g);
	}
};

/* Rectified Linear Unit activation function */
struct relu {
	/**
	 *\brief for executing nodal activation(s):
	  //tex:
	  //$$\sigma(x) = \max(0,x)$$
	 *\param [x] input data with a type of scalar, ::Matrix or ::tensor.
	 */
	template<typename _Ty>
	MATRICE_GLOBAL_INL static _Ty forward(const _Ty& x) {
		using value_t = conditional_t<is_scalar_v<_Ty>,
			_Ty, typename _Ty::value_t>;
		return max(zero<value_t>, x);
	}
	/**
	 *\brief evaluating gradient(s) w.r.t. the input:
	  //tex:
	  //$$\sigma'(x) = \left\{
     //\begin{array}{lr}1,&x>0\\0,&x\leq0\end{array}\right.$$
	 *\param [y] data produced by relu::forward().
	 */
	template<typename _Ty>
	MATRICE_GLOBAL_INL static _Ty grad(const _Ty& y) {
		using value_t = conditional_t<is_scalar_v<_Ty>,
			_Ty, typename _Ty::value_t>;
		if constexpr (is_scalar_v<_Ty>) return y > 0 ? 1 : 0;
		else {
			_Ty _Ret = y;
			_Ret.each([](auto& _val) { _val = _val > 0 ? 1 : 0; });
			return std::forward<_Ty>(_Ret);
		}
	}
	/**
	 *\brief backward update gradient w.r.t. input params in-place.
	 *\param [g] input gradients, [y] tensor produced by forward().
	 */
	template<typename _Ty>
	MATRICE_GLOBAL_INL static _Ty& backward(_Ty& g, const _Ty& y) {
		const auto jac = relu::grad(y);
		return (g);
	}
};

};
}
DGE_MATRICE_END