/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#include "../../../core/matrix.h"

DGE_MATRICE_BEGIN 
namespace dnn {

///<function> implementations for DNNs </function>
struct function {

template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)> 
constexpr static _Ty zero = static_cast<_Ty>(0);
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)> 
constexpr static _Ty one = static_cast<_Ty>(1);
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
constexpr static _Ty two = static_cast<_Ty>(2);

/* Sigmoid activation function */
struct sigmoid {
	/**
	 *\brief for executing nodal activation(s):
	  //tex:
	  //$$\sigma(x) = \dfrac{1}{1+e^{-x}}$$
	 *\param [x] input data with a type of scalar, ::Matrix or ::tensor
	 */
	template<typename _Ty> 
	MATRICE_GLOBAL_INL static _Ty forward(const _Ty& x) {
		using value_t = conditional_t<is_scalar_v<_Ty>, _Ty, 
			typename _Ty::value_t>;
		return one<value_t>/(one<value_t> + exp(value_t(-1)*x));
	}
	/**
	 *\brief for evaluating gradient(s) w.r.t. the input:
	  //tex:
	  //$$\sigma'(x) = \sigma(x)(1-\sigma(x))$$
	 *\param [x] input data with a type of scalar, ::Matrix or ::tensor
	 */
	template<typename _Ty>
	MATRICE_GLOBAL_INL static _Ty backward(const _Ty& x) {
		using value_t = conditional_t<is_scalar_v<_Ty>, _Ty, 
			typename _Ty::value_t>;
		const auto y = forward(x);
		return y * (one<value_t> - y);
	}
};
};

/* Tanh activation function */
struct tanh {
	/**
	 *\brief for executing nodal activation(s):
	  //tex:
	  //$$\sigma(x) = \dfrac{2}{1+e^{-2x}}-1$$
	 *\param [x] input data with a type of scalar, ::Matrix or ::tensor
	 */
	template<typename _Ty>
	MATRICE_GLOBAL_INL static _Ty forward(const _Ty& x) {
		using value_t = conditional_t<is_scalar_v<_Ty>, _Ty,
			typename _Ty::value_t>;
		return two<value_t> / (one<value_t> + exp(value_t(-2)*x));
	}
	/**
	 *\brief evaluating gradient(s) w.r.t. the input:
	  //tex:
	  //$$\sigma'(x) = 1-\sigma^2(x)$$
	 *\param [x] input data with a type of scalar, ::Matrix or ::tensor
	 */
	template<typename _Ty>
	MATRICE_GLOBAL_INL static _Ty backward(const _Ty& x) {
		using value_t = conditional_t<is_scalar_v<_Ty>, _Ty,
			typename _Ty::value_t>;
		const auto y = forward(x);
		return (one<value_t> - sqr(y));
	}
};

}
DGE_MATRICE_END