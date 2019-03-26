/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#include "../../../core/matrix.h"

DGE_MATRICE_BEGIN 
namespace dnn {

struct function {
struct sigmoid {
	/**
	 *\brief execute sigmoid activation function:
	  //tex:
	  //$$\sigma(x) = \dfrac{1}{1+e^{-x}}$$
	 
	 *\param [x] input data with a type of scalar, dgelom::Matrix or dgelom::tensor
	 */
	template<typename _Ty> 
	MATRICE_GLOBAL_INL static _Ty forward(const _Ty& x) {
		using value_t = conditional_t<is_scalar_v<_Ty>, _Ty, typename _Ty::value_t>;
			return value_t(1) / (value_t(1) + exp(value_t(-1)*x));
	}
	/**
	 *\brief evaluate the gradients of sigmoid activation function:
	  //tex:
	  //$$\sigma'(x) = \sigma(x)(1-\sigma(x))$$

	 *\param [x] input data with a type of scalar, dgelom::Matrix or dgelom::tensor
	 */
	template<typename _Ty>
	MATRICE_GLOBAL_INL static _Ty backward(const _Ty& x) {
		using value_t = conditional_t<is_scalar_v<_Ty>, _Ty, typename _Ty::value_t>;
		const auto f = forward(x);
		return f * (value_(1) - f);
	}
};
};

}
DGE_MATRICE_END