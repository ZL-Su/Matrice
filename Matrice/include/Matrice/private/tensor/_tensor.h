/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include "../_matrix_base.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
namespace common = types;

template<typename _Ty>
class _Tensor : public common::Base_<_Tensor<_Ty>>
{
	using _Myt = _Tensor;
	using _Mybase = common::Base_<_Myt>;
	using _Mytraits = matrix_traits<_Myt>;
public:
	using typename _Mybase::value_t;
	using typename _Mybase::value_type;
	using typename _Mybase::const_init_list;
	using tensor_shape = basic_shape<size_t>;

	/**
	 *\brief Empty constructor
	 */
	_Tensor()
		: _Mybase(),_Myshape({0, 0, 0, 0}) {}
	/**
	 *\brief Construct a tensor with shape [1,[1,_Shape]]
	 */
	_Tensor(shape_t<size_t>&& _Shape)
		: _Mybase(_Shape), _Myshape(std::move(_Shape)) {}
	/**
	 *\brief Construct a tensor with shape [1,[1,_Shape]] and fill with _Val
	 */
	_Tensor(shape_t<size_t>&& _Shape, value_t _Val)
		: _Mybase(_Shape, _Val), _Myshape(std::move(_Shape)) {}
	/**
	 *\brief Construct a tensor with shape [1,[_Shape]]
	 */
	_Tensor(shape3_t<size_t>&& _Shape)
		: _Mybase(std::get<2>(_Shape), std::get<1>(_Shape)*std::get<0>(_Shape)), _Myshape(std::move(_Shape)) {}
	/**
	 *\brief Construct a tensor with shape [1,[_Shape]] and fill with _Val
	 */
	_Tensor(shape3_t<size_t>&& _Shape, value_t _Val)
		: _Tensor(std::move(_Shape)) { _Mybase::operator= ({ _Val }); }
	/**
	 *\brief Construct a tensor with shape _Shape
	 */
	_Tensor(shape4_t<size_t>&& _Shape) 
		: _Mybase(std::get<0>(_Shape)*std::get<2>(_Shape), std::get<1>(_Shape)*std::get<3>(_Shape)), _Myshape(std::move(_Shape)) {}
	/**
	 *\brief Construct a tensor with shape _Shape and fill with _Val
	 */
	_Tensor(shape4_t<size_t>&& _Shape, value_t _Val)
		: _Tensor(std::move(_Shape)) {_Mybase::operator= ({_Val });}
	/**
	 *\brief Move constructor
	 */
	_Tensor(_Myt&& _Other) noexcept
		: _Mybase(std::move(_Other)), _Myshape(std::move(_Other.shape()())) {}

	/**
	 *\brief Create a tensor
	 *\param [_Shape] any shape type compatible with tensor_shape 
	 */
	MATRICE_HOST_INL auto& create(tensor_shape&& _Shape) {
		_Myshape = std::move(_Shape);
		create(_Myshape.rows(), _Myshape.cols());
		return (*this);
	}

	/**
	 *\brief Get and set tensor shape
	 */
	MATRICE_HOST_INL auto& shape() {
		return (_Myshape);
	}
	MATRICE_HOST_INL const auto& shape() const {
		return (_Myshape);
	}
	/**
	 *\brief Generate a random filled tensor
	 *\param [_Shape] any shape type compatible with tensor_shape
	 */
	static MATRICE_HOST_INL auto rand(tensor_shape&& _Shape) {
		auto _Ret = _Mybase::rand(_Shape.rows(), _Shape.cols());
		_Ret._Myshape = std::move(_Shape);
		return (std::move(_Ret));
	}

	/**
	 *\brief Internal used function
	 *\param [_1, _2] dummy
	 */
	MATRICE_HOST_INL void create(size_t _1, size_t _2) {
		_Mybase::operator= (std::move(_Mybase(_1, _2)));
	}
private:
	tensor_shape _Myshape; //[batchs, [channels, [height, width]]]
};

template<typename _Ty>
struct matrix_traits<_Tensor<_Ty>> {
	using type = _Ty;
	static constexpr auto _M = 0, _N = 0;
	static constexpr auto Is_base = std::true_type::value;
	struct size {
		struct rows { static constexpr auto value = _M; };
		struct cols { static constexpr auto value = _N; };
		static constexpr auto rows_v = rows::value;
		static constexpr auto cols_v = cols::value;
	};
};
#undef MATRICE_EXPAND_SHAPE
_DETAIL_END
DGE_MATRICE_END
