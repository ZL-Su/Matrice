/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include "../_matrix_base.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
namespace common = types;

using shape3 = tuple<size_t, size_t, size_t>;
using shape4 = tuple<size_t, size_t, size_t, size_t>;

namespace tensor_impl {
struct _Tensor_shape {
#define _TSHAPE_BATCHS std::get<0>(data)
#define _TSHAPE_CHANNELS std::get<0>(std::get<1>(data))
#define _TSHAPE_SHAPE(I) std::get<I>(std::get<1>(std::get<1>(data)))
	using _Myt = _Tensor_shape;
	using type = tuple<size_t, tuple<size_t, shape>>;
	MATRICE_HOST_INL _Tensor_shape(const type& _Shape)
		:data(_Shape) {}
	MATRICE_HOST_INL _Tensor_shape(type&& _Shape)
		: data(std::move(_Shape)) {}
	MATRICE_HOST_INL _Tensor_shape(const shape & _Shape)
		:data{1,{1,_Shape}} {}
	MATRICE_HOST_INL _Tensor_shape(const shape3& _Shape)
		:data{1,{std::get<0>(_Shape), {std::get<1>(_Shape),std::get<2>(_Shape)}}} {}
	MATRICE_HOST_INL _Tensor_shape(const shape4& _Shape)
		:data{ std::get<0>(_Shape), {std::get<1>(_Shape), {std::get<2>(_Shape),std::get<3>(_Shape)}}}{}
	template<typename _Ity>
	MATRICE_HOST_INL _Tensor_shape(std::initializer_list<_Ity> _Shape){
		if (_Shape.size() == 2) {
			data = {1,{1,{(size_t)*_Shape.begin(), (size_t)*(_Shape.begin()+1)}}};
		}
		if (_Shape.size() == 3) {
			data = { 1,{(size_t)*_Shape.begin(),{(size_t)*(_Shape.begin() + 1), (size_t)*(_Shape.begin() + 2)}} };
		}
		if (_Shape.size() == 4) {
			data = { (size_t)*_Shape.begin(),{(size_t)*(_Shape.begin() + 1),{(size_t)*(_Shape.begin() + 2), (size_t)*(_Shape.begin() + 3)}} };
		}
	}
	MATRICE_HOST_INL auto& operator= (const _Myt& _Oth) {
		data = (_Oth.data); return (*this);
	}
	MATRICE_HOST_INL auto& operator= (_Myt&& _Oth) {
		data = std::move(_Oth.data); return (*this);
	}
	MATRICE_HOST_INL constexpr auto operator()() const {
		return shape4{_TSHAPE_BATCHS, _TSHAPE_CHANNELS,
			_TSHAPE_SHAPE(0), _TSHAPE_SHAPE(1) };
	}
	MATRICE_HOST_INL constexpr auto rows() const {
		return _TSHAPE_BATCHS * _TSHAPE_SHAPE(0);
	}
	MATRICE_HOST_INL constexpr auto cols() const {
		return _TSHAPE_CHANNELS * _TSHAPE_SHAPE(1);
	}

	type data;
#undef _TSHAPE_BATCHS
#undef _TSHAPE_CHANNELS
#undef _TSHAPE_SHAPE
};
}

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
	using tensor_shape = tensor_impl::_Tensor_shape;

	/**
	 *\brief Empty constructor
	 */
	_Tensor()
		: _Mybase(),_Myshape({0, 0, 0, 0}) {}
	/**
	 *\brief Construct a tensor with shape [1,[1,_Shape]]
	 */
	_Tensor(const shape& _Shape)
		: _Mybase(_Shape), _Myshape(_Shape) {}
	/**
	 *\brief Construct a tensor with shape [1,[1,_Shape]] and fill with _Val
	 */
	_Tensor(const shape& _Shape, value_t _Val)
		: _Mybase(_Shape, _Val), _Myshape(_Shape) {}
	/**
	 *\brief Construct a tensor with shape [1,[_Shape]]
	 */
	_Tensor(const shape3& _Shape)
		: _Mybase(std::get<2>(_Shape), std::get<1>(_Shape)*std::get<0>(_Shape)), _Myshape(_Shape) {}
	/**
	 *\brief Construct a tensor with shape [1,[_Shape]] and fill with _Val
	 */
	_Tensor(const shape3& _Shape, value_t _Val)
		: _Tensor(_Shape) { _Mybase::operator= ({ _Val }); }
	/**
	 *\brief Construct a tensor with shape _Shape
	 */
	_Tensor(const shape4& _Shape) 
		: _Mybase(std::get<0>(_Shape)*std::get<2>(_Shape), std::get<1>(_Shape)*std::get<3>(_Shape)), _Myshape(_Shape) {}
	/**
	 *\brief Construct a tensor with shape _Shape and fill with _Val
	 */
	_Tensor(const shape4& _Shape, value_t _Val)
		: _Tensor(_Shape) {_Mybase::operator= ({_Val });}
	/**
	 *\brief Move constructor
	 */
	_Tensor(_Myt&& _Other) noexcept
		: _Mybase(std::move(_Other)), _Myshape(_Other.shape().data) {}

	/**
	 *\brief Create a tensor
	 *\param [_Shape] any shape type compatible with tensor_shape 
	 */
	MATRICE_HOST_INL auto& create(const tensor_shape& _Shape) {
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
	static MATRICE_HOST_INL auto rand(const tensor_shape& _Shape) {
		auto _Ret = _Mybase::rand(_Shape.rows(), _Shape.cols());
		_Ret.shape() = std::move(_Shape);
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
