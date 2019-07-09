/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include "../_matrix_base.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
namespace common = types;

template<typename _Ty>
class _Tensor<_Ty, 0> : public common::Base_<_Tensor<_Ty, 0>>
{
	using _Myt = _Tensor;
	using _Mybase = common::Base_<_Myt>;
	using _Mytraits = matrix_traits<_Myt>;
public:
	enum { Size = 0, CompileTimeRows = 0, CompileTimeCols = 0, };
	using typename _Mybase::value_t;
	using typename _Mybase::value_type;
	using typename _Mybase::const_initlist;
	using tensor_shape = basic_shape<size_t>;
	using _Mybase::dims;

	/**
	 *\brief Empty constructor
	 */
	_Tensor() noexcept
		: _Mybase() {}
	/**
	 *\brief Construct a tensor with shape [1,[1,_Shape]]
	 */
	_Tensor(shape_t<size_t>&& _Shape) noexcept
		: _Mybase(_Shape) {
		_Myshape = move(_Shape); _Mybase::_Flush_view_buf();
	}
	/**
	 *\brief Construct a tensor with shape [1,[1,_Shape]] and fill with _Val
	 */
	_Tensor(shape_t<size_t>&& _Shape, value_t _Val) noexcept
		: _Tensor(move(_Shape)) { _Mybase::operator= ({ _Val }); }
	/**
	 *\brief Construct a tensor with shape [1,[_Shape]]
	 */
	_Tensor(shape3_t<size_t>&& _Shape) noexcept
		: _Mybase(get<2>(_Shape), get<1>(_Shape)*get<0>(_Shape)) { 
		_Myshape=move(_Shape); _Mybase::_Flush_view_buf();
	}
	/**
	 *\brief Construct a tensor with shape [1,[_Shape]] and fill with _Val
	 */
	_Tensor(shape3_t<size_t>&& _Shape, value_t _Val) noexcept
		: _Tensor(move(_Shape)) { _Mybase::operator= ({ _Val }); }
	/**
	 *\brief Construct a tensor with shape _Shape
	 */
	_Tensor(shape4_t<size_t>&& _Shape) noexcept
		: _Mybase(get<0>(_Shape)*get<2>(_Shape), get<1>(_Shape)*get<3>(_Shape)) { 
		_Myshape=move(_Shape); _Mybase::_Flush_view_buf();
	}
	/**
	 *\brief Construct a tensor with shape _Shape and fill with _Val
	 */
	_Tensor(shape4_t<size_t>&& _Shape, value_t _Val) noexcept
		: _Tensor(move(_Shape)) {_Mybase::operator= ({_Val });}
	/**
	 *\brief Move constructor
	 */
	_Tensor(_Myt&& _Other) noexcept
		: _Mybase(move(_Other)) {}

	/**
	 *\brief Create a tensor
	 *\param [_Shape] any shape type compatible with tensor_shape 
	 */
	MATRICE_HOST_INL auto& create(tensor_shape&& _Shape) {
		_Myshape = move(_Shape);
		create(_Myshape.rows(), _Myshape.cols());
		return (*this);
	}
	/**
	 *\brief Get element at _Idx
	 *\param [_Idx] input linear index
	 */
	MATRICE_HOST_INL auto& operator()(size_t _Idx) {
		auto[n, c, h, w] = _Myshape.parse(_Idx);

	}
	MATRICE_HOST_INL const auto& operator()(size_t _Idx) const {
		auto[n, c, h, w] = _Myshape.parse(_Idx);
	}
	/**
	 *\brief Get a sub-tensor at this[n,:,:,:]
	 *\param [n] the index of the first dim
	 */
	MATRICE_HOST_INL auto at(size_t n) const {
		_Myt _Ret({ 1,dims().get(1),dims().get(2),dims().get(3) });
		auto _Begin = this->begin()+n*_Ret.size();
		std::copy(_Begin, _Begin + _Ret.size(), _Ret.begin());
		return forward<_Myt>(_Ret);
	}
	/**
	 *\brief Get a sub-tensor at this[n,c,:,:]
	 *\param [n,c] the indices of the first and sencond dim
	 */
	MATRICE_HOST_INL auto at(size_t n, size_t c) const {
		_Myt _Ret({ 1, 1, _Myshape.get(2), _Myshape.get(3) });
		auto _Begin = this->begin()+n*_Ret.size()*_Myshape.get(1)+c*_Ret.size();
		std::copy(_Begin, _Begin + _Ret.size(), _Ret.begin());
		return forward<_Myt>(_Ret);
	}

	/**
	 *\brief Get the view at this[n,c,:,:]
	 *\param [n, c] the indices of the first and sencond dim
	 */
	MATRICE_HOST_INL auto view(size_t n, size_t c) const {
		auto _Off_y = n * dims().get(2), _Off_x = c * dims().get(3);
		return (this->block(_Off_x, _Off_x+ dims().get(3), _Off_y, _Off_y+ dims().get(2)));
	}

	/**
	 *\brief Generate a zero-value filled tensor
	 *\param [_Shape] any shape type compatible with tensor_shape
	 */
	static MATRICE_HOST_INL auto zero(tensor_shape&& _Shape) {
		return (move(_Myt({_Shape.get(0),_Shape.get(1),_Shape.get(2),_Shape.get(3)}, 0)));
	}
	/**
	 *\brief Generate a random filled tensor
	 *\param [_Shape] any shape type compatible with tensor_shape
	 */
	static MATRICE_HOST_INL auto rand(tensor_shape&& _Shape) {
		auto _Ret = _Mybase::rand(_Shape.rows(), _Shape.cols());
		_Ret._Myshape = move(_Shape);
		return (move(_Ret));
	}

	/**
	 *\brief Internal used function
	 *\param [_1, _2] dummy
	 */
	MATRICE_HOST_INL void create(size_t _1, size_t _2) {
		_Mybase::m_storage.create(_1, _2);
		_Mybase::_Flush_view_buf();
	}

	/**
	 *\brief Deleted methods
	 */
	template<ttag _Ltag = ttag::N, ttag _Rtag = ttag::N, typename _Rhs = _Myt>
	MATRICE_GLOBAL_FINL auto inplace_mul(const _Rhs& _Right) = delete;
private:
	using _Mybase::_Myshape; //[batchs, [channels, [height, width]]]
};
#undef MATRICE_EXPAND_SHAPE
_DETAIL_END
DGE_MATRICE_END