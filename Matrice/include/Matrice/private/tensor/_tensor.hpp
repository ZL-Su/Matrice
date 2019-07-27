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
#include "../_plain_base.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ty>
struct index<_Ty, tag::_Tensor_tag> {
	static_assert(is_integral_v<_Ty>, "_Ty must be an integer type");

	MATRICE_GLOBAL_INL index(_Ty d, _Ty h, _Ty w) noexcept {
		data[0] = d, data[1] = h, data[2] = w;
	}
	MATRICE_GLOBAL_INL _Ty value() noexcept {
		return 0;
	}

	std::array<_Ty, 3> data;
};

/**
 *\brief CLASS TEMPLATE dgelom::tensor prototype
 *\param <_Ty> data type
 *\param <_Depth> tensor depth which grows vertically
 */
template<typename _Ty, size_t _Depth>
class _Tensor 
	: public types::Base_<_Tensor<_Ty, _Depth>, tensor_traits<_Tensor<_Ty, _Depth>>>
{
	using _Myt = _Tensor;
	using _Mytraits = tensor_traits<_Myt>;
	using _Mybase = types::Base_<_Myt, _Mytraits>;
public:
	enum { Size = 0, CompileTimeRows = 0, CompileTimeCols = 0 };
	static constexpr auto depth = _Mytraits::depth;
	using typename _Mybase::value_type;
	//using typename _Mybase::scalar_type;
	using _Mybase::operator=;
	using _Mybase::operator();

	/**
	 *\brief default constructor
	 */
	_Tensor() noexcept 
		: _Mybase() {
	}
	/**
	 *\brief constructor
	 *\param [h, w] rows and cols of each tensor cell
	 */
	_Tensor(size_t h, size_t w) noexcept
		:_Mybase(basic_shape_t{ depth, h, w }) {
	}
	/**
	 *\brief constructor with a value initialization
	 *\param [h, w] rows and cols of each tensor cell
	 *\param [_Val] initial value 
	 */
	_Tensor(size_t h, size_t w, value_type _Val) noexcept
		:_Tensor(h, w) {
		_Mybase::operator=(_Val);
	}
	/**
	 *\brief copy constructor
	 *\param [oth] an other tensor
	 */
	_Tensor(const _Myt& oth) noexcept
		:_Mybase(oth) {
	}
	/**
	 *\brief move constructor
	 *\param [oth] an other tensor
	 */
	_Tensor(_Myt&& oth) noexcept
		:_Mybase(move(oth)) {
	}
	/**
	 *\brief template constructor
	 *\param [args...] argument(s) with any supported type(s) 
	 */
	template<typename _Arg>
	_Tensor(_Arg&& arg) noexcept
		:_Mybase(forward<_Arg>(arg)) {
	}

	/**
	 *\brief ewise accessor
	 *\param [d, r, c] the indices in depth, height and width axis respectively 
	 */
	MATRICE_HOST_INL const value_type& operator()(size_t d, size_t r, size_t c) const noexcept {
#if defined _DEBUG || MATRICE_DEBUG
		DGELOM_CHECK(d < depth, "depth index over range.");
#endif
		const auto inner_size = m_shape.hw();
		return (_Mybase::operator[](d*inner_size+r*m_width)[c]);
	}
	MATRICE_HOST_INL value_type& operator()(size_t d, size_t r, size_t c) noexcept {
#if defined _DEBUG || MATRICE_DEBUG
		DGELOM_CHECK(d < _Depth, "depth index over range.");
#endif
		const auto inner_size = m_shape.hw();
		return (_Mybase::operator[](d*inner_size + r * m_width)[c]);
	}

	/**
	 *\brief tensor slice operation along the depth axis, return value is a matrix view.
	 *\param [d] - input index of depth
	 */
	MATRICE_HOST_INL decltype(auto) array(size_t d) noexcept {
#if defined _DEBUG || MATRICE_DEBUG
		DGELOM_CHECK(d < _Depth, "depth index over range.");
#endif
		const auto w = m_shape.w();
		const auto h = m_shape.h();
		const auto r0 = d*h, r1 = (d+1) * h;
		return _Mybase::block(0, w, r0, r1);
	}
	/**
	 *\brief returns a matrix view of the full tensor.
	 */
	MATRICE_HOST_INL decltype(auto) array() noexcept {
		return _Mybase::block(0, this->m_cols, 0, this->m_rows);
	}

	/**
	 *\brief matrix multiplication is disabled for a tensor.
	 */
	template<typename _Rhs>
	MATRICE_GLOBAL_INL auto mul(const _Rhs& _Right) const = delete;

	/**
	 *\brief internal used method within CRTP. 
	 */
	MATRICE_HOST_INL void __create_impl(size_t h, size_t w) {
		this->_Reset({ depth, h, w });
		m_width = m_shape.w();
		m_height = m_shape.h();
	}

private:
	using _Mybase::m_shape;
	size_t m_width = m_shape.w();
	size_t m_height = m_shape.h();
};

/**
 *\brief CLASS TEMPLATE dgelom::tensor prototype with dynamic depth
 *\param <_Ty> data type
 */
template<typename _Ty>
class _Tensor<_Ty, 0> : public types::Base_<_Tensor<_Ty, 0>>
{
	using _Myt = _Tensor;
	using _Mybase = types::Base_<_Myt>;
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
		m_shape = move(_Shape); _Mybase::_Flush_view_buf();
	}
	/**
	 *\brief Construct a tensor with shape [1,[1,_Shape]] and fill with _Val
	 */
	_Tensor(shape_t<size_t>&& _Shape, value_t _Val) noexcept
		: _Tensor(move(_Shape)) {
		_Mybase::operator= ({ _Val });
	}
	/**
	 *\brief Construct a tensor with shape [1,[_Shape]]
	 */
	_Tensor(shape3_t<size_t>&& _Shape) noexcept
		: _Mybase(get<2>(_Shape), get<1>(_Shape)*get<0>(_Shape)) {
		m_shape = move(_Shape); _Mybase::_Flush_view_buf();
	}
	/**
	 *\brief Construct a tensor with shape [1,[_Shape]] and fill with _Val
	 */
	_Tensor(shape3_t<size_t>&& _Shape, value_t _Val) noexcept
		: _Tensor(move(_Shape)) {
		_Mybase::operator= ({ _Val });
	}
	/**
	 *\brief Construct a tensor with shape _Shape
	 */
	_Tensor(shape4_t<size_t>&& _Shape) noexcept
		: _Mybase(get<0>(_Shape)*get<2>(_Shape), get<1>(_Shape)*get<3>(_Shape)) {
		m_shape = move(_Shape); _Mybase::_Flush_view_buf();
	}
	/**
	 *\brief Construct a tensor with shape _Shape and fill with _Val
	 */
	_Tensor(shape4_t<size_t>&& _Shape, value_t _Val) noexcept
		: _Tensor(move(_Shape)) {
		_Mybase::operator= ({ _Val });
	}
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
		m_shape = move(_Shape);
		create(m_shape.rows(), m_shape.cols());
		return (*this);
	}
	/**
	 *\brief Get element at _Idx
	 *\param [_Idx] input linear index
	 */
	MATRICE_HOST_INL auto& operator()(size_t _Idx) {
		auto[n, c, h, w] = m_shape.parse(_Idx);

	}
	MATRICE_HOST_INL const auto& operator()(size_t _Idx) const {
		auto[n, c, h, w] = m_shape.parse(_Idx);
	}
	/**
	 *\brief Get a sub-tensor at this[n,:,:,:]
	 *\param [n] the index of the first dim
	 */
	MATRICE_HOST_INL auto at(size_t n) const {
		_Myt _Ret({ 1,dims().get(1),dims().get(2),dims().get(3) });
		auto _Begin = this->begin() + n * _Ret.size();
		std::copy(_Begin, _Begin + _Ret.size(), _Ret.begin());
		return forward<_Myt>(_Ret);
	}
	/**
	 *\brief Get a sub-tensor at this[n,c,:,:]
	 *\param [n,c] the indices of the first and sencond dim
	 */
	MATRICE_HOST_INL auto at(size_t n, size_t c) const {
		_Myt _Ret({ 1, 1, m_shape.get(2), m_shape.get(3) });
		auto _Begin = this->begin() + n * _Ret.size()*m_shape.get(1) + c * _Ret.size();
		std::copy(_Begin, _Begin + _Ret.size(), _Ret.begin());
		return forward<_Myt>(_Ret);
	}

	/**
	 *\brief Get the view at this[k,d,:,:]
	 *\param [k, d] the indices of the first and sencond dim
	 */
	MATRICE_HOST_INL auto view(size_t k, size_t d) const {
		auto _Off_y = k * dims().get(2), _Off_x = d * dims().get(3);
		return (this->block(_Off_x, _Off_x + dims().get(3), _Off_y, _Off_y + dims().get(2)));
	}

	/**
	 *\brief Generate a zero-value filled tensor
	 *\param [_Shape] any shape type compatible with tensor_shape
	 */
	static MATRICE_HOST_INL auto zero(tensor_shape&& _Shape) {
		return (move(_Myt({ _Shape.get(0),_Shape.get(1),_Shape.get(2),_Shape.get(3) }, 0)));
	}
	/**
	 *\brief Generate a random filled tensor
	 *\param [_Shape] any shape type compatible with tensor_shape
	 */
	static MATRICE_HOST_INL auto rand(tensor_shape&& _Shape) {
		auto _Ret = _Mybase::rand(_Shape.rows(), _Shape.cols());
		_Ret.m_shape = move(_Shape);
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
	using _Mybase::m_shape; //[extent, [depth, [height, width]]]
};
#undef MATRICE_EXPAND_SHAPE


template<typename _Ty, size_t _Depth>
struct tensor_traits<_Tensor<_Ty, _Depth>> {
	using type = _Ty;
	using category = tag::_Tensor_tag;
	static constexpr auto depth = _Depth;
	static constexpr auto _M = 0, _N = 0;
	static constexpr bool Is_base = std::false_type::value;
};
_DETAIL_END
DGE_MATRICE_END