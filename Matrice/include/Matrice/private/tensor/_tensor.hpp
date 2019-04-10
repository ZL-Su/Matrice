/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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

#include <tuple>
#include <valarray>
#include "../_range.h"
#include "_tensor_exp.hpp"

DGE_MATRICE_BEGIN
namespace types {
template<typename _Ty, int _M, int _N> class Matrix_;
}

namespace detail {

template<typename _Ty, int _M = 0, int _N = 0, size_t _K = 0,
	typename matrix_type = types::Matrix_<_Ty, _M, _N>>
class _Tensor_impl MATRICE_NONHERITABLE : public std::valarray<matrix_type> {
	using _Myt = _Tensor_impl;
	using _Mybase = std::valarray<matrix_type>;
	using _My_element_traits = matrix_traits<matrix_type>;
	using _Mytp = _Tensor_impl<_Ty, (_M > 0 ? _N:_M), (_N > 0 ? _M:_N)> ;
public:
	enum{CompileTimeRows = 0, CompileTimeCols = 0};
	using element_type = matrix_type;
	using value_type = typename _My_element_traits::type;
	using value_t = value_type;

	_Tensor_impl() 
		: _Mybase(), m_rows(0), m_cols(0) {}
	_Tensor_impl(size_t _Rows)
		: _Mybase(m_size = _Rows), m_rows(_Rows), m_cols(1) {}
	_Tensor_impl(size_t _Rows, const matrix_type& _Mat)
		: _Mybase(_Mat, m_size = _Rows), m_rows(_Rows), m_cols(1) {}
	_Tensor_impl(size_t _Rows, size_t _Cols) 
		: _Mybase(m_size = _Rows*_Cols), m_rows(_Rows), m_cols(_Cols) {}
	_Tensor_impl(size_t _Rows, size_t _Cols, const matrix_type& _Mat) : _Mybase(_Mat, m_size=(_Rows * _Cols)), m_rows(_Rows), m_cols(_Cols) {}
	_Tensor_impl(const _Myt& _Other)
		: _Mybase(static_cast<_Mybase>(_Other)), m_rows(_Other.m_rows), m_cols(_Other.m_cols), m_size(m_rows*m_cols) {}
	_Tensor_impl(_Myt&& _Other)
		: _Mybase(move(static_cast<_Mybase>(_Other))), m_rows(_Other.m_rows), m_cols(_Other.m_cols), m_size(m_rows*m_cols) {}
	_Tensor_impl(const _Mybase& _Other)
		: _Mybase(_Other), m_rows(m_size=(_Other.size())), m_cols(1) {}

	MATRICE_HOST_INL auto& create(size_t _Rows, size_t _Cols) {
		m_rows = _Rows, m_cols = _Cols;
		_Mybase::resize(m_size = m_rows * m_cols);
		m_data = &(*this)[0];
		return (*this);
	}

	MATRICE_HOST_INL auto& operator()(size_t _R, size_t _C) {
		return m_data[_R*m_cols + _C];
	}
	MATRICE_HOST_INL const auto& operator()(size_t _R, size_t _C) const {
		return m_data[_R*m_cols + _C];
	}
	MATRICE_HOST_INL auto& operator()(size_t _Idx) {
		return m_data[_Idx];
	}
	MATRICE_HOST_INL const auto& operator()(size_t _Idx) const {
		return m_data[_Idx];
	}
	MATRICE_HOST_INL auto& operator= (const _Mybase& _Other) {
		_Mybase::operator= (_Other);
		m_data = &(*this)[0];
		return (*this);
	}
	MATRICE_HOST_INL auto& operator= (_Mybase&& _Other) {
		_Mybase::operator= (move(_Other));
		m_data = &(*this)[0];
		return (*this);
	}
	MATRICE_HOST_INL auto& operator= (const _Myt& _Other) {
		m_rows = _Other.m_rows, m_cols = _Other.m_cols;
		m_size = _Other.m_size;
		_Mybase::operator= (static_cast<_Mybase>(_Other));
		m_data = &(*this)[0];
		return (*this);
	}
	MATRICE_HOST_INL auto& operator= (_Myt&& _Other) {
		m_rows = _Other.m_rows, m_cols = _Other.m_cols;
		m_size = _Other.m_size;
		_Mybase::operator= (move(static_cast<_Mybase>(_Other)));
		m_data = &(*this)[0];
		return (*this);
	}
	MATRICE_HOST_INL auto rows() const { return m_rows; }
	MATRICE_HOST_INL auto cols() const { return m_cols; }
	MATRICE_HOST_INL auto shape() const { return std::tie(m_rows, m_cols); }
	MATRICE_HOST_INL auto& reshape(size_t _Rows) {
		m_rows = _Rows, m_cols = m_size / _Rows;
		return (*this);
	}
	MATRICE_HOST_INL auto reduce() const {
		auto _Ret = m_data[0];
		if (m_size > 1) for (const auto _Idx : range(1, m_size))
				_Ret = _Ret + m_data[_Idx];
		return (_Ret);
	}
	MATRICE_HOST_INL auto t() const {
		_Mytp _Ret(m_cols, m_rows);
		for (const auto _Idx : range(0, m_size))
			_Ret[_Idx] = m_data[_Idx].t();
		return forward<decltype(_Ret)>(_Ret);
	}
	MATRICE_HOST_INL auto mul(const _Mytp& _Rhs) const {
		using _Op = detail::_Tensor_exp_op::_Ewise_mmul<element_type, typename _Mytp::element_type>;
		return detail::_Tensor_exp<_Myt, _Mytp, _Op>(*this, _Rhs);
	}

	/**
	 * \element-wise addition.
	 */
	friend MATRICE_HOST_INL _Myt operator+(const _Myt& _Left, const _Myt& _Right) {
		return ((_Myt)((_Mybase)(_Left)+(_Mybase)(_Right))).reshape(_Left.rows());
	}
	/**
	 * \element-wise subtraction.
	 */
	friend MATRICE_HOST_INL _Myt operator-(const _Myt& _Left, const _Myt& _Right) {
		return ((_Myt)((_Mybase)(_Left)-(_Mybase)(_Right))).reshape(_Left.rows());
	}
	/**
	 * \element-wise multiplication.
	 */
	friend MATRICE_HOST_INL _Myt operator*(const _Myt& _Left, const _Myt& _Right) {
		return ((_Myt)((_Mybase)(_Left)*(_Mybase)(_Right))).reshape(_Left.rows());
	}
	/**
	 * \element-wise division.
	 */
	friend MATRICE_HOST_INL _Myt operator/(const _Myt& _Left, const _Myt& _Right) {
		return ((_Myt)((_Mybase)(_Left)/(_Mybase)(_Right))).reshape(_Left.rows());
	}
	/**
	 * \multiplies each element with the _Left scalar.
	 */
	friend MATRICE_HOST_INL _Myt operator*(value_t _Left, const _Myt& _Right) {
		_Myt _Ret(_Right.rows(), _Right.cols());
		for (const auto& _Idx : range(0, _Ret.size()))
			_Ret(_Idx) = _Left*_Right(_Idx);
		return std::forward<_Myt>(_Ret);
	}
	/**
	 * \element-wise multiplies with the _Left matrix convertible type.
	 */
	template<typename _Lhs, typename = std::enable_if_t<is_matrix_convertible_v<_Lhs>>>
	friend MATRICE_HOST_INL _Myt operator*(const _Lhs& _Left, const _Myt& _Right) {
#ifdef _DEBUG
		if (_Left.shape() != _Right.shape()) throw
			std::runtime_error("_Left and _Right must have a uniform shape.");
#endif // _DEBUG

		_Myt _Ret(_Right.rows(), _Right.cols());
		for (const auto& _Idx : range(0, _Ret.size()))
			_Ret(_Idx) = _Left(_Idx)*_Right(_Idx);
		return std::forward<_Myt>(_Ret);
	}

private:
	std::size_t m_rows, m_cols, m_size;
	std::add_pointer_t<element_type> m_data = &(*this)[0];
};
template<typename _Ty, int _M, int _N, int _K>
struct is_tensor<_Tensor_impl<_Ty, _M, _N, _K>> : std::true_type {};
template<typename _Ty, int _M, int _N, int _K>
struct tensor_traits< _Tensor_impl<_Ty, _M, _N, _K>> {
	using value_type = _Ty;
	using element_type = Matrix_<value_type, _M, _N>;
};

template<typename _Ty, int _M = 0, int _N = _M>
class _Multi_matrix MATRICE_NONHERITABLE
	: public std::vector<types::Matrix_<_Ty, _M, _N>>
{
	using _Mybase = std::vector<types::Matrix_<_Ty, _M, _N>>;
	using _Myt = _Multi_matrix;
public:
	using matrix_type = typename _Mybase::value_type;
	using value_type = typename matrix_type::value_type;
	using value_t = value_type;

	MATRICE_HOST_FINL explicit _Multi_matrix(std::size_t _Count) 
		: _Mybase(_Count) {
	}
	MATRICE_HOST_FINL explicit _Multi_matrix(const matrix_type& _Mat, std::size_t _Count)
		: _Mybase(_Count, _Mat) {
	}
	MATRICE_HOST_FINL _Multi_matrix(const std::initializer_list<matrix_type>& _L)
		: _Mybase(_L.size()) {
		for (auto _Idx = 0; _Idx < this->size(); ++_Idx) {
			this->operator[](_Idx) = *(_L.begin() + _Idx);
		}
	}
	MATRICE_HOST_FINL _Myt& operator= (const std::initializer_list<matrix_type>& _L) {
		if (this->size() < _L.size()) this->resize(_L.size());
		for (auto _Idx = 0; _Idx < this->size(); ++_Idx) {
			this->operator[](_Idx) = *(_L.begin() + _Idx);
		}
	}
	/**
	 * \gets the _N block views of the multi-matrix after pos of _Off-set.
	 */
	template<std::size_t _N, std::size_t _Off = 0, typename _Ity = std::size_t> 
	MATRICE_HOST_FINL auto view_n (const std::tuple<_Ity, _Ity, _Ity, _Ity>& _R) const {
		return tuple_n<_N-1>::_(this->data()+ _Off, [&](const matrix_type& _Mat) {
			return _Mat.block(_R); 
		});
	}
	template<std::size_t _N, std::size_t _Off = 0, typename _Ity = std::size_t>
	MATRICE_HOST_FINL auto view_n(_Ity _L, _Ity _R, _Ity _U, _Ity _D) const {
		return tuple_n<_N - 1>::_(this->data()+ _Off, [&](const matrix_type& _Mat) {
			return _Mat.block(_L, _R, _U, _D);
		});
	}
	template<std::size_t _N, std::size_t _Off = 0, typename _Ity = std::size_t>
	MATRICE_HOST_FINL auto view_n(const range<_Ity>& _Rx, const range<_Ity>& _Ry) const {
		return tuple_n<_N - 1>::_(this->data() + _Off, [&](const matrix_type& _Mat) {
			return _Mat.block(_Rx.begin(), _Rx.end(), _Ry.begin(), _Ry.end());
		});
	}
};
}

DGE_MATRICE_END