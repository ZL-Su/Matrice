/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

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

MATRICE_NAMESPACE_BEGIN(dnn)
namespace detail {

template<typename _Ity> class _Input_layer {
	static_assert(sizeof(_Ity) != sizeof(_Ity),
		"Unsupported type in _Input_layer<_Ity>, \
		only supports dgelom::Matrix_, \
		dgelom::Tensor or array of dgelom::Matrix_.");
};

template<typename _Ty, int _N>
class _Input_layer<Vec_<_Ty, _N>> {
	using _Myt = _Input_layer;
public:
	using value_type = _Ty;
	using input_t = Vec_<value_type, _N>;
	using category = typename _Layer_tag::input;

	_Input_layer(const _Myt& _Other) noexcept
		: _Mydata(_Other._Mydata) {
	}
	_Input_layer(const input_t& _Input) noexcept
		: _Mydata(_Input) {
	}

	/**
	 *\brief Get rows of the input feature matrix.
	 */
	MATRICE_HOST_INL constexpr auto rows() const noexcept {
		return _Mydata.rows();
	}
	/**
	 *\brief Get rows of the input feature matrix.
	 */
	MATRICE_HOST_INL constexpr auto cols() const noexcept {
		return _Mydata.cols();
	}

	MATRICE_HOST_INL decltype(auto) data() const noexcept {
		return (_Mydata);
	}
	MATRICE_HOST_INL decltype(auto) data() noexcept {
		return (_Mydata);
	}
private:
	input_t _Mydata;
};
template<typename _Ty, int _N>
struct _Layer_traits<_Input_layer<Vec_<_Ty, _N>>>
{
	static constexpr auto depth = 1;
	static constexpr auto insize = _N;
	static constexpr auto outsize = _N;
	static constexpr auto has_bias = false;
	using value_type = _Ty;
};

template<typename _Ty, int _M, int _N>
class _Input_layer<Matrix_<_Ty, _M, _N>> {
	using _Myt = _Input_layer;
public:
	using value_type = _Ty;
	using input_t = Matrix_<value_type, _M, _N>;
	using category = typename _Layer_tag::input;

	_Input_layer(const _Myt& _Other) noexcept
		: _Mydata(_Other._Mydata) {
	}
	_Input_layer(const input_t& _Input) noexcept 
		: _Mydata(_Input) {
	}

	/**
	 *\brief Get rows of the input feature matrix. 
	 */
	MATRICE_HOST_INL constexpr auto rows() const noexcept {
		return _Mydata.rows();
	}
	/**
	 *\brief Get rows of the input feature matrix.
	 */
	MATRICE_HOST_INL constexpr auto cols() const noexcept {
		return _Mydata.cols();
	}

	MATRICE_HOST_INL decltype(auto) data() const noexcept {
		return (_Mydata);
	}
	MATRICE_HOST_INL decltype(auto) data() noexcept {
		return (_Mydata);
	}
private:
	input_t _Mydata;
};
template<typename _Ty, int _N>
struct _Layer_traits<_Input_layer<Matrix_<_Ty, 1, _N>>>
{
	static constexpr auto depth = 1;
	static constexpr auto insize = _N;
	static constexpr auto outsize = _N;
	static constexpr auto has_bias = false;
	using value_type = _Ty;
};


template<typename _Ty, int _M, int _N, size_t _D>
class _Input_layer<std::array<Matrix_<_Ty, _M, _N>, _D>> {
	using _Myt = _Input_layer;
public:
	using input_t = std::array<Matrix_<_Ty, _M, _N>, _D>;
	using category = typename _Layer_tag::input;

	_Input_layer() noexcept {}
	_Input_layer(const _Myt&) noexcept {}

	template<typename _Inty>
	_Input_layer(const _Inty&) noexcept {}

private:
	input_t _Mydata;
};
} // namespace detail

template<typename _Input>
using input_layer = detail::_Input_layer<_Input>;

MATRICE_NAMESPACE_END(dnn)