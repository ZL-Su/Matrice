/**********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2022, Zhilong(Dgelom) Su, all rights reserved.

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
**********************************************************************/
#pragma once

#include "expr_base.hpp"
#include "forward.hpp"

namespace dgelom {
namespace xpr {

template<typename _Ty, size_t _Rows=0, size_t _Cols = _Rows>
class identity_matrix_exp : public __xpr__ {

public:
	static constexpr auto rows_at_compiletime = _Rows;
	static constexpr auto cols_at_compiletime = _Cols;
	using value_t = _Ty;
	using value_type = value_t;

	MATRICE_GLOBAL_FINL
	identity_matrix_exp(size_t rows = _Rows) noexcept {};

	MATRICE_GLOBAL_FINL 
	constexpr auto operator()(size_t) const noexcept {
		return value_type(1);
	}

	MATRICE_GLOBAL_FINL
	constexpr const auto rows() const noexcept {
		return rows_at_compiletime;
	}
	MATRICE_GLOBAL_FINL
	constexpr const auto cols() const noexcept {
		return cols_at_compiletime;
	}
	MATRICE_GLOBAL_FINL
	constexpr auto size()const noexcept {
		return rows() * cols();
	}

	MATRICE_GLOBAL_FINL
	auto eval() const noexcept {
		Matrix_<value_type,
			rows_at_compiletime,
			cols_at_compiletime> _Ret(rows(), cols());
		for (auto _Idx = 0; _Idx < _Ret.rows(); ++_Idx) {
			_Ret[_Idx][_Idx] = (*this)(_Idx);
		}
		return _Ret;
	}
};

/// <summary>
/// \brief Specialization for dynamic identity matrix.
/// </summary>
/// <typeparam name="_Ty"></typeparam>
template<typename _Ty>
class identity_matrix_exp<_Ty, 0, 0> : public __xpr__ {

public:
	static constexpr auto rows_at_compiletime = 0;
	static constexpr auto cols_at_compiletime = 0;
	using value_t = _Ty;
	using value_type = value_t;

	MATRICE_GLOBAL_FINL
	identity_matrix_exp(size_t rows) noexcept
		:_Mycols(rows), _Myrows(rows) {
	}

	MATRICE_GLOBAL_FINL
	constexpr auto operator()(size_t) const noexcept {
		return value_type(1);
	}

	MATRICE_GLOBAL_FINL
	const auto rows() const noexcept {
		return _Myrows;
	}
	MATRICE_GLOBAL_FINL
	const auto cols() const noexcept {
		return _Mycols;
	}
	MATRICE_GLOBAL_FINL
	auto size()const noexcept {
		return rows() * cols();
	}

	MATRICE_GLOBAL_FINL
	auto eval() const noexcept {
		Matrix_<value_type, 
			rows_at_compiletime, 
			cols_at_compiletime> _Ret(rows(), cols());
		for (auto _Idx = 0; _Idx < _Ret.rows(); ++_Idx) {
			_Ret[_Idx][_Idx] = (*this)(_Idx);
		}
		return _Ret;
	}

private:
	size_t _Myrows{0}, _Mycols{0};
};

}
}