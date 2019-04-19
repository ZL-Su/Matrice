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

MATRICE_ALGS_BEGIN

template<typename _Derived> MATRICE_HOST_INL
auto _Interpolation_base<_Derived>::_Value_at(const point_type& _Pos) const {
	const auto _Ix = floor<int>(_Pos.x), _Iy = floor<int>(_Pos.y);
	const auto _Dx = _Pos.x - _Ix, _Dy = _Pos.y - _Iy;

	const auto _Diff_x_n = _Mydt_this->_Val_dx_n(_Dx);
	const auto _Diff_y_n = _Mydt_this->_Val_dy_n(_Dy);

	constexpr auto Ldv = decltype(_Diff_x_n)::CompileTimeRows;
	constexpr auto _L = ~-(Ldv >> 1), _R = -~(Ldv >> 1);
	const auto _Coeff = _Mycoeff(_Ix-_L, _Ix+_R, _Iy-_L, _Iy+_R).eval<Ldv, Ldv>();

	const auto& _Kov = _Mydt_this->_Kernel_of_value();
	remove_all_t<decltype(_Kov)> _KtCK;

	auto _Temp = _Kov.t().mul(_Coeff.mul(_Kov).eval()).eval();

	return (_Diff_y_n.mul(_Temp.mul(_Diff_x_n).eval()))(0);
}
template<typename _Derived> MATRICE_HOST_INL
auto _Interpolation_base<_Derived>::_Gradx_at(const point_type& _Pos) const {
	const auto _Ix = floor<int>(_Pos.x), _Iy = floor<int>(_Pos.y);
	const auto _Dx = _Pos.x - _Ix, _Dy = _Pos.y - _Iy;

	const auto _Diff_x_n = _Mydt_this->_Grad_dx_n(_Dx);
	const auto _Diff_y_n = _Mydt_this->_Val_dy_n(_Dy);

	constexpr auto Ldv = decltype(_Diff_y_n)::CompileTimeCols;
	constexpr auto _L = ~-(Ldv >> 1), _R = -~(Ldv >> 1);
	const auto _Coeff = _Mycoeff(_Ix-_L, _Ix+_R, _Iy-_L, _Iy+_R).eval<Ldv, Ldv>();

	const auto& _Kov = _Mydt_this->_Kernel_of_value();
	const auto& _Kog = _Mydt_this->_Kernel_of_grad();
	auto _Temp = _Kov.t().mul(_Coeff.mul(_Kog).eval()).eval();

	return (_Diff_y_n.mul(_Temp.mul(_Diff_x_n).eval()))(0);
}
template<typename _Derived> MATRICE_HOST_INL
auto _Interpolation_base<_Derived>::_Grady_at(const point_type& _Pos) const {
	const auto _Ix = floor<int>(_Pos.x), _Iy = floor<int>(_Pos.y);
	const auto _Dx = _Pos.x - _Ix, _Dy = _Pos.y - _Iy;

	const auto _Diff_x_n = _Mydt_this->_Val_dx_n(_Dx);
	const auto _Diff_y_n = _Mydt_this->_Grad_dy_n(_Dy);

	constexpr auto Ldv = decltype(_Diff_x_n)::CompileTimeRows;
	constexpr auto _L = ~-(Ldv >> 1), _R = -~(Ldv >> 1);
	const auto _Coeff = _Mycoeff(_Ix-_L, _Ix+_R, _Iy-_L, _Iy+_R).eval<Ldv, Ldv>();

	const auto& _Kov = _Mydt_this->_Kernel_of_value();
	const auto& _Kog = _Mydt_this->_Kernel_of_grad();
	auto _Temp = _Kog.t().mul(_Coeff.mul(_Kov).eval()).eval();

	return (_Diff_y_n.mul(_Temp.mul(_Diff_x_n).eval()))(0);
}
MATRICE_ALGS_END