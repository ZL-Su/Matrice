/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for
more detail.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once

#include "../_storage.hpp"

MATRICE_NAMESPACE_BEGIN_ _DETAIL_BEGIN

template<typename _Ty>
template<Location _Loc, size_t _Opt> MATRICE_GLOBAL_FINL
Storage_<_Ty>::DenseBase<_Loc, _Opt>::DenseBase() noexcept
 : my_rows(0), my_cols(0), my_size(0), my_data(nullptr) {}

template<typename _Ty>
template<Location _Loc, size_t _Opt> MATRICE_GLOBAL_FINL
Storage_<_Ty>::DenseBase<_Loc, _Opt>::DenseBase(int_t _rows, int_t _cols, pointer _data)
 : my_rows(_rows), my_cols(_cols), my_size(my_rows*my_cols), my_data(_data) {
	if constexpr (location == OnStack) my_owner = Owner;
	else my_owner = Proxy;
}

template<typename _Ty>
template<Location _Loc, size_t _Opt>
template<Location _From, size_t _Option> MATRICE_GLOBAL_FINL
Storage_<_Ty>::DenseBase<_Loc, _Opt>::DenseBase(const DenseBase<_From, _Option>& _other, pointer _data)
 : DenseBase(_other.rows(), _other.cols(), _data) {
	privt::unified_sync<value_t, _From, _Loc, _Option>::
		op(my_data, _other.data(), my_rows, my_cols, _other.pitch());
}

template<typename _Ty>
template<Location _Loc, size_t _Opt>
template<Location _From, size_t _Option> MATRICE_GLOBAL_FINL
Storage_<_Ty>::DenseBase<_Loc, _Opt>::DenseBase(const DenseBase<_From, _Option>& _other)
 : DenseBase(_other.rows(), _other.cols()) {
	privt::unified_sync<value_t, _From, _Loc, _Loc==OnDevice?_Opt: _Option>::
		op(my_data, _other.data(), my_rows, my_cols, _Loc == OnDevice ? my_pitch : _other.pitch());
}

template<typename _Ty> 
template<Location _Loc, size_t _Opt> MATRICE_GLOBAL_FINL
Storage_<_Ty>::DenseBase<_Loc, _Opt>::DenseBase(
	int_t _rows, int_t _cols, pointer _data, initlist<value_t> _list)
 : DenseBase(_rows, _cols, _data) {
	my_owner = (Owner);
	if (_list.size() == 1)
		privt::unified_fill<value_t, location + option>::op(my_data, *_list.begin(), my_rows, my_cols, my_pitch);
	else 
		std::memcpy((void*)my_data, (void*)&(*_list.begin()), _list.size() * type_bytes_v<value_t>);
}

template<typename _Ty>
template<Location _Loc, size_t _Opt> MATRICE_GLOBAL_FINL
Storage_<_Ty>::DenseBase<_Loc, _Opt>::DenseBase(initlist<value_t> _list)
	: DenseBase(_list.size(), 1) {
	my_owner = Ownership::Owner;
	if (_list.size() == 1)
		privt::unified_fill<value_t, location + option>::op(my_data, *_list.begin(), my_rows, my_cols, my_pitch);
	else
		std::memcpy((void*)my_data, (void*)&(*_list.begin()), _list.size() * type_bytes_v<value_t>);
}

template<typename _Ty>
template<Location _Loc, size_t _Opt> MATRICE_GLOBAL_FINL
Storage_<_Ty>::DenseBase<_Loc, _Opt>::~DenseBase() noexcept {
#if MATRICE_SHARED_STORAGE == 0
	if (my_data && my_owner == Owner){
		internal::_Memory::free<location>(my_data);
	}
#endif
}

template<typename _Ty>
template<Location _Loc, size_t _Opt> MATRICE_GLOBAL_FINL
decltype(auto) Storage_<_Ty>::DenseBase<_Loc, _Opt>::operator=(value_t _val) noexcept {
	privt::unified_fill<value_t, location + option>::
		op(my_data, _val, my_rows, my_cols, my_pitch);

	return (*this);
}

template<typename _Ty> 
template<Location _Loc, size_t _Opt> MATRICE_GLOBAL_FINL 
decltype(auto) Storage_<_Ty>::DenseBase<_Loc, _Opt>::operator=(initlist<value_t> _list) noexcept {
	if (_list.size() == 1) this->operator=(*_list.begin());
	else
		privt::unified_sync<value_t, OnStack, Location(location), LINEAR+COPY>::
			op(my_data, pointer(_list.begin()), my_rows, my_cols, my_pitch);

	return (*this);
}

template<typename _Ty>
template<Location _Loc, size_t _Opt>
template<int _M, int _N> MATRICE_GLOBAL_FINL
Storage_<_Ty>::DenseBase<_Loc, _Opt>::DenseBase(const Allocator<_M, _N, allocator_traits<_M, _N>::value>& _al) noexcept
	: DenseBase(_al.rows(), _al.cols()) {
	constexpr auto _From = Location(remove_all_t<decltype(_al)>::location);
	constexpr auto _Option = remove_all_t<decltype(_al)>::option;
	privt::unified_sync<value_t, _From, _Loc, _Loc==OnDevice?_Opt:_Option>::
		op(my_data, (pointer)_al.data(), my_rows, my_cols, _Loc == OnDevice ? my_pitch : _al.pitch());
}
template<typename _Ty>
template<Location _Loc, size_t _Opt> 
template<int _M, int _N> MATRICE_GLOBAL_FINL
decltype(auto) Storage_<_Ty>::DenseBase<_Loc, _Opt>::operator=(const Allocator<_M, _N, allocator_traits<_M, _N>::value>& _al) noexcept {
	constexpr auto _From = Location(remove_all_t<decltype(_al)>::location);
	constexpr auto _Option = remove_all_t<decltype(_al)>::option;

	privt::unified_sync<value_t, _From, _Loc, _Loc==OnDevice?_Opt:_Option>::
		op(my_data, (pointer)_al.data(), my_rows, my_cols, _Loc == OnDevice ? my_pitch : _al.pitch());
	return (*this);
}
_DETAIL_END
_MATRICE_NAMESPACE_END