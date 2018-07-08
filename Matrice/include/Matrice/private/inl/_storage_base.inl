#pragma once
#include "../_storage.hpp"

MATRICE_NAMESPACE_BEGIN_ namespace details {

template<typename _Ty>
template<Location _Loc, size_t _Opt>
MATRICE_GLOBAL_FINL Storage_<_Ty>::DenseBase<_Loc, _Opt>::DenseBase()
 : my_rows(0), my_cols(0), my_size(0), my_data(0), my_owner(Dummy) 
{
}

template<typename _Ty>
template<Location _Loc, size_t _Opt>
MATRICE_GLOBAL_FINL Storage_<_Ty>::DenseBase<_Loc, _Opt>::DenseBase(int_t _rows, int_t _cols, pointer _data)
 : my_rows(_rows), my_cols(_cols), my_size(my_rows*my_cols), my_data(_data), my_owner(location == OnStack ? Owner : Proxy) 
{
}

template<typename _Ty>
template<Location _Loc, size_t _Opt>
template<Location _From, size_t _Option>
MATRICE_GLOBAL_FINL Storage_<_Ty>::DenseBase<_Loc, _Opt>::DenseBase(const DenseBase<_From, _Option>& _other, pointer _data)
 : DenseBase<_Loc, _Opt>(_other.rows(), _other.cols(), _data)
{
	privt::unified_sync<value_t, _From, _Loc, _Option>::op(my_data, _other.data(), my_rows, my_cols, _other.pitch());
}

template<typename _Ty>
template<Location _Loc, size_t _Opt>
template<Location _From, size_t _Option>
MATRICE_GLOBAL_FINL Storage_<_Ty>::DenseBase<_Loc, _Opt>::DenseBase(const DenseBase<_From, _Option>& _other)
 : DenseBase<_Loc, _Opt>(_other.rows(), _other.cols())
{
	privt::unified_sync<value_t, _From, _Loc, _Loc == OnDevice ? option : _Option>::op(my_data, _other.data(), my_rows, my_cols, _Loc == OnDevice ? my_pitch : _other.pitch());
}

template<typename _Ty> 
template<Location _Loc, size_t _Opt>
MATRICE_GLOBAL_FINL Storage_<_Ty>::DenseBase<_Loc, _Opt>::DenseBase(
	int_t _rows, int_t _cols, pointer _data,
	std::initializer_list<value_t> _list)
 : my_rows(_rows), my_cols(_cols), my_data(_data), 
	my_size(my_rows*my_cols), my_owner(Owner)
{
	if (_list.size() == 1)
		privt::unified_fill<value_t, location + option>::op(my_data, *_list.begin(), my_rows, my_cols, my_pitch);
	else 
		std::memcpy((void*)my_data, (void*)&(*_list.begin()), _list.size() * type_bytes<value_t>::value);
}

template<typename _Ty>
template<Location _Loc, size_t _Opt>
MATRICE_GLOBAL_FINL Storage_<_Ty>::DenseBase<_Loc, _Opt>::~DenseBase()
{
#ifndef __CXX11_SHARED__
	if (my_data && my_owner = Owner){
		privt::unified_free<value_t, location> _Tidy;
		_Tidy(my_data);
	}
#endif
}

template<typename _Ty> 
template<Location _Loc, size_t _Opt> MATRICE_GLOBAL_FINL 
auto& Storage_<_Ty>::DenseBase<_Loc, _Opt>::operator=(std::initializer_list<value_t> _list)
{
	if (_list.size() == 1)
		privt::unified_fill<value_t, location + option>::op(my_data, *_list.begin(), my_rows, my_cols, my_pitch);
	else
		privt::unified_sync<value_t, OnStack, Location(location), LINEAR+COPY>::op(my_data, pointer(_list.begin()), my_rows, my_cols, my_pitch);

	return (*this);
}

} _MATRICE_NAMESPACE_END