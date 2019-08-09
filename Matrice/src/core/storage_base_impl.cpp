/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
#include <complex>
#include "../../include/Matrice/private/_storage.hpp"
#include "../../include/Matrice/private/_type_traits.h"
#include "../../include/Matrice/private/_unified_memory.h"
#include "../../include/Matrice/private/_decl_dev_funcs.h"
#include "../../include/Matrice/util/_exception.h"

#pragma warning(disable: 4715 4661 4224 4267 4244 4819 4199)

#define MATRICE_DENSEBASE_SIG \
template<typename _Ty> \
template<Location _Loc, size_t _Opt> \
Storage_<_Ty>::DenseBase<_Loc, _Opt>

DGE_MATRICE_BEGIN _DETAIL_BEGIN

MATRICE_DENSEBASE_SIG::DenseBase(int_t _Rows, int_t _Cols)
	: my_rows(_Rows), my_cols(_Cols), my_size(_Rows*_Cols), my_owner(Ownership::Owner) {
#if MATRICE_SHARED_STORAGE == 1
	if constexpr (location == Location::OnStack) {
		// donot use shared pointer for managed memory
	}
	else if constexpr (location == Location::OnHeap) {
		my_shared = SharedPtr(privt::aligned_malloc<value_t>(my_size), 
			[=](pointer ptr) { privt::aligned_free(ptr); });
		my_data = my_shared.get();
	}
#ifdef MATRICE_ENABLE_CUDA
	else if constexpr (location == Location::OnGlobal) {
		my_data = privt::global_malloc<value_t>(size_t(my_size));
		my_shared = SharedPtr(my_data, 
			[=](pointer ptr) { privt::device_free(ptr); });
	}
	else if constexpr (location == Location::OnDevice) {
		size_t w = my_cols, h = my_rows;
		if (w == 1) std::swap(w, h);
		my_pitch = w;
		my_data = privt::device_malloc<value_t>(my_pitch, h);
		my_shared = SharedPtr(my_data, 
			[=](pointer ptr) { privt::device_free(ptr); });
		my_pitch = my_pitch == w ? 1 : my_pitch;
	}
#endif
	else DGELOM_ERROR("Oops, unknown device is used.");
#else
	if constexpr (location == Location::OnStack) {
		// donot do anything
	}
	else if constexpr (location == Location::OnHeap) {
		my_data = privt::aligned_malloc<value_t>(my_size);
	} 
#ifdef MATRICE_ENABLE_CUDA
	else if constexpr (location == Location::OnGlobal) {
		my_data = privt::global_malloc<value_t>(size_t(my_size));
	}
	else if constexpr (location == Location::OnDevice) {
		size_t w = my_cols, h = my_rows;
		if (w == 1) std::swap(w, h);
		my_pitch = w;
		my_data = privt::device_malloc<value_t>(my_pitch, h);
		my_pitch = my_pitch == w ? 1 : my_pitch;
	}
#endif
	else DGELOM_ERROR("Oops, unknown device is used.");
#endif
}

MATRICE_DENSEBASE_SIG::DenseBase(int_t rows, int_t cols, const value_t val)
	: DenseBase(rows, cols) {
	if constexpr (_Loc != OnDevice) {
		for (int_t i = 0; i < my_size; ++i) {
			my_data[i] = val;
		}
	}
}	

MATRICE_DENSEBASE_SIG::DenseBase(const DenseBase& _other)
#if MATRICE_SHARED_STORAGE == 1
	: my_rows(_other.my_rows), my_cols(_other.my_cols),
	my_size(_other.my_size), my_pitch(_other.my_pitch) {
	if (this != std::addressof(_other)) {
		if constexpr (_Loc == Location::OnStack) {
			privt::fill_mem(_other.my_data, my_data, my_size);
			my_owner = _other.my_owner;
		}
		else my_shared = _other.my_shared, my_data = _other.my_data, my_owner = Refer;
	}
}
#else
{
	if (this != std::addressof(_other)) {
		this->create(_other.rows(), _other.cols());
		privt::unified_sync<value_t, _other.location, location, option>::op(my_data, _other.my_data, my_rows, my_cols, my_pitch);
	}
}
#endif

MATRICE_DENSEBASE_SIG::DenseBase(DenseBase&& _other) noexcept
	: my_rows(_other.my_rows), my_cols(_other.my_cols), my_size(_other.my_size), my_owner(_other.my_owner), my_pitch(_other.my_pitch), my_data(_other.my_data)
#if MATRICE_SHARED_STORAGE == 1
	 , my_shared(std::move(_other.my_shared))
{
	_other.my_shared = nullptr;
#else
{
#endif
	_other.my_rows = 0, _other.my_cols = 0, _other.my_size = 0;
	_other.my_data = 0, _other.my_owner = Empty;
}

MATRICE_DENSEBASE_SIG&
Storage_<_Ty>::DenseBase<_Loc,_Opt>::create(int_t _Rows, int_t _Cols){
	my_rows = _Rows, my_cols = _Cols;
	my_size = _Rows*_Cols, my_owner = Owner;
#if MATRICE_SHARED_STORAGE == 1
		if constexpr (location == Location::OnStack) {
			// donot use shared pointer for managed memory
		}
		else if constexpr (location == Location::OnHeap) {
			my_shared = SharedPtr(privt::aligned_malloc<value_t>(my_size),
				[=](pointer ptr) { privt::aligned_free(ptr); });
			my_data = my_shared.get();
		}
#ifdef MATRICE_ENABLE_CUDA
		else if constexpr (location == Location::OnGlobal) {
			my_data = privt::global_malloc<value_t>(size_t(my_size));
			my_shared = SharedPtr(my_data,
				[=](pointer ptr) { privt::device_free(ptr); });
		}
		else if constexpr (location == Location::OnDevice) {
			size_t w = my_cols, h = my_rows;
			if (w == 1) std::swap(w, h);
			my_pitch = w;
			my_data = privt::device_malloc<value_t>(my_pitch, h);
			my_shared = SharedPtr(my_data,
				[=](pointer ptr) { privt::device_free(ptr); });
			my_pitch = my_pitch == w ? 1 : my_pitch;
		}
#endif
		else DGELOM_ERROR("Oops, unknown device is used.");
#else
		if constexpr (location == Location::OnStack) {
			// donot do anything
		}
		else if constexpr (location == Location::OnHeap) {
			my_data = privt::aligned_malloc<value_t>(my_size);
		}
#ifdef MATRICE_ENABLE_CUDA
		else if constexpr (location == Location::OnGlobal) {
			my_data = privt::global_malloc<value_t>(size_t(my_size));
		}
		else if constexpr (location == Location::OnDevice) {
			size_t w = my_cols, h = my_rows;
			if (w == 1) std::swap(w, h);
			my_pitch = w;
			my_data = privt::device_malloc<value_t>(my_pitch, h);
			my_pitch = my_pitch == w ? 1 : my_pitch;
		}
#endif
		else DGELOM_ERROR("Oops, unknown device is used.");
#endif
		return (*this);
}

MATRICE_DENSEBASE_SIG&
Storage_<_Ty>::DenseBase<_Loc, _Opt>::operator=(const DenseBase& _other){
	if (this != std::addressof(_other)) {
		my_rows = _other.my_rows, my_cols = _other.my_cols;
		my_size = _other.my_size, my_pitch = _other.my_pitch;
#if MATRICE_SHARED_STORAGE == 1
		if constexpr (_Loc == OnStack) {
			privt::fill_mem(_other.my_data, my_data, my_size);
			my_owner = _other.my_owner;
		}
		else my_shared = _other.my_shared, my_data = _other.my_data, my_owner = Refer;
#else
		if (my_owner == Empty) this->create(my_rows, my_cols);
		privt::unified_sync<value_t, _other.location, location, option>::op(my_data, _other.my_data, my_rows, my_cols, my_pitch);
#endif
	}
	return (*this);
}

MATRICE_DENSEBASE_SIG& 
Storage_<_Ty>::DenseBase<_Loc, _Opt>::operator=(DenseBase&& _other) noexcept{
	if (this != std::addressof(_other)) {
		my_owner = _other.my_owner, my_size = _other.my_size;
		my_rows = _other.my_rows, my_cols = _other.my_cols;
		my_data = _other.my_data;
#if MATRICE_SHARED_STORAGE == 1
		my_shared = std::move(_other.my_shared);
		_other.my_shared = nullptr;
#endif
		_other.my_rows = 0, _other.my_cols = 0, _other.my_size = 0;
		_other.my_data = 0, _other.my_owner = Owner;
	}
	return (*this);
}

///<brief> Explicit Specialization </brief>
template class Storage_<int>::DenseBase<OnStack, LINEAR + COPY>;
template class Storage_<int>::DenseBase<OnHeap, LINEAR + COPY>;
template class Storage_<int>::DenseBase<UnSpecified>;
template class Storage_<char>::DenseBase<OnStack, LINEAR + COPY>;
template class Storage_<char>::DenseBase<OnHeap, LINEAR + COPY>;
template class Storage_<char>::DenseBase<UnSpecified>;
template class Storage_<bool>::DenseBase<OnStack, LINEAR + COPY>;
template class Storage_<bool>::DenseBase<OnHeap, LINEAR + COPY>;
template class Storage_<bool>::DenseBase<UnSpecified>;
template class Storage_<size_t>::DenseBase<OnStack, LINEAR + COPY>;
template class Storage_<size_t>::DenseBase<OnHeap, LINEAR + COPY>;
template class Storage_<size_t>::DenseBase<UnSpecified>;
template class Storage_<float>::DenseBase<OnStack, LINEAR + COPY>;
template class Storage_<float>::DenseBase<OnHeap, LINEAR + COPY>;
template class Storage_<float>::DenseBase<UnSpecified>;
template class Storage_<double>::DenseBase<OnStack, LINEAR + COPY>;
template class Storage_<double>::DenseBase<OnHeap, LINEAR + COPY>;
template class Storage_<double>::DenseBase<UnSpecified>;
template class Storage_<unsigned char>::DenseBase<OnStack, LINEAR + COPY>;
template class Storage_<unsigned char>::DenseBase<OnHeap, LINEAR + COPY>;
template class Storage_<unsigned char>::DenseBase<UnSpecified>;
template class Storage_<int>::DenseBase<OnStack, LINEAR+SHARED>;
template class Storage_<int>::DenseBase<OnHeap, LINEAR+SHARED>;
template class Storage_<char>::DenseBase<OnStack, LINEAR+SHARED>;
template class Storage_<char>::DenseBase<OnHeap, LINEAR+SHARED>;
template class Storage_<bool>::DenseBase<OnStack, LINEAR+SHARED>;
template class Storage_<bool>::DenseBase<OnHeap, LINEAR+SHARED>;
template class Storage_<size_t>::DenseBase<OnStack, LINEAR+SHARED>;
template class Storage_<size_t>::DenseBase<OnHeap, LINEAR+SHARED>;
template class Storage_<float>::DenseBase<OnStack, LINEAR+SHARED>;
template class Storage_<float>::DenseBase<OnHeap, LINEAR+SHARED>;
template class Storage_<double>::DenseBase<OnStack, LINEAR+SHARED>;
template class Storage_<double>::DenseBase<OnHeap, LINEAR+SHARED>;
template class Storage_<unsigned char>::DenseBase<OnStack, LINEAR+SHARED>;
template class Storage_<unsigned char>::DenseBase<OnHeap, LINEAR+SHARED>;
#ifdef MATRICE_ENABLE_CUDA
template class Storage_<int>::DenseBase<OnGlobal, LINEAR + COPY>;
template class Storage_<char>::DenseBase<OnGlobal, LINEAR + COPY>;
template class Storage_<size_t>::DenseBase<OnGlobal, LINEAR + COPY>;
template class Storage_<bool>::DenseBase<OnGlobal, LINEAR + COPY>;
template class Storage_<float>::DenseBase<OnGlobal, LINEAR + COPY>;
template class Storage_<double>::DenseBase<OnGlobal, LINEAR + COPY>;
template class Storage_<unsigned char>::DenseBase<OnGlobal, LINEAR + COPY>;
template class Storage_<int>::DenseBase<OnGlobal, LINEAR>;
template class Storage_<char>::DenseBase<OnGlobal, LINEAR>;
template class Storage_<bool>::DenseBase<OnGlobal, LINEAR>;
template class Storage_<size_t>::DenseBase<OnGlobal, LINEAR>;
template class Storage_<float>::DenseBase<OnGlobal, LINEAR>;
template class Storage_<double>::DenseBase<OnGlobal, LINEAR>;
template class Storage_<unsigned char>::DenseBase<OnGlobal, LINEAR>;
template class Storage_<int>::DenseBase<OnDevice, LINEAR>;
template class Storage_<char>::DenseBase<OnDevice, LINEAR>;
template class Storage_<bool>::DenseBase<OnDevice, LINEAR>;
template class Storage_<float>::DenseBase<OnDevice, LINEAR>;
template class Storage_<double>::DenseBase<OnDevice, LINEAR>;
template class Storage_<unsigned char>::DenseBase<OnDevice, LINEAR>;
template class Storage_<size_t>::DenseBase<OnDevice, LINEAR>;
template class Storage_<int>::DenseBase<OnDevice, PITCHED>;
template class Storage_<char>::DenseBase<OnDevice, PITCHED>;
template class Storage_<bool>::DenseBase<OnDevice, PITCHED>;
template class Storage_<float>::DenseBase<OnDevice, PITCHED>;
template class Storage_<double>::DenseBase<OnDevice, PITCHED>;
template class Storage_<unsigned char>::DenseBase<OnDevice, PITCHED>;
template class Storage_<size_t>::DenseBase<OnDevice, PITCHED>;
#endif
_DETAIL_END DGE_MATRICE_END
#undef MATRICE_DENSEBASE_SIG