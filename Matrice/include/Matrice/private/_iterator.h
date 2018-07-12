
#pragma once
#include <type_traits>
#include <iterator>
#include <vector>
#include "../util/_macros.h"

MATRICE_NAMESPACE_BEGIN_
template<typename _InIt> MATRICE_GLOBAL_FINL
_InIt _End(const _InIt _Begin, size_t _Size, size_t _Stride = 1) {
	return (_Begin + (_Size)*_Stride);
}
//Iterator for forward range [_My_ptr, _My_end), which is compatible with STD::ITERATOR
template<typename _Ty, typename = std::enable_if_t<std::is_arithmetic_v<_Ty>>>
class iterator_base
{
public:
	using iterator_category = std::random_access_iterator_tag;
	using value_type = _Ty;
	using pointer = std::add_pointer_t<value_type>;
	using reference = std::add_lvalue_reference_t<typename std::pointer_traits<pointer>::element_type>;
	using difference_type = std::ptrdiff_t;
	MATRICE_GLOBAL_FINL iterator_base(pointer _Ptr) noexcept
		:_My_ptr(_Ptr), _My_size(0), _My_step(0) {}
	MATRICE_GLOBAL_FINL iterator_base(pointer _Ptr, size_t _Size, size_t _Step =  1) noexcept
		:_My_ptr(_Ptr), _My_size(_Size), _My_step(_Step) {}
	
	MATRICE_GLOBAL_FINL reference operator*() const { 
		return ((reference)(*_My_ptr));
	}
	MATRICE_GLOBAL_FINL pointer operator->() const { 
		return (std::pointer_traits<pointer>::pointer_to(**this)); 
	}
	MATRICE_GLOBAL_FINL iterator_base operator++() { //preincrement
		_My_ptr += _My_step;
		return (*this);
	}
	MATRICE_GLOBAL_FINL iterator_base operator++(int) { //postincrement
		auto _Tmp = *this;
		*this += _My_step;
		return (_Tmp);
	}
	MATRICE_GLOBAL_FINL iterator_base operator--() { //preincrement
		_My_ptr -= _My_step;
		return (*this);
	}
	MATRICE_GLOBAL_FINL iterator_base operator--(int) { //postincrement
		auto _Tmp = *this;
		*this += _My_step;
		return (_Tmp);
	}
	MATRICE_GLOBAL_FINL iterator_base& operator+=(difference_type _Offset) {
		_Offset *= _My_step;
#if _ITERATOR_DEBUG_LEVEL == 2
		if (_Offset != 0) {
			if (_My_ptr + _Offset < _My_ptr || _My_end < _My_ptr + _Offset) {
				_DEBUG_ERROR("iterator + offset out of range");
			}
		}
#endif
		_My_ptr += _Offset;
		return (*this);
	}
	MATRICE_GLOBAL_FINL iterator_base operator+(difference_type _Offset) const {
		auto _Tmp = *this;
		return (_Tmp += _Offset*_My_step);
	}
	MATRICE_GLOBAL_FINL iterator_base operator-=(difference_type _Offset) {
		return (*this += -(_Offset * _My_step));
	}
	MATRICE_GLOBAL_FINL iterator_base operator-(difference_type _Offset) const {
		auto _Tmp = *this;
		return (_Tmp -= (_Offset * _My_step));
	}
	MATRICE_GLOBAL_FINL iterator_base operator-(const iterator_base& _Right) const {
		return (_My_ptr - _Right._My_ptr);
	}
	MATRICE_GLOBAL_FINL reference operator[](difference_type _Offset) const {
		return (*(*this + _Offset * _My_step));
	}
	MATRICE_GLOBAL_FINL bool operator==(const iterator_base& _Right) const {
		return (_My_ptr == _Right._My_ptr);
	}
	MATRICE_GLOBAL_FINL bool operator!=(const iterator_base& _Right) const {
		return (!(*this == _Right));
	}
	MATRICE_GLOBAL_FINL bool operator<(const iterator_base& _Right) const {
		return (_My_ptr < _Right._My_ptr);
	}
	MATRICE_GLOBAL_FINL bool operator>(const iterator_base& _Right) const {
		return (_Right < *this);
	}
	MATRICE_GLOBAL_FINL bool operator<=(const iterator_base& _Right) const {
		return (!(_Right < *this));
	}
	MATRICE_GLOBAL_FINL bool operator>=(const iterator_base& _Right) const {
		return (!(_Right > *this));
	}

	//test for iterator end condition
	MATRICE_GLOBAL_FINL operator bool() { return (_My_ptr != _My_end); }

	//return pointer to current object
	MATRICE_GLOBAL_FINL operator pointer() { return (_My_ptr); }

	//forward range iteration methods for [this->_My_ptr, this->_My_end)
	MATRICE_GLOBAL_FINL auto& begin() { return (*this); }
	MATRICE_GLOBAL_FINL auto end() { auto _Tmp = *this; return (_Tmp += _My_size); }
	MATRICE_GLOBAL_FINL const auto& begin() const { return (*this); }
	MATRICE_GLOBAL_FINL const auto end() const { auto _Tmp = *this; return (_Tmp += _My_size); }

protected:
	size_t _My_size;
	size_t _My_step;
	pointer _My_ptr;
	pointer _My_end = _End(_My_ptr, _My_size, _My_step);
	pointer _My_last = _My_end - _My_step;

};
template<typename _Ty>
MATRICE_GLOBAL_FINL iterator_base<_Ty> operator+ (typename iterator_base<_Ty>::difference_type _Offset, iterator_base<_Ty> _Next) {
	return (_Next += _Offset);
}

template<typename _Ty, typename = std::enable_if_t<std::is_arithmetic_v<_Ty>>>
class _Matrix_forward_iterator : public iterator_base<_Ty>
{
	using base_type = iterator_base<_Ty>;
public:
	template<typename... _Args> 
	MATRICE_GLOBAL_FINL _Matrix_forward_iterator(_Args... _args)
		: base_type(_args...) {}
};
_MATRICE_NAMESPACE_END
