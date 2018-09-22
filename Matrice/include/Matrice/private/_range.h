#pragma once
#include "../core/vector.h"

MATRICE_PRIVATE_BEGIN

// \class range base : [begin, end)
template<typename _Ty, 
	typename  = std::enable_if<std::is_arithmetic_v<_Ty>>>
class _Range_base {
public:
	using value_t = _Ty;
	using const_value_t = std::add_const_t<value_t>;
	_Range_base(const_value_t& _Begin, const_value_t& _End) noexcept
		: m_begin(_Begin), m_end(_End) {}
	_Range_base(const_value_t& _Begin, const_value_t& _End, const_value_t& _Stride) noexcept
		: m_begin(_Begin), m_end(_End), m_stride(_Stride) {}
	_Range_base(const _Range_base& _other) noexcept
		:m_begin(_other.m_begin), m_end(_other.m_end), m_stride(_other.m_stride) {}

	_Range_base& operator= (const _Range_base& _other) {
		m_begin = _other.m_begin, m_end = _other.m_end;
		m_stride = _other.m_stride;
		return (*this);
	}
	_Range_base& operator= (const std::initializer_list<value_t> _L) {
		m_begin = *_L.begin(), m_end = *(_L.begin()+1);
		if (_L.size() == 2) m_stride = value_t(1);
		if (_L.size() == 3) m_stride = *(_L.begin() + 2);
		return (*this);
	}

	value_t operator[](index_t i) { return (m_pos = m_begin + i * m_stride); }
	const_value_t operator[](index_t i) const { return (m_pos = m_begin + i * m_stride); }
	value_t operator()(index_t i) { return operator[](i); }
	const_value_t operator()(index_t i) const { return operator[](i); }

	value_t& value() { return m_pos; }
	const_value_t& value() const { return m_pos; }
	std::size_t size() const { return (m_end - m_begin) / m_stride; }

	MATRICE_GLOBAL_FINL operator bool() const { return (m_pos + m_stride < m_end); }

private:
	value_t m_begin, m_end, m_stride = value_t(1);
	mutable value_t m_pos = m_begin;
};

MATRICE_PRIVATE_END

MATRICE_NAMESPACE_BEGIN_

// \TEMPLATE CLASS range : [begin, end[, stride])
template<typename _Ty>
class range MATRICE_NONHERITABLE : public privt::_Range_base<_Ty> 
{
	using base_type = privt::_Range_base<_Ty>;
public:
	template<typename... _Args>
	range(const _Args&... _args) : base_type(_args...) {}
};
_MATRICE_NAMESPACE_END
