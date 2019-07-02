#pragma once
#include "../core/vector.h"

DGE_MATRICE_BEGIN namespace detail {

// \class range base : [begin, end)
template<typename _Ty, 
	typename _Iy = conditional_t<is_iterator_v<_Ty>, int64_t, _Ty>>
class _Range_base {
	using _Myt = _Range_base;
public:
	using value_t = _Ty;
	using stride_t = _Iy;
	using const_value_t = std::add_const_t<value_t>;
	using const_stride_t = std::add_const_t<stride_t>;
	using reference = std::add_lvalue_reference_t<value_t>;

	class _My_iterator {
	public:
		MATRICE_GLOBAL_INL _My_iterator(const_value_t& _Pos, const_stride_t&_Inc = 1) noexcept
			: _My_pos(_Pos), _My_inc(_Inc) {}

		MATRICE_GLOBAL_INL reference operator++() noexcept {
			return _My_pos += _My_inc; 
		}
		MATRICE_GLOBAL_INL const reference operator++() const noexcept {
			return _My_pos += _My_inc; 
		}

		MATRICE_GLOBAL_INL reference operator*() noexcept { return _My_pos; }
		MATRICE_GLOBAL_INL const reference operator*() const noexcept { return _My_pos; }

		MATRICE_GLOBAL_INL auto operator != (const _My_iterator& _Other)const noexcept {
			return _My_pos < (_Other._My_pos);
		}

		MATRICE_GLOBAL_INL operator value_t() const noexcept { return _My_pos; }

	private:
		value_t _My_pos;
		stride_t _My_inc;
	};

	_Range_base(const_value_t& _Begin, const_value_t& _End, const_stride_t& _Stride = 1) noexcept
		: m_begin(_Begin), m_end(_End), m_stride(_Stride) {}
	_Range_base(const _Range_base& _other) noexcept
		:m_begin(_other.m_begin), m_end(_other.m_end), m_stride(_other.m_stride) {}

	MATRICE_GLOBAL_INL _Range_base& operator= (const _Range_base& _other) noexcept {
		m_begin = _other.m_begin, m_end = _other.m_end;
		m_stride = _other.m_stride;
		return (*this);
	}
	MATRICE_GLOBAL_INL _Range_base& operator= (const initlist<value_t> _L) noexcept {
		m_begin = *_L.begin(), m_end = *(_L.begin()+1);
		if (_L.size() == 2) m_stride = value_t(1);
		if (_L.size() == 3) m_stride = *(_L.begin() + 2);
		return (*this);
	}

	MATRICE_GLOBAL_INL value_t operator[](index_t i) noexcept {
		return (m_pos = m_begin + i * m_stride); 
	}
	MATRICE_GLOBAL_INL const_value_t operator[](index_t i) const noexcept {
		return (m_pos = m_begin + i * m_stride); 
	}
	MATRICE_GLOBAL_INL value_t operator()(index_t i) noexcept {
		return operator[](i); 
	}
	MATRICE_GLOBAL_INL const_value_t operator()(index_t i) const noexcept {
		return operator[](i); 
	}

	MATRICE_GLOBAL_INL reference value() noexcept {
		return m_pos; 
	}
	MATRICE_GLOBAL_INL const reference value() const noexcept {
		return m_pos; 
	}
	MATRICE_GLOBAL_INL size_t size() const noexcept { 
		return static_cast<size_t>((m_end - m_begin) / m_stride); 
	}

	MATRICE_GLOBAL_INL reference front() noexcept {
		return m_begin; 
	}
	MATRICE_GLOBAL_INL const reference front() const noexcept {
		return m_begin; 
	}

	MATRICE_GLOBAL_INL _My_iterator begin() { 
		return _My_iterator(m_pos, m_stride); 
	}
	MATRICE_GLOBAL_INL _My_iterator end() { 
		return _My_iterator(m_end, zero<stride_t>);
	}
	MATRICE_GLOBAL_INL const _My_iterator begin() const {
		return _My_iterator(m_pos, m_stride);
	}
	MATRICE_GLOBAL_INL const _My_iterator end() const {
		return _My_iterator(m_end, zero<stride_t>);
	}

	/**
	 *\brief set the starting position.
	 */
	MATRICE_GLOBAL_INL _Myt& begin(value_t _Beg) noexcept {
		m_begin = _Beg; return (*this);
	}
	/**
	 *\brief set the last position.
	 */
	MATRICE_GLOBAL_INL _Myt& end(value_t _End) noexcept {
		m_end = _End; return (*this);
	}
	/**
	 *\brief set the stride of a range.
	 */
	MATRICE_GLOBAL_INL _Myt& stride(stride_t _Std) noexcept {
		m_stride = _Std; return (*this);
	}
	/**
	 *\brief check if the pos of the range meets the end. 
	 */
	MATRICE_GLOBAL_INL operator bool() const noexcept {
		return (m_pos + m_stride < m_end); 
	}

private:
	value_t m_begin, m_end;
	stride_t m_stride = one<stride_t>;
	mutable value_t m_pos = m_begin;
};

template<typename _Ty> 
class _Rect_impl MATRICE_NONHERITABLE {
public:
	using value_type = _Ty;
	using point_type = Vec_<value_type, 2>;

	_Rect_impl() noexcept {}
	template<typename _U1, typename _U2>
	_Rect_impl(const Vec_<_U1, 2>& _X, const _U2& _W, const _U2& _H)
		:_Mybegin(_X.x, _X.y), _Mywidth(_W), _Myheight(_H) { _My_end(); }
	template<typename _U1, typename _U2>
	_Rect_impl(const _U1& _X, const _U1& _Y, const _U2& _W, const _U2& _H)
		:_Mybegin(_X, _Y), _Mywidth(_W), _Myheight(_H) { _My_end(); }

	MATRICE_HOST_INL auto& operator()(const point_type& _X) {
		_Mybegin = _X; _My_end(); return (*this);
	}

	MATRICE_HOST_INL point_type& begin() {
		return (_Mybegin); 
	}
	MATRICE_HOST_INL const point_type& begin() const {
		return (_Mybegin); 
	}
	MATRICE_HOST_INL point_type& end() {
		return (_Myend); 
	}
	MATRICE_HOST_INL const point_type& end() const {
		return (_Myend); 
	}

	MATRICE_HOST_INL void set(const point_type& p, value_type w, value_type h) noexcept {
		_Mybegin = p;
		_Mywidth = w, _Myheight = h;
		_My_end();
	}

	MATRICE_HOST_INL void set_size(value_type w, value_type h) noexcept {
		_Mywidth = w, _Myheight = h;
		_My_end();
	}

private:
	MATRICE_HOST_INL auto _My_end() {
		_Myend.x = _Mybegin.x + _Mywidth;
		_Myend.y = _Mybegin.y + _Myheight;
	}
	point_type _Mybegin, _Myend;
	value_type _Mywidth, _Myheight;
};
}

// \CLASS TEMPLATE linear range : [begin, end[, stride])
template<typename _Ty, typename _Uy = _Ty>
class range MATRICE_NONHERITABLE 
	: public detail::_Range_base<common_type_t<_Ty, _Uy>> {
	using _Mybase = detail::_Range_base<common_type_t<_Ty, _Uy>>;
public:
	using typename _Mybase::stride_t;
	MATRICE_GLOBAL_INL explicit 
		range(const _Ty& _First, const _Uy& _Last)
		: _Mybase(_First, _Last) {}
	MATRICE_GLOBAL_INL explicit 
		range(const _Ty& _First, const _Uy& _Last, const stride_t _Inc)
		: _Mybase(_First, _Last, _Inc) {}
	MATRICE_GLOBAL_INL explicit 
		range(const _Ty& _First, long _Size)
		: _Mybase(_First, _First+_Size) {}
	MATRICE_GLOBAL_INL explicit 
		range(const _Ty& _First, long _Size, const stride_t _Inc)
		: _Mybase(_First, _First+_Size, _Inc) {}
};

template<typename _Ty> using rect = detail::_Rect_impl<_Ty>;

DGE_MATRICE_END