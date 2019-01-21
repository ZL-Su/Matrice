#pragma once

#include <type_traits>
#include "_macros.h"
#include "_std_wrapper.h"
#include "_exception.h"

#if (defined __enable_cuda__ && !defined __disable_cuda__)
#include <device_functions.h>
#include <thrust\complex.h>
#endif

DGE_MATRICE_BEGIN

static_assert(sizeof(void *) == 8, "MATRICE supports 64 bit only");

template<typename _Ty = long double>
MATRICE_GLOBAL_INL constexpr _Ty ¦Ð = static_cast<_Ty>(3.14159265358979323846264338327950288419716939937510582097494459);

template<typename T1, typename T2, typename _Ret = common_type_t<T1, T2>, typename = enable_if_t<is_arithmetic_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret add(const T1& a, const T2& b) { return a + b; }
template<typename T1, typename T2, typename _Ret = common_type_t<T1, T2>, typename = enable_if_t<is_arithmetic_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret sub(const T1& a, const T2& b) { return a - b; }
template<typename T1, typename T2, typename _Ret = common_type_t<T1, T2>, typename = enable_if_t<is_arithmetic_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret mul(const T1& a, const T2& b) { return a * b; }
template<typename T1, typename T2, typename _Ret = common_type_t<T1, T2>, typename = enable_if_t<is_arithmetic_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret div(const T1& a, const T2& b) { return a / b; }
template<typename T1, typename T2, typename _Ret = common_type_t<T1, T2>, typename = enable_if_t<is_arithmetic_v<_Ret>>> 
MATRICE_GLOBAL_FINL constexpr _Ret max(const T1& a, const T2& b) { return a < b ? b : a; }
template<typename T1, typename T2, typename _Ret = common_type_t<T1, T2>, typename = enable_if_t<is_arithmetic_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret min(const T1& a, const T2& b) { return a < b ? a : b; }
template<typename T, typename _Ret = T, typename = enable_if_t<is_arithmetic_v<T>>> 
MATRICE_GLOBAL_FINL constexpr _Ret sqrt(const T& x) { return std::sqrt(x); }
template<typename T, typename _Ret = T, typename = enable_if_t<is_arithmetic_v<T>>> 
MATRICE_HOST_FINL constexpr _Ret abs(const T& x) { return std::abs(x); }
template<typename T, typename _Ret = T, typename = enable_if_t<is_arithmetic_v<T>>>
MATRICE_HOST_FINL constexpr _Ret exp(const T& x) { return std::exp(x); }
template<typename T, typename _Ret = T, typename = enable_if_t<is_arithmetic_v<T>>>
MATRICE_HOST_FINL constexpr _Ret log(const T& x) { return std::log(x); }
template<typename T, typename _Ret = T, typename = enable_if_t<is_arithmetic_v<T>>>
MATRICE_HOST_FINL constexpr _Ret log2(const T& x) { return std::log2(x); }
template<typename T, typename _Ret = T, typename = enable_if_t<is_arithmetic_v<T>>>
MATRICE_HOST_FINL constexpr _Ret log10(const T& x) { return std::log10(x); }
template<typename T, typename _Ret = T, typename = enable_if_t<is_arithmetic_v<T>>>
MATRICE_HOST_FINL constexpr _Ret floor(const T& x) { return static_cast<_Ret>(std::floor(x)); }
template<typename T, typename _Ret = T, typename = enable_if_t<is_arithmetic_v<T>>> 
MATRICE_HOST_FINL constexpr _Ret ceil(const T& x) { return static_cast<_Ret>(std::ceil(x)); }
template<typename T1, typename T2, typename _Ret = common_type_t<T1, T2>, typename = enable_if_t<is_arithmetic_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret pow(const T1& x, const T2& y) { return std::pow(x,y); }

namespace detail {
	template<typename _Ty> struct string_to_numval {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return (_Str); }
	};
	template<> struct string_to_numval<int> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stoi(_Str); }
	};
	template<> struct string_to_numval<long> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stol(_Str); }
	};
	template<> struct string_to_numval<float> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stof(_Str); }
	};
	template<> struct string_to_numval<double> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stod(_Str); }
	};
	template<> struct string_to_numval<long double> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stold(_Str); }
	};
	template<> struct string_to_numval<long long> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stoll(_Str); }
	};
	template<> struct string_to_numval<unsigned long> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stoul(_Str); }
	};
	template<> struct string_to_numval<unsigned long long> {
		static MATRICE_HOST_FINL auto value(const std::string& _Str) { return std::stoull(_Str); }
	};
}

/**
 * \Convert a string to a user specified (numerical) value.
 * Example: auto _Ret = dgelom::stonv<float>("1.0");
 */
template<typename T = std::string> MATRICE_HOST_FINL
constexpr T stonv(const std::string& _Str) { return detail::string_to_numval<T>::value(_Str); }

/**
 * \get T-typed zero value
 */
template<typename T, typename = enable_if_t<is_arithmetic_v<T>>>
struct zero { static const constexpr T value = T(0);};
template<typename T> MATRICE_HOST_INL constexpr auto zero_v = zero<T>::value;

/**
 * \append a T-typed element into tuple _Tpl
 */
template<typename T, typename... U> MATRICE_HOST_FINL
std::tuple<U..., T> tuple_append(const std::tuple<U...>& _Tpl, const T& _Val) {
	return std::tuple_cat(_Tpl, std::make_tuple(_Val));
}

/**
 * \packs the first _N element from _E into a tuple
 */
template<size_t _N> struct tuple_n {
	template<typename U> MATRICE_HOST_FINL static auto _(const U& _E) {
		return tuple_append(tuple_n<_N - 1>::_(_E), _E);
	}
	template<typename U> MATRICE_HOST_FINL static auto _(const U* _E) {
		return tuple_append(tuple_n<_N - 1>::_(_E), _E[_N]);
	}
	template<typename U, typename F> MATRICE_HOST_FINL static auto _(const U* _E, F&& _Op) {
		return tuple_append(tuple_n<_N - 1>::_(_E, _Op), _Op(_E[_N]));
	}
};
template<> struct tuple_n<0> {
	template<typename U> MATRICE_HOST_FINL static auto _(const U& _E) {
		return std::make_tuple(_E);
	}
	template<typename U> MATRICE_HOST_FINL static auto _(const U* _E) {
		return std::make_tuple(_E[0]);
	}
	template<typename U, typename F> MATRICE_HOST_FINL static auto _(const U* _E, F&& _Op) {
		return std::make_tuple(_Op(_E[0]));
	}
};

/**
 * \2D shape type, auto [width, height] = shape(width, height)
 */
using shape = tuple<size_t, size_t>;

template<typename _Ity = size_t>
using shape_t = tuple<_Ity,_Ity>;
template<typename _Ity = size_t>
using shape3_t = tuple<_Ity, _Ity, _Ity>;
template<typename _Ity = size_t>
using shape4_t = tuple<_Ity, _Ity, _Ity, _Ity>;

/**
 * \transform functor definitions
 */
struct transforms {
	template<typename _Ty> struct scale {
		using value_type = _Ty;
		template<typename _Uy>
		MATRICE_GLOBAL_INL scale(const _Uy& _Scale = _Uy(1)) : _Myscale(_Scale) {}
		MATRICE_GLOBAL_INL auto operator()(const value_type& _Val)const { return (_Myscale*_Val); }
		value_type _Myscale = 1.;
	};
	template<typename _Ty> struct clamp {
		using value_type = _Ty;
		template<typename _Uy>
		MATRICE_GLOBAL_INL clamp(const _Uy& _Lower, const _Uy& _Upper) : _Mylower(_Lower),_Myupper(_Upper) {}
		MATRICE_GLOBAL_INL auto operator()(const value_type& _Val)const { return min(max(_Val,_Mylower),_Myupper); }
		value_type _Mylower = std::numeric_limits<value_type>::min();
		value_type _Myupper = std::numeric_limits<value_type>::max();
	};
	template<typename _Ty> struct relu {
		using value_type = _Ty;
		MATRICE_GLOBAL_INL relu() {}
		MATRICE_GLOBAL_INL auto operator()(const value_type& _Val)const { return max(_Myzero,_Val); }
		value_type _Myzero = zero_v<value_type>;
	};
};


/**
 *\brief CLASS TEMPLATE basic shape type
 *\param <_Ity> an integral template type
 */
template<typename _Ity = size_t,
	typename = enable_if_t<is_integral_v<_Ity>>>
class basic_shape {
	using _Myt = basic_shape;
public:
	using value_type = _Ity;
	using view_shape = tuple<value_type,value_type>;
	using hist_shape = tuple<value_type,view_shape>;
	using full_shape = tuple<value_type,hist_shape>;

	MATRICE_GLOBAL_INL basic_shape(shape_t<value_type>&& _Shape) noexcept
		: _Data{ 1,{1,_Shape} } {}
	MATRICE_GLOBAL_INL basic_shape(shape3_t<value_type>&& _Shape) noexcept
		: _Data{ 1,{std::get<0>(_Shape), {std::get<1>(_Shape),std::get<2>(_Shape)}} } {}
	MATRICE_GLOBAL_INL basic_shape(shape4_t<value_type>&& _Shape) noexcept
		: _Data{ std::get<0>(_Shape), {std::get<1>(_Shape), {std::get<2>(_Shape),std::get<3>(_Shape)}} } {}
	template<typename _Jty>
	MATRICE_GLOBAL_INL basic_shape(std::initializer_list<_Jty> _Shape) noexcept {
		if (_Shape.size() == 2) {
			_Data = { 1,{1,{(value_type)*_Shape.begin(), (value_type)*(_Shape.begin() + 1)}} };
		}
		if (_Shape.size() == 3) {
			_Data = { 1,{(value_type)*_Shape.begin(),{(value_type)*(_Shape.begin() + 1), (value_type)*(_Shape.begin() + 2)}} };
		}
		if (_Shape.size() == 4) {
			_Data = { (value_type)*_Shape.begin(),{(value_type)*(_Shape.begin() + 1),{(value_type)*(_Shape.begin() + 2), (value_type)*(_Shape.begin() + 3)}} };
		}
	}
	MATRICE_GLOBAL_INL basic_shape(const _Myt& _Other) noexcept
		: _Data(_Other._Data) {}
	MATRICE_GLOBAL_INL basic_shape(_Myt&& _Other) noexcept
		: _Data(std::move(_Other._Data)) {}

	MATRICE_GLOBAL_INL auto& operator= (const _Myt& _Oth) {
		_Data = (_Oth._Data); return (*this);
	}
	MATRICE_GLOBAL_INL auto& operator= (_Myt&& _Oth) {
		_Data = std::move(_Oth._Data); return (*this);
	}
	/**
	 *\brief Get unrolled shape data
	 */
	MATRICE_GLOBAL_INL constexpr auto&& operator()() const {
		return shape4_t<value_type>{ get(0), get(1), get(2), get(3) };
	}
	/**
	 *\brief Get dim value at _Dim
	 */
	MATRICE_GLOBAL_INL constexpr auto get(uint8_t _Dim) const {
		if (_Dim == 0) return std::get<0>(_Data);
		if (_Dim == 1) return std::get<0>(std::get<1>(_Data));
		if (_Dim == 2) return std::get<0>(std::get<1>(std::get<1>(_Data)));
		if (_Dim == 3) return std::get<1>(std::get<1>(std::get<1>(_Data)));
		DGELOM_CHECK(_Dim > 3, "_Dim over range of _Data.");
	}
	/**
	 *\brief Get full rows
	 */
	MATRICE_GLOBAL_INL constexpr auto rows() const {
		return get(0) * get(2);
	}
	/**
	 *\brief Get full cols
	 */
	MATRICE_GLOBAL_INL constexpr auto cols() const {
		return get(1) * get(3);
	}

private:
	full_shape _Data = { value_type(0),{value_type(0),{value_type(0),value_type(0)}} };
};

using basic_shape_t = basic_shape<>;
DGE_MATRICE_END
