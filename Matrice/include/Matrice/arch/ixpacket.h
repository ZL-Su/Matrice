#pragma once
#include "_ixbase.h"

#ifdef __AVX__
MATRICE_ARCH_BEGIN

template<int _Elems, typename _Scalar> inline constexpr
auto _Set_packet(const std::initializer_list<_Scalar> _list) {};
template<> inline
auto _Set_packet<8>(const std::initializer_list<float> _list) 
{
	auto _Data = _list.begin();
	return _mm256_set_ps(*_Data, *(_Data++), *(_Data++), *(_Data++), *(_Data++), *(_Data++), *(_Data++), *(_Data++));
};
template<> inline
auto _Set_packet<4>(const std::initializer_list<double> _list)
{
	auto _Data = _list.begin();
	return _mm256_set_pd(*_Data, *(_Data++), *(_Data++), *(_Data++));
};
template<int _Elems, typename _Scalar> inline constexpr
auto _Set_packet(const _Scalar* _Data) {};
template<> inline
auto _Set_packet<8>(const float* _Data)
{
	return _mm256_set_ps(*_Data, *(_Data++), *(_Data++), *(_Data++), *(_Data++), *(_Data++), *(_Data++), *(_Data++));
};
template<> inline
auto _Set_packet<4>(const double* _Data)
{
	return _mm256_set_pd(*_Data, *(_Data++), *(_Data++), *(_Data++));
};

template<typename T, int _Elems>
class Packet_ MATRICE_NONHERITABLE : public simd::simd_base_<T, _Elems>
{
	using xbase_t = simd::simd_base_<T, _Elems>;
	using internal_t = typename xbase_t::internal_t;
	using initlist_t = typename xbase_t::initlist_t;
	using xbase_t::m_data;
public:
	using value_t = typename simd::simd_base_<T, _Elems>::value_t;
	using pointer = value_t*;
	MATRICE_HOST_INL Packet(const value_t _value) noexcept {
		xbase_t::_Set(_value); }
	MATRICE_HOST_INL Packet(const pointer _data) noexcept {
		xbase_t::_Set(_data); }
	MATRICE_HOST_INL Packet(const initlist_t _list) noexcept { 
		xbase_t::_Set(_list.begin()); }


};
#include "./inl/_ixpacket.inl"
MATRICE_ARCH_END
#endif
