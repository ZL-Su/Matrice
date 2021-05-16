#pragma once
#include "../transform.h"
#ifdef MATRICE_SIMD_ARCH
#include "arch/ixpacket.h"
#endif

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
/// <summary>
/// \brief X axis represetation
/// </summary>
template<> struct _Axis_type<0, 3> {
	MATRICE_GLOBAL_FINL auto operator()() const {
		return Vec3_<bool>(1, 0, 0);
	}
};
/// <summary>
/// \brief Y axis represetation
/// </summary>
template<> struct _Axis_type<1, 3> {
	MATRICE_GLOBAL_FINL auto operator()() const {
		return Vec3_<bool>(0, 1, 0);
	}
};
/// <summary>
/// \brief Z axis represetation
/// </summary>
template<> struct _Axis_type<2, 3> {
	MATRICE_GLOBAL_FINL auto operator()() const {
		return Vec3_<bool>(0, 0, 1);
	}
};

template<typename _Ty> 
class _Axis_angle_rep<_Ty, 3> {
public:
	using angle_type = _Ty;
	using axis_type = Vec3_<angle_type>;

	_Axis_angle_rep(const axis_type& _Axis, const angle_type& _Angl)
		:_Mydata(_Axis, _Angl) {}
	template<uint8_t... _Dims>
	_Axis_angle_rep(const _Axis_type<_Dims...>& _Axis, const angle_type& _Angl)
		: _Mydata(_Axis(), _Angl) {}

	MATRICE_GLOBAL_INL constexpr auto& axis() const noexcept { 
		return MATRICE_STD(get)<0>(_Mydata);
	}
	MATRICE_GLOBAL_INL constexpr auto& angle() const noexcept { 
		return MATRICE_STD(get) <1>(_Mydata);
	}

	/**
	 *\brief Converts to an equivalent rotation matrix
	 *\cite{https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation}
	 */
	template<size_t _Dim = 3, MATRICE_ENABLE_IF(_Dim==3||_Dim==4)>
	MATRICE_GLOBAL_INL auto matrix() const {
		using sqr_matrix_t = Matrix_<angle_type, _Dim, _Dim>;

		const auto s = sin(angle());
		const auto c = cos(angle());
		const auto K = cross_prod_matrix<_Dim>(axis());
		
		return sqr_matrix_t(sqr_matrix_t::diag(1) + s*K + c*K*K);
	}
private:
	tuple<axis_type, angle_type> _Mydata;
};

#ifdef MATRICE_SIMD_ARCH
template<typename _Ty> struct _Avxmat4_ {
	using type = simd::Packet_<_Ty, 4>;
	using iterator = typename type::pointer;
	MATRICE_HOST_INL _Avxmat4_(const iterator _Data)
		:_pa(_Data),
		_pb(_Data + type::size),
		_pc(_Data + type::size * 2),
		_pd(_Data + type::size * 3) {}
	MATRICE_HOST_INL _Avxmat4_(const iterator _Data, size_t)
		: _pa(_Data), _pb(_Data + 3), _pc(_Data + 6), _pd({ 0,0,0,1 }) {}

	template<size_t _N>
	MATRICE_HOST_INL auto mul(iterator _Data) {
		type pd(_Data);
		if constexpr (_N == ~- type::size) {
			pd[3] = typename type::value_t(1);
		}
		_Data[0] = simd::reduce(_pa*pd);
		_Data[1] = simd::reduce(_pb*pd);
		_Data[2] = simd::reduce(_pc*pd);
		if constexpr (_N == type::size) {
			_Data[3] = simd::reduce(_pd*pd);
		}
	}

	type _pa, _pb, _pc, _pd;
};
#endif

template<typename _Ty> struct _Geotf_base {
public:
	using value_type = _Ty;
	using vec4_type = Vec4_<value_type>;
	using vec3_type = Vec3_<value_type>;
	using matrix_type = Matrix_<value_type, 4, 4>;

#ifdef MATRICE_SIMD_ARCH
	/**
	 *\brief Constructor from 4-by-4 transformation matrix
	 *\param [m] 4-by-4 transformation matrix
	 */
	_Geotf_base(const matrix_type& m)
		: _Mymat(m.data()) {
	}
	/**
	 *\brief Constructor from axis-angle representation
	 *\param [a] axis-angle data
	 */
	_Geotf_base(const _Axis_angle_rep<value_type, 3>& a)
		: _Geotf_base(a.matrix<4>()) {
	}
	/**
	 *\brief Constructor from a rotation matrix and a translation vector
	 *\param [R] 3-by-3 rotation matrix
	 *\param [T] 3-translation vector
	 */
	_Geotf_base(const Matrix_<value_type, 3>& R, const vec3_type& T)
		: _Mymat(R.data(), R.rows()) {
		_Mymat._pa[3] = T[0], _Mymat._pb[3] = T[1], _Mymat._pc[3] = T[2];
	}

	/**
	 *\brief In-place transformation
	 *\param [p] input point, overwritten by the transformed data
	 */
	template<typename _Vty, MATRICE_ENABLE_IF(is_fxdvector_v<_Vty>)>
	MATRICE_HOST_INL decltype(auto) operator()(const _Vty& p) const {
		constexpr auto _Size = p.Size;
		static_assert(_Size >= 3, 
			"Size of the input 'p' must be not less than 3.");
		_Mymat.mul<_Size>(p.data());
		return (p);
	}
protected:
	mutable _Avxmat4_<value_type> _Mymat;
#else
	/**
	 * \brief Ctor from 4-by-4 transformation matrix.
	 * \param 'm' 4-by-4 transformation matrix.
	 */
	_Geotf_base(const matrix_type& m) noexcept
		: _Mymat(m) {
	}
	/**
	 * \brief Ctor from axis-angle representation.
	 * \param 'a' axis-angle data.
	 */
	_Geotf_base(const _Axis_angle_rep<value_type, 3>& a) noexcept
		: _Geotf_base(a.matrix<matrix_type::rows_at_compiletime>()) {
	}
	/**
	 * \brief Ctor from a rotation matrix and a translation vector.
	 * \param 'R' 3-by-3 rotation matrix.
	 * \param 'T' 3-translation vector.
	 */
	_Geotf_base(const Matrix_<value_type, 3>& R, const vec3_type& T) noexcept {
		_Mymat.block(0, 3, 0, 3) = R;
		_Mymat[0][3] = T[0], _Mymat[1][3] = T[1], _Mymat[2][3] = T[2];
		_Mymat.rview(3) = { 0, 0, 0, 1 };
	}

	/**
	 *\brief In-place transformation.
	 *\param 'p' input point, overwritten by the transformed data.
	 */
	template<typename _Vty, MATRICE_ENABLE_IF(is_fxdvector_v<_Vty>)>
	MATRICE_HOST_INL decltype(auto) operator()(const _Vty& p) const {
		constexpr auto _Size = p.Size;
		static_assert(_Size >= 3,
			"Size of the input 'p' must be not less than 3.");
		vec4_type tmp(p);
		if constexpr (_Size == 3) {
			tmp(3) = 1;
		}
		p = _Mymat.mul(tmp);
		return (p);
	}
protected:
	matrix_type _Mymat;
#endif

public:
	/**
	 *\brief Compute rotation matrix between two 3d vectors
	 *\param [v1, v2] the given two vectors
	*/
	static MATRICE_HOST_INL auto rotation(const vec3_type& v1, const vec3_type& v2) {
		return (_Rotation_between(v1, v2));
	}
};
_DETAIL_END

template<>
class detail::_Geotf_isometry<float> 
	: public detail::_Geotf_base<float> {
	using _Mybase = detail::_Geotf_base<float>;
public:
	using _Mybase::_Geotf_base;
	using _Mybase::operator();
private:
};

template<>
class detail::_Geotf_isometry<double> 
	: public detail::_Geotf_base<double> {
	using _Mybase = detail::_Geotf_base<double>;
public:
	using _Mybase::_Geotf_base;
	using _Mybase::operator();
private:
};

DGE_MATRICE_END