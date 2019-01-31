#pragma once
#include "../transform.h"
#include "../../../arch/ixpacket.h"

DGE_MATRICE_BEGIN
namespace detail {
template<> struct _Axis_type<0, 3> {
	MATRICE_GLOBAL_FINL auto operator()() const {
		return Vec3_<bool>(1, 0, 0);
	}
};
template<> struct _Axis_type<1, 3> {
	MATRICE_GLOBAL_FINL auto operator()() const {
		return Vec3_<bool>(0, 1, 0);
	}
};
template<> struct _Axis_type<2, 3> {
	MATRICE_GLOBAL_FINL auto operator()() const {
		return Vec3_<bool>(0, 0, 1);
	}
};

template<typename _Ty> struct _Avxmat4_ {
	using type = simd::Packet_<_Ty, 4>;
	using iterator = typename type::pointer;
	MATRICE_HOST_INL _Avxmat4_(const iterator _Data)
		:_pa(_Data), 
		_pb(_Data + type::size), 
		_pc(_Data + type::size*2),
		_pd(_Data + type::size*3){}
	MATRICE_HOST_INL _Avxmat4_(const iterator _Data, size_t)
		: _pa(_Data), _pb(_Data + 3), _pc(_Data + 6), _pd({ 0,0,0,1 }) {}

	template<size_t _N>
	MATRICE_HOST_INL auto mul(iterator _Data) {
		type pd(_Data);
		if constexpr (_N == ~-type::size) {
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

template<typename _Ty> class _Axis_angle_rep<_Ty, 3> {
public:
	using angle_type = _Ty;
	using axis_type = Vec3_<angle_type>;

	_Axis_angle_rep(const axis_type& _Axis, const angle_type& _Angl)
		:_Data(_Axis, _Angl) {}
	template<uint8_t... _Dims>
	_Axis_angle_rep(const _Axis_type<_Dims...>& _Axis, const angle_type& _Angl)
		:_Data(_Axis(), _Angl) {}

	MATRICE_GLOBAL_INL constexpr auto& axis() const { return std::get<0>(_Data); }
	MATRICE_GLOBAL_INL constexpr auto& angle() const { return std::get<1>(_Data); }

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
	tuple<axis_type, angle_type> _Data;
};

template<typename _Ty> struct _Geotf_base {
public:
	using value_type = _Ty;
	using vec4_type = Vec4_<value_type>;
	using vec3_type = Vec3_<value_type>;
	using matrix_type = Matrix_<value_type, 4, 4>;

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
	_Geotf_base(const Matrix_<value_type, 3, 3>& R, const vec3_type& T)
		: _Mymat(R.data(), R.rows()) {
		_Mymat._pa[3] = T[0], _Mymat._pb[3] = T[1], _Mymat._pc[3] = T[2];
	}

	/**
	 *\brief In-place transformation
	 *\param [p] input point, overwritten by the transformed data
	 */
	template<typename _Vty, typename = enable_if_t<is_fxdvector_v<_Vty>>>
	MATRICE_HOST_INL auto& operator()(const _Vty& p) const {
		_Mymat.mul<p.Size>(p.data());
		return (p);
	}

	/**
	 *\brief Compute rotation matrix between two 3d vectors
	 *\param [v1, v2] the given two vectors
	*/
	static MATRICE_HOST_INL auto rotation(const vec3_type& v1, const vec3_type& v2) {
		return std::move(_Rot_from(v1, v2));
	}

protected:
	mutable _Avxmat4_<value_type> _Mymat;
};
}

template<>
class detail::_GeoTransform_isometry<float> 
	: public detail::_Geotf_base<float> {
	using _Mybase = detail::_Geotf_base<float>;
public:
	using _Mybase::_Geotf_base;
	using _Mybase::operator();
private:
};

template<>
class detail::_GeoTransform_isometry<double> 
	: public detail::_Geotf_base<double> {
	using _Mybase = detail::_Geotf_base<double>;
public:
	using _Mybase::_Geotf_base;
	using _Mybase::operator();
private:
};

DGE_MATRICE_END
