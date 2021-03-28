/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library for
3D vision and photo-mechanics.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.
*********************************************************************/
#pragma once
#include "../core.hpp"

MATRICE_ALG_BEGIN(vision)

template<typename _Ty>
concept real_type = is_floating_point_v<_Ty>;

template<typename _Ty>
using vec3d_t = auto_vector_t<_Ty, 3>;
template<typename _Ty>
using vec4d_t = auto_vector_t<_Ty, 4>;

namespace detail {

// \brief Cross product of two 3-vectors.
template<typename _Ty> requires real_type<_Ty>
vec3d_t<_Ty> _Cross(vec3d_t<_Ty> left, vec3d_t<_Ty> right) noexcept {
	const auto _X = left[1] * right[2] - left[2] * right[1];
	const auto _Y = -left[0] * right[2] + left[2] * right[0];
	const auto _Z = left[0] * right[1] - left[1] * right[0];
	return {_X, _Y, _Z};
}

// \brief Compute normalization of a 3-vector.
template<typename _Ty> requires real_type<_Ty>
vec3d_t<_Ty> _Normalize(_Ty x, _Ty y, _Ty z) noexcept {
	MATRICE_USE_STD(pow);
	const auto a = pow(sqsum(x, y, z), 0.5);
	return{x / a, y / a, z / a};
}

// \brief Compute intersection of a line and a plane.
// \param 'point' end point of the line.
// \param 'normal' direction of the line.
// \param 'plane' coefs A, B, C, and D of the plane: Ax + By + Cz + D = 0.
template<typename _Ty> requires real_type<_Ty>
vec3d_t<_Ty> _Line_and_plane_intersection(vec3d_t<_Ty> point, vec3d_t<_Ty> normal, vec4d_t<_Ty> plane)
{
	const auto x = -((plane[3] * normal[0] - plane[1] * normal[1] * point[0] - plane[2] * normal[2] * point[0] + plane[1] * normal[0] * point[1] + plane[2] * normal[0] * point[2]) / (plane[0] * normal[0] + plane[1] * normal[1] + plane[2] * normal[2]));
	const auto y = -((plane[3] * normal[1] + plane[0] * normal[1] * point[0] - plane[0] * normal[0] * point[1] - plane[2] * normal[2] * point[1] + plane[2] * normal[1] * point[2]) / (plane[0] * normal[0] + plane[1] * normal[1] + plane[2] * normal[2]));
	const auto z = -((plane[3] * normal[2] + plane[0] * normal[2] * point[0] + plane[1] * normal[2] * point[1] - plane[0] * normal[0] * point[2] - plane[1] * normal[1] * point[2]) / (plane[0] * normal[0] + plane[1] * normal[1] + plane[2] * normal[2]));
	
	return { x, y, z };
}

// \brief Compute the included angle of two 3-vectors.
template<typename _Ty> requires real_type<_Ty>
_Ty _Vector_angle(vec3d_t<_Ty> vec1, vec3d_t<_Ty> vec2) noexcept {
	MATRICE_USE_STD(acos);

	const auto vec11 = _Normalize(vec1[0], vec1[1], vec1[2]);
	const auto vec22 = _Normalize(vec2[0], vec2[1], vec2[2]);
	return acos(vec11[0] * vec22[0] + vec11[1] * vec22[1] + vec11[2] * vec22[2]);
}

/// <summary>
/// \brief Compute the unit direction of the refracted ray. 
/// </summary>
/// <param name="'vec1'">Incident vector</param>
/// <param name="'N1'">Normal of the interface</param>
/// <param name="'n1'">Refractive index of incident medium</param>
/// <param name="'n2'">Refractive index of emergent medium</param>
/// <param name="'theta1'">Incident angle</param>
/// <param name="'theta2'">Emergent angle</param>
/// <returns>Unit direction vector</returns>
template<typename _Ty> requires real_type<_Ty>
vec3d_t<_Ty> _Get_refracted_vector(vec3d_t<_Ty> vec1, vec3d_t<_Ty> N1, 
	_Ty n1, _Ty n2, _Ty theta1, _Ty theta2) noexcept {
	MATRICE_USE_STD(cos);

	const auto _X = n1 / n2 * vec1[0] - (n1 / n2 * cos(theta1) - cos(theta2)) * N1[0];
	const auto _Y = n1 / n2 * vec1[1] - (n1 / n2 * cos(theta1) - cos(theta2)) * N1[1];
	const auto _Z = n1 / n2 * vec1[2] - (n1 / n2 * cos(theta1) - cos(theta2)) * N1[2];

	return { _X, _Y, _Z };
}

/// <summary>
/// \brief Compute the intersection of two lines.
/// </summary>
/// <param name="'point1'">End point of line 1;</param>
/// <param name="'vec1'">Direction of line 1;</param>
/// <param name="'point2'">End point of line 2;</param>
/// <param name="'vec2'">Direction of line 2.</param>
/// <returns>Coordinates of the intersection</returns>
template<typename _Ty> requires real_type<_Ty>
vec3d_t<_Ty> _Intersection_of_lines(vec3d_t<_Ty> point1, vec3d_t<_Ty> vec1,
	vec3d_t<_Ty> point2, vec3d_t<_Ty> vec2) noexcept {
	const auto vec3 = _Cross(vec1, vec2);
	const auto x1 = -((point1[2] * vec1[0] * vec2[1] * vec3[0] - point2[2] * vec1[0] * vec2[1] * vec3[0] -
		point1[0] * vec1[2] * vec2[1] * vec3[0] - point1[1] * vec1[0] * vec2[2] * vec3[0] +
		point2[1] * vec1[0] * vec2[2] * vec3[0] + point1[0] * vec1[1] * vec2[2] * vec3[0] -
		point1[2] * vec1[0] * vec2[0] * vec3[1] + point2[2] * vec1[0] * vec2[0] * vec3[1] +
		point1[0] * vec1[2] * vec2[0] * vec3[1] - point2[0] * vec1[0] * vec2[2] * vec3[1] +
		point1[1] * vec1[0] * vec2[0] * vec3[2] - point2[1] * vec1[0] * vec2[0] * vec3[2] -
		point1[0] * vec1[1] * vec2[0] * vec3[2] + point2[0] * vec1[0] * vec2[1] * vec3[2]) /
		(vec1[2] * vec2[1] * vec3[0] - vec1[1] * vec2[2] * vec3[0] - vec1[2] * vec2[0] * vec3[1] +
			vec1[0] * vec2[2] * vec3[1] + vec1[1] * vec2[0] * vec3[2] - vec1[0] * vec2[1] * vec3[2]));
	const auto y1 = -((-point1[2] * vec1[1] * vec2[1] * vec3[0] + point2[2] * vec1[1] * vec2[1] * vec3[0] +
		point1[1] * vec1[2] * vec2[1] * vec3[0] - point2[1] * vec1[1] * vec2[2] * vec3[0] +
		point1[2] * vec1[1] * vec2[0] * vec3[1] - point2[2] * vec1[1] * vec2[0] * vec3[1] -
		point1[1] * vec1[2] * vec2[0] * vec3[1] + point1[1] * vec1[0] * vec2[2] * vec3[1] -
		point1[0] * vec1[1] * vec2[2] * vec3[1] + point2[0] * vec1[1] * vec2[2] * vec3[1] +
		point2[1] * vec1[1] * vec2[0] * vec3[2] - point1[1] * vec1[0] * vec2[1] * vec3[2] +
		point1[0] * vec1[1] * vec2[1] * vec3[2] - point2[0] * vec1[1] * vec2[1] * vec3[2]) /
		(-vec1[2] * vec2[1] * vec3[0] + vec1[1] * vec2[2] * vec3[0] + vec1[2] * vec2[0] * vec3[1] -
			vec1[0] * vec2[2] * vec3[1] - vec1[1] * vec2[0] * vec3[2] + vec1[0] * vec2[1] * vec3[2]));
	const auto z1 = -((point2[2] * vec1[2] * vec2[1] * vec3[0] - point1[2] * vec1[1] * vec2[2] * vec3[0] +
		point1[1] * vec1[2] * vec2[2] * vec3[0] - point2[1] * vec1[2] * vec2[2] * vec3[0] -
		point2[2] * vec1[2] * vec2[0] * vec3[1] + point1[2] * vec1[0] * vec2[2] * vec3[1] -
		point1[0] * vec1[2] * vec2[2] * vec3[1] + point2[0] * vec1[2] * vec2[2] * vec3[1] +
		point1[2] * vec1[1] * vec2[0] * vec3[2] - point1[1] * vec1[2] * vec2[0] * vec3[2] +
		point2[1] * vec1[2] * vec2[0] * vec3[2] - point1[2] * vec1[0] * vec2[1] * vec3[2] +
		point1[0] * vec1[2] * vec2[1] * vec3[2] - point2[0] * vec1[2] * vec2[1] * vec3[2]) /
		(-vec1[2] * vec2[1] * vec3[0] + vec1[1] * vec2[2] * vec3[0] + vec1[2] * vec2[0] * vec3[1] -
			vec1[0] * vec2[2] * vec3[1] - vec1[1] * vec2[0] * vec3[2] + vec1[0] * vec2[1] * vec3[2]));

	const auto x2 = -((-point1[2] * vec1[1] * vec2[0] * vec3[0] + point2[2] * vec1[1] * vec2[0] * vec3[0] +
		point1[1] * vec1[2] * vec2[0] * vec3[0] - point2[1] * vec1[2] * vec2[0] * vec3[0] +
		point2[0] * vec1[2] * vec2[1] * vec3[0] - point2[0] * vec1[1] * vec2[2] * vec3[0] +
		point1[2] * vec1[0] * vec2[0] * vec3[1] - point2[2] * vec1[0] * vec2[0] * vec3[1] -
		point1[0] * vec1[2] * vec2[0] * vec3[1] + point2[0] * vec1[0] * vec2[2] * vec3[1] -
		point1[1] * vec1[0] * vec2[0] * vec3[2] + point2[1] * vec1[0] * vec2[0] * vec3[2] +
		point1[0] * vec1[1] * vec2[0] * vec3[2] - point2[0] * vec1[0] * vec2[1] * vec3[2]) /
		(-vec1[2] * vec2[1] * vec3[0] + vec1[1] * vec2[2] * vec3[0] + vec1[2] * vec2[0] * vec3[1] -
			vec1[0] * vec2[2] * vec3[1] - vec1[1] * vec2[0] * vec3[2] + vec1[0] * vec2[1] * vec3[2]));
	const auto y2 = -((-point1[2] * vec1[1] * vec2[1] * vec3[0] + point2[2] * vec1[1] * vec2[1] * vec3[0] +
		point1[1] * vec1[2] * vec2[1] * vec3[0] - point2[1] * vec1[1] * vec2[2] * vec3[0] -
		point2[1] * vec1[2] * vec2[0] * vec3[1] + point1[2] * vec1[0] * vec2[1] * vec3[1] -
		point2[2] * vec1[0] * vec2[1] * vec3[1] - point1[0] * vec1[2] * vec2[1] * vec3[1] +
		point2[0] * vec1[2] * vec2[1] * vec3[1] + point2[1] * vec1[0] * vec2[2] * vec3[1] +
		point2[1] * vec1[1] * vec2[0] * vec3[2] - point1[1] * vec1[0] * vec2[1] * vec3[2] +
		point1[0] * vec1[1] * vec2[1] * vec3[2] - point2[0] * vec1[1] * vec2[1] * vec3[2]) /
		(-vec1[2] * vec2[1] * vec3[0] + vec1[1] * vec2[2] * vec3[0] +
			vec1[2] * vec2[0] * vec3[1] - vec1[0] * vec2[2] * vec3[1] -
			vec1[1] * vec2[0] * vec3[2] + vec1[0] * vec2[1] * vec3[2]));
	const auto z2 = -((point2[2] * vec1[2] * vec2[1] * vec3[0] - point1[2] * vec1[1] * vec2[2] * vec3[0] +
		point1[1] * vec1[2] * vec2[2] * vec3[0] - point2[1] * vec1[2] * vec2[2] * vec3[0] -
		point2[2] * vec1[2] * vec2[0] * vec3[1] + point1[2] * vec1[0] * vec2[2] * vec3[1] -
		point1[0] * vec1[2] * vec2[2] * vec3[1] + point2[0] * vec1[2] * vec2[2] * vec3[1] +
		point2[2] * vec1[1] * vec2[0] * vec3[2] - point2[2] * vec1[0] * vec2[1] * vec3[2] -
		point1[1] * vec1[0] * vec2[2] * vec3[2] + point2[1] * vec1[0] * vec2[2] * vec3[2] +
		point1[0] * vec1[1] * vec2[2] * vec3[2] - point2[0] * vec1[1] * vec2[2] * vec3[2]) /
		(-vec1[2] * vec2[1] * vec3[0] + vec1[1] * vec2[2] * vec3[0] +
			vec1[2] * vec2[0] * vec3[1] - vec1[0] * vec2[2] * vec3[1] -
			vec1[1] * vec2[0] * vec3[2] + vec1[0] * vec2[1] * vec3[2]));

	return { (x1 + x2) / 2,	(y1 + y2) / 2, (z1 + z2) / 2 };
}
}

/********************************************************************************************
 * \brief Driver for refractive 3D reconstruction. The function is thread-safe.
 * \param 'raw_point' distorted object point P.
 * \param 'translation' translation vector of the right camera relative to the frame O-XYZ.
 * \param 'interf1' refracting interface Ax + By + Cz +D = 0, with N = (A. B. C).
 * \param 'thickness' distance between the two interfaces.
 * \param 'n1, n2, n3' refractive indices.
 * \returns The true 3D coordinates of the object point.
 ********************************************************************************************/
template<typename _Ty> requires real_type<_Ty>
vec3d_t<_Ty> refractive_3d_reconstruction(
	vec3d_t<_Ty> raw_point, 
	vec3d_t<_Ty> translation, 
	vec4d_t<_Ty> interf1, 
	_Ty thickness, _Ty n1, _Ty n2, _Ty n3)
{
	MATRICE_USE_STD(pow);
	MATRICE_USE_STD(sin);
	MATRICE_USE_STD(cos);

	const auto d = thickness;
	const auto t = translation;
	
	// Build the interface 2
	const auto D2 = interf1[3] - d * pow(sqsum(interf1[0], interf1[1], interf1[2]), 0.5);
	decltype(interf1) interf2 = { interf1[0], interf1[1], interf1[2], D2 };

	// Compute the normalized normal of the interfaces 1 and 2
	const auto N1 = detail::_Normalize(interf1[0], interf1[1], interf1[2]);
	const auto N2 = detail::_Normalize(interf2[0], interf2[1], interf2[2]);

	// Direction vectors of the ray L1 and L1'
	auto nlineL1 = detail::_Normalize(raw_point[0], raw_point[1], raw_point[2]);
	auto nlineR1 = detail::_Normalize(raw_point[0] + t[0], raw_point[1] + t[1], raw_point[2] + t[2]);

	// Incident points P1 and P1' of L1 and L1' at the interface 1
	auto pointL1 = detail::_Line_and_plane_intersection({ 0.0, 0.0, 0.0 }, nlineL1, interf1);
	auto pointR1 = detail::_Line_and_plane_intersection({ -t[0], -t[1], -t[2] }, nlineR1, interf1);

	// Incident angles of rays L1 and L1'
	const auto angleL1 = detail::_Vector_angle(nlineL1, N1);
	const auto angleR1 = detail::_Vector_angle(nlineR1, N1);

	// Refracted angles of rays L1 and L1' 
	const auto angleL2 = asin(n1 * sin(angleL1) / n2);
	const auto angleR2 = asin(n1 * sin(angleR1) / n2);

	// Direction vectors of the ray L2 and L2'
	auto nlineL2 = detail::_Get_refracted_vector(nlineL1, N1, n1, n2, angleL1, angleL2);
	auto nlineR2 = detail::_Get_refracted_vector(nlineR1, N1, n1, n2, angleR1, angleR2);
	nlineL2 = detail::_Normalize(nlineL2[0], nlineL2[1], nlineL2[2]);
	nlineR2 = detail::_Normalize(nlineR2[0], nlineR2[1], nlineR2[2]);

	// Incident points P2 and P2' of L2 and L2' at the interface 2
	auto pointL2 = detail::_Line_and_plane_intersection(pointL1, nlineL2, interf2);
	auto pointR2 = detail::_Line_and_plane_intersection(pointR1, nlineR2, interf2);

	// Incident angles of rays L2 and L2'
	const auto angleL3 = angleL2;
	const auto angleR3 = angleR2;

	// Refracted angles of rays L2 and L2'
	auto angleL4 = asin(n1 * sin(angleL3) / n2);
	auto angleR4 = asin(n1 * sin(angleR3) / n2);

	// Direction vectors of the ray L3 and L3'
	const auto nlineL3 = detail::_Get_refracted_vector(nlineL2, N2, n2, n3, angleL3, angleL4);
	const auto nlineR3 = detail::_Get_refracted_vector(nlineR2, N2, n2, n3, angleR3, angleR4);

	//auto nCommonVerticalLine3 = detail::_Cross(nlineL3, nlineR3);
	
	// Compute the true object point Q
	return detail::_Intersection_of_lines(pointL2, nlineL3, pointR2, nlineR3);
}

MATRICE_ALG_END(vision)