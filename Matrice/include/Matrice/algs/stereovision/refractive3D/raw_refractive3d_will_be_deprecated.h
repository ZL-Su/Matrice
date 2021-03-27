/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library for
3D vision and photo-mechanics.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.
*********************************************************************/
#pragma once
#include "../core.hpp"

MATRICE_ALG_BEGIN(vision)

template<typename _Ty> 
using point3d_t = std::array<_Ty, 3>;

namespace detail{
//�������
template<typename _Ty>
point3d_t<_Ty> _Cross(point3d_t<_Ty> vector1, point3d_t<_Ty> vector2) {
	point3d_t<_Ty> result;
	result[0] = vector1[1] * vector2[2] - vector1[2] * vector2[1];
	result[1] = -vector1[0] * vector2[2] + vector1[2] * vector2[0];
	result[2] = vector1[0] * vector2[1] - vector1[1] * vector2[0];
	return result;
}

//��һ��
template<typename _Ty>
point3d_t<_Ty> _Normalize(_Ty x, _Ty y, _Ty z) {
	MATRICE_USE_STD(pow);

	point3d_t<_Ty> result;
	auto a = pow((pow(x, 2) + pow(y, 2) + pow(z, 2)), 0.5);
	result[0] = x / a;
	result[1] = y / a;
	result[2] = z / a;
	return result;
}

//��㷨ʽֱ�ߺ�һ��ʽƽ��Ľ���point-ֱ�ߵĵ㣬normalVectorֱ�ߵķ�����,plane-ƽ���abcd
template<typename _Ty>
point3d_t<_Ty> _Line_and_plane_intersection(point3d_t<_Ty> point, point3d_t<_Ty> normalVector, std::array<_Ty, 4> plane)
{
	const auto x = -((plane[3] * normalVector[0] - plane[1] * normalVector[1] * point[0] - plane[2] * normalVector[2] * point[0] + plane[1] * normalVector[0] * point[1] + plane[2] * normalVector[0] * point[2]) / (plane[0] * normalVector[0] + plane[1] * normalVector[1] + plane[2] * normalVector[2]));
	const auto y = -((plane[3] * normalVector[1] + plane[0] * normalVector[1] * point[0] - plane[0] * normalVector[0] * point[1] - plane[2] * normalVector[2] * point[1] + plane[2] * normalVector[1] * point[2]) / (plane[0] * normalVector[0] + plane[1] * normalVector[1] + plane[2] * normalVector[2]));
	const auto z = -((plane[3] * normalVector[2] + plane[0] * normalVector[2] * point[0] + plane[1] * normalVector[2] * point[1] - plane[0] * normalVector[0] * point[2] - plane[1] * normalVector[1] * point[2]) / (plane[0] * normalVector[0] + plane[1] * normalVector[1] + plane[2] * normalVector[2]));
	point3d_t<_Ty> result = { x, y, z };
	return result;
}


//�������н�
template<typename _Ty>
_Ty _Vector_angle(point3d_t<_Ty> vector1, point3d_t<_Ty> vector2) {
	MATRICE_USE_STD(acos);

	const auto vector11 = _Normalize(vector1[0], vector1[1], vector1[2]);
	const auto vector22 = _Normalize(vector2[0], vector2[1], vector2[2]);
	return acos(vector11[0] * vector22[0] + vector11[1] * vector22[1] + vector11[2] * vector22[2]);
}


//��������߷����������������߷�������.vector1 ֱ�߷���������N1 ƽ�淨������n1 n2�����������ʵ������ʣ�angel1 angel2 ����������
template<typename _Ty>
point3d_t<_Ty> _Get_refracted_vector(point3d_t<_Ty> vector1, point3d_t<_Ty> N1, 
	double n1, double n2, double angle1, double angle2) {
	MATRICE_USE_STD(cos);

	point3d_t<_Ty> result;
	result[0] = n1 / n2 * vector1[0] - (n1 / n2 * cos(angle1) - cos(angle2)) * N1[0];
	result[1] = n1 / n2 * vector1[1] - (n1 / n2 * cos(angle1) - cos(angle2)) * N1[1];
	result[2] = n1 / n2 * vector1[2] - (n1 / n2 * cos(angle1) - cos(angle2)) * N1[2];
	return result;
}


// ��ռ������㷨ʽֱ�ߵĽ������ꡣpoint1��vector1��ֱ��1�ĵ�ͷ���point2��vector2��ֱ��2�ĵ�ͷ�.
template<typename _Ty>
point3d_t<_Ty> _Intersection_of_lines(point3d_t<_Ty> point1, point3d_t<_Ty> vector1,
	point3d_t<_Ty> point2, point3d_t<_Ty> vector2)
{
	point3d_t<_Ty> result;
	const auto vector3 = _Cross(vector1, vector2);
	const auto x1 = -((point1[2] * vector1[0] * vector2[1] * vector3[0] - point2[2] * vector1[0] * vector2[1] * vector3[0] -
		point1[0] * vector1[2] * vector2[1] * vector3[0] - point1[1] * vector1[0] * vector2[2] * vector3[0] +
		point2[1] * vector1[0] * vector2[2] * vector3[0] + point1[0] * vector1[1] * vector2[2] * vector3[0] -
		point1[2] * vector1[0] * vector2[0] * vector3[1] + point2[2] * vector1[0] * vector2[0] * vector3[1] +
		point1[0] * vector1[2] * vector2[0] * vector3[1] - point2[0] * vector1[0] * vector2[2] * vector3[1] +
		point1[1] * vector1[0] * vector2[0] * vector3[2] - point2[1] * vector1[0] * vector2[0] * vector3[2] -
		point1[0] * vector1[1] * vector2[0] * vector3[2] + point2[0] * vector1[0] * vector2[1] * vector3[2]) /
		(vector1[2] * vector2[1] * vector3[0] - vector1[1] * vector2[2] * vector3[0] - vector1[2] * vector2[0] * vector3[1] +
			vector1[0] * vector2[2] * vector3[1] + vector1[1] * vector2[0] * vector3[2] - vector1[0] * vector2[1] * vector3[2]));
	const auto y1 = -((-point1[2] * vector1[1] * vector2[1] * vector3[0] + point2[2] * vector1[1] * vector2[1] * vector3[0] +
		point1[1] * vector1[2] * vector2[1] * vector3[0] - point2[1] * vector1[1] * vector2[2] * vector3[0] +
		point1[2] * vector1[1] * vector2[0] * vector3[1] - point2[2] * vector1[1] * vector2[0] * vector3[1] -
		point1[1] * vector1[2] * vector2[0] * vector3[1] + point1[1] * vector1[0] * vector2[2] * vector3[1] -
		point1[0] * vector1[1] * vector2[2] * vector3[1] + point2[0] * vector1[1] * vector2[2] * vector3[1] +
		point2[1] * vector1[1] * vector2[0] * vector3[2] - point1[1] * vector1[0] * vector2[1] * vector3[2] +
		point1[0] * vector1[1] * vector2[1] * vector3[2] - point2[0] * vector1[1] * vector2[1] * vector3[2]) /
		(-vector1[2] * vector2[1] * vector3[0] + vector1[1] * vector2[2] * vector3[0] + vector1[2] * vector2[0] * vector3[1] -
			vector1[0] * vector2[2] * vector3[1] - vector1[1] * vector2[0] * vector3[2] + vector1[0] * vector2[1] * vector3[2]));
	const auto z1 = -((point2[2] * vector1[2] * vector2[1] * vector3[0] - point1[2] * vector1[1] * vector2[2] * vector3[0] +
		point1[1] * vector1[2] * vector2[2] * vector3[0] - point2[1] * vector1[2] * vector2[2] * vector3[0] -
		point2[2] * vector1[2] * vector2[0] * vector3[1] + point1[2] * vector1[0] * vector2[2] * vector3[1] -
		point1[0] * vector1[2] * vector2[2] * vector3[1] + point2[0] * vector1[2] * vector2[2] * vector3[1] +
		point1[2] * vector1[1] * vector2[0] * vector3[2] - point1[1] * vector1[2] * vector2[0] * vector3[2] +
		point2[1] * vector1[2] * vector2[0] * vector3[2] - point1[2] * vector1[0] * vector2[1] * vector3[2] +
		point1[0] * vector1[2] * vector2[1] * vector3[2] - point2[0] * vector1[2] * vector2[1] * vector3[2]) /
		(-vector1[2] * vector2[1] * vector3[0] + vector1[1] * vector2[2] * vector3[0] + vector1[2] * vector2[0] * vector3[1] -
			vector1[0] * vector2[2] * vector3[1] - vector1[1] * vector2[0] * vector3[2] + vector1[0] * vector2[1] * vector3[2]));

	const auto x2 = -((-point1[2] * vector1[1] * vector2[0] * vector3[0] + point2[2] * vector1[1] * vector2[0] * vector3[0] +
		point1[1] * vector1[2] * vector2[0] * vector3[0] - point2[1] * vector1[2] * vector2[0] * vector3[0] +
		point2[0] * vector1[2] * vector2[1] * vector3[0] - point2[0] * vector1[1] * vector2[2] * vector3[0] +
		point1[2] * vector1[0] * vector2[0] * vector3[1] - point2[2] * vector1[0] * vector2[0] * vector3[1] -
		point1[0] * vector1[2] * vector2[0] * vector3[1] + point2[0] * vector1[0] * vector2[2] * vector3[1] -
		point1[1] * vector1[0] * vector2[0] * vector3[2] + point2[1] * vector1[0] * vector2[0] * vector3[2] +
		point1[0] * vector1[1] * vector2[0] * vector3[2] - point2[0] * vector1[0] * vector2[1] * vector3[2]) /
		(-vector1[2] * vector2[1] * vector3[0] + vector1[1] * vector2[2] * vector3[0] + vector1[2] * vector2[0] * vector3[1] -
			vector1[0] * vector2[2] * vector3[1] - vector1[1] * vector2[0] * vector3[2] + vector1[0] * vector2[1] * vector3[2]));
	const auto y2 = -((-point1[2] * vector1[1] * vector2[1] * vector3[0] + point2[2] * vector1[1] * vector2[1] * vector3[0] +
		point1[1] * vector1[2] * vector2[1] * vector3[0] - point2[1] * vector1[1] * vector2[2] * vector3[0] -
		point2[1] * vector1[2] * vector2[0] * vector3[1] + point1[2] * vector1[0] * vector2[1] * vector3[1] -
		point2[2] * vector1[0] * vector2[1] * vector3[1] - point1[0] * vector1[2] * vector2[1] * vector3[1] +
		point2[0] * vector1[2] * vector2[1] * vector3[1] + point2[1] * vector1[0] * vector2[2] * vector3[1] +
		point2[1] * vector1[1] * vector2[0] * vector3[2] - point1[1] * vector1[0] * vector2[1] * vector3[2] +
		point1[0] * vector1[1] * vector2[1] * vector3[2] - point2[0] * vector1[1] * vector2[1] * vector3[2]) /
		(-vector1[2] * vector2[1] * vector3[0] + vector1[1] * vector2[2] * vector3[0] +
			vector1[2] * vector2[0] * vector3[1] - vector1[0] * vector2[2] * vector3[1] -
			vector1[1] * vector2[0] * vector3[2] + vector1[0] * vector2[1] * vector3[2]));
	const auto z2 = -((point2[2] * vector1[2] * vector2[1] * vector3[0] - point1[2] * vector1[1] * vector2[2] * vector3[0] +
		point1[1] * vector1[2] * vector2[2] * vector3[0] - point2[1] * vector1[2] * vector2[2] * vector3[0] -
		point2[2] * vector1[2] * vector2[0] * vector3[1] + point1[2] * vector1[0] * vector2[2] * vector3[1] -
		point1[0] * vector1[2] * vector2[2] * vector3[1] + point2[0] * vector1[2] * vector2[2] * vector3[1] +
		point2[2] * vector1[1] * vector2[0] * vector3[2] - point2[2] * vector1[0] * vector2[1] * vector3[2] -
		point1[1] * vector1[0] * vector2[2] * vector3[2] + point2[1] * vector1[0] * vector2[2] * vector3[2] +
		point1[0] * vector1[1] * vector2[2] * vector3[2] - point2[0] * vector1[1] * vector2[2] * vector3[2]) /
		(-vector1[2] * vector2[1] * vector3[0] + vector1[1] * vector2[2] * vector3[0] +
			vector1[2] * vector2[0] * vector3[1] - vector1[0] * vector2[2] * vector3[1] -
			vector1[1] * vector2[0] * vector3[2] + vector1[0] * vector2[1] * vector3[2]));

	result[0] = (x1 + x2) / 2;	result[1] = (y1 + y2) / 2; result[2] = (z1 + z2) / 2;
	return result;
}
}

//����������
//originPoint-δ��������ά��
//translationVector-ƽ������
//plane1-������һ�ķ���ax+by+cz+d=0��abcd�ĸ�����
//thickness-�����������
//n1,n2,n3-������.
//�������
template<typename _Ty>
point3d_t<_Ty> refractive_3d_reconstruction(
	point3d_t<_Ty> originPoint, 
	point3d_t<_Ty> translationVector, 
	std::array<_Ty, 4> plane1, 
	_Ty thickness, _Ty n1, _Ty n2, _Ty n3)
{
	MATRICE_USE_STD(pow);
	MATRICE_USE_STD(sin);
	MATRICE_USE_STD(cos);

	////����������1�ķ�����������2�ķ��̡�
	auto d2 = plane1[3] - thickness * pow(plane1[0] * plane1[0] + plane1[1] * plane1[1] + plane1[2] * plane1[2], 0.5);
	decltype(plane1) plane2 = { plane1[0], plane1[1], plane1[2], d2 };

	////������1��2�ķ���������һ�����
	auto N1 = detail::_Normalize(plane1[0], plane1[1], plane1[2]);
	auto N2 = detail::_Normalize(plane2[0], plane2[1], plane2[2]);

	////����1�ķ�������
	auto nlineL1 = detail::_Normalize(originPoint[0], originPoint[1], originPoint[2]);
	auto nlineR1 = detail::_Normalize(originPoint[0] + translationVector[0], originPoint[1] + translationVector[1], originPoint[2] + translationVector[2]);

	////�����1��������1�Ľ�������
	auto pointL1 = detail::_Line_and_plane_intersection({ 0.0, 0.0, 0.0 }, nlineL1, plane1);
	auto pointR1 = detail::_Line_and_plane_intersection({ -translationVector[0], -translationVector[1], -translationVector[2] }, nlineR1, plane1);

	//��������1�����
	auto angleL1 = detail::_Vector_angle(nlineL1, N1);
	auto angleR1 = detail::_Vector_angle(nlineR1, N1);

	//��������1�����
	auto angleL2 = asin(n1 * sin(angleL1) / n2);
	auto angleR2 = asin(n1 * sin(angleR1) / n2);

	//�����2�ķ�������
	auto nlineL2 = detail::_Get_refracted_vector(nlineL1, N1, n1, n2, angleL1, angleL2);
	auto nlineR2 = detail::_Get_refracted_vector(nlineR1, N1, n1, n2, angleR1, angleR2);
	nlineL2 = detail::_Normalize(nlineL2[0], nlineL2[1], nlineL2[2]);
	nlineR2 = detail::_Normalize(nlineR2[0], nlineR2[1], nlineR2[2]);

	//�����2��������2�Ľ�������
	auto pointL2 = detail::_Line_and_plane_intersection(pointL1, nlineL2, plane2);
	auto pointR2 = detail::_Line_and_plane_intersection(pointR1, nlineR2, plane2);

	//��������2�����
	const auto angleL3 = angleL2;
	const auto angleR3 = angleR2;

	//��������2�����
	auto angleL4 = asin(n1 * sin(angleL3) / n2);
	auto angleR4 = asin(n1 * sin(angleR3) / n2);

	//�����3�ķ�������
	const auto nlineL3 = detail::_Get_refracted_vector(nlineL2, N2, n2, n3, angleL3, angleL4);
	const auto nlineR3 = detail::_Get_refracted_vector(nlineR2, N2, n2, n3, angleR3, angleR4);

	////�����3�Ĺ����߷�������
	//auto nCommonVerticalLine3 = detail::_Cross(nlineL3, nlineR3);
	
	//�����3�Ľ�������
	return detail::_Intersection_of_lines(pointL2, nlineL3, pointR2, nlineR3);
}

MATRICE_ALG_END(vision)