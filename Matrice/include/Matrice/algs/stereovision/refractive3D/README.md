# CODE NOTICE
The first release of "Refractive 3D Reconstruction" is in header file "_refractive3d_funcs.h", where the function template:
```
template<typename _Ty>
vec3d_t<_Ty> refractive_3d_reconstruction(vec3d_t<_Ty>, vec3d_t<_Ty>, vec4d_t<_Ty>, _Ty...)
```
is the driver for performing the 3D estimation in step by step. However, it is not cache-friendly.

A more efficient and portable version is being released in the file "_refractive3d.h", with the following class template:
```
template<typename _Ty> requires is_floating_point_v<_Ty>
class _Refractive_reconstruction;
```
Of course, user friendly interface is provided with the alias template as:
```
template<typename _Ty>
using refractive_reconstruction = detail::_Refractive_reconstruction;
```
In near future, the tutoriul of this module will be published here.
