#pragma once

#include "interpolation\_interpolation.h"

MATRICE_NAMESPACE_BEGIN_
enum {
	bilinear = algs::INTERP | algs::BILINEAR,
	bcspline = algs::INTERP | algs::BSPLINE | algs::BICUBIC,
	bqspline = algs::INTERP | algs::BSPLINE | algs::BIQUINTIC,
	bsspline = algs::INTERP | algs::BSPLINE | algs::BISEPTIC,
};
template<typename _Ty, size_t _Options = bcspline, typename = std::enable_if_t<std::is_scalar_v<_Ty>>>
using interpolation = algs::Interpolation<_Ty, _Options>;
_MATRICE_NAMESPACE_END
