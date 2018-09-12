#include "../interpolation.h"
#pragma once

MATRICE_ALGS_BEGIN
namespace details { namespace internals {
template<size_t _Values> struct static_size
{
	enum { value = 
		(_Values & bilinear == bilinear) ? 3 :
		(_Values & bcspline == bcspline) ? 4 :
		(_Values & bqspline == bqspline) ? 6 :
		(_Values & bsspline == bsspline) ? 8 : _Values
	};
};

template<typename _Solver_type>
struct gn_solver_traits { 
	using type = typename _Solver_type::value_t;
	enum { options = _Solver_type::options };
};

}}
MATRICE_ALGS_END
