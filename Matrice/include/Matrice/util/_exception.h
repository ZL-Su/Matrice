/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include <exception>
#include <string>
#include <vector>
#include "_macros.h"

#if defined(_MSC_VER) && _MSC_VER <= 1900
#define __func__ __FUNCTION__
#endif

DGE_MATRICE_BEGIN _DETAIL_BEGIN
/**
 * \record source code location
 * \Example: _Source_location _Loc{__func__, __FILE__, __LINE__}
 */
struct _Source_location {
	const char* _Func = nullptr;
	const char* _File = nullptr;
	long        _Line;
};
#define MATRICE_EXCLOC {__func__, __FILE__, __LINE__}

/**
 * \exception process
 */
struct _Exception_wrapper
{
	class error : public std::exception {
		using _Myt = error;
		using _Mybase = std::exception;
		using _Myloc_type = _Source_location;
	public:
		using msg_type = std::string;
		using msg_list = std::vector<msg_type>;

		error(const _Myloc_type& _Loc, const msg_type& _Msg="None")
			:_Mymsg(_Msg), _Myloc(_Loc) {}

		MATRICE_HOST_INL auto message() const {

		}

	private:
		const void* _Mycaller = nullptr;
		msg_type _Mymsg;
		msg_list _Mymsgstack;
		_Myloc_type _Myloc;
	};
};
_DETAIL_END
using exception = detail::_Exception_wrapper;
DGE_MATRICE_END