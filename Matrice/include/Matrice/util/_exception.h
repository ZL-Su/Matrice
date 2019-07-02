/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include <functional>
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
	const char* _func = nullptr;
	const char* _file = nullptr;
	long        _line;
};
#define __exceploc__ {__func__, __FILE__, __LINE__}

/**
 * \exception process
 */
struct _Exception_wrapper
{
	using msg_type = std::string;
	using msg_list = std::vector<msg_type>;
	using loc_type = _Source_location;

	class error : public std::exception {
		using _Myt = error;
		using _Mybase = std::exception;
	public:
		using _Mybase::exception;
		error(const loc_type& _Loc, const msg_type& _Msg="None") noexcept
			:_Myloc(_Loc), _Mybase(_Msg.c_str()) {}

		/**
		 * \return exception location
		 */
		MATRICE_HOST_INL loc_type location() const noexcept {
			return (_Myloc);
		}

	private:
		const void* _Mycaller = nullptr;
		msg_type _Mymsg = this->what();
		msg_list _Mymsgstack;
		loc_type _Myloc;
	};

	using handler = std::function<void(const msg_type&, const char*)>;
	static MATRICE_HOST_INL auto warning(loc_type _Loc, msg_type _Msg);
};
_DETAIL_END
using exception = detail::_Exception_wrapper;
DGE_MATRICE_END

#define DGELOM_ERROR(...) \
	throw ::dgelom::exception::error(__exceploc__, \
			::dgelom::exception::msg_type(__VA_ARGS__))

#define DGELOM_CHECK(_Cond, ...) \
	if(!(_Cond)) { \
		DGELOM_ERROR(::dgelom::exception::msg_type(__VA_ARGS__)); \
	}