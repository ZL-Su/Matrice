/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#include "../io.hpp"

DGE_MATRICE_BEGIN namespace io { _DETAIL_BEGIN
/**
 * \Specified to collect folder(s) in a given path 
 */
template<> struct _Collector<folder_tag> {
	MATRICE_HOST_FINL static auto get(std::string&& path) {
		std::vector<std::string> _Ret;
		if (fs::path _Path(path); fs::exists(_Path)) {
			fs::directory_iterator _End;
			for (decltype(_End)_Begin(_Path);_Begin != _End; ++_Begin) {
				if (fs::is_directory(_Begin->status())) {
					_Ret.emplace_back(_Begin->path().string().substr(path.size()));
				}
			}
			_Ret.shrink_to_fit();
		}
		return forward<decltype(_Ret)>(_Ret);
	}
};
/**
 * \Specified to collect file name(s) in a given path
 */
template<> struct _Collector<file_tag> {
	MATRICE_HOST_FINL static auto get(std::string&& path) {
		std::vector<std::string> _Ret;
		if (fs::path _Path(path); fs::exists(_Path)) {
			fs::directory_iterator _End;
			for (decltype(_End) _Begin(_Path); _Begin != _End; ++_Begin) {
				if (fs::is_regular_file(_Begin->status())) {
					_Ret.emplace_back(_Begin->path().string().substr(path.size()));
				}
			}
			_Ret.shrink_to_fit();
		}
		return std::forward<decltype(_Ret)>(_Ret);
	}
};

/**
 * \Specified to construct directory for a given _Root
 */
template<> class _Dir_impl<folder_tag> {
public:
	using category = folder_tag;
	using value_type = std::string;
	using path_type = fs::path;
	using container = std::vector<value_type>;
	/**
	 * \_Root must has form "./folder_name/.../folder_name" or "/folder_name/.../folder_name"
	 */
	_Dir_impl(value_type _Root, path_type&& _Workdir = fs::current_path())
		:_Mypath(_Workdir) {
		if (_Root[0] == '.') _Mypath.concat(_Root.begin() + 1, _Root.end());
		else _Mypath.concat(_Root.begin(), _Root.end());

		_Mysubfolders = _Collector<folder_tag>::get(_Mypath.string());
	}

	/**
	 * \return full path of _Idx-th folder
	 */
	MATRICE_HOST_INL auto operator[] (size_t _Idx) const {
		return (_Mypath.string() + folder(_Idx));
	}

	/**
	 * \return _Idx-th folder name under work-path _Mypath
	 */
	MATRICE_HOST_INL value_type folder(size_t _Idx) const {
		if (_Mysubfolders.size() == 0) return value_type();
		else {
#ifdef _DEBUG
			DGELOM_CHECK(_Idx < _Mysubfolders.size(), "_Idx over range of _Mysubfolders")
#endif
			return (_Mysubfolders[_Idx]);
		}
	}

	/**
	 * \return path of current directory
	 */
	MATRICE_HOST_INL value_type path() const {
		return forward<value_type>(_Mypath.string());
	}

	/**
	 * \return size (number of subfolders) of current directory
	 */
	MATRICE_HOST_INL const auto size() const {
		return (_Mysubfolders.size());
	}

private:
	path_type _Mypath;
	container _Mysubfolders;
};

_DETAIL_END } DGE_MATRICE_END