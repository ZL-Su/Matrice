/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#include "../io.hpp"
#include "../../util/_exception.h"
#include "../../util/_std_wrapper.h"

DGE_MATRICE_BEGIN namespace io { _DETAIL_BEGIN

/**
 *\brief <class> _Data_loader_base </class>
 */
template<typename _Derived> class _Data_loader_base {
	using _Myt = _Data_loader_base;
	using _Mydt = _Derived;
public:
	using directory_type = _Dir_impl<folder_tag>;
protected:
	struct _Loader_iterator {
		_Loader_iterator(const std::add_pointer_t<_Myt> _This)
			:_Mythis(_This), _Mybatchs(_Mythis->batch_size()),
			_Mypos(_Mythis->pos()) {}

		MATRICE_HOST_FINL auto operator*() const {
			directory_type::container _Ret;
			for (auto _Idx = 0; _Idx < _Mybatchs; ++_Idx) {
				_Ret.emplace_back(_Mythis->directory()[_Idx] + 
					_Mythis->file_names(_Idx)[_Mypos]);
			}
			return forward<decltype(_Ret)>(_Ret);
		}
		MATRICE_HOST_FINL _Loader_iterator operator++() {
			_Mythis->pos() += 1;
			return (*this);
		}
		MATRICE_HOST_FINL _Loader_iterator operator++(int) {
			auto _Tmp = *this;
			_Mythis->pos() += 1;
			return (_Tmp);
		}
	private:
		size_t& _Mypos;
		size_t  _Mybatchs = 0;
		std::add_pointer_t<_Myt> _Mythis = nullptr;
		directory_type::container::iterator _Myitr;
	};
public:
	using iterator = _Loader_iterator;

	/**
	 * \set loader iterator to begin (zero) pos
	 */
	MATRICE_HOST_INL bool begin() const {
		return (_Mypos = -1);
	}
	/**
	 * \set loader iterator to reverse begin (_Mydepth) pos
	 */
	MATRICE_HOST_INL bool rbegin() const {
		return (_Mypos = _Mydepth);
	}
	/**
	 * \check if loader iterator meets the upper bound
	 */
	MATRICE_HOST_INL bool end() const {
		return (_Mypos >= _Mydepth);
	}
	/**
	 * \check if loader iterator meets the lower bound
	 */
	MATRICE_HOST_INL bool rend() const {
		return (_Mypos < 0);
	}
	/**
	 * \brief Move loader iterator _Off steps and return the loader
	 */
	MATRICE_HOST_INL decltype(auto) shift(index_t _Off) const {
		_Mypos += _Off;
#ifdef _DEBUG
		_COND_EXCEPTION(end() || rend(),
			"_Off over range of loader depth")
#endif
			return (*this);
	}

	/**
	 * \return iterator pos
	 */
	MATRICE_HOST_FINL auto& pos() { return (_Mypos); }
	MATRICE_HOST_FINL const auto& pos() const { return (_Mypos); }

	MATRICE_HOST_INL directory_type& directory() {
		return (_Mydir);
	}
	MATRICE_HOST_INL const directory_type& directory() const {
		return (_Mydir);
	}
	MATRICE_HOST_FINL size_t batch_size() const {
		return (_Mydir.size());
	}

	/**
	 * \return all file names in currenct work path
	 */
	MATRICE_HOST_INL std::vector<directory_type::container>& file_names() {
		return (_Mynames);
	}
	MATRICE_HOST_INL const std::vector<directory_type::container>& file_names() const {
		return (_Mynames);
	}

	/**
	 * \return all file names in _Idx-th subfolder for currenct work path
	 */
	MATRICE_HOST_FINL directory_type::container& file_names(size_t _Idx) {
#ifdef _DEBUG
		DGELOM_CHECK(_Idx < _Mynames.size(), "_Idx over range of field ::_Mynames.");
#endif // _DEBUG

		return (_Mynames)[_Idx];
	}
	MATRICE_HOST_FINL const directory_type::container& file_names(size_t _Idx) const {
#ifdef _DEBUG
		DGELOM_CHECK(_Idx < _Mynames.size(), "_Idx over range of field ::_Mynames.");
#endif // _DEBUG
		return (_Mynames)[_Idx];
	}

	/**
	 * \return depth (the number of files) of the loader
	 */
	MATRICE_HOST_FINL index_t depth() const {
		return (_Mydepth);
	}

	/**
	 *\brief Operator to check if all files are loaded
	 */
	MATRICE_HOST_INL operator bool() const {
		return (_Mypos < _Mydepth);
	}
protected:
	MATRICE_HOST_INL void _Collect_fnames() {
		using collector_t = _Collector<file_tag>;
		if (_Mydir.size() == 0) {//no subfolder(s)
			_Mynames.emplace_back(collector_t::get(_Mydir[0]));
			_Mydepth = _Mynames.front().size();
			return;
		}
		_Mynames.resize(_Mydir.size());
		for (const auto _Idx : range(0, _Mynames.size())) {
			_Mynames[_Idx] = collector_t::get(_Mydir[_Idx]);
			if (auto _Cnt = _Mynames[_Idx].size(); _Cnt < _Mydepth)
				_Mydepth = _Cnt;
		}
	}

	directory_type _Mydir;
	mutable index_t _Mypos = 0;
	index_t _Mydepth = std::numeric_limits<index_t>::max();
	std::vector<directory_type::container> _Mynames;
	uint32_t m_nworkes = 1; //indicates the number of work threads being used
};

/**
 * \Specified to load data from a given directory
 */
template<typename _Ty>
class _Data_loader_impl<_Ty, loader_tag> {
	using _Mydir_type = _Dir_impl<folder_tag>;
	using _Myt = _Data_loader_impl;
	struct _Loader_iterator {
		_Loader_iterator(const std::add_pointer_t<_Myt> _This)
			:_Mythis(_This), _Mybatchs(_Mythis->batch_size()),
			_Mypos(_Mythis->pos()) {}

		MATRICE_HOST_FINL auto operator*() const {
			_Mydir_type::container _Ret;
			for (auto _Idx = 0; _Idx < _Mybatchs; ++_Idx) {
				_Ret.emplace_back(_Mythis->directory()[_Idx] + _Mythis->file_names(_Idx)[_Mypos]);
			}
			return forward<decltype(_Ret)>(_Ret);
		}
		MATRICE_HOST_FINL _Loader_iterator operator++() {
			_Mythis->pos() += 1;
			return (*this);
		}
		MATRICE_HOST_FINL _Loader_iterator operator++(int) {
			auto _Tmp = *this;
			_Mythis->pos() += 1;
			return (_Tmp);
		}
	private:
		size_t& _Mypos;
		size_t _Mybatchs = 0;
		std::add_pointer_t<_Myt> _Mythis = nullptr;
		_Mydir_type::container::iterator _Myitr;
	};
public:
	using category = loader_tag;
	using dir_type = _Mydir_type;
	using iterator = _Loader_iterator;
	using data_type = Matrix<_Ty>;

	_Data_loader_impl(const _Mydir_type& _Dir)
		: _Mydir(_Dir) {
		_Collect_fnames();
	}
	template<typename _Fn>
	_Data_loader_impl(const _Mydir_type& _Dir, _Fn&& _Op)
		: _Mydir(_Dir), _Myloader(std::forward<_Fn>(_Op)) {
		_Collect_fnames();
	}
	_Data_loader_impl(_Mydir_type&& _Dir)
		: _Mydir(std::forward<_Mydir_type>(_Dir)) {
		_Collect_fnames();
	}
	template<typename _Fn>
	_Data_loader_impl(_Mydir_type&& _Dir, _Fn&& _Op)
		: _Mydir(std::forward<_Mydir_type>(_Dir)),
		_Myloader(std::forward<_Fn>(_Op)) {
		_Collect_fnames();
	}

	/**
	 * \set loader iterator to begin (zero) pos
	 */
	MATRICE_HOST_INL bool begin() const {
		return (_Mypos = -1);
	}
	/**
	 * \set loader iterator to reverse begin (_Mydepth) pos
	 */
	MATRICE_HOST_INL bool rbegin() const {
		return (_Mypos = _Mydepth);
	}
	/**
	 * \check if loader iterator meets the upper bound
	 */
	MATRICE_HOST_INL bool end() const {
		return (_Mypos >= _Mydepth);
	}
	/**
	 * \check if loader iterator meets the lower bound
	 */
	MATRICE_HOST_INL bool rend() const {
		return (_Mypos < 0);
	}
	/**
	 * \brief Move loader iterator _Off steps and return the loader
	 */
	MATRICE_HOST_INL decltype(auto) shift(index_t _Off) const {
		_Mypos += _Off;
#ifdef _DEBUG
		_COND_EXCEPTION(end()||rend(),
			"_Off over range of loader depth")
#endif
			return (*this);
	}

	/**
	 * \forward iterate to retrieve data paths with common loader.
	 */
	MATRICE_HOST_INL auto forward() const {
		std::vector<data_type> _Data;
		const auto _Size = _Mydir.size() == 0 ? 1 : _Mydir.size();
		for (const auto& _Idx : range(0, _Size)) {
			const auto& _Names = _Mynames[_Idx];
#ifdef _DEBUG
			DGELOM_CHECK(_Mypos<_Names.size(), "file list subscript out of range.");
#endif
			_Data.emplace_back(_Myloader(_Mydir[_Idx] + _Names[_Mypos]));
		}
		_Mypos++;
		return std::forward<decltype(_Data)>(_Data);
	}
	/**
	 * \forward iterate to retrieve data paths with given _Loader. The data type depends on the return-type of _Loader.
	 */
	template<typename _Fn>
	MATRICE_HOST_INL auto forward(_Fn&& _Loader) const {
		auto _Op = [&](auto i){return _Loader(_Mydir[i]+_Mynames[i][_Mypos]); };
		const auto _First = _Op(0);
		std::vector<remove_all_t<decltype(_First)>> _Data;
		_Data.emplace_back(_First);
		for (const auto& _Idx : range(1, _Mydir.size()-1)) {
			_Data.emplace_back(_Op(_Idx));
		}
		_Mypos++;
		return std::forward<decltype(_Data)>(_Data);
	}
	/**
	 * \reverse iterate to retrieve data paths
	 */
	MATRICE_HOST_INL auto reverse() const {
		_Mypos--;
		std::vector<data_type> _Data;

		return std::forward<decltype(_Data)>(_Data);
	}

	/**
	 *\brief load the first file
	 */
	MATRICE_HOST_INL auto front() const {
		std::vector<data_type> _Data;
		const auto _Size = _Mydir.size() == 0 ? 1 : _Mydir.size();
		for (const auto& _Idx : range(0, _Size)) {
			const auto& _Name = _Mynames[_Idx].front();
			_Data.emplace_back(_Myloader(_Mydir[_Idx] + _Name));
		}
		return std::forward<decltype(_Data)>(_Data);
	}
	/**
	 *\brief load the last file
	 */
	MATRICE_HOST_INL auto back() const {
		std::vector<data_type> _Data;
		const auto _Size = _Mydir.size() == 0 ? 1 : _Mydir.size();
		for (const auto& _Idx : range(0, _Size)) {
			const auto& _Name = _Mynames[_Idx].back();
			_Data.emplace_back(_Myloader(_Mydir[_Idx] + _Name));
		}
		return std::forward<decltype(_Data)>(_Data);
	}

	/**
	 *\brief load file at i-th position
	 *\param [i] index of a file name
	 */
	MATRICE_HOST_INL auto at(size_t i) const {
		std::vector<data_type> _Data;
		const auto _Size = _Mydir.size()==0?1:_Mydir.size();
		for (const auto& _Idx : range(0, _Size)) {
			const auto& _Names = _Mynames[_Idx];
#ifdef _DEBUG
			DGELOM_CHECK(i < _Names.size(), "file list subscript out of range.");
#endif
			_Data.emplace_back(_Myloader(_Mydir[_Idx]+_Names[i]));
		}
		return std::forward<decltype(_Data)>(_Data);
	}

	/**
	 * \return iterator pos
	 */
	MATRICE_HOST_FINL auto& pos() { return (_Mypos); }
	MATRICE_HOST_FINL const auto& pos() const { return (_Mypos); }

	MATRICE_HOST_INL dir_type& directory() {
		return (_Mydir);
	}
	MATRICE_HOST_INL const dir_type& directory() const {
		return (_Mydir);
	}
	MATRICE_HOST_FINL size_t batch_size() const {
		return (_Mydir.size());
	}

	/**
	 * \return all file names in currenct work path
	 */
	MATRICE_HOST_INL std::vector<_Mydir_type::container>& file_names() {
		return (_Mynames);
	}
	MATRICE_HOST_INL const std::vector<_Mydir_type::container>& file_names() const {
		return (_Mynames);
	}

	/**
	 * \return all file names in _Idx-th subfolder for currenct work path
	 */
	MATRICE_HOST_FINL _Mydir_type::container& file_names(size_t _Idx) {
#ifdef _DEBUG
		DGELOM_CHECK(_Idx<_Mynames.size(), "_Idx over range of field ::_Mynames.");
#endif // _DEBUG

		return (_Mynames)[_Idx];
	}
	MATRICE_HOST_FINL const _Mydir_type::container& file_names(size_t _Idx) const {
#ifdef _DEBUG
		DGELOM_CHECK(_Idx<_Mynames.size(), "_Idx over range of field ::_Mynames.");
#endif // _DEBUG
		return (_Mynames)[_Idx];
	}

	/**
	 * \return depth (the number of files) of the loader
	 */
	MATRICE_HOST_FINL auto depth() const {
		return (_Mydepth);
	}

	/**
	 *\brief Operator to check if all files are loaded
	 */
	MATRICE_HOST_INL operator bool() const {
		return (_Mypos < _Mydepth);
	}

private:
	MATRICE_HOST_INL void _Collect_fnames() {
		using collector_t = _Collector<file_tag>;
		if (_Mydir.size() == 0) {//no subfolder(s)
			_Mynames.emplace_back(collector_t::get(_Mydir[0]));
			_Mydepth = _Mynames.front().size();
			return;
		}
		_Mynames.resize(_Mydir.size());
		for (const auto _Idx : range(0, _Mynames.size())) {
			_Mynames[_Idx] = collector_t::get(_Mydir[_Idx]);
			if (auto _Cnt = _Mynames[_Idx].size(); _Cnt < _Mydepth)
				_Mydepth = _Cnt;
		}
	}

	_Mydir_type _Mydir;
	mutable index_t _Mypos = 0;
	index_t _Mydepth = std::numeric_limits<index_t>::max();
	std::vector<_Mydir_type::container> _Mynames;
	std::function<data_type(_Mydir_type::value_type&&)> _Myloader;
};

template<typename _Ty>
struct _Loader_impl<_Ty, loader_tag::tiff> {
	using value_type = _Ty;
	using category = loader_tag::tiff;

	MATRICE_HOST_INL decltype(auto) operator()(std::string path) {
		return read_tiff_file<value_type>(path.c_str());
	}
};
_DETAIL_END } DGE_MATRICE_END