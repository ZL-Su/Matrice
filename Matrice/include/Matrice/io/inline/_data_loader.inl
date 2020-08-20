/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include "../io.hpp"
#include "util/_exception.h"
#include "util/_std_wrapper.h"

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
#ifdef MATRICE_DEBUG
		MATRICE_COND_EXCEPTION(end() || rend(),
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
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_Idx < _Mynames.size(), "_Idx over range of field ::_Mynames.");
#endif // _DEBUG

		return (_Mynames)[_Idx];
	}
	MATRICE_HOST_FINL const directory_type::container& file_names(size_t _Idx) const {
#ifdef MATRICE_DEBUG
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
		: _Mydir(_Dir), _Myloader(_Op) {
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
	MATRICE_HOST_INL bool begin() const noexcept {
		return (_Mypos = -1);
	}
	/**
	 * \set loader iterator to reverse begin (_Mydepth) pos
	 */
	MATRICE_HOST_INL bool rbegin() const noexcept {
		return (_Mypos = _Mydepth);
	}
	/**
	 * \check if loader iterator meets the upper bound
	 */
	MATRICE_HOST_INL bool end() const noexcept {
		return (_Mypos >= _Mydepth);
	}
	/**
	 * \check if loader iterator meets the lower bound
	 */
	MATRICE_HOST_INL bool rend() const noexcept {
		return (_Mypos < 0);
	}
	/**
	 * \brief Move loader iterator _Off steps and return the loader
	 */
	MATRICE_HOST_INL decltype(auto) shift(index_t _Off) const {
		_Mypos += _Off;
#ifdef MATRICE_DEBUG
		MATRICE_COND_EXCEPTION(end()||rend(),
			"_Off over range of loader depth")
#endif
			return (*this);
	}

	/**
	 * \brief Load current file in each of subfolders in one folder
	 */
	MATRICE_HOST_INL auto forward() const {
		std::vector<data_type> _Data;
		const auto _Size = _Mydir.size() == 0 ? 1 : _Mydir.size();
		for (const auto& _Idx : range(0, _Size)) {
			const auto& _Names = _Mynames[_Idx];
#ifdef MATRICE_DEBUG
			DGELOM_CHECK(_Mypos<_Names.size(), "file list subscript out of range.");
#endif
			auto _File = _Myloader(_Mydir[_Idx] + _Names[_Mypos]);
			_Cond_push_file(_Data, _File);
		}
		_Mypos++;
		return std::forward<decltype(_Data)>(_Data);
	}
	/**
	 *\forward Iterate to retrieve data paths with given _Loader. 
		The data type depends on the return-type of _Loader.
	 *\note	Load the current file in each of subfolders in one folder
	 */
	template<typename _Fn>
	MATRICE_HOST_INL auto forward(_Fn&& _Loader) const {
		auto _Op = [&](auto i){return _Loader(_Mydir[i]+_Mynames[i][_Mypos]); };
		const auto _First = _Op(0);
		std::vector<remove_all_t<decltype(_First)>> _Data;
		_Data.emplace_back(_First);
		for (const auto& _Idx : range(1, _Mydir.size()-1)) {
			auto _File = _Op(_Idx);
			static_assert(is_same_v<data_type, remove_all_t<decltype(_File)>>, "Returned type of the user defined loader [_Fn] is not same with the [data_type].");
			if(_File.size()!=0) _Data.push_back(_File);
			else continue;
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
	 *\brief Load the first file in each of subfolders in one folder
	 */
	MATRICE_HOST_INL auto front() const {
		std::vector<data_type> _Data;
		const auto _Size = _Mydir.size() == 0 ? 1 : _Mydir.size();
		for (const auto& _Idx : range(0, _Size)) {
			const auto& _Name = _Mynames[_Idx].front();
			auto _File = _Myloader(_Mydir[_Idx] + _Name);
			_Cond_push_file(_Data, _File);
		}
		return std::forward<decltype(_Data)>(_Data);
	}
	/**
	 *\brief Load the last file in each of subfolders in one folder.
	 */
	MATRICE_HOST_INL auto back() const {
		std::vector<data_type> _Data;
		const auto _Size = _Mydir.size() == 0 ? 1 : _Mydir.size();
		for (const auto& _Idx : range(0, _Size)) {
			const auto& _Name = _Mynames[_Idx].back();
			auto _File = _Myloader(_Mydir[_Idx] + _Name);
			_Cond_push_file(_Data, _File);
		}
		return std::forward<decltype(_Data)>(_Data);
	}

	/**
	 *\brief Load the i-th file in each of subfolders in one folder
	 *\param [i] index of the file to be loaded
	 */
	MATRICE_HOST_INL auto at(size_t i) const {
		std::vector<data_type> _Data;
		const auto _Size = _Mydir.size()==0?1:_Mydir.size();
		for (const auto& _Idx : range(0, _Size)) {
			const auto& _Names = _Mynames[_Idx];
#ifdef MATRICE_DEBUG
			DGELOM_CHECK(i < _Names.size(), "file list subscript out of range.");
#endif
			auto _File = _Myloader(_Mydir[_Idx] + _Names[i]);
			_Cond_push_file(_Data, _File);
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
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_Idx < _Mynames.size(), "_Idx over range of field ::_Mynames.");
#endif

		return (_Mynames)[_Idx];
	}
	MATRICE_HOST_FINL const _Mydir_type::container& file_names(size_t _Idx) const {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_Idx < _Mynames.size(), "_Idx over range of field ::_Mynames.");
#endif
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

	template<class _Dty, class _Fty>
	MATRICE_HOST_INL void _Cond_push_file(_Dty& _Data, const _Fty& _File)const {
		using file_type = remove_all_t<decltype(_File)>;
		if (_File.size()) {
			if constexpr (is_same_v<file_type, image_instance>) {
				auto _Mat = _File.matrix<data_type::value_t>();
				_Data.push_back(_Mat);
			}
			else {
				_Data.push_back(_File);
			}
		}
	}

	_Mydir_type _Mydir;
	mutable index_t _Mypos = 0;
	index_t _Mydepth = std::numeric_limits<index_t>::max();
	std::vector<_Mydir_type::container> _Mynames;
	std::function<image_instance(_Mydir_type::value_type&&)> _Myloader;
};

template<typename _Ty>
struct _Loader_impl<_Ty, loader_tag::tiff> {
	using value_type = _Ty;
	using category = loader_tag::tiff;

	MATRICE_HOST_INL image_instance operator()(std::string path) {
		if( string_helper::split(path, '.').back() == "tif" ||
			string_helper::split(path, '.').back() == "tiff" )
			return read_tiff_file(path.c_str());
		else return image_instance{};
	}
};
_DETAIL_END } DGE_MATRICE_END