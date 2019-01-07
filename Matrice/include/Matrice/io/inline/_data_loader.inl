/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#include "../io.h"
#include "../../util/_exception.h"

DGE_MATRICE_BEGIN namespace io { _DETAIL_BEGIN

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
			return std::forward<decltype(_Ret)>(_Ret);
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
		std::size_t& _Mypos;
		std::size_t _Mybatchs = 0;
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
	 * \move loader iterator _Off steps
	 */
	MATRICE_HOST_INL auto shift(index_t _Off) const {
		_Mypos += _Off;
#ifdef _DEBUG
		_COND_EXCEPTION(end() || rend(), "_Off over range of loader depth")
#endif
			return (_Mypos);
	}

	/**
	 * \forward iterate to retrieve data paths
	 */
	MATRICE_HOST_INL auto forward() const {
		_Mypos++;
		std::vector<data_type> _Data;
		for (const auto& _Idx : range(0, _Mydir.size())) {
			_Data.emplace_back(_Myloader(_Mydir[_Idx] + _Mynames[_Idx][_Mypos]));
		}
		return std::forward<decltype(_Data)>(_Data);
	}
	/**
	 * \forward iterate to retrieve data paths with given _Loader
	 */
	template<typename _Fn>
	MATRICE_HOST_INL auto forward(_Fn&& _Loader) const {
		_Mypos++;
		std::vector<data_type> _Data;
		for (const auto& _Idx : range(0, _Mydir.size())) {
			_Data.emplace_back(_Loader(_Mydir[_Idx] + _Mynames[_Idx][_Mypos]));
		}
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
	MATRICE_HOST_FINL std::size_t batch_size() const {
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
	MATRICE_HOST_FINL _Mydir_type::container& file_names(std::size_t _Idx) {
#ifdef _DEBUG
		DGELOM_CHECK(_Idx < _Mynames.size(), "_Idx over range of field ::_Mynames.");
#endif // _DEBUG

		return (_Mynames)[_Idx];
	}
	MATRICE_HOST_FINL const _Mydir_type::container& file_names(std::size_t _Idx) const {
#ifdef _DEBUG
		DGELOM_CHECK(_Idx < _Mynames.size(), "_Idx over range of field ::_Mynames.");
#endif // _DEBUG
		return (_Mynames)[_Idx];
	}

	/**
	 * \return depth (the number of files) of the loader
	 */
	MATRICE_HOST_FINL auto depth() const {
		return (_Mydepth);
	}

private:
	MATRICE_HOST_INL void _Collect_fnames() {
		_Mynames.resize(_Mydir.size());
		for (const auto _Idx : range(0, _Mynames.size())) {
			_Mynames[_Idx] = _Collector<file_tag>::get(_Mydir[_Idx]);
			if (auto _Cnt = _Mynames[_Idx].size(); _Cnt < _Mydepth)
				std::swap(_Mydepth, _Cnt);
		}
	}

	_Mydir_type _Mydir;
	mutable index_t _Mypos = -1;
	std::vector<_Mydir_type::container> _Mynames;
	std::function<data_type(_Mydir_type::value_type&&)> _Myloader;
	std::size_t _Mydepth = std::numeric_limits<std::size_t>::max();
};

_DETAIL_END } DGE_MATRICE_END