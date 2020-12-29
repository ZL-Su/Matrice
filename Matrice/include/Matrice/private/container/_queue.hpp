/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2020, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/
#pragma once
#include "util/_macros.h"
#include <vector>
#include <algorithm>

DGE_MATRICE_BEGIN
template<typename _Ty, typename _Pre = std::less<_Ty>>
class priority_queue {
public:
    using value_type = _Ty;
    using container = std::vector<value_type>;
    using iterator = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    priority_queue() : _Mymax_size(0) {}
    priority_queue(size_t max_size) : _Mymax_size(max_size) {}

    const_iterator begin() const { return _Myc.begin(); }
    iterator begin() { return _Myc.begin(); }
    const_iterator end() const { return _Myc.end(); }
    iterator end() { return _Myc.end(); }

    MATRICE_HOST_INL void push(const value_type& x) {
        MATRICE_USE_STD(make_heap);
        if (_Myc.size() == _Mymax_size) {
            MATRICE_USE_STD(min_element);
            auto it = min_element(_Myc.begin(), _Myc.end(), _Mycmp);
            if (*it < x) {
                *it = x;
                make_heap(_Myc.begin(), _Myc.end(), _Mycmp);
            }
        }
        else {
            _Myc.push_back(x);
            make_heap(_Myc.begin(), _Myc.end(), _Mycmp);
        }
    }

    MATRICE_HOST_INL void pop() {
        if (_Myc.empty()) return;

        MATRICE_USE_STD(pop_heap);
        pop_heap(_Myc.begin(), _Myc.end(), _Mycmp);
        _Myc.pop_back();
    }

    MATRICE_HOST_INL const value_type& top() const {
        return _Myc.front();
    }

    MATRICE_HOST_INL const bool empty() const noexcept {
        return _Myc.empty();
    }

    MATRICE_HOST_INL const size_t size() const noexcept {
        return _Myc.size();
    }

    MATRICE_HOST_INL void enlarge_max_size(size_t max_size) {
        if (_Mymax_size < max_size)
            _Mymax_size = max_size;
    }

protected:
    container _Myc;
    size_t _Mymax_size;
    _Pre _Mycmp;

private:
    // heap allocation is not allowed
    void* operator new (size_t) = delete;
    void* operator new[] (size_t) = delete;
    void operator delete (void*) = delete;
    void operator delete[] (void*) = delete;
};
DGE_MATRICE_END