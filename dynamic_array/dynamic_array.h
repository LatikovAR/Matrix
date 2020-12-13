#pragma once

#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <new>

#include "mem_storage/mem_storage.h"

using namespace mem_storage;

namespace dyn_arr {

template <typename T>
class dynamic_array final:
        protected Vector_Storage<T>
{
private:
    using Vector_Storage<T>::size_;
    using Vector_Storage<T>::used_;
    using Vector_Storage<T>::data_;
    using Vector_Storage<T>::swap;
public:
    dynamic_array(size_t size = 0);

    dynamic_array(const std::vector<T>& rhs);

    dynamic_array(const dynamic_array& rhs);

    dynamic_array<T>& operator= (const dynamic_array& rhs)&;

    size_t size() const { return size_; }

    T& operator[](size_t i);

    const T& operator[](size_t i) const;
};

template <typename T>
dynamic_array<T>::dynamic_array(size_t size):
    Vector_Storage<T>(size)
{
    for(size_t i = 0; i < size; ++i) {
        new (&(data_[i])) T();
        ++used_;
    }
}

template <typename T>
dynamic_array<T>::dynamic_array(const std::vector<T>& rhs):
    Vector_Storage<T>(rhs.size())
{
    for(size_t i = 0; i < size_; ++i) {
        new (&(data_[i])) T(rhs[i]);
        ++used_;
    }
}

template <typename T>
dynamic_array<T>::dynamic_array(const dynamic_array& rhs):
    Vector_Storage<T> (rhs.size())
{
    for(size_t i = 0; i < size_; ++i) {
        new (&(data_[i])) T(rhs[i]);
        ++used_;
    }
}

template <typename T>
dynamic_array<T>& dynamic_array<T>::operator= (const dynamic_array& rhs)& {
    dynamic_array<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T>
T& dynamic_array<T>::operator[](size_t i) {
    if(i >= size_) {
        throw std::out_of_range("invalid dyn_array iterator");
    }
    return data_[i];
}

template <typename T>
const T& dynamic_array<T>::operator[](size_t i) const {
    if(i >= size_) {
        throw std::out_of_range("invalid dyn_array iterator");
    }
    return data_[i];
}
}
