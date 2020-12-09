#pragma once

#include <cstdlib>
#include <new>
#include <utility>

namespace mem_storage {

template<typename T> class Matrix_Storage;

template <typename T>
class Vector_Storage {
protected:
    T *data_;
    size_t size_;

    Vector_Storage(size_t size): size_(size) {
        if(size > 0) {
            data_ = static_cast<T*>(::operator new(sizeof(T) * size));
        }
    }

    virtual ~Vector_Storage() {
        if(size_ == 0) return;
        for(size_t i = 0; i < size_; ++i) {
            data_[i].~T();
        }
        delete data_;
    }

    Vector_Storage(const Vector_Storage&) = delete;
    Vector_Storage& operator=(const Vector_Storage&) = delete;

    virtual void swap(Vector_Storage& rhs) noexcept {
        std::swap(size_, rhs.size_);
        std::swap(data_, rhs.data_);
    }
};


template <typename T>
class Matrix_Storage: protected Vector_Storage<T*>, protected Vector_Storage<T> {
protected:
    const size_t& mem_size_ = Vector_Storage<T>::size_;
    const size_t& column_size_ = Vector_Storage<T*>::size_;
    size_t row_size_;
    using Vector_Storage<T*>::data_;

    Matrix_Storage(size_t column_size, size_t row_size):
        Vector_Storage<T*>(column_size),
        Vector_Storage<T>(column_size * row_size),
        row_size_(row_size)
    {
        for(size_t i = 0; i < column_size; ++i) {
            data_[i] = &(Vector_Storage<T>::data_[i * row_size]);
        }
    }

    virtual ~Matrix_Storage() {}

    Matrix_Storage(const Matrix_Storage<T>&) = delete;
    Matrix_Storage& operator=(const Matrix_Storage<T>&) = delete;

    virtual void swap(Matrix_Storage<T>& rhs) noexcept {
        Vector_Storage<T>::swap(rhs);
        Vector_Storage<T*>::swap(rhs);
        std::swap(row_size_, rhs.row_size_);
    }
};
}
