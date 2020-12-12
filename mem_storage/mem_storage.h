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
    size_t used_;

    Vector_Storage(size_t size): size_(size), used_(0) {
        if(size > 0) {
            data_ = static_cast<T*>(::operator new(sizeof(T) * size));
        }
    }

    virtual ~Vector_Storage() {
        if(size_ == 0) return;
        for(size_t i = 0; i < used_; ++i) {
            data_[i].~T();
        }
        ::operator delete(data_);
    }

    Vector_Storage(const Vector_Storage&) = delete;
    Vector_Storage& operator=(const Vector_Storage&) = delete;

    virtual void swap(Vector_Storage& rhs) noexcept {
        std::swap(used_, rhs.used_);
        std::swap(size_, rhs.size_);
        std::swap(data_, rhs.data_);
    }
};


template <typename T>
class Matrix_Storage:
        protected Vector_Storage<T*>,
        protected Vector_Storage<T>
{
private:
    const size_t& mem_size_ = Vector_Storage<T>::size_;
protected:
    const size_t& column_size_ = Vector_Storage<T*>::size_;
    size_t row_size_;
    using Vector_Storage<T*>::data_;
    using Vector_Storage<T>::used_;

    Matrix_Storage(size_t column_size, size_t row_size):
        Vector_Storage<T*>(column_size),
        Vector_Storage<T>(column_size * row_size),
        row_size_(row_size)
    {
        for(size_t i = 0; i < column_size; ++i) {
            data_[i] = &(Vector_Storage<T>::data_[i * row_size]);
            Vector_Storage<T*>::used_ += 1;
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


template <typename T>
class Symmetric_Matrix_Storage:
        protected Vector_Storage<T*>,
        protected Vector_Storage<T>
{
private:
    static size_t calculate_mem_size(size_t str_size) {
        return ((1 + str_size) * str_size) / 2;
    }

    const size_t& mem_size_ = Vector_Storage<T>::size_;
protected:
    using Vector_Storage<T*>::size_;
    using Vector_Storage<T*>::data_;
    using Vector_Storage<T>::used_;

    Symmetric_Matrix_Storage(size_t size):
        Vector_Storage<T*>(size),
        Vector_Storage<T>(calculate_mem_size(size))
    {
        size_t p_pos = 0;
        for(size_t i = 0; i < size; ++i) {
            data_[i] = &(Vector_Storage<T>::data_[p_pos]);
            Vector_Storage<T*>::used_ += 1;

            p_pos += (i + 1);
        }
    }

    virtual ~Symmetric_Matrix_Storage() {}

    Symmetric_Matrix_Storage(const Symmetric_Matrix_Storage<T>&) = delete;
    Symmetric_Matrix_Storage& operator=(const Symmetric_Matrix_Storage<T>&) = delete;

    virtual void swap(Symmetric_Matrix_Storage<T>& rhs) noexcept {
        Vector_Storage<T>::swap(rhs);
        Vector_Storage<T*>::swap(rhs);
    }
};
}
