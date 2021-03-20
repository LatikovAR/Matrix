#pragma once

#include <cstdlib>
#include <new>
#include <utility>

namespace mem_storage {

template<typename T> class Matrix_Storage;

template <typename T>
class Vector_Storage final {
public:
    T *data_;
    size_t size_;
    size_t used_;

    Vector_Storage(size_t size): size_(size), used_(0) {
        if(size > 0) {
            data_ = static_cast<T*>(::operator new(sizeof(T) * size));
        }
    }

    ~Vector_Storage() {
        if(size_ == 0) return;
        for(size_t i = 0; i < used_; ++i) {
            data_[i].~T();
        }
        ::operator delete(data_);
    }

    Vector_Storage(const Vector_Storage&) = delete;
    Vector_Storage& operator=(const Vector_Storage&) = delete;

    Vector_Storage(Vector_Storage&& rhs) noexcept {
        data_ = rhs.data_;
        rhs.data_ = nullptr;
        size_ = rhs.size_;
        rhs.size_ = 0;
        used_ = rhs.used_;
        rhs.used_ = 0;
    }

    T& operator()(size_t i)& {
        return data_[i];
    }

    const T& operator()(size_t i) const& {
        return data_[i];
    }

    void swap(Vector_Storage& rhs) noexcept {
        std::swap(used_, rhs.used_);
        std::swap(size_, rhs.size_);
        std::swap(data_, rhs.data_);
    }
};


template <typename T>
class Matrix_Storage final {
private:
    Vector_Storage<T*> rows_;
    Vector_Storage<T> data_;
    size_t row_size_;

public:
    size_t mem_size() const { return data_.size_; }
    size_t column_size() const { return rows_.size_; }
    size_t row_size() const { return row_size_; }
    size_t& used()& { return data_.used_; }
    T& operator()(size_t row_num, size_t column_num)& {
        return *(rows_(row_num) + column_num);
    }
    const T& operator()(size_t row_num, size_t column_num) const& {
        return *(rows_(row_num) + column_num);
    }
    T*& operator()(size_t row_num)& {
        return rows_(row_num);
    }
    const T*& operator()(size_t row_num) const& {
        return rows_(row_num);
    }

    Matrix_Storage(size_t column_size, size_t row_size):
        rows_(column_size),
        data_(column_size * row_size),
        row_size_(row_size)
    {
        for(size_t i = 0; i < column_size; ++i) {
            rows_(i) = &(data_(i * row_size));
            rows_.used_ += 1;
        }
    }

    ~Matrix_Storage() = default;

    Matrix_Storage(const Matrix_Storage<T>&) = delete;
    Matrix_Storage& operator=(const Matrix_Storage<T>&) = delete;

    Matrix_Storage(Matrix_Storage&& rhs) noexcept:
        rows_{std::move(rhs.rows_)},
        data_{std::move(rhs.data_)},
        row_size_(rhs.row_size_)
    {
        rhs.row_size_ = 0;
    }

    void swap(Matrix_Storage<T>& rhs) noexcept {
        data_.swap(rhs.data_);
        rows_.swap(rhs.rows_);
        std::swap(row_size_, rhs.row_size_);
    }
};


template <typename T>
class Symmetric_Matrix_Storage final {
private:
    Vector_Storage<T*> rows_;
    Vector_Storage<T> data_;

    static size_t calculate_mem_size(size_t str_size) {
        return ((1 + str_size) * str_size) / 2;
    }

public:
    size_t mem_size() const { return data_.size_; }
    size_t size() const { return rows_.size_; }
    size_t& used()& { return data_.used_; }
    T& operator()(size_t row_num, size_t column_num)& {
        if(row_num >= column_num) {
            return *(rows_(row_num) + column_num);
        }
        else {
            return *(rows_(column_num) + row_num);
        }
    }
    const T& operator()(size_t row_num, size_t column_num) const& {
        if(row_num >= column_num) {
            return *(rows_(row_num) + column_num);
        }
        else {
            return *(rows_(column_num) + row_num);
        }
    }

    Symmetric_Matrix_Storage(size_t size):
        rows_(size),
        data_(calculate_mem_size(size))
    {
        size_t p_pos = 0;
        for(size_t i = 0; i < size; ++i) {
            rows_(i) = &(data_(p_pos));
            rows_.used_ += 1;

            p_pos += (i + 1);
        }
    }

    ~Symmetric_Matrix_Storage() = default;

    Symmetric_Matrix_Storage(const Symmetric_Matrix_Storage<T>&) = delete;
    Symmetric_Matrix_Storage& operator=(const Symmetric_Matrix_Storage<T>&) = delete;

    Symmetric_Matrix_Storage(Symmetric_Matrix_Storage&& rhs) noexcept:
        rows_{std::move(rhs.rows_)},
        data_{std::move(rhs.data_)} {}

    void swap(Symmetric_Matrix_Storage<T>& rhs) noexcept {
        data_.swap(rhs.data_);
        rows_.swap(rhs.rows_);
    }
};

} //namespace mem_storage
