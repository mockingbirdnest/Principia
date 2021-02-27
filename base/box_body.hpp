#pragma once

#include "glog/logging.h"

namespace principia {
namespace base {

template <typename T, typename Allocator>
Box<T, Allocator> Box<T, Allocator>::UnsafeMakeUninitialized() {
  return Box<T, Allocator>(UninitTag());
}

template <typename T, typename Allocator>
Box<T, Allocator> Box<T, Allocator>::FromNullableChecked(
    NullableBox<T, Allocator>&& nullable) {
  CHECK(nullable.get() != nullptr);
  Box<T, Allocator> box{NullTag()};
  box.box_.swap(nullable);
  return box;
}

template <typename T, typename Allocator>
template <typename... Args>
Box<T, Allocator>::Box(Args&&... args) : box_() {
  Allocator alloc;
  allocator_traits::construct(alloc, box_.get(), std::forward<Args>(args)...);
}

template <typename T, typename Allocator>
Box<T, Allocator>::Box(Box&& other) : box_(std::move(other.box_)) {}

template <typename T, typename Allocator>
template <typename U, typename E>
Box<T, Allocator>::Box(Box<U, E>&& other) : box_(std::move(other.box_)) {}

template <typename T, typename Allocator>
Box<T, Allocator>::Box(Box<T, Allocator>::UninitTag) : box_() {}

template <typename T, typename Allocator>
Box<T, Allocator>::Box(Box<T, Allocator>::NullTag) : box_(nullptr) {}

template <typename T, typename Allocator>
Box<T, Allocator>& Box<T, Allocator>::operator=(Box<T, Allocator>&& other) {
  box_.swap(other.box_);
  return *this;
}
template <typename T, typename Allocator>
NullableBox<T, Allocator> Box<T, Allocator>::IntoNullable() && {
  return std::move(box_);
}

template <typename T, typename Allocator>
void Box<T, Allocator>::swap(Box& other) {
  box_.swap(other.box_);
}

template <typename T, typename Allocator>
T* Box<T, Allocator>::get() const {
  return box_.get();
}

template <typename T, typename Allocator>
typename std::add_lvalue_reference<T>::type Box<T, Allocator>::operator*()
    const {
  return *box_;
}
template <typename T, typename Allocator>
T* Box<T, Allocator>::operator->() const {
  return box_.operator->();
}

template <typename T, typename Allocator>
bool Box<T, Allocator>::operator==(T* other) const {
  return box_ == other;
}
template <typename T, typename Allocator>
bool Box<T, Allocator>::operator==(Box const& other) const {
  return box_ == other.box_;
}
template <typename T, typename Allocator>
bool Box<T, Allocator>::operator!=(T* other) const {
  return box_ != other;
}
template <typename T, typename Allocator>
bool Box<T, Allocator>::operator!=(Box const& other) const {
  return box_ != other.box_;
}

template <typename T, typename Allocator>
bool Box<T, Allocator>::operator<(Box const& other) const {
  return box_ < other.box_;
}
template <typename T, typename Allocator>
bool Box<T, Allocator>::operator<=(Box const& other) const {
  return box_ <= other.box_;
}
template <typename T, typename Allocator>
bool Box<T, Allocator>::operator>=(Box const& other) const {
  return box_ >= other.box_;
}
template <typename T, typename Allocator>
bool Box<T, Allocator>::operator>(Box const& other) const {
  return box_ > other.box_;
}

template <typename T, typename Allocator>
NullableBox<T, Allocator>::NullableBox(Box<T>&& box) : ptr_(nullptr) {
  *this = std::move(box).IntoNullable();
}

template <typename T, typename Allocator>
NullableBox<T, Allocator>::NullableBox(std::nullptr_t) : ptr_(nullptr) {}

template <typename T, typename Allocator>
NullableBox<T, Allocator>::NullableBox(NullableBox&& other) : ptr_(nullptr) {
  swap(other);
}

template <typename T, typename Allocator>
template <typename U, typename E, typename>
NullableBox<T, Allocator>::NullableBox(NullableBox<U, E>&& other) {
  ptr_ = other.ptr_;
  other.ptr_ = nullptr;
}

template <typename T, typename Allocator>
NullableBox<T, Allocator>::NullableBox() {
  Allocator alloc;
  ptr_ = allocator_traits::allocate(alloc, 1);
}

template <typename T, typename Allocator>
NullableBox<T, Allocator>& NullableBox<T, Allocator>::operator=(
    NullableBox<T, Allocator>&& other) {
  this->swap(other);
  return *this;
}

template <typename T, typename Allocator>
NullableBox<T, Allocator>::~NullableBox() {
  if (ptr_ != nullptr) {
    Allocator alloc;
    allocator_traits::destroy(alloc, ptr_);
    allocator_traits::deallocate(alloc, ptr_, 1);
  }
}

template <typename T, typename Allocator>

Box<T, Allocator> NullableBox<T, Allocator>::AssertNotNull() && {
  return Box<T, Allocator>::FromNullableChecked(std::move(*this));
}

template <typename T, typename Allocator>
void NullableBox<T, Allocator>::swap(NullableBox& other) {
  std::swap(ptr_, other.ptr_);
}

template <typename T, typename Allocator>
T* NullableBox<T, Allocator>::get() const {
  return ptr_;
}

template <typename T, typename Allocator>
typename std::add_lvalue_reference<T>::type
NullableBox<T, Allocator>::operator*() const {
  return *ptr_;
}
template <typename T, typename Allocator>
T* NullableBox<T, Allocator>::operator->() const {
  return ptr_;
}

template <typename T, typename Allocator>
bool NullableBox<T, Allocator>::operator==(T* other) const {
  return ptr_ == other;
}
template <typename T, typename Allocator>
bool NullableBox<T, Allocator>::operator==(NullableBox const& other) const {
  return ptr_ == other.ptr_;
}
template <typename T, typename Allocator>
bool NullableBox<T, Allocator>::operator!=(T* other) const {
  return ptr_ != other;
}
template <typename T, typename Allocator>
bool NullableBox<T, Allocator>::operator!=(NullableBox const& other) const {
  return ptr_ != other.ptr_;
}

template <typename T, typename Allocator>
bool NullableBox<T, Allocator>::operator<(NullableBox const& other) const {
  return ptr_ < other.ptr_;
}
template <typename T, typename Allocator>
bool NullableBox<T, Allocator>::operator<=(NullableBox const& other) const {
  return ptr_ <= other.ptr_;
}
template <typename T, typename Allocator>
bool NullableBox<T, Allocator>::operator>=(NullableBox const& other) const {
  return ptr_ >= other.ptr_;
}
template <typename T, typename Allocator>
bool NullableBox<T, Allocator>::operator>(NullableBox const& other) const {
  return ptr_ > other.ptr_;
}

template <typename T, typename Allocator>
void swap(Box<T, Allocator>& lhs, Box<T, Allocator>& rhs) {
  lhs.swap(rhs);
}

template <typename T, typename Allocator>
void swap(NullableBox<T, Allocator>& lhs, NullableBox<T, Allocator>& rhs) {
  lhs.swap(rhs);
}

}  // namespace base
}  // namespace principia
