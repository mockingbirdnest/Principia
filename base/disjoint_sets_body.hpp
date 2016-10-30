#pragma once

#include "base/disjoint_sets.hpp"

#include <utility>

namespace principia {
namespace base {

template<typename T>
Subset<T> Unite(Subset<T> left, Subset<T> right) {
  not_null<SubsetNode<T>*> const left_root = left.node_;
  not_null<SubsetNode<T>*> const right_root = right.node_;
  if (left_root == right_root) {
    return left;
  } else if (left_root->rank_ < right_root->rank_) {
    left_root->parent_ = right_root;
    right_root->properties_.value.MergeWith(left_root->properties_.value);
    return Subset<T>(right_root);
  } else if (right_root->rank_ < left_root->rank_) {
    right_root->parent_ = left_root;
    left_root->properties_.value.MergeWith(right_root->properties_.value);
    return Subset<T>(left_root);
  } else {
    right_root->parent_ = left_root;
    ++left_root->rank_;
    left_root->properties_.value.MergeWith(right_root->properties_.value);
    return Subset<T>(left_root);
  }
}

template<typename T>
SubsetProperties<T> const& Subset<T>::properties() {
  return node_->properties_.value;
}

template<typename T>
Subset<T>::Subset(not_null<SubsetNode<T>*> node) : node_(node) {}

template<typename T>
template<typename... SubsetPropertiesArgs>
Subset<T> Subset<T>::MakeSingleton(
    T& element,
    SubsetPropertiesArgs... subset_properties_args) {
  not_null<SubsetNode<T>*> node = GetSubsetNode(element);
  node->parent_ = node;
  node->rank_ = 0;
  node->properties_.value = SubsetProperties<T>(
      std::forward<SubsetPropertiesArgs>(subset_properties_args)...);
  return Subset<T>(node);
}

template<typename T>
Subset<T> Subset<T>::Find(T& element) {
  return Subset<T>(GetSubsetNode(element)->Root());
}

template<typename T>
SubsetNode<T>::SubsetNode() : parent_(this), properties_({0xDB}) {}

template<typename T>
not_null<SubsetNode<T>*> SubsetNode<T>::Root() {
  if (parent_ != this) {
    // Path compression.
    parent_ = parent_->Root();
  }
  return parent_;
}

}  // namespace base
}  // namespace principia
