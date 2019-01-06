
#pragma once

#include "base/disjoint_sets.hpp"

#include <utility>

namespace principia {
namespace base {

template<typename T>
typename Subset<T>::Properties const& Subset<T>::properties() const {
  return *node_->properties_;
}

template<typename T>
typename Subset<T>::Properties& Subset<T>::mutable_properties() {
  return *node_->properties_;
}

template<typename T>
Subset<T>::Subset(not_null<Node*> node) : node_(std::move(node)) {}

template<typename T>
template<typename... SubsetPropertiesArgs>
Subset<T> Subset<T>::MakeSingleton(
    T& element,
    SubsetPropertiesArgs... subset_properties_args) {
  not_null<Node*> node = Node::Get(element);
  node->parent_ = node;
  node->rank_ = 0;
  node->properties_ =
      Properties(std::forward<SubsetPropertiesArgs>(subset_properties_args)...);
  return Subset(node);
}

template<typename T>
Subset<T> Subset<T>::Unite(Subset left, Subset right) {
  not_null<Node*> const left_root = left.node_;
  not_null<Node*> const right_root = right.node_;
  if (left_root == right_root) {
    return left;
  } else if (left_root->rank_ < right_root->rank_) {
    left_root->parent_ = right_root;
    right_root->properties_->MergeWith(*left_root->properties_);
    return Subset(right_root);
  } else if (right_root->rank_ < left_root->rank_) {
    right_root->parent_ = left_root;
    left_root->properties_->MergeWith(*right_root->properties_);
    return Subset(left_root);
  } else {
    right_root->parent_ = left_root;
    ++left_root->rank_;
    left_root->properties_->MergeWith(*right_root->properties_);
    return Subset(left_root);
  }
}

template<typename T>
Subset<T> Subset<T>::Find(T& element) {
  return Subset(Node::Get(element)->Root());
}

template<typename T>
Subset<T>::Node::Node() : parent_(this) {}

template<typename T>
not_null<typename Subset<T>::Node*> Subset<T>::Node::Root() {
  if (parent_ != this) {
    // Path compression.
    parent_ = parent_->Root();
  }
  return parent_;
}

}  // namespace base
}  // namespace principia
