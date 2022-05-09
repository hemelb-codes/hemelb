// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef ENUMERATE_ITERATOR_HPP
#define ENUMERATE_ITERATOR_HPP

template <typename Container>
class EnumIter {
public:
  using iterator = decltype(std::begin(std::declval<Container>()));

  EnumIter(iterator begin) : iter_{begin}, index_{0} {
  }

  EnumIter& operator++() {
    iter_++;
    index_++;
    return *this;
  }

  bool operator!=(const EnumIter& rhs) {
    return iter_ != rhs.iter_; // or self.index_ != rhs.index_;
  }

  using iter_ref = typename std::iterator_traits<iterator>::reference;

  std::pair<int, iter_ref> operator*() const {
    return { index_, *iter_ };
  }

private:
  iterator iter_;
  int index_ = 0;
};

template <typename Container>
class EnumerationAdaptor
{
public:
  EnumerationAdaptor(Container& container) : container_(container) {}
  EnumIter<Container> begin() const { return std::begin(container_); }
  EnumIter<Container> end() const { return std::end(container_); }

private:
  Container& container_;
};

template <typename Container>
EnumerationAdaptor<Container> 
enumerate(Container& container) {
  return container;
}

template <typename Container>
EnumerationAdaptor<const Container> 
const_enumerate(const Container& container) {
  return container;
}

#endif
