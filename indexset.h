/*
 * indexset.h
 *
 * Implement an efficient set if the universe is fixed from 0..capacity-1.
 * Created on: Jul 30, 2009
 * Author: William Deng
 */

#ifndef INDEX_SET_H
#define INDEX_SET_H
#include <cassert>
#include "indexmap.h"

class IndexSet {
public:
  IndexSet(int s);
  IndexSet(const IndexSet& source);
  void operator=(const IndexSet& source);
  ~IndexSet() {
    delete map;
  }
  bool contains(int elem) const {
    return map->containsKey(elem);
  }
  void insert(int elem);
  void remove(int elem);
  int operator [](int position) const;
  int getIthElement(int position) const {
    return (*this)[position];
  }
  int size() const {
    return map->size();
  }
  bool hasIntersetion(const IndexSet& other) const;
  bool isSubsetOf(const IndexSet& superSet) const;
  bool isEmpty() const {
    return map->size() == 0;
  }
  //make the set to be a universe, that is, contain all numbers from [0..capacity-1]
  void makeUniverse();
  void clear() {
    map->clear();
  }
  int getCapacity() const {
    return map->getCapacity();
  }
  //invariant of the class
  bool repOK() const {
    return map->repOK();
  }

private:
  const static int VALUE = 1;
  IndexMap* map;
};
bool operator==(const IndexSet& set1, const IndexSet& set2);
bool operator<(const IndexSet& set1, const IndexSet& set2);

#endif /* INDEX_SET_H */
