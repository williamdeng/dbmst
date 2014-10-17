/*
 * indexset.cpp
 *
 *  Created on: Aug 6, 2009
 *      Author: deng
 */
#include "indexset.h"

IndexSet::IndexSet(int s) {
  map = new IndexMap(s);
}
IndexSet::IndexSet(const IndexSet& source) {
  map = new IndexMap(source.map->getCapacity());
  *map = *(source.map);
}
void IndexSet::operator=(const IndexSet& source) {
  if (this == &source)
    return;
  *map = *(source.map);
}

void IndexSet::insert(int elem) {
  assert(!contains(elem));
  map->put(elem, VALUE);
}
void IndexSet::remove(int elem) {
  assert(contains(elem));
  int ret = map->remove(elem);
  assert (ret == VALUE);
}
int IndexSet::operator [](int position) const {
  return map->getIthKey(position);
}
bool IndexSet::hasIntersetion(const IndexSet& other) const {
  for (int i = 0; i < other.size(); ++i) {
    if (contains(other[i])) {
      return true;
    }
  }
  return false;
}
bool IndexSet::isSubsetOf(const IndexSet& superSet) const {
  for(int i=0; i<this->size(); ++i) {
    if(!superSet.contains((*this)[i]))
      return false;
  }
  return true;
}
void IndexSet::makeUniverse() {
  map->clear();
  for(int i=0; i<map->getCapacity(); ++i) {
    map->put(i, VALUE);
  }
}
bool operator==(const IndexSet& set1, const IndexSet& set2) {
  if (&set1 == &set2)
    return true;
  if (set1.getCapacity() != set2.getCapacity() || set1.size() != set2.size())
    return false;
  for (int i = 0; i < set1.getCapacity(); ++i) {
    if (set1.contains(i) && !set2.contains(i) || !set1.contains(i)
        && set2.contains(i))
      return false;
  }
  return true;
}
bool operator<(const IndexSet& set1, const IndexSet& set2) {
  if (&set1 == &set2)
    return false;
  assert(set1.getCapacity() == set2.getCapacity());
  for (int i = 0; i < set1.getCapacity(); ++i) {
    if (set1.contains(i) && !set2.contains(i))
      return false;
    else if(!set1.contains(i) && set2.contains(i))
      return true;
  }
  return false;
}
