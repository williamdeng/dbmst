/*
 * A map that supports the nth key operation. The key ranges from
 * [0, capacity) where capacity is fixed.
 *
 *  Created on: Aug 7, 2009
 *  Author: William Deng
 */

#ifndef INDEXMAP_H_
#define INDEXMAP_H_
#include <cassert>
class IndexMap {
public:
  IndexMap(int capacity);
  virtual ~IndexMap();
  IndexMap(const IndexMap& source);
  void operator=(const IndexMap& source);
  bool containsKey(int key) const {
    assert(inRangeKey(key));
    return (keys[key]!= -1);
  }
  void put(int key, int value);
  int remove(int key);
  int getIthKey(int index) const {
    assert(index <numEntries);
    return keyList[index];
  }
  //get the value which has is mapped by key
  int operator[](int key) const {
    assert (containsKey(key));
    return values[key];
  }
  int size() const {
    return numEntries;
  }
  bool isEmpty() const {
    return numEntries == 0;
  }
  void clear() ;
  //invariant of the class
  bool repOK() const;
  int getCapacity() const {
    return capacity;
  }

private:
  bool inRangeKey(int key) const {
    return 0 <= key && key < capacity;
  }
  //invariant: the universe of the key set is [0..capacity) and current number of entries of the map is numEntries
  // and keys[0..size-1] are the indicators of whether a key is presented and if key is a key of the map then
  //keyList[keys[key]] == key otherwise keys[key] == -1.
  int capacity;
  int numEntries;
  int* values; //key -> value
  int* keyList; //key list index -> key
  int* keys;//key->key list index
};
bool operator==(const IndexMap& set1, const IndexMap& set2);

#endif /* INDEXMAP_H_ */
