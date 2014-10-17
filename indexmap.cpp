/*
 * IndexMap.cpp
 * Implementation of a map where the keys are from [0..size-1] and
 * values are integers too.
 *
 *  Created on: Aug 7, 2009
 *      Author: deng
 */

#include <cassert>
#include "indexmap.h"

IndexMap::IndexMap(int capacity) {
  assert (capacity> 0);
  this->capacity = capacity;
  numEntries = 0;
  values = new int[capacity];
  keyList = new int[capacity];
  keys = new int[capacity];
  for (int i = 0; i < capacity; i++)
    keys[i] = -1;
}

IndexMap::~IndexMap() {
  delete[] keys;
  delete[] keyList;
  delete[] values;
}

IndexMap::IndexMap(const IndexMap& source) {
  capacity = source.capacity;
  numEntries = source.numEntries;
  keyList = new int[capacity];
  keys = new int[capacity];
  values = new int[capacity];
  for (int i = 0; i < capacity; ++i) {
    keyList[i] = source.keyList[i];
    keys[i] = source.keys[i];
    values[i] = source.values[i];
  }
}
void IndexMap::operator=(const IndexMap& source) {
  assert(capacity == source.capacity);
  numEntries = source.numEntries;
  for (int i = 0; i < capacity; ++i) {
    keyList[i] = source.keyList[i];
    keys[i] = source.keys[i];
    values[i] = source.values[i];
  }
}

void IndexMap::put(int key, int value) {
  assert(inRangeKey(key));
  if (containsKey(key)) {
    values[key] = value;
  } else {
    keyList[numEntries] = key;
    keys[key] = numEntries;
    ++numEntries;
    values[key] = value;
  }
}
int IndexMap::remove(int key) {
  assert(containsKey(key));
  keyList[keys[key]] = keyList[numEntries - 1];
  keys[keyList[numEntries - 1]] = keys[key];
  --numEntries;
  keys[key] = -1;
  return values[key];
}
void IndexMap::clear() {
  numEntries = 0;
  for (int i = 0; i < capacity; ++i)
    keys[i] = -1;
}
//invariant of the class
bool IndexMap::repOK() const {
  for (int i = 0; i < numEntries; ++i) {
    if (!inRangeKey(keyList[i])) {
      return false;
    }
    if (keys[keyList[i]] != i) {
      return false;
    }
  }
  for (int i = 0; i < capacity; ++i) {
    bool found = false;
    for (int j = 0; j < numEntries; ++j) {
      if (i == keyList[j]) {
        found = true;
        break;
      }
    }
    if (found) {
      if (keys[i] < 0 || keys[i] >= numEntries) {
        return false;
      }
      if (keyList[keys[i]] != i) {
        return false;
      }
    }
    if (!found && keys[i] != -1) {
      return false;
    }
  }
  return true;

}
bool operator==(const IndexMap& map1, const IndexMap& map2) {
  if (map1.getCapacity() != map2.getCapacity())
    return false;
  if (map1.size() != map2.size())
    return false;
  for (int i = 0; i < map1.size(); i++) {
    int key = map1.getIthKey(i);
    if (!map2.containsKey(key) || map2[key] != map1[key])
      return false;
  }
  return true;
}
