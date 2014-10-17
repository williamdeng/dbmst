#include "indexmap.h"
#include "gtest/gtest.h"

TEST(isEmpty, Empty) {
  IndexMap map(5);
  EXPECT_TRUE(map.isEmpty());
}

TEST(isEmpty, NotEmpty) {
  IndexMap map(5);
  map.put(1, 5);
  EXPECT_FALSE(map.isEmpty());
}
