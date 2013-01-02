#include <limits.h>
#include "main/common.hpp"
#include "gtest/gtest.h"


TEST(CommonTest, Trim) {
	char a[]="  test   ";
	ASSERT_STREQ("test", trim(a));
}
