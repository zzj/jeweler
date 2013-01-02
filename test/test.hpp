// A header only file to predefine some macros for unit testing
// Copied from
// http://reference-man.blogspot.com/2009/07/checking-vectors-of-doubles-with-google.html
// EXPECT_ITERABLE_DOUBLE_EQ(std::vector, expectedTemp, radTemperature);

#define EXPECT_ITERABLE_DOUBLE_EQ(TYPE, ref, target) \
    { \
    const TYPE& _ref(ref); \
    const TYPE& _target(target); \
    TYPE::const_iterator tarIter   = _target.begin(); \
    TYPE::const_iterator refIter = _ref.begin(); \
    unsigned int i = 0; \
    while(refIter != _ref.end()) { \
    if (tarIter == _target.end()) { \
    ADD_FAILURE() << #target \
    " has a smaller length than " #ref ; \
    break; \
    } \
        EXPECT_DOUBLE_EQ(* refIter, * tarIter) \
        << "Vectors " #ref  " (refIter) " \
           "and " #target " (tarIter) " \
            "differ at index " << i; \
            ++refIter; ++tarIter; ++i; \
    } \
    EXPECT_TRUE(tarIter == _target.end()) \
        << #ref " has a smaller length than " \
        #target ; \
    }


#define EXPECT_ITERABLE_EQ(TYPE, ref, target) \
    { \
    const TYPE& _ref(ref); \
    const TYPE& _target(target); \
    TYPE::const_iterator tarIter   = _target.begin(); \
    TYPE::const_iterator refIter = _ref.begin(); \
    unsigned int i = 0; \
    while(refIter != _ref.end()) { \
    if (tarIter == _target.end()) { \
    ADD_FAILURE() << #target \
    " has a smaller length than " #ref ; \
    break; \
    } \
        EXPECT_EQ(* refIter, * tarIter) \
        << "Vectors " #ref  " (refIter) " \
           "and " #target " (tarIter) " \
            "differ at index " << i; \
            ++refIter; ++tarIter; ++i; \
    } \
    EXPECT_TRUE(tarIter == _target.end()) \
        << #ref " has a smaller length than " \
        #target ; \
    }
