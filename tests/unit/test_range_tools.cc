#include <gtest/gtest.h>

#include <utility/range_tools.hpp>

#include <cstdint>
#include <cstring>
#include <string>
#include <tuple>
#include <vector>

#define TEST_DEBUG_VERBOSE_

TEST(RangeTools, FindRange1) {
    std::vector<int32_t> data = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0};

    auto result = istl::FindRanges<int32_t>(data);

    std::vector<std::pair<size_t, size_t>> expected = {
                                            {0, 1},
                                            {1, 2},
                                            {2, 3},
                                            {3, 4},
                                            {4, 5},
                                            {5, 6},
                                            {6, 7},
                                            {7, 8},
                                            {8, 9},
                                            {9, 10},
                                                        };

    ASSERT_EQ(result, expected);
}

TEST(RangeTools, FindRange2) {
    std::vector<int32_t> data = {1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 0};

    auto result= istl::FindRanges<int32_t>(data);

    std::vector<std::pair<size_t, size_t>> expected = {
                                            {0, 1},
                                            {1, 3},
                                            {3, 4},
                                            {4, 5},
                                            {5, 6},
                                            {6, 7},
                                            {7, 8},
                                            {8, 9},
                                            {9, 10},
                                            {10, 11},
                                                        };

    ASSERT_EQ(result, expected);
}

TEST(RangeTools, FindRange3) {
    std::vector<int32_t> data = {1, 2, 2, 3, 3, 3, 3, 3, 4};

    auto result= istl::FindRanges<int32_t>(data);

    std::vector<std::pair<size_t, size_t>> expected = {
                                            {0, 1},
                                            {1, 3},
                                            {3, 8},
                                            {8, 9},
                                                        };

    ASSERT_EQ(result, expected);
}

TEST(RangeTools, FindRange4) {
    std::vector<int32_t> data = {1, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4};

    auto result= istl::FindRanges<int32_t>(data);

    std::vector<std::pair<size_t, size_t>> expected = {
                                            {0, 1},
                                            {1, 3},
                                            {3, 8},
                                            {8, 13},
                                                        };

    ASSERT_EQ(result, expected);
}

TEST(RangeTools, FindRange5) {
    std::vector<int32_t> data = {};

    auto result= istl::FindRanges<int32_t>(data);

    std::vector<std::pair<size_t, size_t>> expected = { };

    ASSERT_EQ(result, expected);
}

TEST(RangeTools, FindRange6) {
    std::vector<int32_t> data = {1};

    auto result= istl::FindRanges<int32_t>(data);

    std::vector<std::pair<size_t, size_t>> expected = {
                                            {0, 1},
                                                        };

    ASSERT_EQ(result, expected);
}

class TestDataClass {
    public:
    TestDataClass(int64_t _target_id, bool _target_rev)
                : target_id(_target_id), target_rev(_target_rev) { }
    int64_t target_id;
    bool target_rev;
};

TEST(RangeTools, FindRangeInCustomClass1) {
    std::vector<TestDataClass> data = {
                                    TestDataClass(0, false),
                                    TestDataClass(0, false),
                                    TestDataClass(0, false),
                                    TestDataClass(0, true),
                                    TestDataClass(0, true),
                                    TestDataClass(1, true),
                                    TestDataClass(1, true),
                                    TestDataClass(1, true),
                                    TestDataClass(1, true),
                                    TestDataClass(1, true),
                                    TestDataClass(17, false),
                                };

    auto result = istl::FindRanges<TestDataClass>(data,
                                    [](const TestDataClass& a, const TestDataClass& b) {
                                        return (a.target_id == b.target_id && a.target_rev == b.target_rev);
                                    }
                                );

    std::vector<std::pair<size_t, size_t>> expected = {
                                            {0, 3},
                                            {3, 5},
                                            {5, 10},
                                            {10, 11},
                                                        };

    ASSERT_EQ(result, expected);
}

TEST(RangeTools, FindRangeInCustomClass2) {
		std::string data = "122233455555";

		auto result = istl::FindRanges<char, std::string>(data, [](char a, char b) { return a == b; });

    std::vector<std::pair<size_t, size_t>> expected = {
                                            {0, 1},
                                            {1, 4},
                                            {4, 6},
                                            {6, 7},
                                            {7, 12},
                                                        };

    ASSERT_EQ(result, expected);
}

