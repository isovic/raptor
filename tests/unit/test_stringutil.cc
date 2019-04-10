#include <gtest/gtest.h>

#include <utility/stringutil.h>

TEST(StringUtil, TokenizeToWhitespaces1) {
    std::string input("This is  just     a test!");
    std::vector<std::string> expected = {"This", "is", "just", "a", "test!"};
    std::vector<std::string> result = raptor::TokenizeToWhitespaces(input);
    ASSERT_EQ(result, expected);
}

TEST(StringUtil, TokenizeToWhitespaces2) {
    std::string input("This is  just\na test!");
    std::vector<std::string> expected = {"This", "is", "just", "a", "test!"};
    std::vector<std::string> result = raptor::TokenizeToWhitespaces(input);
    ASSERT_EQ(result, expected);
}

TEST(StringUtil, TokenizeToWhitespaces3) {
    std::string input("This is  just\na test!");
    std::vector<std::string> expected = {"This", "is", "just", "a", "test!"};
    std::vector<std::string> result = raptor::TokenizeToWhitespaces(input);
    ASSERT_EQ(result, expected);
}

TEST(StringUtil, TokenizeToWhitespaces4) {
    std::string input("");
    std::vector<std::string> expected = {};
    std::vector<std::string> result = raptor::TokenizeToWhitespaces(input);
    ASSERT_EQ(result, expected);
}

TEST(StringUtil, TrimToFirstWhiteSpace1) {
    std::string input("This is  just     a test!");
    std::string expected("This");
    std::string result = raptor::TrimToFirstWhiteSpace(input);
    ASSERT_EQ(result, expected);
}

TEST(StringUtil, TrimToFirstWhiteSpace2) {
    std::string input("");
    std::string expected("");
    std::string result = raptor::TrimToFirstWhiteSpace(input);
    ASSERT_EQ(result, expected);
}

TEST(StringUtil, TrimToFirstWhiteSpace3) {
    std::string input("Test");
    std::string expected("Test");
    std::string result = raptor::TrimToFirstWhiteSpace(input);
    ASSERT_EQ(result, expected);
}
