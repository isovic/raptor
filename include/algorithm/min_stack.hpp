/*
 * min_stack.hpp
 *
 *  Created on: Dec 08, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_MIN_STACK_H_
#define SRC_MIN_STACK_H_

#include <cstdint>
#include <stack>
#include <functional>
#include <array>

namespace raptor {

template<class T>
class StackElement {
public:
    StackElement(const T& _val, const T&_min_val) : val{_val}, min_val{_min_val} { }

    T val;
    T min_val;
};

/*
 * @brief MinStack is a stack implementation with minimum lookup in O(1) time.
*/
template<class T>
class MinStack {
public:
    MinStack(std::function<bool(const T& a, const T& b)> comp =
                                      [](const T& a, const T& b)
                                      { return a < b; } )
                                      : comp_{comp} {
    }

    ~MinStack() = default;

    void Push(const T& val) {
        if (stack_.size() == 0) {
            stack_.emplace(StackElement<T>(val, val));
            return;
        }

        const auto& prev_min = stack_.top().min_val;
        auto min_val = (comp_(val, prev_min)) ? val : prev_min;
        stack_.emplace(StackElement<T>(val, min_val));
    }

    void Pop() {
        if (stack_.size() == 0) {
            throw std::logic_error("The MinStack is empty. Cannot Pop().");
        }
        stack_.pop();
    }

    T Min() const {
        if (stack_.size() == 0) {
            throw std::logic_error("The MinStack is empty. Cannot Pop().");
        }
        return stack_.top().min_val;
    }

    StackElement<T> Top() {
        if (stack_.size() == 0) {
            throw std::logic_error("The MinStack is empty. Cannot Pop().");
        }
        return stack_.top();
    }

    size_t Size() const {
        return stack_.size();
    }

private:
    std::stack<StackElement<T>> stack_;
    std::function<bool(const T& a, const T& b)> comp_;
};

}

#endif
