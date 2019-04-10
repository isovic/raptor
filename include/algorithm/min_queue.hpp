/*
 * min_stack.hpp
 *
 *  Created on: Dec 08, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_MIN_QUEUE_H_
#define SRC_MIN_QUEUE_H_

#include <cstdint>
#include <functional>
#include <vector>
#include <algorithm/min_stack.hpp>

namespace raptor {

/*
 * @brief MinStack is a stack implementation with minimum lookup in O(1) time.
*/
template<class T>
class MinQueue {
public:
    MinQueue(std::function<bool(const T& a, const T& b)> comp =
                                      [](const T& a, const T& b)
                                      { return a < b; } )
                                      : stacks_{raptor::MinStack<T>(comp), raptor::MinStack<T>(comp)},
                                        comp_{comp} {
    }

    ~MinQueue() = default;

    void Enqueue(const T& val) {
        stacks_[1].Push(val);
    }

    raptor::StackElement<T> Dequeue() {
        if (stacks_[0].Size() == 0) {
            while (stacks_[1].Size()) {
                stacks_[0].Push(stacks_[1].Top().val);
                stacks_[1].Pop();
            }
        }
        auto ret = stacks_[0].Top();
        stacks_[0].Pop();
        return ret;
    }

    T Min() const {
        if (stacks_[0].Size()) {
            if (stacks_[1].Size()) {
                auto min_0 = stacks_[0].Min();
                auto min_1 = stacks_[1].Min();
                return ((comp_(min_0, min_1)) ? (min_0) : (min_1));
            }
            return stacks_[0].Min();
        }
        return stacks_[1].Min();
    }

private:
    std::vector<raptor::MinStack<T>> stacks_;
    std::function<bool(const T& a, const T& b)> comp_;
};

}

#endif
