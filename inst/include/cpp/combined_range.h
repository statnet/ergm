/*  File inst/include/cpp/combined_range.h in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#pragma once
#include <optional>
#include <iterator>

namespace ergm {
inline namespace v1 {

// Generic combined iterator for two ranges of the same iterator type
// Templated on Iterator type
// Assumes Iterator supports operator*, operator++, operator!=
template <typename Range>
class CombinedIterator {
public:
    using Iterator = typename Range::iterator;
    using iterator_category = std::forward_iterator_tag;
    using value_type = typename Iterator::value_type;
    using difference_type = typename Iterator::difference_type;
    using pointer = typename Iterator::pointer;
    using reference = typename Iterator::reference;

    CombinedIterator(const Range& range1,
                    std::optional<Range> range2 = std::nullopt,
                    bool end = false)
        : range2_(range2), in1_(!end || !range2_.has_value()) {
        if (in1_ && !end) {
            it_ = range1.begin();
            end_ = range1.end();
            if (!(it_ != end_)) {
                if (range2_.has_value()) {
                    in1_ = false;
                    it_ = range2_->begin();
                    end_ = range2_->end();
                }
            }
        } else {
            if (!range2_.has_value()) {
                it_ = range1.end();
                end_ = range1.end();
            } else {
                it_ = range2_->end();
                end_ = range2_->end();
            }
        }
    }

    value_type operator*() const { return *it_; }

    CombinedIterator& operator++() {
        ++it_;
        if (range2_.has_value() && in1_ && !(it_ != end_)) {
            in1_ = false;
            it_ = range2_->begin();
            end_ = range2_->end();
        }
        return *this;
    }

    bool operator!=(const CombinedIterator& other) const {
        return in1_ != other.in1_ || it_ != other.it_;
    }

private:
    std::optional<Range> range2_;
    Iterator it_, end_;
    bool in1_;
};

// Generic combined range for two ranges of the same iterator type
// Templated on Range type (must have begin() and end() methods)
template <typename Range>
class CombinedRange {
public:
    using Iterator = typename Range::iterator;
    using CombinedIteratorType = CombinedIterator<Range>;

    CombinedRange(Range range1, std::optional<Range> range2 = std::nullopt)
        : range1_(range1), range2_(range2) {}

    CombinedIteratorType begin() const {
        if (range2_.has_value()) {
            return CombinedIteratorType(range1_, range2_, false);
        } else {
            return CombinedIteratorType(range1_, std::nullopt, false);
        }
    }
    CombinedIteratorType end() const {
        if (range2_.has_value()) {
            return CombinedIteratorType(range1_, range2_, true);
        } else {
            return CombinedIteratorType(range1_, std::nullopt, true);
        }
    }
private:
    Range range1_;
    std::optional<Range> range2_;
};

} // namespace v1
} // namespace ergm
