//
// @author: jiahong
// @date  : 24/02/18 5:22 PM
//

#include <iostream>
#include "germline.h"


namespace immulator {
std::ostream
&operator<<(std::ostream &os, const immulator::Germline &germ) {
    return os << germ.ascnum_  << '\t' << germ.name_ << '\n' << germ.seq_;
}

immulator::Germline
operator+(immulator::Germline lhs, const immulator::Germline &rhs) {
    lhs.name_ += "," + rhs.name_;
    lhs.ascnum_ += "," + rhs.ascnum_;
    recombine(lhs, rhs);
    return lhs;
}

void
recombine(immulator::Germline &lhs, const immulator::Germline &rhs) {
    lhs.seq_ += rhs.seq_;
}

}

