//
// @author: jiahong
// @date  : 24/02/18 4:29 PM
//

#ifndef IMMULATOR_GERMLINE_H
#define IMMULATOR_GERMLINE_H

#include <string>

namespace immulator {
struct Germline {
    Germline(const std::string& name, const std::string &seq) :
            name_(name), seq_(seq) {}
    const std::string seq_;
    const std::string name_;
};

}

#endif //IMMULATOR_GERMLINE_H
