//
// @author: jiahong
// @date  : 24/02/18 4:29 PM
//

#ifndef IMMULATOR_GERMLINE_H
#define IMMULATOR_GERMLINE_H

#include <string>
#include <iostream>

namespace immulator {
struct Germline {
    Germline(const std::string& name, const std::string &ascnum, const std::string &seq) :
            name_(name), ascnum_(ascnum), seq_(seq) {}
    const std::string family_name() const {
        return name_.substr(0, name_.find_first_of('-'));
    }
    const std::string gene_name() const {
        return name_.substr(0, name_.find_first_of('*'));
    }
    const std::string seq_;
    const std::string name_;
    const std::string ascnum_;
};

}

std::ostream
&operator<<(std::ostream &os, const immulator::Germline &germ);


#endif //IMMULATOR_GERMLINE_H
