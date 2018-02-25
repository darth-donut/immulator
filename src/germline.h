//
// @author: jiahong
// @date  : 24/02/18 4:29 PM
//

#ifndef IMMULATOR_GERMLINE_H
#define IMMULATOR_GERMLINE_H

#include <string>
#include <iostream>
#include <numeric>
#include <algorithm>
#include "immutils.h"


namespace immulator {

class Germline {
friend std::ostream &operator<<(std::ostream &os, const immulator::Germline &germ);
friend Germline operator+(Germline lhs, const Germline& rhs);
public:
    Germline(const std::string& name, const std::string &ascnum, const std::string &seq) :
            name_(name), ascnum_(ascnum), seq_(seq) {}

    Germline &operator+=(const Germline &other) {
        name_ += immulator::Germline::recombination_delim + other.name_;
        ascnum_ += immulator::Germline::recombination_delim + other.ascnum_;
        seq_ += other.seq_;
        return *this;
    }
public:

    const std::string family_name() const {
        if (name_.find(recombination_delim) == std::string::npos) {
            return name_.substr(0, name_.find_first_of('-'));
        } else {
            auto tokens = immulator::split_string(name_, std::string(1, recombination_delim));
            std::transform(tokens.begin(), tokens.end(), tokens.begin(), [] (auto &tok) {
                return tok.substr(0, tok.find_first_of('-'));
            });
            return immulator::join_string(tokens.cbegin(), tokens.cend(), std::string(1, recombination_delim));
        }
    }
    const std::string gene_name() const {
        if (name_.find(recombination_delim) == std::string::npos) {
            return name_.substr(0, name_.find_first_of('*'));
        } else {
            auto tokens = immulator::split_string(name_, std::string(1, recombination_delim));
            std::transform(tokens.begin(), tokens.end(), tokens.begin(), [] (auto &tok) {
                return tok.substr(0, tok.find_first_of('*'));
            });
            return immulator::join_string(tokens.cbegin(), tokens.cend(), std::string(1, recombination_delim));
        }
    }

    const std::string name() const {
        return name_;
    }

private:
    std::string seq_;
    std::string name_;
    std::string ascnum_;
    static constexpr char recombination_delim = ',';
};

std::ostream
&operator<<(std::ostream &os, const immulator::Germline &germ);

immulator::Germline operator+(Germline lhs, const Germline& rhs);


}



#endif //IMMULATOR_GERMLINE_H
