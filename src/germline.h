//
// @author: jiahong
// @date  : 24/02/18 4:29 PM
//

#ifndef IMMULATOR_GERMLINE_H
#define IMMULATOR_GERMLINE_H

#include <string>
#include <iostream>


namespace immulator {

class Germline {
friend std::ostream &operator<<(std::ostream &os, const immulator::Germline &germ);
friend Germline operator+(Germline lhs, const Germline& rhs);
friend void recombine(immulator::Germline &lhs, const immulator::Germline &rhs);
public:
    Germline(const std::string& name, const std::string &ascnum, const std::string &seq) :
            name_(name), ascnum_(ascnum), seq_(seq) {}
    const std::string family_name() const {
        return name_.substr(0, name_.find_first_of('-'));
    }
    const std::string gene_name() const {
        return name_.substr(0, name_.find_first_of('*'));
    }


private:
    std::string seq_;
    std::string name_;
    std::string ascnum_;
};

std::ostream
&operator<<(std::ostream &os, const immulator::Germline &germ);

immulator::Germline operator+(Germline lhs, const Germline& rhs);


}



#endif //IMMULATOR_GERMLINE_H
