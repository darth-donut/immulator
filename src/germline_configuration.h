//
// @author: jiahong
// @date  : 24/02/18 5:41 PM
//

#ifndef IMMULATOR_GERMLINE_CONFIGURATION_H
#define IMMULATOR_GERMLINE_CONFIGURATION_H

#include <string>
#include <random>
#include <map>
#include <algorithm>
#include <iostream>

namespace immulator {
class GermlineConfiguration {
friend std::ostream &operator<<(std::ostream &os, const immulator::GermlineConfiguration &gcfg);
public:
    GermlineConfiguration() {}
    GermlineConfiguration(const std::string &filename, bool percentage = false,
                          const std::size_t seed = std::random_device{}()) :
            filename_(filename), mersenne_(seed) {
        parse_file(percentage);
    }
    std::string next_roll() const {
        if (filename_.empty()) {
            return "";
        } else {
            auto rand_d = dist_(mersenne_);
            return nearest_key(rand_d);
        }
    }


private:
    void parse_file(bool percentage);
    std::string nearest_key(double roll) const;
private:
    const std::string filename_;
    std::multimap<double, std::string> germline_distribution_;
    mutable std::mt19937 mersenne_;
    mutable std::uniform_real_distribution<double> dist_{0, 1};
};

}

namespace immulator {
std::ostream
&operator<<(std::ostream &os, const immulator::GermlineConfiguration &gcfg);

}

#endif //IMMULATOR_GERMLINE_CONFIGURATION_H
