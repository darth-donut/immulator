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
#include <cassert>

namespace immulator {
class GermlineConfiguration {
friend std::ostream &operator<<(std::ostream &os, const immulator::GermlineConfiguration &gcfg);
public:
    GermlineConfiguration() {}
    GermlineConfiguration(const std::string &filename, bool percentage = false) :
            filename_(filename) {
        parse_file(percentage);
    }

    template<typename T>
    std::string next_roll(T &generator) const {
        if (filename_.empty()) {
            return "";
        } else {
            auto rand_d = dist_(generator);
            return nearest_key(rand_d, generator);
        }
    }


private:
    void parse_file(bool percentage);
    template<typename T>
    std::string nearest_key(double roll, T &generator) const;
private:
    const std::string filename_;
    std::multimap<double, std::string> germline_distribution_;
    mutable std::uniform_real_distribution<double> dist_{0, 1};
};

std::ostream &operator<<(std::ostream &os, const immulator::GermlineConfiguration &gcfg);

template<typename T>
std::string
GermlineConfiguration::nearest_key(double roll, T &generator) const {
    double distance = std::numeric_limits<double>::infinity();
    std::vector<std::string> best_selection;
    for (auto &keypair : germline_distribution_) {
        auto dist = std::abs(roll - keypair.first);
        if (dist <= distance) {
            // dist == distance => multimap has similar keys, need to store them and return a random one
            distance = dist;
            best_selection.push_back(keypair.second);
        } else {
            break;
        }
    }
    assert(!best_selection.empty());
    if (best_selection.size() == 1) {
        return best_selection.front();
    } else {
        std::uniform_int_distribution<std::vector<std::string>::size_type> dst(0, best_selection.size()-1);
        std::vector<std::string>::size_type index = dst(generator);
        return best_selection[index];
    }
}

}   // namespace immulator


#endif //IMMULATOR_GERMLINE_CONFIGURATION_H
