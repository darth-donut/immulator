//
// @author: jiahong
// @date  : 24/02/18 5:41 PM
//

#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include "germline_configuration.h"
#include "immutils.h"

using std::string;
using std::getline;
using std::fstream;

void
immulator::GermlineConfiguration::parse_file(bool percentage) {
    fstream ifs(filename_);
    string buffer;
    if (ifs) {
        while (getline(ifs, buffer)) {
            string germ_name = strip_string(buffer.substr(0, buffer.find_first_of(',')), " ");
            double prob = 1 - (stod(buffer.substr(buffer.find_first_of(',') + 1)) / (percentage ? 100 : 1));
            germline_distribution_.insert(std::make_pair(prob, germ_name));
        }
    }
}

namespace immulator {

std::ostream
&operator<<(std::ostream &os, const immulator::GermlineConfiguration &gcfg) {
    std::string buffer;
    std::for_each(gcfg.germline_distribution_.cbegin(), gcfg.germline_distribution_.cend(), [&buffer] (auto &keypair) {
        buffer += std::to_string(keypair.first) + "\t" + keypair.second + "\n";
    });
    // remove last \n
    buffer.pop_back();
    return os << buffer;
}

std::string
GermlineConfiguration::nearest_key(double roll) const {
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
        std::vector<std::string>::size_type index = dst(mersenne_);
        return best_selection[index];
    }
}

}

