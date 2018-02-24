//
// @author: jiahong
// @date  : 24/02/18 5:52 PM
//

#ifndef IMMULATOR_IMMUTILS_H
#define IMMULATOR_IMMUTILS_H


#include <string>

namespace immulator {
inline std::string strip_ws(const std::string &str);

std::string
strip_ws(const std::string &str) {
    std::string left_strip = str.substr(str.find_first_not_of(' '));
    std::string rstring(left_strip.rbegin(), left_strip.rend());
    if (rstring.find_first_of(' ') != std::string::npos) {
        return left_strip.substr(left_strip.size() - rstring.find_first_of(' '));
    } else {
        return left_strip;
    }
}

}

#endif //IMMULATOR_IMMUTILS_H
