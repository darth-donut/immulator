//
// @author: jiahong
// @date  : 24/02/18 5:52 PM
//

#ifndef IMMULATOR_IMMUTILS_H
#define IMMULATOR_IMMUTILS_H


#include <string>
#include <vector>

namespace immulator {
inline std::string strip_string(const std::string &str, const std::string &delim);
inline std::vector<std::string> split_string(const std::string &str, const std::string &delim);
template<typename In> inline std::string join_string(const In &begin, const In &end, const std::string &delim);

std::string
strip_string(const std::string &str, const std::string &delim) {
    std::string left_strip = str.substr(str.find_first_not_of(delim));
    std::string rstring(left_strip.rbegin(), left_strip.rend());
    if (rstring.find_first_of(' ') != std::string::npos) {
        return left_strip.substr(left_strip.size() - rstring.find_first_of(delim));
    } else {
        return left_strip;
    }
}

std::vector<std::string>
split_string(const std::string &str, const std::string &delim) {
    std::vector<std::string> collector;
    std::string::size_type pos = 0;
    std::string::size_type old_pos = pos;

    while ((pos = str.find_first_of(delim, old_pos)) != std::string::npos) {
        collector.push_back(str.substr(old_pos, pos - old_pos));
        old_pos = pos + 1;
    }
    if (old_pos < str.size()) {
        collector.push_back(str.substr(old_pos));
    }
    return collector;
}

template<typename In>
std::string
join_string(const In &begin, const In &end, const std::string &delim) {
    return std::accumulate(begin, end, std::string(""), [&delim] (auto &acc, auto &tok) {
        return acc + (acc.empty() ? "" : delim) + tok;
    });
}

}

#endif //IMMULATOR_IMMUTILS_H
