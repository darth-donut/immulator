//
// @author: jiahong
// @date  : 24/02/18 5:52 PM
//

#ifndef IMMULATOR_IMMUTILS_H
#define IMMULATOR_IMMUTILS_H


#include <string>
#include <vector>
#include <unordered_map>
#include <cctype>

namespace immulator {
inline std::string strip_string(const std::string &str, const std::string &delim);
inline std::vector<std::string> split_string(const std::string &str, const std::string &delim);
inline std::string translate(const std::string &ntseq);
inline std::string &toupper(std::string &str);
inline std::string toupper(const std::string &str);
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

/// warning: when using this function, it's doing touooer INPLACE, then returning the parameter itself
/// \param str string
/// \return str, convert to uppercase inplace
std::string
&toupper(std::string &str) {
    std::transform(str.begin(), str.end(), str.begin(), [] (char c) { return std::toupper(c); });
    return str;
}

std::string
toupper(const std::string &str) {
    return std::accumulate(str.cbegin(),
                           str.cend(),
                           std::string(""),
                           [] (const std::string &acc,  char c) -> std::string {
                               return acc + std::string(1, std::toupper(c));
                           });
}

template<typename In>
std::string
join_string(const In &begin, const In &end, const std::string &delim) {
    return std::accumulate(begin, end, std::string(""), [&delim] (auto &acc, auto &tok) {
        return acc + (acc.empty() ? "" : delim) + tok;
    });
}

std::string
translate(const std::string &ntseq) {
    static constexpr char STOP_CODON = '*';
    static std::unordered_map<std::string, char> codon_table = {
            {"TTT", 'F'}, {"TCT", 'S'}, {"TAT", 'Y'}, {"TGT", 'C'},
            {"TTC", 'F'}, {"TCC", 'S'}, {"TAC", 'Y'}, {"TGC", 'C'},
            {"TTA", 'L'}, {"TCA", 'S'}, {"TAA", STOP_CODON}, {"TGA", STOP_CODON},
            {"TTG", 'L'}, {"TCG", 'S'}, {"TAG", STOP_CODON}, {"TGG", 'W'},
            {"CTT", 'L'}, {"CCT", 'P'}, {"CAT", 'H'}, {"CGT", 'R'},
            {"CTC", 'L'}, {"CCC", 'P'}, {"CAC", 'H'}, {"CGC", 'R'},
            {"CTA", 'L'}, {"CCA", 'P'}, {"CAA", 'Q'}, {"CGA", 'R'},
            {"CTG", 'L'}, {"CCG", 'P'}, {"CAG", 'Q'}, {"CGG", 'R'},
            {"ATT", 'I'}, {"ACT", 'T'}, {"AAT", 'N'}, {"AGT", 'S'},
            {"ATC", 'I'}, {"ACC", 'T'}, {"AAC", 'N'}, {"AGC", 'S'},
            {"ATA", 'I'}, {"ACA", 'T'}, {"AAA", 'K'}, {"AGA", 'R'},
            {"ATG", 'M'}, {"ACG", 'T'}, {"AAG", 'K'}, {"AGG", 'R'},
            {"GTT", 'V'}, {"GCT", 'A'}, {"GAT", 'D'}, {"GGT", 'G'},
            {"GTC" ,'V'}, {"GCC", 'A'}, {"GAC", 'D'}, {"GGC", 'G'},
            {"GTA", 'V'}, {"GCA", 'A'}, {"GAA", 'E'}, {"GGA", 'G'},
            {"GTG", 'V'}, {"GCG", 'A'}, {"GAG", 'E'}, {"GGG", 'G'}
    };
    std::string nt(ntseq);

    std::string aa;
    for (std::string::size_type i = 0; i < nt.size() / 3; ++i) {
        auto codon = immulator::toupper(nt.substr(i * 3, 3));
        if (codon_table.find(codon) == codon_table.end()) {
            std::cerr << "WARNING: Unknown codon encountered: " << codon << '\n';
        } else {
            aa += codon_table[codon];
        }
    }
    return aa;
}

template<typename T>
class optional {
public:
    optional() = default;
    optional(const T& item) : item_(item), b_(true) {}

    const T &operator*() const {
        return item_;
    }

    const T *operator->() const {
        return &item_;
    }

    T &operator*() {
        return item_;
    }

    T *operator->() {
        return &item_;
    }

    operator bool() const {
        return this->has_value();
    }

    bool has_value() const {
        return b_;
    }

    const T &value() const {
        return **this;
    }

    T &value() {
        return **this;
    }

    const T &value_or(const T& other) const {
        if (this->has_value()) {
            return **this;
        } else {
            return other;
        }
    }

    T &value_or(const T& other) {
        if (this->has_value()) {
            return **this;
        } else {
            return other;
        }
    }

private:
    T item_;
    bool b_ = false;
};

}

#endif //IMMULATOR_IMMUTILS_H


