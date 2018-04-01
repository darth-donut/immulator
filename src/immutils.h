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

template<typename In>
inline std::string join_string(const In &begin, const In &end, const std::string &delim);

inline std::tuple<double, std::string::size_type, std::string::size_type/*, std::string, std::string*/>
local_align(const std::string &string1, const std::string &string2, double ins, double del,
            std::unordered_map<char, std::unordered_map<char, double>> scoring_matrix);

template<typename T>
const T &max(const T &x, const T &x1);

template<typename T, typename... Args>
const T &max(const T &x, const T &x1, Args... xs);


inline bool double_eq(double d1, double d2, double epsilon = 1e-5);

inline void print_alignment(const std::string &string1, const std::string &string2, std::string::size_type start,
                            std::string::size_type end, double score, std::ostream &os = std::cout);

inline std::string allowed_nts(const std::string &rem);


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

/// warning: when using this function, it's doing toupper INPLACE, then returning the parameter itself
/// \param str string
/// \return str, convert to uppercase inplace
std::string
&toupper(std::string &str) {
    std::transform(str.begin(), str.end(), str.begin(), [](char c) { return std::toupper(c); });
    return str;
}

std::string
toupper(const std::string &str) {
    return std::accumulate(str.cbegin(),
                           str.cend(),
                           std::string(""),
                           [](const std::string &acc, char c) -> std::string {
                               return acc + static_cast<char>(std::toupper(c));
                           });
}

template<typename In>
std::string
join_string(const In &begin, const In &end, const std::string &delim) {
    return std::accumulate(begin, end, std::string(""), [&delim](auto &acc, auto &tok) {
        return acc + (acc.empty() ? "" : delim) + tok;
    });
}

std::string
allowed_nts(const std::string &rem) {
    // TAA TGA TAG -> stop codons
    // the only nt we need to becareful of are TA and TG
    // banned is a map of cautionary nucleotide sequences to "permitted" nucleotide suffix
    // where permitted will not cause the translated sequence to contain a stop codon
    static std::unordered_map<std::string, std::string> banned = {
            {"TA", "CT"}, {"TG", "CGT"}
    };
    return banned.find(rem) == banned.end() ? "ACGT" : banned[rem];
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

    std::string aa;
    for (std::string::size_type i = 0; i < ntseq.size() / 3; ++i) {
        auto codon = immulator::toupper(ntseq.substr(i * 3, 3));
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

    optional(const T &item) : item_(item), b_(true) {}

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

    const T &value_or(const T &other) const {
        if (this->has_value()) {
            return **this;
        } else {
            return other;
        }
    }

    T &value_or(const T &other) {
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


std::tuple<double, std::string::size_type, std::string::size_type/*, std::string, std::string*/>
local_align(const std::string &string1, const std::string &string2, double ins, double del,
            std::unordered_map<char, std::unordered_map<char, double>> scoring_matrix) {
    std::vector<std::vector<double>> matrix;

    /*
     *  -s t r i n g 1
     * -+---------------------
     * s|0 0 0 0 0 0 0 0 0 ...
     * t|0
     * r|0
     * i|0
     * n|0
     * g|0
     * 2|0
     *  |.
     *  |.
     *  |.
     */
    matrix.resize(string2.size() + 1);
    for (auto &row : matrix) {
        row.resize(string1.size() + 1);
        row[0] = 0;
    }
    matrix[0] = std::vector<double>(string1.size() + 1, 0);

    double max_val = -std::numeric_limits<double>::infinity();
    std::pair<std::size_t, std::size_t> best_index;

    for (std::string::size_type i = 1; i <= string2.size(); ++i) {
        for (std::string::size_type j = 1; j <= string1.size(); ++j) {
            matrix[i][j] = immulator::max(matrix[i - 1][j] + del,
                                          matrix[i][j - 1] + ins,
                                          matrix[i - 1][j - 1] + scoring_matrix[string1[j - 1]][string2[i - 1]],
                                          0.0
            );
            if (matrix[i][j] > max_val) {
                max_val = matrix[i][j];
                best_index = {i, j};
            }
        }
    }

    auto i = best_index.first;
    auto j = best_index.second;

//    std::string re_string1;
//    std::string re_string2;

    while (i > 0 && j > 0) {
        auto val = matrix[i][j];
        if (double_eq(val, matrix[i - 1][j] + del)) {
//            re_string1 += '-';
//            re_string2 += string2[i-1];
            --i;
        } else if (double_eq(val, matrix[i][j - 1] + ins)) {
//            re_string2 += '-';
//            re_string1 += string1[j-1];
            --j;
        } else if (double_eq(val, matrix[i - 1][j - 1] + scoring_matrix[string1[j - 1]][string2[i - 1]])) {
//            re_string2 += string2[i-1];
//            re_string1 += string1[j-1];
            --i, --j;
        } else {
//            if (j < i) {
//                re_string1 = std::string(i-j, '-') + string1;
//                re_string2 = string2;
//            } else {
//                re_string2 = std::string(j-i, '-') + string2;
//                re_string1 = string1;
//            }
//            if (re_string1.size() < re_string2.size()) {
//                re_string1 += std::string(re_string2.size() - re_string1.size(), '-');
//            } else {
//                re_string2 += std::string(re_string1.size() - re_string2.size(), '-');
//            }
            return std::make_tuple(max_val, std::max(i, j), std::max(best_index.first, best_index.second)
                    /*, re_string1, re_string2*/);
        }
    }
//    std::string trailing1, trailing2;
//    if (i == 0) {
//        re_string2 += std::string(j, '-');
//        re_string1 = std::string(string1.crbegin(), string1.crend());
//    } else {
//        re_string1 += std::string(i, '-');
//        re_string2 = std::string(string2.crbegin(), string2.crend());
//    }
//
//    std::string::size_type diff;
//    if (re_string1.size() > re_string2.size()) {
//        diff = re_string1.size() - re_string2.size();
//    } else {
//        diff = re_string2.size() - re_string1.size();
//    }
//    if (diff > 0) {
//        if (re_string1.size() > re_string2.size()) {
//            trailing2 += std::string(diff, '-');
//        } else {
//            trailing1 += std::string(diff, '-');
//        }
//    }
    return std::make_tuple(max_val, i == 0 ? j : i, best_index.second + i/*,
                           std::string(re_string1.crbegin(), re_string1.crend()) + trailing1,
                           std::string(re_string2.crbegin(), re_string2.crend()) + trailing2*/);


};


// yes.. we can use std::max({1,2,3,... }). whatever.
template<typename T>
const T &max(const T &x, const T &x1) {
    return std::max(x, x1);
}

template<typename T, typename... Args>
const T &max(const T &x, const T &x1, Args... xs) {
    return immulator::max(std::max(x, x1), xs...);
};

bool
double_eq(double d1, double d2, double epsilon) {
    return std::abs(d1 - d2) <= epsilon;
}

void
print_alignment(const std::string &string1, const std::string &string2, std::string::size_type start,
                std::string::size_type end, double score, std::ostream &os) {
    os << "Alignment score: " << score << '\n'
       << string1 << '\n'
       << std::string(start, ' ') << std::string(end - start, '|') << '\n'
       << string2 << '\n';
}

} // namespace immulator


#endif //IMMULATOR_IMMUTILS_H


