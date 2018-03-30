#include <iostream>
#include <string>
#include <random>
#include <cassert>
#include <unordered_map>
#include <cmath>
#include <tuple>
#include <regex>

#include "germline_factory.h"

using std::string;
using immulator::Germline;

// Testing
Germline vdj_recombination(const Germline &vgerm, const Germline &dgerm, const Germline &jgerm, bool prod = true,
                           bool multiple = true);

template<typename Gen>
immulator::optional<std::pair<Germline, immulator::Germline::size_type>> vcutter(Germline vgerm, Gen &generator);

template<typename Gen>
std::tuple<Germline, immulator::Germline::size_type, bool>
dcutter(Germline dgerm, Gen &generator, const std::string &rem, bool check = true);

template<typename Gen>
immulator::optional<std::tuple<Germline, immulator::Germline::size_type, bool>>
jcutter(Germline jgerm, Gen &generator, const std::string &rem, bool check = true, bool multiple = true);

int
main() {
    constexpr int seqs = 1'000;
    std::mt19937 mersenne(std::random_device{}());

    immulator::GermlineConfiguration gcfg("../germ.cfg", true);
    const string title(80, '=');
    std::cout << title << '\n'
              << "\t\t\tConfiguration file found\n" << title << '\n'
              << gcfg << std::endl;

    immulator::GermlineFactory vgermlines("../vgerm.txt"/*, gcfg*/);
    immulator::GermlineFactory dgermlines("../dgerm.txt"/*, gcfg*/);
    immulator::GermlineFactory jgermlines("../jgerm.txt"/*, gcfg*/);

    std::map<std::string, int> counter;
    for (int i = 0; i < seqs; i++) {
        counter[vgermlines(mersenne).family_name()] += 1;
    }

    std::cout << "\nTesting recombination dummy: "
              << vdj_recombination(vgermlines(mersenne), dgermlines(mersenne), jgermlines(mersenne))
              << std::endl;

    // count germline distribution used
    double total = std::accumulate(counter.begin(), counter.end(), 0, [](auto val, auto &keypair) {
        return val + keypair.second;
    });
    std::cout << title << '\n'
              << "\t\t\tCalculating germline distribution ...\n"
              << title << std::endl;

    for (auto &i : counter) {
        std::cout << i.first << "\t" << i.second / total << std::endl;
    }
    return 0;
}


// Testing
Germline
vdj_recombination(const Germline &vgerm, const Germline &dgerm, const Germline &jgerm, bool prod, bool multiple) {
    using size_type = immulator::Germline::size_type;
    static std::mt19937 mersenne(std::random_device{}());
    auto v = vcutter(vgerm, mersenne);
    if (v) {
        // todo:: V - D INSERTION, D-J insertion
        Germline d;
        size_type dsize;
        bool d_prod;
        // starts AFTER Cys (and convert to 1-index)
        size_type cdr3_start_pos = v->second + 3 + 1;
        std::tie(d, dsize, d_prod) = dcutter(dgerm,
                                             mersenne,
                // todo:: insertion here!
                                             v->first.substr(v->first.size() - (v->first.size() % 3)),
                                             prod);
        if (d_prod) {
            Germline j;
            size_type fwgxg_conserved_index;
            bool j_prod;
            // todo:: insertion between D-J instead of just D
            auto jtry = jcutter(jgerm, mersenne, d.substr(d.size() - (v->first.size() + d.size()) % 3), prod, multiple);
            if (jtry) {
                std::tie(j, fwgxg_conserved_index, j_prod) = *jtry;
                size_type cdr3_end_pos = v->first.size() + d.size() + fwgxg_conserved_index;
                std::cerr << "fwgxg index " << fwgxg_conserved_index << std::endl;
                std::cerr << "CDR3 start: " << cdr3_start_pos << " CDR3 END: " << cdr3_end_pos << std::endl;
                return (v->first + d + j);
            } else {
                std::cerr << "J gene failed to find consensus [FW]GXG region\n";
            }
        } else {
            std::cerr << "D gene failed to be found" << std::endl;
        }
    } else {
        // todo
        std::cerr << "V gene failed to find consensus Cys aa\n";
    }
    return vgerm + dgerm + jgerm;
}

///
/// \tparam Gen
/// \param vgerm
/// \param generator
/// \return  std::pair<Trimmed V Germline, NT index of last occurring Cys> if Cys can be found.
template<typename Gen>
immulator::optional<std::pair<Germline, immulator::Germline::size_type>>
vcutter(Germline vgerm, Gen &generator) {
    using size_type = immulator::Germline::size_type;
    auto aa = immulator::translate(vgerm);
    size_type cys;
    if ((cys = aa.find_last_of('C')) == std::string::npos) {
        std::cerr << "WARNING: Cys consensus failed to be located in:\n"
                  << "\t" << vgerm << '\n';
        return {};
    } else {
//        std::cerr << "DEBUG: " << "Cys found at index " << cys << std::endl;
        size_type nuc_index = cys * 3;

        // remaining nucleotides that we can cut
        size_type nt_rem = vgerm.size() - (nuc_index + 3) - 1;

        // we can cut anywhere between 0 - nt_rem nucleotides
        std::uniform_int_distribution<size_type> idist(0, nt_rem);
        size_type final_length = vgerm.size() - idist(generator);
        return std::make_pair(vgerm.trim(0, final_length), nuc_index);
    }
}

template<typename Gen>
std::tuple<Germline, immulator::Germline::size_type, bool>
dcutter(Germline dgerm, Gen &generator, const std::string &rem, bool check) {
    using size_type = immulator::Germline::size_type;
    constexpr std::size_t MAX_ATTEMPTS = 100'000;

    bool productive = true;

    /* ---------------------------------------------------------------------------- *
     *          Determine how to cut the front nt seqs from D Germline              *
     *                                                                              *
     * ---------------------------------------------------------------------------- */
    constexpr double FRONT_CUT_PERC = 30.0 / 100;

    auto max_front_cut_size = static_cast<size_type>(std::ceil(FRONT_CUT_PERC * dgerm.size()));
    std::uniform_int_distribution<size_type> front_idist(0, max_front_cut_size);

    // cut front of D gene by "front_cut" much
    auto front_cut = front_idist(generator);
    if (check) {
        auto aa = immulator::translate(rem + dgerm.substr(front_cut));
        std::size_t attempt = 0;
        for (; attempt < MAX_ATTEMPTS && aa.find('*') != std::string::npos; ++attempt) {
            front_cut = front_idist(generator);
            aa = immulator::translate(rem + dgerm.substr(front_cut));
        }
        if (attempt == MAX_ATTEMPTS) {
            std::cerr << "WARNING: Tried too hard, but in the end, nothing matters.\n";
            productive = false;
        }
    }


    /* ---------------------------------------------------------------------------- *
     *          Determine how to cut the end in nt seqs from D Germline             *
     *                                                                              *
     * ---------------------------------------------------------------------------- */

    constexpr double BACK_CUT_PERC = 30.0 / 100;
    auto max_back_cut_size = static_cast<size_type>(std::ceil(BACK_CUT_PERC * dgerm.size()));
    std::uniform_int_distribution<size_type> back_idist(0, max_back_cut_size);

    // cut back of D gene by "back_cut" much
    auto back_cut = front_idist(generator);
    return std::make_tuple(dgerm.trim(front_cut, dgerm.size() - back_cut - front_cut),
                           dgerm.size() - front_cut - back_cut,
                           productive);
}

template<typename Gen>
immulator::optional<std::tuple<Germline, immulator::Germline::size_type, bool>>
jcutter(Germline jgerm, Gen &generator, const std::string &rem, bool check, bool multiple) {
    using size_type = immulator::Germline::size_type;
    constexpr std::size_t MAX_ATTEMPTS = 100'000;

    bool productive = true;

    /* ------------------------------------------------------------------------------ *
     *                          Find consensus FR4 [FW]G.G                            *
     *                                                                                *
     * ------------------------------------------------------------------------------ */

    static std::unordered_map<std::string, std::unordered_map<std::string, std::string>> FR4_CONSENSUS_AA = {
            {
                    "H.SAPIENS", {
                                         {"hv", "WGQGTXVTVSS"},
                                         {"kv", "FGXGTKLEIK"},
                                         {"lv", "FGXGTKLTVL"}
                                 }
            },
    };
    static std::unordered_map<std::string, std::unordered_map<std::string, std::string>> FR4_CONSENSUS_DNA = {
            {
                    "H.SAPIENS", {
                                         {"hv", "TGGGGCCAGGGCACCNNNGTGACCGTGAGCAGC"},
                                         {"kv", "TTTGGCCAGGGGACCAAGCTGGAGATCAAA"},
                                         {"lv", "TTCGGCGGAGGGACCAAGCTGACCGTCCTA"}
                                 }
            },
    };
    // first, find all the matching positions of this pattern
    auto jaa = immulator::translate(jgerm);

    /*
     *   # https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
     *   blosum62 = """\
     *      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
     *   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
     *   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
     *   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
     *   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
     *   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
     *   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
     *   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
     *   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
     *   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
     *   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
     *   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
     *   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
     *   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
     *   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
     *   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
     *   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
     *   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
     *   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
     *   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
     *   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
     *   B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
     *   Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
     *   X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
     *   * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
     *   """
     *   tokens = {}
     *   rows = blosum62.split('\n')
     *   ref = rows[0].split()
     *   rest = rows[1:-1]
     *   print(ref)
     *   for row in rest:
     *       aa, scores = row.split()[0], row.split()[1:]
     *       tokens[aa] = {}
     *       for i, r in enumerate(ref):
     *           tokens[aa][r] = scores[i]
     *   print("{")
     *   for aa in tokens:
     *       print("{{'{}', {{".format(aa), end='')
     *       for maa, score in tokens[aa].items():
     *           print("{{'{}', {}}}".format(maa, score), end=',')
     *       print("}},")
     *   print("};")
     */
    // use above python program to generate the table below:
    static std::unordered_map<char, std::unordered_map<char, double>> scoring_matrix = {
            {'A', {{'A', 4},  {'R', -1}, {'N', -2}, {'D', -2}, {'C', 0},  {'Q', -1}, {'E', -1}, {'G', 0},  {'H', -2}, {'I', -1}, {'L', -1}, {'K', -1}, {'M', -1}, {'F', -2}, {'P', -1}, {'S', 1},  {'T', 0},  {'W', -3}, {'Y', -2}, {'V', 0},  {'B', -2}, {'Z', -1}, {'X', 0},  {'*', -4},}},
            {'R', {{'A', -1}, {'R', 5},  {'N', 0},  {'D', -2}, {'C', -3}, {'Q', 1},  {'E', 0},  {'G', -2}, {'H', 0},  {'I', -3}, {'L', -2}, {'K', 2},  {'M', -1}, {'F', -3}, {'P', -2}, {'S', -1}, {'T', -1}, {'W', -3}, {'Y', -2}, {'V', -3}, {'B', -1}, {'Z', 0},  {'X', -1}, {'*', -4},}},
            {'N', {{'A', -2}, {'R', 0},  {'N', 6},  {'D', 1},  {'C', -3}, {'Q', 0},  {'E', 0},  {'G', 0},  {'H', 1},  {'I', -3}, {'L', -3}, {'K', 0},  {'M', -2}, {'F', -3}, {'P', -2}, {'S', 1},  {'T', 0},  {'W', -4}, {'Y', -2}, {'V', -3}, {'B', 3},  {'Z', 0},  {'X', -1}, {'*', -4},}},
            {'D', {{'A', -2}, {'R', -2}, {'N', 1},  {'D', 6},  {'C', -3}, {'Q', 0},  {'E', 2},  {'G', -1}, {'H', -1}, {'I', -3}, {'L', -4}, {'K', -1}, {'M', -3}, {'F', -3}, {'P', -1}, {'S', 0},  {'T', -1}, {'W', -4}, {'Y', -3}, {'V', -3}, {'B', 4},  {'Z', 1},  {'X', -1}, {'*', -4},}},
            {'C', {{'A', 0},  {'R', -3}, {'N', -3}, {'D', -3}, {'C', 9},  {'Q', -3}, {'E', -4}, {'G', -3}, {'H', -3}, {'I', -1}, {'L', -1}, {'K', -3}, {'M', -1}, {'F', -2}, {'P', -3}, {'S', -1}, {'T', -1}, {'W', -2}, {'Y', -2}, {'V', -1}, {'B', -3}, {'Z', -3}, {'X', -2}, {'*', -4},}},
            {'Q', {{'A', -1}, {'R', 1},  {'N', 0},  {'D', 0},  {'C', -3}, {'Q', 5},  {'E', 2},  {'G', -2}, {'H', 0},  {'I', -3}, {'L', -2}, {'K', 1},  {'M', 0},  {'F', -3}, {'P', -1}, {'S', 0},  {'T', -1}, {'W', -2}, {'Y', -1}, {'V', -2}, {'B', 0},  {'Z', 3},  {'X', -1}, {'*', -4},}},
            {'E', {{'A', -1}, {'R', 0},  {'N', 0},  {'D', 2},  {'C', -4}, {'Q', 2},  {'E', 5},  {'G', -2}, {'H', 0},  {'I', -3}, {'L', -3}, {'K', 1},  {'M', -2}, {'F', -3}, {'P', -1}, {'S', 0},  {'T', -1}, {'W', -3}, {'Y', -2}, {'V', -2}, {'B', 1},  {'Z', 4},  {'X', -1}, {'*', -4},}},
            {'G', {{'A', 0},  {'R', -2}, {'N', 0},  {'D', -1}, {'C', -3}, {'Q', -2}, {'E', -2}, {'G', 6},  {'H', -2}, {'I', -4}, {'L', -4}, {'K', -2}, {'M', -3}, {'F', -3}, {'P', -2}, {'S', 0},  {'T', -2}, {'W', -2}, {'Y', -3}, {'V', -3}, {'B', -1}, {'Z', -2}, {'X', -1}, {'*', -4},}},
            {'H', {{'A', -2}, {'R', 0},  {'N', 1},  {'D', -1}, {'C', -3}, {'Q', 0},  {'E', 0},  {'G', -2}, {'H', 8},  {'I', -3}, {'L', -3}, {'K', -1}, {'M', -2}, {'F', -1}, {'P', -2}, {'S', -1}, {'T', -2}, {'W', -2}, {'Y', 2},  {'V', -3}, {'B', 0},  {'Z', 0},  {'X', -1}, {'*', -4},}},
            {'I', {{'A', -1}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', -1}, {'Q', -3}, {'E', -3}, {'G', -4}, {'H', -3}, {'I', 4},  {'L', 2},  {'K', -3}, {'M', 1},  {'F', 0},  {'P', -3}, {'S', -2}, {'T', -1}, {'W', -3}, {'Y', -1}, {'V', 3},  {'B', -3}, {'Z', -3}, {'X', -1}, {'*', -4},}},
            {'L', {{'A', -1}, {'R', -2}, {'N', -3}, {'D', -4}, {'C', -1}, {'Q', -2}, {'E', -3}, {'G', -4}, {'H', -3}, {'I', 2},  {'L', 4},  {'K', -2}, {'M', 2},  {'F', 0},  {'P', -3}, {'S', -2}, {'T', -1}, {'W', -2}, {'Y', -1}, {'V', 1},  {'B', -4}, {'Z', -3}, {'X', -1}, {'*', -4},}},
            {'K', {{'A', -1}, {'R', 2},  {'N', 0},  {'D', -1}, {'C', -3}, {'Q', 1},  {'E', 1},  {'G', -2}, {'H', -1}, {'I', -3}, {'L', -2}, {'K', 5},  {'M', -1}, {'F', -3}, {'P', -1}, {'S', 0},  {'T', -1}, {'W', -3}, {'Y', -2}, {'V', -2}, {'B', 0},  {'Z', 1},  {'X', -1}, {'*', -4},}},
            {'M', {{'A', -1}, {'R', -1}, {'N', -2}, {'D', -3}, {'C', -1}, {'Q', 0},  {'E', -2}, {'G', -3}, {'H', -2}, {'I', 1},  {'L', 2},  {'K', -1}, {'M', 5},  {'F', 0},  {'P', -2}, {'S', -1}, {'T', -1}, {'W', -1}, {'Y', -1}, {'V', 1},  {'B', -3}, {'Z', -1}, {'X', -1}, {'*', -4},}},
            {'F', {{'A', -2}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', -2}, {'Q', -3}, {'E', -3}, {'G', -3}, {'H', -1}, {'I', 0},  {'L', 0},  {'K', -3}, {'M', 0},  {'F', 6},  {'P', -4}, {'S', -2}, {'T', -2}, {'W', 1},  {'Y', 3},  {'V', -1}, {'B', -3}, {'Z', -3}, {'X', -1}, {'*', -4},}},
            {'P', {{'A', -1}, {'R', -2}, {'N', -2}, {'D', -1}, {'C', -3}, {'Q', -1}, {'E', -1}, {'G', -2}, {'H', -2}, {'I', -3}, {'L', -3}, {'K', -1}, {'M', -2}, {'F', -4}, {'P', 7},  {'S', -1}, {'T', -1}, {'W', -4}, {'Y', -3}, {'V', -2}, {'B', -2}, {'Z', -1}, {'X', -2}, {'*', -4},}},
            {'S', {{'A', 1},  {'R', -1}, {'N', 1},  {'D', 0},  {'C', -1}, {'Q', 0},  {'E', 0},  {'G', 0},  {'H', -1}, {'I', -2}, {'L', -2}, {'K', 0},  {'M', -1}, {'F', -2}, {'P', -1}, {'S', 4},  {'T', 1},  {'W', -3}, {'Y', -2}, {'V', -2}, {'B', 0},  {'Z', 0},  {'X', 0},  {'*', -4},}},
            {'T', {{'A', 0},  {'R', -1}, {'N', 0},  {'D', -1}, {'C', -1}, {'Q', -1}, {'E', -1}, {'G', -2}, {'H', -2}, {'I', -1}, {'L', -1}, {'K', -1}, {'M', -1}, {'F', -2}, {'P', -1}, {'S', 1},  {'T', 5},  {'W', -2}, {'Y', -2}, {'V', 0},  {'B', -1}, {'Z', -1}, {'X', 0},  {'*', -4},}},
            {'W', {{'A', -3}, {'R', -3}, {'N', -4}, {'D', -4}, {'C', -2}, {'Q', -2}, {'E', -3}, {'G', -2}, {'H', -2}, {'I', -3}, {'L', -2}, {'K', -3}, {'M', -1}, {'F', 1},  {'P', -4}, {'S', -3}, {'T', -2}, {'W', 11}, {'Y', 2},  {'V', -3}, {'B', -4}, {'Z', -3}, {'X', -2}, {'*', -4},}},
            {'Y', {{'A', -2}, {'R', -2}, {'N', -2}, {'D', -3}, {'C', -2}, {'Q', -1}, {'E', -2}, {'G', -3}, {'H', 2},  {'I', -1}, {'L', -1}, {'K', -2}, {'M', -1}, {'F', 3},  {'P', -3}, {'S', -2}, {'T', -2}, {'W', 2},  {'Y', 7},  {'V', -1}, {'B', -3}, {'Z', -2}, {'X', -1}, {'*', -4},}},
            {'V', {{'A', 0},  {'R', -3}, {'N', -3}, {'D', -3}, {'C', -1}, {'Q', -2}, {'E', -2}, {'G', -3}, {'H', -3}, {'I', 3},  {'L', 1},  {'K', -2}, {'M', 1},  {'F', -1}, {'P', -2}, {'S', -2}, {'T', 0},  {'W', -3}, {'Y', -1}, {'V', 4},  {'B', -3}, {'Z', -2}, {'X', -1}, {'*', -4},}},
            {'B', {{'A', -2}, {'R', -1}, {'N', 3},  {'D', 4},  {'C', -3}, {'Q', 0},  {'E', 1},  {'G', -1}, {'H', 0},  {'I', -3}, {'L', -4}, {'K', 0},  {'M', -3}, {'F', -3}, {'P', -2}, {'S', 0},  {'T', -1}, {'W', -4}, {'Y', -3}, {'V', -3}, {'B', 4},  {'Z', 1},  {'X', -1}, {'*', -4},}},
            {'Z', {{'A', -1}, {'R', 0},  {'N', 0},  {'D', 1},  {'C', -3}, {'Q', 3},  {'E', 4},  {'G', -2}, {'H', 0},  {'I', -3}, {'L', -3}, {'K', 1},  {'M', -1}, {'F', -3}, {'P', -1}, {'S', 0},  {'T', -1}, {'W', -3}, {'Y', -2}, {'V', -2}, {'B', 1},  {'Z', 4},  {'X', -1}, {'*', -4},}},
            {'X', {{'A', 0},  {'R', -1}, {'N', -1}, {'D', -1}, {'C', -2}, {'Q', -1}, {'E', -1}, {'G', -1}, {'H', -1}, {'I', -1}, {'L', -1}, {'K', -1}, {'M', -1}, {'F', -1}, {'P', -2}, {'S', 0},  {'T', 0},  {'W', -2}, {'Y', -1}, {'V', -1}, {'B', -1}, {'Z', -1}, {'X', -1}, {'*', -4},}},
            {'*', {{'A', -4}, {'R', -4}, {'N', -4}, {'D', -4}, {'C', -4}, {'Q', -4}, {'E', -4}, {'G', -4}, {'H', -4}, {'I', -4}, {'L', -4}, {'K', -4}, {'M', -4}, {'F', -4}, {'P', -4}, {'S', -4}, {'T', -4}, {'W', -4}, {'Y', -4}, {'V', -4}, {'B', -4}, {'Z', -4}, {'X', -4}, {'*', 1},}},
    };

    // try getting AA position first
    std::string::size_type start, end;
    std::tie(std::ignore, start, end) = immulator::local_align(jaa, FR4_CONSENSUS_AA["H.SAPIENS"]["hv"], -5, -5,
                                                               scoring_matrix);
    if (start < end) {
        // convert to NT start position
        start *= 3;
        end *= 3;
        std::cerr << "AMINO DISCOVERY: " << start << std::endl;
    } else {
        // try AA position
        static constexpr double NT_MATCH = 5;
        static constexpr double NT_MISMATCH = -5;
        static std::unordered_map<char, std::unordered_map<char, double>> nt_scoring_matrix = {
                {'A', {{'A', NT_MATCH},    {'C', NT_MISMATCH}, {'G', NT_MISMATCH}, {'T', NT_MISMATCH}}},
                {'C', {{'A', NT_MISMATCH}, {'C', NT_MATCH},    {'G', NT_MISMATCH}, {'T', NT_MISMATCH}}},
                {'G', {{'A', NT_MISMATCH}, {'C', NT_MISMATCH}, {'G', NT_MATCH},    {'T', NT_MISMATCH}}},
                {'T', {{'A', NT_MISMATCH}, {'C', NT_MISMATCH}, {'G', NT_MISMATCH}, {'T', NT_MATCH}}},
        };
        std::tie(std::ignore, start, end) = immulator::local_align(jgerm,
                                                                   FR4_CONSENSUS_DNA["H.SAPIENS"]["hv"], -5, -5,
                                                                   nt_scoring_matrix);
        if (start < end) {
            std::cerr << "Tried nucleotide consensus FR4 region with no luck\n";
            return {};
        }
        std::cerr << "NT SE DISCOVERY: " << start << std::endl;
    }

    /* -------------------------------------------------------------------------------- *
     *                       Determine how to cut the front nt seqs                     *
     *                                                                                  *
     * -------------------------------------------------------------------------------- */
    constexpr double FRONT_CUT_PERC = 30.0 / 100;
    auto max_front_cut_size = static_cast<size_type>(std::ceil(FRONT_CUT_PERC * jgerm.size()));
    std::uniform_int_distribution<size_type> front_idist(0, std::min(max_front_cut_size, start));

    auto front_cut = front_idist(generator);
    if (check) {
        auto aa = immulator::translate(rem + jgerm.substr(front_cut));
        std::size_t attempt = 0;
        for (; attempt < MAX_ATTEMPTS && aa.find('*') != std::string::npos; ++attempt) {
            front_cut = front_idist(generator);
            aa = immulator::translate(rem + jgerm.substr(front_cut));
        }
        if (attempt == MAX_ATTEMPTS) {
            std::cerr << "WARNING: Tried too hard, but in the end, nothing matters.\n";
            productive = false;
        }
    }

    /* -------------------------------------------------------------------------------- *
     *                       Determine how to cut the front nt seqs                     *
     *                                                                                  *
     * -------------------------------------------------------------------------------- */
    constexpr double BACK_CUT_PERC = 30.0 / 100;
    auto max_back_cut_size = static_cast<size_type>(std::ceil(BACK_CUT_PERC * jgerm.size()));
    std::uniform_int_distribution<size_type> back_idist(0, max_back_cut_size);

    auto back_cut = back_idist(generator);
    if (multiple) {
        if (auto rem_nt = (jgerm.size() - front_cut + rem.size() - back_cut) % 3) {
            back_cut += rem_nt;
        }
    }
    assert(front_cut <= start);
    return std::make_tuple(jgerm.trim(front_cut, jgerm.size() - back_cut - front_cut), start - front_cut, productive);
}


