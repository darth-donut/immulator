#include <iostream>
#include <string>
#include "germline_factory.h"
#include <random>
#include <cmath>
#include <tuple>
#include <regex>

using std::string;
using immulator::Germline;

// Testing
Germline vdj_recombination(const Germline &vgerm, const Germline &dgerm, const Germline &jgerm, bool prod = true);

template<typename Gen>
immulator::optional<std::pair<Germline, immulator::Germline::size_type>> vcutter(Germline vgerm, Gen &generator);

template<typename Gen>
std::tuple<Germline, immulator::Germline::size_type, bool>
dcutter(Germline dgerm, Gen &generator, const std::string &rem, bool check = true);

template<typename Gen>
immulator::optional<std::tuple<Germline, immulator::Germline::size_type, bool>>
jcutter(Germline jgerm, Gen &generator, const std::string &rem, bool check = true);

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

//    std::cout << "\nTesting recombination dummy: " << (vgermlines(mersenne) + vgermlines(mersenne)) << std::endl;
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
vdj_recombination(const Germline &vgerm, const Germline &dgerm, const Germline &jgerm, bool prod) {
    using size_type = immulator::Germline::size_type;
    static std::mt19937 mersenne(std::random_device{}());
    auto v = vcutter(vgerm, mersenne);
    if (v) {
        // todo:: V - D INSERTION, D-J insertion
        Germline d;
        size_type dsize;
        bool d_prod;
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
            auto jtry = jcutter(jgerm, mersenne, d.substr(d.size() - (v->first.size() + d.size()) % 3), prod);
            if (jtry) {
                std::tie(j, fwgxg_conserved_index, j_prod) = *jtry;
                std::cerr << "HERE: " << std::endl;
                std::cerr << v->first + d + j << std::endl;
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
        std::uniform_int_distribution<int> idist(0, nt_rem);
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

    auto max_front_cut_size = static_cast<int>(std::ceil(FRONT_CUT_PERC * dgerm.size()));
    std::uniform_int_distribution<int> front_idist(0, max_front_cut_size);

    // cut front of D gene by "front_cut" much
    int front_cut = front_idist(generator);
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
    auto max_back_cut_size = static_cast<int>(std::ceil(BACK_CUT_PERC * dgerm.size()));
    std::uniform_int_distribution<int> back_idist(0, max_back_cut_size);

    // cut back of D gene by "back_cut" much
    int back_cut = front_idist(generator);
    return std::make_tuple(dgerm.trim(front_cut, dgerm.size() - back_cut - front_cut),
                           dgerm.size() - front_cut - back_cut,
                           productive);
}

template<typename Gen>
immulator::optional<std::tuple<Germline, immulator::Germline::size_type, bool>>
jcutter(Germline jgerm, Gen &generator, const std::string &rem, bool check) {
    using size_type = immulator::Germline::size_type;
    constexpr std::size_t MAX_ATTEMPTS = 100'000;

    bool productive = true;

//    FR4_CONSENSUS = {'hv':"WGQGTXVTVSS", 'kv':'FGGGTQ', 'lv':'FGGGTQ'}
//    FR4_CONSENSUS_DNA = {'hv':"TGGGGCCAGGGCACCNNNGTGACCGTGAGCAGC",
//    'kv':'TTCGGCGGCGGCACCCAG', 'lv':'TTCGGCGGCGGCACCCAG'}

    /* ------------------------------------------------------------------------------ *
     *                          Find consensus FR4 [FW]G.G                            *
     *                                                                                *
     * ------------------------------------------------------------------------------ */

    // first, find all the matching positions of this pattern
    const std::regex CONSERVED_PTN(R"([FW]G.G)");
    std::smatch match;
    std::vector<std::string::size_type> match_indices;
    auto jaa = immulator::translate(jgerm);
    if (std::regex_search(jaa, match, CONSERVED_PTN)) {
        for (auto i = 0; i < match.size(); ++i) {
            match_indices.push_back(match.position(i) * 3);
        }
    } else {
        // can't find a consensus!
        std::cerr << "J gene" << jaa << "\n\t has no consensus fwgxg\n";
        return {};
    }

    // next, get the first of that, and generate a cut value based on a uniform int dist
    std::uniform_int_distribution<> idist(0, match_indices.front());
    int front_cut = idist(generator);
    if (check) {
        auto aa = immulator::translate(rem + jgerm.substr(front_cut));
        std::size_t attempt = 0;
        for (; attempt < MAX_ATTEMPTS && aa.find('*') != std::string::npos; ++attempt) {
            front_cut = idist(generator);
            aa = immulator::translate(rem + jgerm.substr(front_cut));
        }
        if (attempt == MAX_ATTEMPTS) {
            std::cerr << "WARNING: Tried too hard, but in the end, nothing matters.\n";
            productive = false;
        }
    }

    return std::make_tuple(jgerm.trim(front_cut), match_indices.front(), productive);
}



