#include <iostream>
#include <string>
#include <random>
#include <cassert>
#include <unordered_map>
#include <cmath>
#include <tuple>
#include <regex>
#include <fstream>

#include "cxxopts.hpp"
#include "germline_factory.h"

#define VERSION "Immulator v0.0.99"

using std::string;
using immulator::Germline;

// Testing
std::tuple<immulator::optional<Germline>, immulator::Germline::size_type, immulator::Germline::size_type>
vdj_recombination(const Germline &vgerm, const Germline &dgerm, const Germline &jgerm, bool prod = true,
                  bool multiple = true);

template<typename Gen>
immulator::optional<std::pair<Germline, immulator::Germline::size_type>> vcutter(Germline vgerm, Gen &generator);

template<typename Gen>
std::tuple<Germline, immulator::Germline::size_type, bool>
dcutter(Germline dgerm, Gen &generator, const std::string &rem, bool check = true);

template<typename Gen>
immulator::optional<std::tuple<Germline, immulator::Germline::size_type, bool>>
jcutter(Germline jgerm, Gen &generator, const std::string &rem,
        std::string::size_type extras, bool check);

template<typename Gen>
immulator::optional<std::string>
palindromic(std::string::size_type n, Gen &generator, const std::string &rem, bool productive = true);

template<typename Gen>
std::string
random_nts(std::string::size_type n, Gen &generator, const std::string &rem, bool productive = true);

void
write_reference(std::ofstream &os, const std::string &name,
                immulator::Germline::size_type cdr3_start,
                immulator::Germline::size_type cdr3_end);

int
main(int argc, char *argv[]) {
    auto seed = std::random_device{}();
    std::string reference_filename("immulator.csv");
    cxxopts::Options options(argv[0], "Immunoglobulin simulator - simulates V region antibody sequences.");
    options.add_options()
            ("n,num", "number of sequences to simulate", cxxopts::value<std::size_t>())
            ("p,productive", "approx percentage of productive/functional sequences to simulate",
             cxxopts::value<double>())
            ("s,seed", "seed random generator; keep this between the range of"
                       " 0 to 2^32", cxxopts::value<unsigned int>())
            ("v,version", "print immulator version and exits")
            ("h,help", "print this help message and exits")
            ("g,germlinecfg", "comma separated germline configuration file; describes V-(D)-J germline "
                              "distributions", cxxopts::value<std::string>())
            ("r,reference", "germline and CDR3 information will be saved in this file, defaults "
                            "to immulator.csv", cxxopts::value<std::string>())
            ("V,vdb", "path to V germline database (FASTA file)", cxxopts::value<std::string>())
            ("D,ddb", "path to D germline database (FASTA file)", cxxopts::value<std::string>())
            ("J,jdb", "path to J germline database (FASTA file)", cxxopts::value<std::string>());
    auto args = options.parse(argc, argv);
    double prod_perc = .75;
    std::string vdb, ddb, jdb;

    if (args.count("help")) {
        std::cout << options.help() << std::endl;
        return (EXIT_SUCCESS);
    }
    if (args.count("version")) {
        std::cout << VERSION << std::endl;
        return (EXIT_SUCCESS);
    }
    if (!args.count("num")) {
        std::cout << options.help() << std::endl;
        return (EXIT_FAILURE);
    }
    if (args.count("seed")) {
        seed = args["seed"].as<unsigned int>();
    }
    if (args.count("productive")) {
        prod_perc = args["productive"].as<double>() / 100;
    }
    if (args.count("reference")) {
        reference_filename = args["reference"].as<std::string>();
    }
    std::unordered_map<char, std::string> database;
    for (const char c : std::string("vdj")) {
        auto key = std::string(1, c) + "db";
        if (args.count(key)) {
            database[c] = args[key].as<std::string>();
            if (!immulator::exists(database[c])) {
                std::cout << database[c] << " file not found. Aborting." << std::endl;
                return(EXIT_FAILURE);
            }
        } else {
            std::cout << "-" << static_cast<char>(std::toupper(c)) << " <path to "
                      << static_cast<char>(std::toupper(c))
                      << " germline FASTA file> is required" << std::endl;
            return(EXIT_FAILURE);
        }
    }
    std::cerr << "This simulation run is generated with seed " << seed << std::endl;
    std::mt19937 mersenne(seed);
    std::ofstream refos(reference_filename);

    if (args.count("germlinecfg")) {
        const std::size_t seqs = args["num"].as<std::size_t>();
        immulator::GermlineConfiguration gcfg(args["germlinecfg"].as<std::string>(), true);
        const string title(80, '=');
        std::cerr << title << '\n'
                  << "\t\t\tConfiguration file found\n" << title << '\n'
                  << gcfg << std::endl;
        immulator::GermlineFactory vgermlines(database['v'], gcfg, false);
        immulator::GermlineFactory dgermlines(database['d'], gcfg, false);
        immulator::GermlineFactory jgermlines(database['j'], gcfg, false);
        for (auto i = 0; i < seqs; ++i) {
            immulator::optional<Germline> recombined;
            immulator::Germline::size_type cdr3_start, cdr3_end;
            do {
                std::tie(recombined, cdr3_start, cdr3_end) = vdj_recombination(vgermlines(mersenne),
                                                                               dgermlines(mersenne),
                                                                               jgermlines(mersenne),
                                                                               immulator::coin_flip(mersenne, prod_perc)
                                                                               );
            } while (!recombined);
            std::cout << ">" << i << *recombined << std::endl;
            write_reference(refos, recombined->name(), cdr3_start, cdr3_end);
        }

    } else {
        const std::size_t seqs = args["num"].as<std::size_t>();
        immulator::GermlineFactory vgermlines(database['v'], false);
        immulator::GermlineFactory dgermlines(database['d'], false);
        immulator::GermlineFactory jgermlines(database['j'], false);
        for (auto i = 0; i < seqs; ++i) {
            immulator::optional<Germline> recombined;
            immulator::Germline::size_type cdr3_start, cdr3_end;
            do {
                std::tie(recombined, cdr3_start, cdr3_end) = vdj_recombination(vgermlines(mersenne),
                                                                               dgermlines(mersenne),
                                                                               jgermlines(mersenne),
                                                                               immulator::coin_flip(mersenne, prod_perc)
                                                                               );
            } while (!recombined);
            std::cout << ">" << i << *recombined << std::endl;
            write_reference(refos, recombined->name(), cdr3_start, cdr3_end);
        }
    }
    return (EXIT_SUCCESS);
}


// Testing
std::tuple<immulator::optional<Germline>, immulator::Germline::size_type, immulator::Germline::size_type>
vdj_recombination(const Germline &vgerm, const Germline &dgerm, const Germline &jgerm, bool prod, bool multiple) {
    using size_type = immulator::Germline::size_type;
    static std::mt19937 mersenne(std::random_device{}());
    static constexpr std::size_t MAX_ATTEMPTS = 100'000;
    auto v = vcutter(vgerm, mersenne);
    std::size_t attempts_insertion = 0;
    Germline buffer;
    std::uniform_int_distribution<std::string::size_type> palin_rand(0, 8);
    std::uniform_int_distribution<std::string::size_type> ins_rand(0, 5);

    if (v) {
        buffer = v->first;
        bool d_prod = false;
        Germline d;
        size_type dsize;
        // starts AFTER Cys (and convert to 1-index)
        size_type cdr3_start_pos = v->second + 3 + 1;
        std::size_t attempt_d = 0;
        auto p1 = palindromic(palin_rand(mersenne), mersenne, buffer.remainder(), prod);
        while (!p1 && prod && ++attempts_insertion < MAX_ATTEMPTS) {
            p1 = palindromic(palin_rand(mersenne), mersenne, buffer.remainder(), prod);
        }
        buffer += *p1;
        auto n1 = random_nts(ins_rand(mersenne), mersenne, buffer.remainder(), prod);
        buffer += n1;
        auto p2 = palindromic(palin_rand(mersenne), mersenne, buffer.remainder(), prod);
        while (!p2 && prod && ++attempts_insertion < MAX_ATTEMPTS) {
            p2 = palindromic(palin_rand(mersenne), mersenne, buffer.remainder(), prod);
        }
        buffer += *p2;
        do {
            std::tie(d, dsize, d_prod) = dcutter(dgerm,
                                                 mersenne,
                                                 buffer.remainder(),
                                                 prod);
            ++attempt_d;
        } while (!d_prod && prod && attempt_d < MAX_ATTEMPTS);
        buffer += d;
        auto p3 = palindromic(palin_rand(mersenne), mersenne, buffer.remainder(), prod);
        while (!p3 && prod && ++attempts_insertion < MAX_ATTEMPTS) {
            p3 = palindromic(palin_rand(mersenne), mersenne, buffer.remainder(), prod);
        }
        buffer += *p3;
        auto n2 = random_nts(ins_rand(mersenne), mersenne, buffer.remainder(), prod);
        buffer += n2;
        auto p4 = palindromic(palin_rand(mersenne), mersenne, buffer.remainder(), prod);
        while (!p4 && prod && ++attempts_insertion < MAX_ATTEMPTS) {
            p4 = palindromic(ins_rand(mersenne), mersenne, buffer.remainder(), prod);
        }
        buffer += *p4;
        auto current_incomplete_cdr3_length = (v->first.size() - cdr3_start_pos + 1)
                                              + p1->size() + n1.size() + p2->size()
                                              + d.size() + p3->size() + n2.size() + p4->size();
        Germline j;
        std::size_t attempt_j = 0;
        size_type fwgxg_conserved_index;
        auto jtry = jcutter(jgerm, mersenne, buffer.remainder(),
                            (3 - (current_incomplete_cdr3_length % 3)) % 3,
                            prod);

        if (jtry) {
            bool j_prod = false;
            std::tie(j, fwgxg_conserved_index, j_prod) = *jtry;
            while (!j_prod && prod && attempt_j++ < MAX_ATTEMPTS) {
                jtry = jcutter(jgerm, mersenne,
                               buffer.remainder(),
                               (3 - (current_incomplete_cdr3_length % 3)) % 3,
                               prod);
                if (!jtry) {
                    // fail to find J gene anchor - fail immediately
                    return {};
                }
                std::tie(j, fwgxg_conserved_index, j_prod) = *jtry;
            }
            size_type cdr3_end_pos = buffer.size() + fwgxg_conserved_index;
            buffer += j;
            if (multiple && (buffer.size() % 3)) {
                buffer.trim(0, buffer.size() - (buffer.size() % 3));
                assert(buffer.size() % 3 == 0);
            }
            return std::make_tuple(buffer, cdr3_start_pos, cdr3_end_pos);
        } else {
            // fail to find J gene anchor [FW]G.G region
            return {};
        }

    } else {
        // fail to find anchor Cys region
        return {};
    }
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
        std::cerr << "WARNING: Cys anchor failed to be located in:\n"
                  << "\t" << vgerm << '\n';
        return {};
    } else {
        size_type nuc_index = cys * 3;

        // V germlines are usually > 200 (actually, >250)
        if (nuc_index < 200) {
            return {};
        }

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
    auto back_cut = back_idist(generator);
    return std::make_tuple(dgerm.trim(front_cut, dgerm.size() - back_cut - front_cut),
                           dgerm.size() - front_cut - back_cut,
                           productive);
}

template<typename Gen>
immulator::optional<std::tuple<Germline, immulator::Germline::size_type, bool>>
jcutter(Germline jgerm, Gen &generator, const std::string &rem,
        std::string::size_type extras, bool check) {
    assert(extras >= 0 && extras <= 2 && "Extras is expected to be an integer between 0 and 2 inclusive");
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
    size_type orf = 0, used_orf = 0;
    double best_score = 0;

    for (; orf < 3; ++orf) {
        double score;
        std::string::size_type current_start, current_end;
        std::tie(score, current_start, current_end) = immulator::local_align(jaa, FR4_CONSENSUS_AA["H.SAPIENS"]["hv"],
                                                                             -5, -5,
                                                                             scoring_matrix);
        if (score > best_score) {
            start = current_start;
            end = current_end;
            best_score = score;
            used_orf = orf;
        }
        jaa = immulator::translate(jgerm.substr(orf + 1));
    }
    if (start < end) {
        // convert to NT start position
        start = start * 3 + used_orf;
        end = end * 3 + used_orf;
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
    }

    /* -------------------------------------------------------------------------------- *
     *                       Determine how to cut the front nt seqs                     *
     *                                                                                  *
     * -------------------------------------------------------------------------------- */
    constexpr double FRONT_CUT_PERC = 30.0 / 100;
    auto max_front_cut_size = static_cast<size_type>(std::ceil(FRONT_CUT_PERC * jgerm.size()));
    std::uniform_int_distribution<size_type> front_idist(0, std::min(max_front_cut_size, start));

    auto front_cut = front_idist(generator);
    // to maintain the V-J frame, FWGXG index - extras % 3 should be 0
    if (front_cut < start) {
        auto offset = (start - front_cut - extras) % 3;
        auto offset_by = (3 - offset) % 3;
        // if we can afford to trim the front or if we CAN'T extend the back, use the front
        if (offset_by <= front_cut && (immulator::coin_flip(generator) || front_cut + offset > start)) {
            front_cut -= offset_by;
        } else {
            front_cut += offset;
        }
    } else {
        assert(front_cut == start && extras <= start);
        // scale back to allow extras to consume the "scaled" back nt
        front_cut -= extras;
    }

    if (check) {
        auto aa = immulator::translate(rem + jgerm.substr(front_cut));
        std::size_t attempt = 0;
        for (; attempt < MAX_ATTEMPTS && aa.find('*') != std::string::npos; ++attempt) {
            front_cut = front_idist(generator);
            if (front_cut < start) {
                auto offset = (start - front_cut - extras) % 3;
                auto offset_by = (3 - offset) % 3;
                // 50% chance of offsetting either from the front or back, but if adding the offset to the back will
                // cause a trim on the conserved anchor position (start), then force offsetting to happen from the front
                // if offsetting from the front will cause a negative index (offset > front_cut), use the back as offset
                // regardless of whether or not we lose the conserved region
                if (offset_by <= front_cut && (immulator::coin_flip(generator) || front_cut + offset > start)) {
                    front_cut -= offset_by;
                } else {
                    front_cut += offset;
                }
            } else {
                assert(front_cut == start && extras <= start);
                front_cut -= extras;
            }
            aa = immulator::translate(rem + jgerm.substr(front_cut));
        }
        if (attempt == MAX_ATTEMPTS) {
            std::cerr << "WARNING: Tried too hard, but in the end, nothing matters.\n";
            productive = false;
        }
    }
    assert(front_cut > start || (start - front_cut - extras) % 3 == 0);
    /* -------------------------------------------------------------------------------- *
     *                       Determine how to cut the back nt seqs                      *
     *                                                                                  *
     * -------------------------------------------------------------------------------- */
    constexpr double BACK_CUT_PERC = 0 / 100;
    auto max_back_cut_size = static_cast<size_type>(std::ceil(BACK_CUT_PERC * jgerm.size()));
    std::uniform_int_distribution<size_type> back_idist(0, max_back_cut_size);

    auto back_cut = back_idist(generator);
    // when front_cut > start, it means we compensated V-J frame with additional cut INTO the conserved region,
    // so naturally CDR3 starts as early as 0
    return std::make_tuple(jgerm.trim(front_cut, jgerm.size() - back_cut - front_cut),
                           front_cut <= start ? start - front_cut : 0, productive);
}

template<typename Gen>
immulator::optional<std::string>
palindromic(std::string::size_type n, Gen &generator, const std::string &rem, bool productive) {
    assert(rem.size() <= 2);
    constexpr static char NTS[] = {'A', 'C', 'G', 'T'};
    constexpr static std::size_t MAX_ATTEMPTS = 100'000;
    static std::unordered_map<char, char> COMPLEMENT_NT = {
            {'A', 'T'},
            {'T', 'A'},
            {'C', 'G'},
            {'G', 'C'}
    };
    std::string nt_seq;

    if (!productive) {
        std::uniform_int_distribution<std::string::size_type> idist(0, 3);      // 0 to len(NTS) - 1
        // first half
        for (auto i = 0; i < n / 2; ++i) {
            nt_seq.push_back(NTS[idist(generator)]);
        }

        // the middle nucleotide (if n is odd)
        if (n % 2) {
            nt_seq.push_back(NTS[idist(generator)]);
        }

        // the remaining (second) half
        for (auto i = 0; i < n / 2; ++i) {
            nt_seq.push_back(COMPLEMENT_NT[nt_seq[n / 2 - i - 1]]);
        }
        return immulator::join_string(nt_seq.cbegin(), nt_seq.cend(), "");
    } else {
        bool is_productive = false;
        std::size_t attempts = 0;
        std::string aa_seq;
        do {
            std::string current_codon = rem;
            nt_seq.clear();
            // first half
            for (auto i = 0; i < n / 2; ++i) {
                auto allowed_nt = immulator::allowed_nts(current_codon);
                std::uniform_int_distribution<std::string::size_type> idist(0, allowed_nt.size() - 1);
                char nt = allowed_nt[idist(generator)];
                nt_seq.push_back(nt);
                current_codon += nt;
                current_codon = current_codon.size() == 3 ? "" : current_codon;
            }

            // middle nucleotide (when n is odd)
            if (n % 2) {
                auto allowed_nt = immulator::allowed_nts(current_codon);
                std::uniform_int_distribution<std::string::size_type> idist(0, allowed_nt.size() - 1);
                char nt = allowed_nt[idist(generator)];
                nt_seq.push_back(nt);
                current_codon += nt;
                current_codon = current_codon.size() == 3 ? "" : current_codon;
            }

            for (auto i = 0; i < n / 2; ++i) {
                nt_seq.push_back(COMPLEMENT_NT[nt_seq[n / 2 - i - 1]]);
            }
            aa_seq = immulator::join_string(nt_seq.cbegin(), nt_seq.cend(), "");
            is_productive = immulator::translate(aa_seq).find('*') == std::string::npos;
        } while (!is_productive && ++attempts < MAX_ATTEMPTS);
        return is_productive ? aa_seq : immulator::optional<std::string>();
    }

}


template<typename Gen>
std::string
random_nts(std::string::size_type n, Gen &generator, const std::string &rem, bool productive) {
    assert(rem.size() <= 2);
    constexpr static char NTS[] = {'A', 'C', 'G', 'T'};

    std::string nt_seq;

    if (!productive) {
        std::uniform_int_distribution<std::string::size_type> idist(0, 3);
        for (auto i = 0; i < n; ++i) {
            nt_seq += NTS[idist(generator)];
        }
    } else {
        std::string current_codon = rem;
        for (auto i = 0; i < n; ++i) {
            auto allowed_nt = immulator::allowed_nts(current_codon);
            std::uniform_int_distribution<std::string::size_type> idist(0, allowed_nt.size() - 1);
            char nt = allowed_nt[idist(generator)];
            nt_seq += nt;
            current_codon += nt;
            current_codon = current_codon.size() == 3 ? "" : current_codon;
        }
    }
    return nt_seq;
}


void
write_reference(std::ofstream &os, const std::string &name,
                immulator::Germline::size_type cdr3_start,
                immulator::Germline::size_type cdr3_end) {
    static bool first_ = true;
    if (first_) {
        first_ = false;
        // the header need only be written once (on the first call)
        os << "Genes,CDR3.start,CDR3.end\n";
    }
    os << name << "," << cdr3_start << "," << cdr3_end << '\n';
}
