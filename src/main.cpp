#include <iostream>
#include <string>
#include "germline_factory.h"
#include <random>

using std::string;
using immulator::Germline;

// Testing
Germline vdj_recombination(const Germline &vgerm, const Germline &dgerm, const Germline &jgerm);
template<typename Gen> immulator::optional<std::pair<Germline, immulator::Germline::size_type>> vcutter(Germline vgerm, Gen &generator);

int
main() {
    constexpr int seqs = 1'000;
    std::mt19937 mersenne(std::random_device{}());

    immulator::GermlineConfiguration gcfg("../germ.cfg", true);
    const string title(80, '=');
    std::cout << title << '\n'
              << "\t\t\tConfiguration file found\n" << title  << '\n'
              << gcfg << std::endl;
    immulator::GermlineFactory vgermlines("../tmp.txt"/*, gcfg*/);
    std::map<std::string, int> counter;
    for (int i = 0; i < seqs; i++) {
        counter[vgermlines(mersenne).family_name()] += 1;
    }

//    std::cout << "\nTesting recombination dummy: " << (vgermlines(mersenne) + vgermlines(mersenne)) << std::endl;
    std::cout << "\nTesting recombination dummy: "
              << vdj_recombination(vgermlines(mersenne) , vgermlines(mersenne), vgermlines(mersenne))
              << std::endl;

    // count germline distribution used
    double total = std::accumulate(counter.begin(), counter.end(), 0, [] (auto val, auto &keypair) {
        return val + keypair.second;
    });
    std::cout << title << '\n'
              << "\t\t\tCalculating germline distribution ...\n"
              << title << std::endl;

    for (auto &i : counter) {
        std::cout << i.first << "\t" << i.second/total << std::endl;
    }
    return 0;
}



// Testing
Germline
vdj_recombination(const Germline &vgerm, const Germline &dgerm, const Germline &jgerm) {
    static std::mt19937 mersenne(std::random_device{}());
    auto v = vcutter(vgerm, mersenne);
    if (v) {
        // todo
    } else {
        // todo
    }
    return vgerm + dgerm + jgerm;
}

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
        return std::make_pair(vgerm.trim(0, final_length), final_length - 1);
    }
}
