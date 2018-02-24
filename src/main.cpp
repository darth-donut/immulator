#include <iostream>
#include <string>
#include "germline_factory.h"

using std::string;

int
main() {
    constexpr int seqs = 1'000'000;
    immulator::GermlineConfiguration gcfg("../germ.cfg", true);
    const string title(80, '=');
    std::cout << title << '\n'
              << "\t\t\tConfiguration file found\n" << title  << '\n'
              << gcfg << std::endl;
    immulator::GermlineFactory vgermlines("../tmp.txt", gcfg);
    std::map<std::string, int> counter;
    for (int i = 0; i < seqs; i++) {
        counter[vgermlines().family_name()] += 1;
    }

    std::cout << "\nTesting recombination dummy: " << (vgermlines() + vgermlines()) << std::endl;

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