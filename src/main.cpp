#include <iostream>
#include "germline_factory.h"

int
main() {
    immulator::GermlineConfiguration gcfg("../germ.cfg", true);
    std::cout << gcfg << std::endl;
    immulator::GermlineFactory vgermlines("../tmp.txt", gcfg);
    std::map<std::string, int> counter;
    for (int i = 0; i < 1'000'000; i++) {
        counter[vgermlines().family_name()] += 1;
    }
    double total = std::accumulate(counter.begin(), counter.end(), 0, [] (auto val, auto &keypair) { return val + keypair.second; });
    for (auto &i : counter) {
        std::cout << i.first << "\t" << i.second/total << std::endl;
    }
    return 0;
}