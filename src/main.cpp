#include <iostream>
#include "germline_factory.h"

int main() {
    immulator::GermlineFactory vgermlines("../tmp.txt");
    std::cout << vgermlines().family_name() << std::endl;
    std::cout << vgermlines().gene_name() << std::endl;
    return 0;
}