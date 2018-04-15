//
// @author: jiahong
// @date  : 24/02/18 4:33 PM
//

#include <fstream>
#include <iostream>
#include "germline_factory.h"

using immulator::Germline;


void
immulator::GermlineFactory::parse_file() {
    std::ifstream ifs(filename_);
    if (ifs) {
        std::string buffer;
        std::getline(ifs, buffer);
        while (!buffer.empty()) {
            // keep "getting" until we reach the first entry
            while (buffer.find_first_of('>') != 0 && std::getline(ifs, buffer));

//            auto tokens = immulator::split_string(buffer, "|");
//            auto gene_name = tokens[1];
//            auto asc_name = tokens[0].substr(1);
//
//            // keep "getting" until we reach the next entry
//            std::string sequence;
//            while (std::getline(ifs, buffer) && buffer.find_first_of('>') != 0) {
//                sequence += buffer;
//            }
//            germline_collection_.emplace_back(std::move(gene_name), std::move(asc_name), std::move(sequence));
            auto gene_name = buffer.substr(1);
            std::string seq;
            while (std::getline(ifs, buffer) && buffer.find_first_of('>') != 0) {
                seq += buffer;
            }
            germline_collection_.emplace_back(gene_name, seq);
        }
    }
}
