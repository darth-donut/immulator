//
// @author: jiahong
// @date  : 24/02/18 4:33 PM
//

#include <fstream>
#include <iostream>
#include "germline_factory.h"

using immulator::Germline;


Germline
immulator::GermlineFactory::operator()() const {
    return germline_collection_[2];
}

void
immulator::GermlineFactory::parse_file() {
    std::ifstream ifs(filename_);
    if (ifs) {
        std::string buffer;
        std::getline(ifs, buffer);
        while (!buffer.empty()) {
            // keep "getting" until we reach the first entry
            while (buffer.find_first_of('>') != 0 && std::getline(ifs, buffer))
                ;

            auto cache_pipe_pos = buffer.find_first_of('|');
            auto gene_name = buffer.substr(cache_pipe_pos + 1,
                                                  buffer.find_first_of('|', cache_pipe_pos + 1) - cache_pipe_pos - 1);
            auto asc_name = buffer.substr(1, cache_pipe_pos - 1);

            // keep "getting" until we reach the next entry
            std::getline(ifs, buffer);
            std::string sequence;
            while (std::getline(ifs, buffer) && buffer.find_first_of('>') != 0) {
                sequence += buffer;
            }
            germline_collection_.emplace_back(std::move(gene_name), std::move(asc_name), std::move(sequence));
        }
    }
}
