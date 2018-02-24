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
    auto query = gcfg_.next_roll();
    if (!query.empty()) {
        std::vector<Germline> filtered_germlines;
        std::copy_if(germline_collection_.begin(),
                     germline_collection_.end(),
                     std::back_inserter(filtered_germlines),
                     [&query](const Germline &key) {
                         if (query.find('*') != std::string::npos) {
                             // user provided full (family-gene*allele)
                             return key.name() == query;
                         } else if (query.find('-') != std::string::npos) {
                             // user provided gene (family-gene)
                             return key.gene_name() == query;
                         } else {
                             // user provided family
                             return key.family_name() == query;
                         }
                     });
        if (filtered_germlines.empty()) {
            return random_germline();
        } else {
            std::uniform_int_distribution<
                    std::vector<Germline>::size_type> dist(0, filtered_germlines.size() - 1);
            std::vector<Germline>::size_type index = dist(mersenne_);
            return filtered_germlines[index];
        }
    } else {
        return random_germline();
    }
}

void
immulator::GermlineFactory::parse_file() {
    std::ifstream ifs(filename_);
    if (ifs) {
        std::string buffer;
        std::getline(ifs, buffer);
        while (!buffer.empty()) {
            // keep "getting" until we reach the first entry
            while (buffer.find_first_of('>') != 0 && std::getline(ifs, buffer));

            auto tokens = immulator::split_string(buffer, "|");
            auto gene_name = tokens[1];
            auto asc_name = tokens[0].substr(1);

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

immulator::Germline
immulator::GermlineFactory::random_germline() const {
    std::uniform_int_distribution<
            std::vector<Germline>::size_type> dist(0, germline_collection_.size() - 1);
    std::vector<Germline>::size_type index = dist(mersenne_);
    return germline_collection_[index];
}
