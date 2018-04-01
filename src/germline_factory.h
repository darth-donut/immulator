//
// @author: jiahong
// @date  : 24/02/18 4:33 PM
//

#ifndef IMMULATOR_GERMLINE_FACTORY_H
#define IMMULATOR_GERMLINE_FACTORY_H

#include <string>
#include <vector>
#include "germline.h"
#include "germline_configuration.h"

namespace immulator{

class GermlineFactory {
public:
    GermlineFactory(const std::string &filename) :
            filename_(filename) {
        parse_file();
    }
    GermlineFactory(const std::string &filename, const immulator::GermlineConfiguration &gcfg) :
            filename_(filename), gcfg_(gcfg) {
        parse_file();
    }

    template<typename T>
    immulator::Germline operator()(T &) const;
private:
    void parse_file();
    template<typename T>
    immulator::Germline random_germline(T &) const;
private:
    const std::string filename_;
    const immulator::GermlineConfiguration gcfg_;
    std::vector<immulator::Germline> germline_collection_;
};


template<typename T>
Germline
immulator::GermlineFactory::operator()(T &rand) const {
    auto query = gcfg_.next_roll(rand);
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
            return random_germline(rand);
        } else {
            std::uniform_int_distribution<
                    std::vector<Germline>::size_type> dist(0, filtered_germlines.size() - 1);
            std::vector<Germline>::size_type index = dist(rand);
            return filtered_germlines[index];
        }
    } else {
        return random_germline(rand);
    }
}

template<typename T>
immulator::Germline
immulator::GermlineFactory::random_germline(T &rand) const {
    std::uniform_int_distribution<
            std::vector<Germline>::size_type> dist(0, germline_collection_.size() - 1);
    std::vector<Germline>::size_type index = dist(rand);
    return germline_collection_[index];
}

}
#endif //IMMULATOR_GERMLINE_FACTORY_H
