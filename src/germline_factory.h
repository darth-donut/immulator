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
    GermlineFactory(const std::string &filename, std::size_t seed=std::random_device{}()) :
            filename_(filename), mersenne_(seed) {
        parse_file();
    }
    GermlineFactory(const std::string &filename, const immulator::GermlineConfiguration &gcfg) :
            filename_(filename), gcfg_(gcfg) {
        parse_file();
    }

    /// returns a random germline
    /// \return Germline object
    immulator::Germline operator()() const;
private:
    void parse_file();
    immulator::Germline random_germline() const;
private:
    const std::string filename_;
    const immulator::GermlineConfiguration gcfg_;
    std::vector<immulator::Germline> germline_collection_;
    mutable std::mt19937 mersenne_;
};


}
#endif //IMMULATOR_GERMLINE_FACTORY_H
