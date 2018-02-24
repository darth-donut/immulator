//
// @author: jiahong
// @date  : 24/02/18 4:33 PM
//

#ifndef IMMULATOR_GERMLINE_FACTORY_H
#define IMMULATOR_GERMLINE_FACTORY_H

#include <string>
#include <vector>
#include "germline.h"

namespace immulator{

class GermlineFactory {
public:
    GermlineFactory(const std::string &filename) :
            filename_(filename) {
        parse_file();
    }

    /// returns a random germline
    /// \return Germline object
    immulator::Germline operator()() const;
private:
    void parse_file();
private:
    const std::string filename_;
    std::vector<immulator::Germline> germline_collection_;
};


}
#endif //IMMULATOR_GERMLINE_FACTORY_H
