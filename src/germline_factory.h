//
// @author: jiahong
// @date  : 24/02/18 4:33 PM
//

#ifndef IMMULATOR_GERMLINE_FACTORY_H
#define IMMULATOR_GERMLINE_FACTORY_H

#include <string>

namespace immulator{

class GermlineFactory {
public:
    GermlineFactory(const std::string &filename) { }


private:
    const std::string filename_;
};


}
#endif //IMMULATOR_GERMLINE_FACTORY_H
