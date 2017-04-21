#ifndef __INVDET_CORE_H
#define __INVDET_CORE_H

#include <vector>
#include <string>
#include "region.h"

namespace loon
{

class InvDector
{
private:
    std::vector<Region> regions;
public:
    void read(const std::string& fname);
    void gen_graph_edges(const std::string& fname, 
        int min_cvg = 0, double min_cvg_percent = 0.0);
};

}// namespace loon

#endif
