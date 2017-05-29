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
    void gen_graphs(const std::string& fname, 
            int min_cvg = 0, double min_cvg_percent = 0.0,
            int min_overlap=0, const std::string& nucmer_prefix="");
    void report_inversions(const std::string& graph_fname, 
            const std::string& maxcut_fname, 
            const std::string& inversion_fname);
};

}// namespace loon

#endif
