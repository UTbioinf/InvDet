#include <fstream>
#include "invdet_core.h"

namespace loon
{

void InvDector::read(const std::string& fname)
{
    ifstream fin(fname.c_str());
    size_t n;
    fin >> n;
    regions.assign(n, Region());

    long long len;
    std::string name;
    for(size_t i = 0; i < n; ++i)
    {
        fin >> len;
        getline(fin, name);
        regions[i].add_ref_info(len, name);
    }

    size_t ref_id;
    long long r_start, r_end, q_start, q_end;
    short score;
    char dir;
    while(fin >> ref_id)
    {
        fin >> name >> len >> r_start >> r_end >> q_start >> q_end >> score >> dir;
        regions[ ref_id ].push_back( name, len, r_start, r_end,
                q_start, q_end, score, dir);
    }

    fin.close();
}

void InvDector::gen_graphs(const std::string& fname,
        int min_cvg/* = 0*/, double min_cvg_percent/* = 0.0 */)
{
    ofstream fout(fname.c_str());
    for(size_t i = 0; i < regions.size(); ++i)
    {
        regions[i].remove_low_coverage_reads(min_cvg, min_cvg_percent);
        regions[i].gen_vertices();
        regions[i].make_pairs();
        regions[i].write_graph( i, fout );
    }
    fout.close();
}

}// namespace loon
