#include <iostream>
#include <fstream>
#include <cstdlib>
#include "invdet_core.h"

namespace loon
{

void InvDector::read(const std::string& fname)
{
    std::ifstream fin(fname.c_str());
    size_t n = 0;
    if(!(fin >> n)) return;
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
    std::ofstream fout(fname.c_str());
    for(size_t i = 0; i < regions.size(); ++i)
    {
        //regions[i].remove_low_coverage_reads(min_cvg, min_cvg_percent);
        regions[i].gen_vertices();
        regions[i].make_pairs();
        regions[i].write_graph( i, fout );
    }
    fout.close();
}

void InvDector::report_inversions(const std::string& graph_fname,
        const std::string& maxcut_fname,
        const std::string& inversion_fname)
{
    std::ifstream graph_fin( graph_fname.c_str() );
    std::ifstream maxcut_fin( maxcut_fname.c_str() );
    std::ofstream inv_fout( inversion_fname.c_str() );

    size_t n_vertices, n_edges;
    size_t ref_id;
    while(graph_fin >> n_edges >> ref_id)
    {
        size_t rid;
        maxcut_fin >> n_vertices >> rid;
        if(ref_id != rid)
        {
            std::cerr << "[ERROR]: The graph file [" << graph_fname << "] and maxcut file [" << maxcut_fname << "] don't match!" << std::endl;
            exit(-1);
        }
        regions[ref_id].report_inversions( ref_id, n_vertices, n_edges,
                graph_fin, maxcut_fin, inv_fout);
    }

    graph_fin.close();
    maxcut_fin.close();
    inv_fout.close();
}

}// namespace loon
