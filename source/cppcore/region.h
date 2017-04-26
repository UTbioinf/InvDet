#ifndef __INVDET_REGION_H
#define __INVDET_REGION_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <RMQnlogn.hpp>

namespace loon
{

class VertexPair
{
public:
    size_t u, v, w;
public:
    VertexPair();
    VertexPair(size_t u, size_t v, size_t weight = 1);
    bool operator<(const VertexPair& rhs) const;
    void set_uvw(size_t u, size_t v, size_t w = 1);
    void inc_weight() const;
};

class OneAln
{
public:
    //std::string q_name;
    size_t q_id;
    long long q_length;
    long long r_start, r_end;
    long long q_start, q_end;
    short mapping_quality;
    char direction, q_pos;
public:
    OneAln(): q_length(-1) {}
    OneAln(const std::string& qname, long long qlen,
            long long rstart, long long rend,
            long long qstart, long long qend,
            short mapQ, char dir);
    bool is_forward() const;
    bool is_backward() const;
    bool is_5end() const;
    void invalidate();
    bool is_valid() const;
    bool operator<(const OneAln& rhs) const;
};

class SegEndpoint
{
public:
    size_t seg_id;
    size_t loc;
    bool is_start;
public:
    SegEndpoint(size_t s_id = -1, size_t l=-1, bool is_left=false):
        seg_id(s_id), loc(l), is_start(is_left)
    {}
    bool operator<(const SegEndpoint& se) const
    {
        if(loc < se.loc)    return true;
        if(loc > se.loc)    return false;
        if(!is_start && se.is_start)    return true;
        return false;
    }
};

class Region
{
private:
    long long r_length;
    std::string r_name;
    std::vector< OneAln > regional_alns;
    std::vector< size_t > segs_start;
    std::vector< size_t > segs_end;
    std::vector< size_t > seg_ids;
    std::map<size_t, std::vector<size_t> > pair_5, pair_3;
public:
    void add_ref_info(long long len, const std::string& name);
    void push_back(const std::string& q_name, long long q_length, 
            long long r_start, long long r_end,
            long long q_start, long long q_end, 
            short mapping_quality, char direction);
    void clear();
    void clear_name();
    void remove_low_coverage_reads(int min_cvg = 0, double min_cvg_percent = 0.0);
    void gen_vertices();
    void make_pairs();
    void write_graph(size_t r_id, std::ostream& out);
};

}//namespace loon

#endif
