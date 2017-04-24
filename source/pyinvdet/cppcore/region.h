#ifndef __INVDET_REGION_H
#define __INVDET_REGION_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>

namespace loon
{

class OneAln
{
public:
    std::string q_name;
    long long q_length;
    long long r_start, r_end;
    long long q_start, q_end;
    short mapping_quality;
    char direction;
public:
    OneAln(){}
    OneAln(const std::string& qname, long long qlen,
            long long rstart, long long rend,
            long long qstart, long long qend,
            short mapQ, char dir);
    bool is_forward() const;
    bool is_backward() const;
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
public:
    void add_ref_info(long long len, const std::string& name);
    void push_back(const std::string& q_name, long long q_length, 
            long long r_start, long long r_end,
            long long q_start, long long q_end, 
            short mapping_quality, char direction);
    void clear();
    void clear_name();
    void gen_graph_edges(size_t r_id, std::ostream& out);
};

}//namespace loon

#endif
