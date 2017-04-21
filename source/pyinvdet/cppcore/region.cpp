#include "region.h"

namespace loon
{

OneAln::OneAln(const std::string& qname, long long qlen,
        long long rstart, long long rend,
        long long qstart, long long qend,
        short mapQ, char dir):
        q_name(qname), q_length(qlen),
        r_start(rstart), r_end(rend),
        q_start(qstart), q_end(qend),
        mapping_quality(mapQ), direction(dir)
{
}

bool OneAln::is_forward() const
{
    return direction == 'F';
}

bool OneAln::is_backward() const
{
    return direction == 'R';
}







void Region::add_ref_info(long long len, const std::string& name)
{
    r_length = len;
    r_name = name;
}

void Region::push_back(const std::string& q_name, long long q_length, 
        long long r_start, long long r_end,
        long long q_start, long long q_end,
        short mapping_quality, char direction)
{
    regional_alns.push_back( OneAln(q_name, q_length,
                                r_start, r_end,
                                q_start, q_end,
                                mapping_quality, direction) );
}

void Region::clear()
{
    regional_alns.clear();
}

void Region::clear_name()
{
    r_name.clear();
}

void Region::gen_graph_edges(size_t r_id, std::ostream& out)
{
    // sort
    // 
}


}// namespace loon
