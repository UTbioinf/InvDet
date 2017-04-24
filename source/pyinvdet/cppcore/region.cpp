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
    // convert to segment endpoints and sort
    std::vector< SegEndpoint > segEndpoints;
    size_t nn = regional_alns.size();
    segEndpoints.reserve(nn << 1);
    for(size_t i = 0; i < nn; ++i)
    {
        segEndpoints.push_back( SegEndpoint(i, regional_alns[i].r_start, true) );
        segEndpoints.push_back( SegEndpoint(i, regional_alns[i].r_end, false) );
    }
    std::sort(segEndpoints.begin(), segEndpoints.end());
    // convert reference to segments and small segments w/i each segment
    std::vector<size_t > ref_segments;
    ref_segments.reserve( nn );
    std::vector<size_t> q_segStart( nn );
    std::vector<size_t> q_segEnd( nn );
    size_t last_pos = 0;
    ref_segments.push_back(0);
    for(size_t i = 0; i < (nn << 1); ++i)
    {
        if(segEndpoints[i].is_start)
        {
            if(segEndpoints[i].loc == last_pos)
            {
                q_segStart[ segEndpoints[i].seg_id ] = ref_segments.size() - 1;
                ++ref_segments.back();
            }
            else
            {
                q_segStart[ segEndpoints[i].seg_id ] = ref_segments.size();
                ref_segments.push_back( ref_segments.back() + 1 );
                last_pos = segEndpoints[i].loc;
            }
        }
        else
        {
            if(last_pos == segEndpoints[i].loc)
            {
                q_segEnd[ segEndpoints[i].seg_id ] = ref_segments.size() - 2;
                --ref_segments.back();
            }
            else
            {
                q_segEnd[ segEndpoints[i].seg_id ] = ref_segments.size() - 1;
                ref_segments.push_back( ref_segments.back() - 1 );
                last_pos = segEndpoints[i].loc;
            }
        }
    }
    // debug
    std::cout << "coverage of each segment" << std::endl;
    for(size_t i = 0; i < ref_segments.size(); ++i)
        std::cout << std::setw(5) << i;
    std::cout << std::endl;
    for(size_t i = 0; i < ref_segments.size(); ++i)
        std::cout << std::setw(5) << ref_segments[i];
    std::cout << std::endl << "range of each read" << std::endl;
    for(size_t i = 0; i < nn; ++i)
        std::cout << "[" << regional_alns[i].r_start << ", " << regional_alns[i].r_end << "): "
                << "[" << q_segStart[i] << ", " << q_segEnd[i] << "]" << std::endl;
    // Range maximum query and remove low coverage reads
}


}// namespace loon
