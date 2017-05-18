#include <relabel.h>
#include "region.h"

namespace loon
{

VertexPair::VertexPair():
    u(0), v(0), w(0)
{}

VertexPair::VertexPair(size_t u, size_t v, size_t weight/*=1*/):
    w(weight)
{
    if(u < v)
    {
        this->u = u;
        this->v = v;
    }
    else
    {
        this->u = v;
        this->v = u;
    }
}

void VertexPair::set_uvw(size_t u, size_t v, size_t w/* = 1*/)
{
    this->w = w;
    if(u < v)
    {
        this->u = u;
        this->v = v;
    }
    else
    {
        this->u = v;
        this->v = u;
    }
}

bool VertexPair::operator<(const VertexPair& rhs) const
{
    return (u < rhs.u || (u == rhs.u && v < rhs.v));
}

void VertexPair::inc_weight() const
{
    ++const_cast<VertexPair*>(this)->w;
}

OneAln::OneAln(const std::string& qname, long long qlen,
        long long rstart, long long rend,
        long long qstart, long long qend,
        short mapQ, char dir):
        q_length(qlen),
        r_start(rstart), r_end(rend),
        q_start(qstart), q_end(qend),
        mapping_quality(mapQ), direction(dir)
{
    std::istringstream iss( qname.substr(qname.find("afun") + std::string("afun").length()) );
    iss >> q_id >> q_pos >> q_pos;
}

bool OneAln::is_forward() const
{
    return direction == 'F';
}

bool OneAln::is_backward() const
{
    return direction == 'R';
}

bool OneAln::is_5end() const
{
    return q_pos == '5';
}

void OneAln::invalidate()
{
    q_length = -1;
}

bool OneAln::is_valid() const
{
    return q_length > 0;
}

bool OneAln::operator<(const OneAln& rhs) const
{
    if(is_valid() && (r_start < rhs.r_start || 
            (r_start == rhs.r_start && r_end > rhs.r_end)))   return true;
    return false;
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
    segs_start.clear();
    segs_end.clear();
    seg_ids.clear();
    pair_5.clear();
    pair_3.clear();
}

void Region::clear_name()
{
    r_name.clear();
}

void Region::remove_low_coverage_reads(int min_cvg/* = 0*/, double min_cvg_percent/* = 0.0*/)
{
    // convert to segment endpoints and sort
    std::vector< SegEndpoint > segEndpoints;
    size_t nn = regional_alns.size();
    segEndpoints.reserve(nn << 1);
    
    double total_len = 0;
    for(size_t i = 0; i < nn; ++i)
    {
        segEndpoints.push_back( SegEndpoint(i, regional_alns[i].r_start, true) );
        segEndpoints.push_back( SegEndpoint(i, regional_alns[i].r_end, false) );
        total_len += regional_alns[i].r_end - regional_alns[i].r_start;
    }
    min_cvg = std::max(min_cvg, int(total_len * min_cvg_percent / r_length));
    std::sort(segEndpoints.begin(), segEndpoints.end());
    // convert reference to segments and small segments w/i each segment
    //std::vector<size_t > ref_segments;
    //ref_segments.reserve( nn );
    std::vector<size_t> q_segStart( nn );
    std::vector<size_t> q_segEnd( nn );
    size_t last_pos = 0;
    //ref_segments.push_back(0);
    RMQnlogn<size_t> rmq;
    size_t cur_cvg = 0, cur_index = 1;
    for(size_t i = 0; i < (nn << 1); ++i)
    {
        if(segEndpoints[i].is_start)
        {
            if(segEndpoints[i].loc == last_pos)
            {
                q_segStart[ segEndpoints[i].seg_id ] = cur_index - 1;//ref_segments.size() - 1;
                //++ref_segments.back();
            }
            else
            {
                q_segStart[ segEndpoints[i].seg_id ] = cur_index;//ref_segments.size();
                rmq.push_back( nn - cur_cvg );
                ++cur_index;
                //ref_segments.push_back( ref_segments.back() + 1 );
                last_pos = segEndpoints[i].loc;
            }
            ++cur_cvg;
        }
        else
        {
            if(last_pos == segEndpoints[i].loc)
            {
                q_segEnd[ segEndpoints[i].seg_id ] = cur_index - 2;//ref_segments.size() - 2;
                //--ref_segments.back();
            }
            else
            {
                q_segEnd[ segEndpoints[i].seg_id ] = cur_index - 1;//ref_segments.size() - 1;
                rmq.push_back( nn - cur_cvg );
                ++cur_index;
                //ref_segments.push_back( ref_segments.back() - 1 );
                last_pos = segEndpoints[i].loc;
            }
            --cur_cvg;
        }
    }
    rmq.push_back( nn-cur_cvg );
    // Range maximum query and remove low coverage reads
    rmq.preprocess();
    for(size_t i = 0; i < nn; ++i)
    {
        size_t max_cvg = nn - rmq.query(q_segStart[i], q_segEnd[i]);
        if(max_cvg < min_cvg)   regional_alns[i].invalidate();
    }
}

#include <iostream>
void Region::gen_vertices()
{
    const int min_overlap = 20;
    std::sort( regional_alns.begin(), regional_alns.end() );
    for(size_t i = 0; i < regional_alns.size(); ++i)
    {
        if(regional_alns[i].is_valid())
        {
            if(segs_start.empty() || segs_end.back() <= regional_alns[i].r_start + min_overlap)
            {// new segment
                segs_start.push_back( regional_alns[i].r_start );
                segs_end.push_back( regional_alns[i].r_end );
            }
            else if(segs_end.back() < regional_alns[i].r_end)
            {// update right-end
                segs_end.back() = regional_alns[i].r_end;
            }
            seg_ids.push_back( segs_start.size() - 1 );
        }
        else
        {
            regional_alns.resize( i );
            break;
        }
    }
    std::cerr << "segments start" << std::endl;
    for(size_t i = 0; i < segs_start.size(); ++i)
        std::cerr << segs_start[i] << ' ' << segs_end[i] << std::endl;
    std::cerr << "segments end" << std::endl << std::endl;
}

void Region::make_pairs()
{
    std::map<size_t, std::vector<size_t> >::iterator it;
    for(size_t i = 0; i < regional_alns.size(); ++i)
    {
        if(regional_alns[i].is_5end())
        {
            it = pair_5.find( regional_alns[i].q_id );
            if(it == pair_5.end())
                pair_5[ regional_alns[i].q_id ] = std::vector<size_t>(1, i);
            else
                it->second.push_back( i );
        }
        else
        {
            it = pair_3.find( regional_alns[i].q_id );
            if(it == pair_3.end())
                pair_3[ regional_alns[i].q_id ] = std::vector<size_t>(1, i);
            else
                it->second.push_back( i );
        }
    }
}

void Region::write_graph(size_t r_id, std::ostream& out)
{
    std::set<VertexPair> edges;
    std::set<VertexPair>::iterator eit;
    VertexPair e;
    std::map<size_t, std::vector<size_t> >::iterator it5, it3;
    for(it5 = pair_5.begin(); it5 != pair_5.end(); ++it5)
    {
        it3 = pair_3.find( it5->first );
        if(it3 == pair_3.end()) continue;
        for(std::vector<size_t>::iterator uit = it5->second.begin();
                uit != it5->second.end(); ++uit)
            for(std::vector<size_t>::iterator vit = it3->second.begin();
                    vit != it3->second.end(); ++vit)
            {
                if(seg_ids[ *uit ] == seg_ids[ *vit ] || regional_alns[ *uit ].is_forward() == regional_alns[ *vit ].is_forward())
                {
                    // if the validated region still have inversions inside, ignore it
                    // or the two reads are in the same direction, ignore it
                    continue;
                }
                e.set_uvw( seg_ids[ *uit ], seg_ids[ *vit ] );
                eit = edges.find( e );
                if(eit == edges.end())
                    edges.insert( e );
                else
                    (*eit).inc_weight();
            }
    }

    if(not edges.empty())
    {
        out << edges.size() << ' ' << r_id << std::endl;
        for(eit = edges.begin(); eit != edges.end(); ++eit)
            out << (*eit).u << ' ' << (*eit).v << ' ' << (*eit).w << std::endl;
    }
}

void Region::report_inversions(size_t r_id, size_t n_vertices, size_t n_edges,
        std::istream& graph_in, std::istream& maxcut_in, std::ostream& inv_out)
{
    RelabelSmallPosInt<int, int> node_relabel;
    std::vector<bool> max_cut(n_vertices, false);

    int u_name, v_name, e_weight;
    bool u_party;
    for(size_t i = 0; i < n_vertices; ++i)
    {
        maxcut_in >> u_name >> u_party;
        max_cut[ node_relabel.add_raw_id( u_name ) ] = u_party;
    }

    std::vector< std::vector<int> > graph( n_vertices );
    for(size_t i = 0; i < n_edges; ++i)
    {
        graph_in >> u_name >> v_name >> e_weight;
        graph[ node_relabel.get_new_id( u_name ) ].push_back( node_relabel.get_new_id( v_name ) );
        graph[ node_relabel.get_new_id( v_name ) ].push_back( node_relabel.get_new_id( u_name ) );
    }

    std::vector<bool> visited( n_vertices, false );
    std::vector<int> flip_id;
    int flip_left;
    size_t weight_false = 0, weight_true = 0;
    size_t tmp_weight;
    for(int i = 0; i < n_vertices; ++i)
    {
        if(! visited[i])
        {
            weight_true = weight_false = 0;
            flip_id.clear();
            flip_left = 0;

            visited[ i ] = true;
            flip_id.push_back( i );
            int raw_id = node_relabel.get_raw_id(i);
            tmp_weight = segs_end[ raw_id ] - segs_start[ raw_id ];
            if( max_cut[i] )    weight_true += tmp_weight;
            else    weight_false += tmp_weight;

            while(flip_left < flip_id.size())
            {
                int uid = flip_id[ flip_left++ ];
                for(std::vector<int>::iterator vit = graph[ uid ].begin();
                        vit != graph[ uid ].end(); ++vit)
                    if(! visited[ *vit ])
                    {
                        visited[ *vit ] = true;
                        flip_id.push_back( *vit );
                        raw_id = node_relabel.get_raw_id( *vit );
                        tmp_weight = segs_end[ raw_id ] - segs_start[ raw_id ];
                        if( max_cut[*vit] ) weight_true += tmp_weight;
                        else    weight_false += tmp_weight;
                    }
            }
            if(weight_true < weight_false)
            {
                for(std::vector<int>::iterator it = flip_id.begin(); it != flip_id.end(); ++it)
                    max_cut[ *it ] = !max_cut[ *it ];
            }
        }
    }

    int false_cnt = std::count(max_cut.begin(), max_cut.end(), false);
    if(false_cnt == 0)  return;
    inv_out << false_cnt << ' ' << r_id << std::endl;
    for(int i = 0; i < n_vertices; ++i)
    {
        if(!max_cut[i])
        {
            int raw_id = node_relabel.get_raw_id( i );
            inv_out << segs_start[ raw_id ] << ' ' << segs_end[ raw_id ] << std::endl;
        }
    }
}

}// namespace loon
