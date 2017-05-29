#include <relabel.h>
#include <stdexcept>
#include "region.h"

#define DEBUG_SEG_EXTRACTION
#define DEBUG_B_GRAPH_PRINT
#define DEBUG_EDGE_WEIGHT_FILTER


namespace loon
{

InvertedRepeats::InvertedRepeats(const std::string& filedir):
    prefix(filedir)
{}

void InvertedRepeats::parse_two_ints(const std::string& line)
{
    size_t i = 0;
    size_t num = 0;
    while(line[i] >= '0' && line[i] <= '9')
        num = num * 10 + line[i++] - '0';
    r_starts.push_back( num - 1 );

    num = 0;    ++i;
    while(line[i] >= '0' && line[i] <= '9')
        num = num * 10 + line[i++] - '0';
    r_ends.push_back( num );
}

void InvertedRepeats::quicksort(long long L, long long R)
{
    long long i = L, j = R;
    size_t s_mid = r_starts[ (i + j) >> 1], e_mid = r_ends[ (i + j) >> 1 ];
    while(i < j)
    {
        while(r_starts[i] < s_mid || (r_starts[i] == s_mid && r_ends[i] < e_mid))    ++i;
        while(r_starts[j] > s_mid || (r_starts[j] == s_mid && r_ends[j] > e_mid))    --j;
        if(i <= j)
        {
            std::swap(r_starts[i], r_starts[j]);
            std::swap(r_ends[i], r_ends[j]);
            ++i; --j;
        }
    }
    if(i < R)   quicksort(i, R);
    if(L < j)   quicksort(L, j);
}

void InvertedRepeats::open(size_t file_id)
{
    std::ostringstream oss;
#if defined(_WIN32) || defined(__CYGWIN__)
    oss << "\\nc_aln." << file_id << ".delta";
#else
    oss << "/nc_aln." << file_id << ".delta";
#endif
    std::string fname = prefix + oss.str();
    fin.open( fname.c_str() );
    if(!fin.is_open())
        throw std::runtime_error("region: cannot open file [" + fname + "]");
}

void InvertedRepeats::close()
{
    fin.close();
}

void InvertedRepeats::read()
{
    clear();
    std::string line;
    getline(fin, line);
    getline(fin, line);
    while(getline(fin, line))
    {
        if(line[0] == '>')
            continue;
        parse_two_ints( line );
        for(getline(fin, line); line[0] != '0'; getline(fin, line))
            ;
    }
    if(! r_starts.empty())
        quicksort(0, r_starts.size() - 1);
}

void InvertedRepeats::clear()
{
    r_starts.clear();
    r_ends.clear();
}

size_t InvertedRepeats::size() const
{
    return r_starts.size();
}

bool InvertedRepeats::contained(size_t i, size_t j) const
{
    // to be improved
    double uncovered_ratio = 1.0, tmp_ratio;
    size_t len = j - i;
    size_t left_ext, right_ext;
    for(size_t ii = 0; ii < r_starts.size(); ++ii)
    {
        left_ext = right_ext = 0;
        if(i < r_starts[ii])    left_ext = r_starts[ii] - i;
        if(j > r_ends[ii])      right_ext = j - r_ends[ii];
        tmp_ratio = double(left_ext + right_ext) / len;
        if(tmp_ratio < uncovered_ratio) uncovered_ratio = tmp_ratio;
    }
    return (uncovered_ratio < 0.001);
}

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

#if defined(DEBUG_SEG_EXTRACTION) || defined(DEBUG_EDGE_WEIGHT_FILTER) || defined(DEBUG_B_GRAPH_PRINT)
#include <iostream>
#endif
void Region::gen_vertices(int min_overlap)
{
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
#ifdef DEBUG_SEG_EXTRACTION
    std::cerr << "segments start" << std::endl;
    for(size_t i = 0; i < segs_start.size(); ++i)
        std::cerr << segs_start[i] << ' ' << segs_end[i] << std::endl;
    std::cerr << "segments end" << std::endl << std::endl;
#endif
}

void Region::make_pairs(InvertedRepeats* inv_repeats/* = NULL */)
{
    std::map<size_t, std::vector<size_t> >::iterator it;
    for(size_t i = 0; i < regional_alns.size(); ++i)
    {
        if(inv_repeats && inv_repeats->contained( regional_alns[i].r_start, regional_alns[i].r_end ))
            continue;
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
        {
            for(std::vector<size_t>::iterator vit = it3->second.begin();
                    vit != it3->second.end(); ++vit)
            {
                if(seg_ids[ *uit ] == seg_ids[ *vit ] || 
                        regional_alns[ *uit ].is_forward() == regional_alns[ *vit ].is_forward())
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
    }

    if(not edges.empty())
    {
    #ifdef DEBUG_EDGE_WEIGHT_FILTER
        size_t edge_count = 0;
        const int min_weight = 5;
        for(eit = edges.begin(); eit != edges.end(); ++eit)
            if((*eit).w > min_weight)
                ++edge_count;
        if(edge_count > 0)
        {
            out << edge_count << ' ' << r_id << std::endl;
            for(eit = edges.begin(); eit != edges.end(); ++eit)
                if((*eit).w > min_weight)
                    out << (*eit).u << ' ' << (*eit).v << ' ' << (*eit).w << std::endl;
        }
    #else
        out << edges.size() << ' ' << r_id << std::endl;
        for(eit = edges.begin(); eit != edges.end(); ++eit)
            out << (*eit).u << ' ' << (*eit).v << ' ' << (*eit).w << std::endl;
    #endif
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
    int u_id, v_id;
#ifdef DEBUG_B_GRAPH_PRINT
    std::cout << "Print graph" << std::endl;
#endif
    for(size_t i = 0; i < n_edges; ++i)
    {
        graph_in >> u_name >> v_name >> e_weight;
        u_id = node_relabel.get_new_id( u_name );
        v_id = node_relabel.get_new_id( v_name );
        if( max_cut[ u_id ] == max_cut[ v_id ])
            continue;
        graph[ u_id ].push_back( v_id );
        graph[ v_id ].push_back( u_id );
    #ifdef DEBUG_B_GRAPH_PRINT
        std::cout << "(" << u_name << ',' << v_name <<"; " << e_weight << ")" << std::endl;
    #endif
    }
#ifdef DEBUG_B_GRAPH_PRINT
    std::cout << "Print graph finished" << std::endl;
#endif

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
#ifdef DEBUG_B_GRAPH_PRINT
    std::cout << "solution: " << std::endl;
#endif
    for(int i = 0; i < n_vertices; ++i)
    {
        if(!max_cut[i])
        {
            int raw_id = node_relabel.get_raw_id( i );
            inv_out << segs_start[ raw_id ] << ' ' << segs_end[ raw_id ] << std::endl;
        }
    #ifdef DEBUG_B_GRAPH_PRINT
        if(max_cut[i])
        {
            int raw_id = node_relabel.get_raw_id( i );
            std::cout << "(1, " << raw_id << "): " << segs_start[ raw_id ] << ' ' << segs_end[ raw_id ] << "; weight=" << (segs_end[ raw_id] - segs_start[ raw_id ])<< std::endl;
        }
        else
        {
            int raw_id = node_relabel.get_raw_id( i );
            std::cout << "(0, " << raw_id << "): " << segs_start[ raw_id ] << ' ' << segs_end[ raw_id ] << "; weight=" << (segs_end[raw_id] - segs_start[raw_id]) << std::endl;
        }
    #endif
    }
#ifdef DEBUG_B_GRAPH_PRINT
    std::cout << "solution finished" << std::endl;
#endif
}

}// namespace loon
