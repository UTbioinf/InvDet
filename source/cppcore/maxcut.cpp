#include "maxcut.h"

namespace loon
{

MaxCut::MaxCut(int threshold/* = 15 */):
    small_threshold(threshold)
{}

MaxCut::~MaxCut()
{
    clear();
}

void MaxCut::clear()
{
    for(std::vector<std::vector<size_t>* >::iterator it = graph.begin(); it != graph.end(); ++it)
        delete (*it);
    node_relabel.clear();
    edge_list.clear();
    graph.clear();
    solution.clear();
}

void MaxCut::set_small_threshold(int threshold)
{
    small_threshold = threshold;
}

void MaxCut::add_edge(int u, int v, int w)
{
    int uu = node_relabel.add_raw_id(u);
    int vv = node_relabel.add_raw_id(v);
    size_t edge_size = edge_list.size();

    if(graph.size() <= uu)  graph.push_back( new std::vector<size_t>(1, edge_size) );
    else    graph[ uu ]->push_back( edge_size );
    if(graph.size() <= vv)  graph.push_back( new std::vector<size_t>(1, edge_size));
    else    graph[ vv ]->push_back( edge_size );

    edge_list.push_back( OneEdge(uu, vv, w));
}

size_t MaxCut::number_of_nodes() const
{
    return node_relabel.size();
}

size_t MaxCut::number_of_edges() const
{
    return edge_list.size();
}

int MaxCut::get_edge_u(size_t i) const
{
    if(i >= edge_list.size())   return -1;
    return edge_list[ i ].u;
}

int MaxCut::get_edge_v(size_t i) const
{
    if(i >= edge_list.size())   return -1;
    return edge_list[ i ].v;
}

int MaxCut::get_edge_w(size_t i) const
{
    if(i >= edge_list.size())   return -1;
    return edge_list[ i ].w;
}

int MaxCut::get_node_rawid(size_t i) const
{
    return node_relabel.get_raw_id( i );
}

const std::vector<bool>& MaxCut::get_solution() const
{
    return solution;
}

double MaxCut::get_value() const
{
    return value;
}

bool MaxCut::solve()
{
    if( is_bipartite() )    return true;
    if(exact_algorithm())   return true;
    max_spanning_tree();
    return false;
}

bool MaxCut::is_bipartite()
{
    int n = node_relabel.size();
    solution.assign(n, false);
    std::vector<bool> node_visisted(n, false);
    std::vector<int> node_queue;
    for(int i = 0; i < n; ++i)
    {
        if(node_visisted[i])    continue;
        node_visisted[i] = true;
        node_queue.push_back( i );
        while(!node_queue.empty())
        {
            int cur = node_queue.back();
            node_queue.pop_back();

            for(std::vector<size_t>::iterator eit = graph[cur]->begin();
                    eit != graph[cur]->end(); ++eit)
            {
                int t = edge_list[*eit].the_other_node( cur );
                if(node_visisted[t])
                {
                    if(solution[cur] == solution[t])    return false;
                }
                else
                {
                    node_visisted[t] = true;
                    node_queue.push_back( t );
                    solution[t] = !solution[cur];
                }
            }
        }
    }
    return true;
}

bool MaxCut::exact_algorithm()
{
    int N = node_relabel.size();
    if(N > small_threshold) return false;
    solution.assign(N, false);
    if(N == 1)  return true;
    if(N == 2)
    {
        solution[1] = true;
        return true;
    }

    int n = 1 << (N - 1);
    value = -1e10;
    int max_sol = 0;
    for(int i = 0; i < n; ++i)
    {
        double tmp_value = 0;
        for(size_t j = 0; j < edge_list.size(); ++j)
            if ( (!( (1<<edge_list[j].u) & (i<<1) )) ^ (!( (1 << edge_list[j].v) & (i<<1))))
                tmp_value += edge_list[j].w;
        if(tmp_value > value)
        {
            value = tmp_value;
            max_sol = i;
        }
    }
    max_sol <<= 1;
    for(int i = 0; i < N; ++i)
        solution[i] = (1 << i) & max_sol;
    return true;
}

void MaxCut::max_spanning_tree()
{
    int N = node_relabel.size();

    sort(edge_list.begin(), edge_list.end());
    std::vector<int> union_set(N, -1);
    for(std::vector<std::vector<size_t>* >::iterator it = graph.begin(); it != graph.end(); ++it)
        (*it)->clear();
    for(size_t i = 0; i < edge_list.size(); ++i)
    {
        int root_u = find_root(union_set, edge_list[i].u);
        int root_v = find_root(union_set, edge_list[i].v);
        if(root_u != root_v)
        {
            union_set[root_u] = root_v;
            graph[ edge_list[i].u ]->push_back( i );
            graph[ edge_list[i].v ]->push_back( i );
        }
    }
    is_bipartite();
    value = 0;
    for(size_t i = 0; i < edge_list.size(); ++i)
    {
        if(solution[ edge_list[i].u ] != solution[ edge_list[i].v ])
            value += edge_list[i].w;
    }
}

int MaxCut::find_root(std::vector<int>& union_set, int i)
{
    if(union_set[i] == -1)  return i;
    return (union_set[ i ] = find_root( union_set, union_set[i] ));
}

}// namespace loon
