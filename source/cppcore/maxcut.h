#ifndef __CORE_MAXCUT_H
#define __CORE_MAXCUT_H

/*
if( the graph is bipartite )
    return the (-1, 1) labeling of the graph
else
    if( the graph is small )
        run the exact algorithm
    else
        run maximum spanning tree
*/

#include <vector>
#include <algorithm>
#include <relabel.h>

namespace loon
{
class MaxCut
{
private:
    class OneEdge{
    public:
        int u, v, w;
    public:
        OneEdge(){}
        OneEdge(int uu, int vv, int ww):
            u(uu), v(vv), w(ww)
        {}
        int the_other_node(int node) const
        {
            return (node == u ? v : u);
        }
        bool operator<(const OneEdge& rhs) const
        {
            return w > rhs.w;
        }
    };
    RelabelSmallPosInt<int, int> node_relabel;
    std::vector<OneEdge> edge_list;
    std::vector<std::vector<size_t>* > graph;
    std::vector<bool> solution;
    int small_threshold;
    double value;

    int find_root(std::vector<int>& union_set, int i);
public:
    MaxCut(int threshold = 15); // threshold should be <= 30
    ~MaxCut();
    void clear();
    void set_small_threshold(int threshold);// threshold should be <= 30
    void add_edge(int u, int v, int w);
    size_t number_of_nodes() const;
    size_t number_of_edges() const;
    int get_edge_u(size_t i) const;
    int get_edge_v(size_t i) const;
    int get_edge_w(size_t i) const;
    int get_node_rawid(size_t i) const;
    const std::vector<bool>& get_solution() const;
    double get_value() const;
    bool solve(); // return: true if optimal, false if not
    bool is_bipartite();
    bool exact_algorithm();
    void max_spanning_tree();
};

}// namespace loon

#endif
