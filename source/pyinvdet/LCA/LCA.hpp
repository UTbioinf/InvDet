/*
LCA: Michael A. Bender & Martin Farach-Colton's algorithm
<O(n), O(1)>

author: zijuexiansheng
*/

/*
In order to this this LCA class, the TreeNodeType should satisfy the following conditions
    * There is an extra member in the TreeNodeType class so that it could be modified
        * a member function `set_nodeId( id )` should be provided to modify its value
        * a member function `get_nodeId()` should be provided to read its value
        * The type should be `size_t`
    * `get_child()`, `get_child(i)` and `get_sibling()` should be provided
    * If the tree is of the form left-child-right-sibling:
        * `get_child(i)` must be provided, but can do anything
        * `get_child()` is used to return the left child
        * `get_sibling()` is used to return the sibling
    * If the tree is just of a generic representation
        * `get_child(i)` is used to return the i-th child
            * if it returns NULL, it means that the i-th, (i+1)-th, ... children are all NULL
        * `get_child()` and `get_sibling` must be provided but can do anything
    * The prototype of these member variable/functions are as follows:
        * `size_t nodeId;`: the variable name `nodeId` is not an requirement
        * `void set_nodeId(size_t id);`: the `void` return type is not an requirement
        * `size_t get_nodeId() const`
        * `TreeNodeType* get_child() const`
        * `TreeNodeType* get_child(size_t i) const`: This one can be combined with the one above
        * `TreeNodeType* get_sibling() const`
*/

#ifndef __LCA_H
#define __LCA_H

#include <vector>

namespace loon
{

template<class TreeNodeType>
class LCA
{
private:
    class ArrayElement
    {
    public:
        size_t key;
        TreeNodeType* value;
    };
    class RMQ_Restricted
    {
    public:
        std::vector<size_t> log_table;
        std::vector<ArrayElement> A;
        std::vector<std::vector<std::vector<size_t> > > blocks;
        std::vector<size_t> block_ids;
        size_t block_size;

        std::vector<std::vector<ArrayElement*> > RMQ_log_minima;
    public:
        RMQ_Restricted();
        void reserve(size_t n);
        size_t size() const;
        void push_back(const ArrayElement& e);
        void clear();
        void preprocess();
        ArrayElement* query(size_t i, size_t j);
        ArrayElement* short_query(size_t i, size_t j);
        size_t get_block_id(size_t i, size_t j);

        // return minimum value
        void RMQ_log_Preprocess(const std::vector<ArrayElement*>& X);
        ArrayElement* RMQ_log_query(size_t i, size_t j);
        // return minmum index within the block
        void RMQ_block_Preprocess(std::vector<std::vector<size_t> >& m, size_t i, size_t j);
        size_t RMQ_block_query(const std::vector<std::vector<size_t> >& m, size_t i, size_t j);

        size_t log2(size_t n);
    };

    size_t nodeId;
    ArrayElement one_ele;
    std::vector<size_t> treeArray;
    RMQ_Restricted rmq;
    bool is_leftchild_rightsibling;
    
    // is a fat tree
    void generate_tree_array_fattree(TreeNodeType* node, size_t level);
    // is left child right sibling
    void generate_tree_array_lcrs(TreeNodeType* node, size_t level);
public:
    void preprocess(TreeNodeType* root, bool is_leftchild_rightsibling = false);// this will clear data automatically
    TreeNodeType* query(TreeNodeType* node1, TreeNodeType* node2);

};

/*LCA<TreeNodeType>::RMQ_Restricted*/
template<class TreeNodeType>
LCA<TreeNodeType>::RMQ_Restricted::RMQ_Restricted(): log_table(2, 0)
{}

template<class TreeNodeType>
void LCA<TreeNodeType>::RMQ_Restricted::reserve(size_t n)
{
    A.reserve(n);
}

template<class TreeNodeType>
size_t LCA<TreeNodeType>::RMQ_Restricted::size() const
{
    return A.size();
}

template<class TreeNodeType>
void LCA<TreeNodeType>::RMQ_Restricted::push_back(const LCA<TreeNodeType>::ArrayElement& e)
{
    A.push_back( e );
}

template<class TreeNodeType>
void LCA<TreeNodeType>::RMQ_Restricted::clear()
{
    A.clear();
    blocks.clear();
    block_ids.clear();
}

template<class TreeNodeType>
void LCA<TreeNodeType>::RMQ_Restricted::preprocess()
{
    size_t n = A.size();
    if(n < 4)   return;
    block_size = log2(A.size()) >> 1;
    blocks.resize( (1 << block_size) - 1);

    std::vector<LCA<TreeNodeType>::ArrayElement*> blockmin;
    for(size_t i = 0; i < n; i += block_size)
    {
        size_t tmp_end = std::min(i + block_size - 1, n - 1);
        size_t blockid = get_block_id(i, tmp_end);
        if(blocks[ blockid ].empty())
            RMQ_block_Preprocess( blocks[blockid], i, tmp_end );
        size_t block_min_pos = RMQ_block_query( blocks[blockid], 0, tmp_end - i ) + i;
        block_ids.push_back(blockid);
        blockmin.push_back( &A[ block_min_pos ] );
    }
    RMQ_log_Preprocess(blockmin);
}

template<class TreeNodeType>
typename LCA<TreeNodeType>::ArrayElement* LCA<TreeNodeType>::RMQ_Restricted::query(size_t i, size_t j)
{
    if(i == j)  return &A[i];
    else if(i > j)  std::swap(i, j);
    if(i + 3 > j)   return short_query(i, j);
    size_t firstblock = i / block_size;
    size_t lastblock  = j / block_size;
    if(firstblock == lastblock)
    {
        size_t block_start = firstblock * block_size;
        return &A[ RMQ_block_query( blocks[block_ids[ firstblock ]], i-block_start, j-block_start ) + block_start ];
    }
    else
    {
        size_t first_start = firstblock * block_size;
        size_t last_start = lastblock * block_size;

        size_t first_block_min = RMQ_block_query( blocks[block_ids[ firstblock ]], i-first_start, block_size-1 ) + first_start;
        size_t last_block_min = RMQ_block_query( blocks[block_ids[ lastblock ]], 0, j-last_start) + last_start;

        LCA<TreeNodeType>::ArrayElement* ret = A[first_block_min].key < A[last_block_min].key ? &A[ first_block_min ] : &A[ last_block_min ];
        if(firstblock + 1 < lastblock)
        {
            LCA<TreeNodeType>::ArrayElement* min_arr = RMQ_log_query( firstblock + 1, lastblock-1 );
            if(min_arr->key < ret->key)
                ret = min_arr;
        }
        return ret;
    }
}

template<class TreeNodeType>
typename LCA<TreeNodeType>::ArrayElement* LCA<TreeNodeType>::RMQ_Restricted::short_query(size_t i, size_t j)
{
    LCA<TreeNodeType>::ArrayElement* ret = &A[i];
    for(++i; i<=j; ++i)
        if(ret->key > A[i].key)
            ret = &A[i];
    return ret;
}

template<class TreeNodeType>
size_t LCA<TreeNodeType>::RMQ_Restricted::get_block_id(size_t i, size_t j)
{
    size_t blockid = 0;
    for(size_t k=i; k<j; ++k)
        blockid |= (A[k+1].key > A[k].key) << (k-i);
    return blockid;
}

template<class TreeNodeType>
void LCA<TreeNodeType>::RMQ_Restricted::RMQ_log_Preprocess(const std::vector<LCA<TreeNodeType>::ArrayElement*>& X)
{
    size_t n = X.size();
    size_t log2_n = log2(n);
    RMQ_log_minima.assign( log2_n + 1, X );
    for(size_t i = 1; i <= log2_n; ++i)
        for(size_t j = 0; j < n; ++j)
        {
            RMQ_log_minima[i][j] = RMQ_log_minima[i-1][std::min(n-1, j+(1<<(i-1)))];
            if(RMQ_log_minima[i-1][j]->key < RMQ_log_minima[i][j]->key)
                RMQ_log_minima[i][j] = RMQ_log_minima[i-1][j];
        }
}

template<class TreeNodeType>
typename LCA<TreeNodeType>::ArrayElement* LCA<TreeNodeType>::RMQ_Restricted::RMQ_log_query(size_t i, size_t j)
{
    if(i == j)  return RMQ_log_minima[0][i];
    size_t k = log2(j-i);

    LCA<TreeNodeType>::ArrayElement* ret = RMQ_log_minima[k][j+1-(1<<k)];
    if(ret->key > RMQ_log_minima[k][i]->key)
        ret = RMQ_log_minima[k][i];
    return ret;
}

template<class TreeNodeType>
void LCA<TreeNodeType>::RMQ_Restricted::RMQ_block_Preprocess(std::vector<std::vector<size_t> >& m, size_t i, size_t j)
{
    size_t n = j + 1 - i;
    m.assign(n, std::vector<size_t>(n));
    for(size_t ii = 0; ii < n; ++ii)
    {
        m[ii][ii] = ii;
        for(size_t jj = ii + 1; jj < n; ++jj)
        {
            m[ii][jj] = jj;
            if(A[i + m[ii][jj]].key > A[i + m[ii][jj-1]].key)
                m[ii][jj] = m[ii][jj-1];
        }
    }
}

template<class TreeNodeType>
size_t LCA<TreeNodeType>::RMQ_Restricted::RMQ_block_query(const std::vector<std::vector<size_t> >& m, size_t i, size_t j)
{
    if(i > j)   std::swap(i, j);
    return m[i][j];
}

template<class TreeNodeType>
size_t LCA<TreeNodeType>::RMQ_Restricted::log2(size_t n)
{
    while(log_table.size() <= n)
        log_table.insert( log_table.end(), log_table.size(), 1 + log_table.back() );
    return log_table[n];
}

/*LCA<TreeNodeType>*/
template<class TreeNodeType>
void LCA<TreeNodeType>::generate_tree_array_fattree(TreeNodeType* node, size_t level)
{
    node->set_nodeId( ++nodeId );
    treeArray.push_back( rmq.size() );

    one_ele.key = level;
    one_ele.value = node;
    rmq.push_back( one_ele );

    TreeNodeType* child;
    for(size_t i = 0; (child = node->get_child(i))!=NULL; ++i)
    {
        generate_tree_array_fattree(child, level+1);
        one_ele.key = level;
        one_ele.value = node;
        rmq.push_back( one_ele );
    }
}

template<class TreeNodeType>
void LCA<TreeNodeType>::generate_tree_array_lcrs(TreeNodeType* node, size_t level)
{
    node->set_nodeId( ++nodeId );
    treeArray.push_back( rmq.size() );

    one_ele.key = level;
    one_ele.value = node;
    rmq.push_back( one_ele );
    TreeNodeType* child;
    for(child = node->get_child(); child != NULL; child=child->get_sibling())
    {
        generate_tree_array_lcrs(child, level + 1);
        one_ele.key = level;
        one_ele.value = node;
        rmq.push_back( one_ele );
    }
}

template<class TreeNodeType>
void LCA<TreeNodeType>::preprocess(TreeNodeType* root, bool is_leftchild_rightsibling/*=false*/)
{
    this->is_leftchild_rightsibling = is_leftchild_rightsibling;
    nodeId = -1;
    treeArray.clear();
    rmq.clear();
    if( is_leftchild_rightsibling )
        generate_tree_array_lcrs(root, 0);
    else
        generate_tree_array_fattree(root, 0);
    rmq.preprocess();
}

template<class TreeNodeType>
TreeNodeType* LCA<TreeNodeType>::query(TreeNodeType* node1, TreeNodeType* node2)
{
    LCA<TreeNodeType>::ArrayElement* ret = rmq.query( treeArray[ node1->get_nodeId() ], treeArray[ node2->get_nodeId() ]);
    return ret->value;
}

}// namespace loon



#endif
