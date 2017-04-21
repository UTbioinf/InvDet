/*
RMQ: Michael A. Bender & Martin Farach-Colton's algorithm
<O(n), O(1)>

author: zijuexiansheng
*/

#ifndef __RMQ_LINEAR_H
#define __RMQ_LINEAR_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include "LCA.hpp"

namespace loon
{

template<class ValType>
class RMQLinear
{
private:
    class TreeNode;
    std::vector<TreeNode> a;
    size_t root;

    class TreeNode
    {
    public:
        size_t nodeId;
        ValType val;
        size_t left_child, right_child;
    public:
        TreeNode();
        TreeNode(const ValType& value);

        void set_nodeId(size_t id);
        size_t get_nodeId() const;
        TreeNode* get_child(size_t i) const;
        TreeNode* get_sibling() const;
    };

    LCA<TreeNode> lca;
public:
    void reserve(size_t n);
    void push_back(const ValType& value);
    void clear();
    void preprocess();
    ValType query(size_t i, size_t j) const;
};

template<class ValType>
RMQLinear<ValType>::TreeNode::TreeNode(): 
    left_child(-1), right_child(-1), root(-1)
{}

template<class ValType>
RMQLinear<ValType>::TreeNode::TreeNode(const ValType& value):
    val(value), left_child(-1), right_child(-1), root(-1)
{}

template<class ValType>
void RMQLinear<ValType>::TreeNode::set_nodeId(size_t id)
{
    nodeId = id;
}

template<class ValType>
size_t RMQLinear<ValType>::TreeNode::get_nodeId() const
{
    return nodeId;
}

template<class ValType>
typename RMQLinear<ValType>::TreeNode* RMQLinear<ValType>::TreeNode::get_child(size_t i) const
{
    size_t child_idx = (i == 0 ? left_child : (i == 1 ? right_child : -1));
    if(child_idx == -1)   return NULL;
    else    return &a[ left_child ];
}

template<class ValType>
typename RMQLinear<ValType>::TreeNode* RMQLinear<ValType>::TreeNode::get_sibling() const
{
    return NULL;
}

template<class ValType>
void RMQLinear<ValType>::reserve(size_t n)
{
    a.reserve(n);
}

template<class ValType>
void RMQLinear<ValType>::push_back(const ValType& value)
{
    a.push_back( value );
    // build the tree
}

template<class ValType>
void RMQLinear<ValType>::clear()
{
    root = -1;
    a.clear();
}

template<class ValType>
void RMQLinear<ValType>::preprocess()
{
    if(root == -1)  return;
    lca.preprocess(&a[root]);
}

template<class ValType>
ValType RMQLinear<ValType>::query(size_t i, size_t j) const
{
    if(root == -1)
    {
        std::cerr << "[ERROR] [RMQLinear]: the data has not been processed yet!!!" << std::endl;
        exit(-1);
    }
    return lca.query( a[i], a[j] )->val;
}

}//namespace loon

#endif
