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
    std::vector<TreeNode*> a;
    TreeNode* root;

    class TreeNode
    {
    public:
        size_t nodeId;
        ValType val;
        TreeNode *left_child, *right_child, *parent;
    public:
        TreeNode();
        TreeNode(const ValType& value);

        void set_nodeId(size_t id);
        size_t get_nodeId() const;
        TreeNode* get_child(size_t i = 0) const;
        TreeNode* get_sibling() const;
    };

    LCA<TreeNode> lca;
public:
    ~RMQLinear();
    void reserve(size_t n);
    void push_back(const ValType& value);
    ///void print_tree() const;
    ///void print_tree(TreeNode* node, int level = 0) const;
    void clear();
    void preprocess();
    ValType query(size_t i, size_t j);
};

template<class ValType>
RMQLinear<ValType>::TreeNode::TreeNode(): 
    left_child(NULL), right_child(NULL), parent(NULL), nodeId(-1)
{}

template<class ValType>
RMQLinear<ValType>::TreeNode::TreeNode(const ValType& value):
    val(value), left_child(NULL), right_child(NULL), parent(NULL), nodeId(-1)
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
typename RMQLinear<ValType>::TreeNode* RMQLinear<ValType>::TreeNode::get_child(size_t i/* = 0*/) const
{
    //return (i == 0 ? (left_child ? left_child : right_child) : 
    //                 (i == 1 ? (left_child ? right_child : NULL) : NULL));
    return (i == 0 ? left_child : (i==1? right_child : NULL));
}

template<class ValType>
typename RMQLinear<ValType>::TreeNode* RMQLinear<ValType>::TreeNode::get_sibling() const
{
    return NULL;
}

template<class ValType>
RMQLinear<ValType>::~RMQLinear()
{
    for(size_t i = 0; i < a.size(); ++i)
        delete a[i];
}

template<class ValType>
void RMQLinear<ValType>::reserve(size_t n)
{
    a.reserve(n);
}

template<class ValType>
void RMQLinear<ValType>::push_back(const ValType& value)
{
    if(a.empty())
    {
        a.push_back( new typename RMQLinear<ValType>::TreeNode(value) );
        root = a.back();
        return;
    }
    typename RMQLinear<ValType>::TreeNode* cur = a.back();
    while(cur != NULL && cur->val > value)
        cur = cur->parent;
    a.push_back( new typename RMQLinear<ValType>::TreeNode(value) );
    if(cur == NULL)
    {
        a.back()->left_child = root;
        root->parent = a.back();
        root =a.back();
    }
    else
    {
        a.back()->left_child = cur->right_child;
        if(cur->right_child != NULL)
            cur->right_child->parent = a.back();
        cur->right_child = a.back();
        a.back()->parent = cur;
    }
}

///template<class ValType>
///void RMQLinear<ValType>::print_tree() const
///{
///    std::cout << "tree array" << std::endl;
///    for(size_t i = 0; i < a.size(); ++i)
///        std::cout << a[i]->val << ' ';
///    std::cout << std::endl << "tree:" << std::endl;
///    print_tree(root);
///    std::cout << "=================================" << std::endl << std::endl;
///}
///
///template<class ValType>
///void RMQLinear<ValType>::print_tree(RMQLinear<ValType>::TreeNode* node, int level) const
///{
///    for(int i = 1; i < level; ++i)
///        std::cout << "|  ";
///    if(level)   std::cout << "└--";
///    if(node == NULL)   std::cout << std::endl;
///    else
///    {
///        std::cout << node->val << "(" << node << ")" << std::endl;
///        print_tree(node->left_child, level + 1);
///        print_tree(node->right_child, level + 1);
///    }
///}

template<class ValType>
void RMQLinear<ValType>::clear()
{
    root = NULL;
    for(size_t i = 0; i < a.size(); ++i)
        delete a[i];
    a.clear();
}

template<class ValType>
void RMQLinear<ValType>::preprocess()
{
    if(root == NULL)  return;
    for(size_t i = 0; i < a.size(); ++i)
    {
        if(a[i]->left_child == NULL)
            std::swap(a[i]->left_child, a[i]->right_child);
    }
    lca.preprocess(root);
}

template<class ValType>
ValType RMQLinear<ValType>::query(size_t i, size_t j)
{
    if(root == NULL)
    {
        std::cerr << "[ERROR] [RMQLinear]: the data has not been processed yet!!!" << std::endl;
        exit(-1);
    }
    return lca.query( a[i], a[j] )->val;
}

}//namespace loon

#endif
