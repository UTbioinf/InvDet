/*A sample LCATree
author: zijuexiansheng
*/
#ifndef __LCATREE_H
#define __LCATREE_H

#include <vector>

namespace loon
{

class LCATreeNode;
class LCATreeNode
{
private:
    size_t nodeId; // nodeId is required for LCA
    std::vector<LCATreeNode* > children;
public:
    LCATreeNode();
    LCATreeNode* append_child();
    void reserve_child(size_t n);

    /* The following are the required member functions for LCA*/
    void set_nodeId(size_t id);
    size_t get_nodeId() const;
    LCATreeNode* get_child(size_t i = -1) const;
    LCATreeNode* get_sibling() const;
};

LCATreeNode::LCATreeNode(): nodeId(-1)
{}

LCATreeNode* LCATreeNode::append_child()
{
    children.push_back( new LCATreeNode() );
    return children.back();
}

void LCATreeNode::reserve_child(size_t n)
{
    children.reserve(n);
}

void LCATreeNode::set_nodeId(size_t id)
{
    nodeId = id;
}

size_t LCATreeNode::get_nodeId() const
{
    return nodeId;
}

LCATreeNode* LCATreeNode::get_child(size_t i/* = 0*/) const
{
    if(i >= children.size())
        return NULL;
    else
        return children[i];
}

LCATreeNode* LCATreeNode::get_sibling() const
{
    return NULL;
}

}

#endif
