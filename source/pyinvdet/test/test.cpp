#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <ctime>

#include <RMQ.hpp>
#include <LCATree.hpp>
#include <LCA.hpp>

using namespace std;

int naive_find_min(const vector<int>& A, int i, int j)
{
    if(i > j)   swap(i, j);
    int ret = A[i];
    for(++i; i <= j; ++i)
        ret = min(ret, A[i]);
    return ret;
}

void test_RMQ_nlogn()
{
    const int n = 10000;
    vector<int> a;
    loon::RMQ<int> rmq;
    
    srand(time(NULL));
    cout << "The list of integers: " << endl;
    for(int i = 0; i < n; ++i)
    {
        cout << setw(4) << i;
        a.push_back( rand() % (n * 100) + 1);
    }
    cout << endl;

    for(int i = 0; i < n; ++i)
    {
        cout << setw(4) << a[i];
        rmq.push_back(a[i]);
    }
    cout << endl;
    rmq.preprocess();

    cout << "Queries" << endl;
    for(int i = 0; i < 30 * n; ++i)
    {
        int ii = rand() % n;
        int jj = rand() % n;
        int correct = naive_find_min(a, ii, jj);
        int rmq_ret = rmq.query(ii, jj);
        cout << "[" << ii << ", " << jj << "]: " << correct << endl;
        if(correct != rmq_ret)
        {
            cerr << "ERROR" << endl;
            exit(1);
        }
    }
    cout << "Done!" << endl;
}

void generate_rand_tree(size_t n, size_t max_branch, loon::LCATreeNode& root, std::vector<loon::LCATreeNode*>& nodeQueue)
{
    nodeQueue.push_back( &root );
    loon::LCATreeNode* cur_node;
    size_t node_cnt = 1;
    size_t cur_pos = 0;

    while(node_cnt < n)
    {
        int n_children = rand() % max_branch + 1;
        cur_node = nodeQueue[ cur_pos++ ];

        for(int i = 0; i < n_children; ++i)
        {
            nodeQueue.push_back( cur_node->append_child() );
            ++node_cnt;
            if(node_cnt >= n)   return;
        }
    }
}

void print_rand_tree(loon::LCATreeNode* node, size_t level, size_t child_id)
{
    for(int i = 1; i < level; ++i)
        cout << "|  ";
    if(level > 0)   cout << "|--";
    cout << "(" << level << ", " << child_id << ")" << endl;

    loon::LCATreeNode* child;
    for(int i = 0; (child = node->get_child(i)); ++i)
        print_rand_tree(child, level + 1, i + 1);
}

void print_rand_tree(loon::LCATreeNode* node, size_t level = 0)
{
    for(int i = 1; i < level; ++i)
        cout << "|  ";
    if(level > 0)   cout << "└--";
    cout << (node->get_nodeId()) << endl;

    loon::LCATreeNode* child;
    for(int i = 0; (child = node->get_child(i)); ++i)
        print_rand_tree(child, level + 1);
}

void test_LCA()
{
    const int n = 5;
    loon::LCATreeNode root;
    loon::LCA<loon::LCATreeNode> lca;
    std::vector<loon::LCATreeNode*> tree_nodes;

    srand(3);
    generate_rand_tree(n, 3, root, tree_nodes);
    lca.preprocess(&root);

    cout << "print tree:" << endl;
    print_rand_tree(&root);

    cout << "do query" << endl;
    srand(time(NULL));
    for(int i = 0; i < 20; ++i)
    {
        int ii = rand() % n;
        int jj = rand() % n;
        loon::LCATreeNode* node1 = tree_nodes[ii];
        loon::LCATreeNode* node2 = tree_nodes[jj];
        cout << "(" << (node1->get_nodeId()) << ", " << (node2->get_nodeId()) << "): ";
        cout << flush;
        loon::LCATreeNode* res = lca.query(node1, node2);
        cout << (res->get_nodeId()) << endl;
    }
    cout << "Done!" << endl;
}

int main()
{
    //test_RMQ_nlogn();
    test_LCA();
    return 0;
}
