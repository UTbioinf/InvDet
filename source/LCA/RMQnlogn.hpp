/*

RMP: range minimum query
Michael A. Bender & Martin Farach-Colton's algorithm
<O(nlogn), O(1)>

author: zijuexiansheng
*/

#ifndef __RMQNLOGN_H
#define __RMQNLOGN_H

#include <vector>

namespace loon
{

template<class ValType>
class RMQnlogn
{
protected:
    static std::vector<size_t> log_table;
    std::vector<ValType> A;
    std::vector<std::vector<ValType> > M;

    static size_t log2(size_t n);
public:
    void reserve(size_t n);
    void clear();
    void push_back(const ValType& val);
    void preprocess();
    ValType query(size_t i, size_t j);
};

template<class ValType>
/*static*/ std::vector<size_t> RMQnlogn<ValType>::log_table(2, 0);

template<class ValType>
/*static*/ size_t RMQnlogn<ValType>::log2(size_t n)
{
    while(log_table.size() <= n)
        log_table.insert(log_table.end(), log_table.size(), 1 + log_table.back());
    return log_table[n];
}

template<class ValType>
void RMQnlogn<ValType>::reserve(size_t n)
{
    A.reserve(n);
}

template<class ValType>
void RMQnlogn<ValType>::clear()
{
    A.clear();
}

template<class ValType>
void RMQnlogn<ValType>::push_back(const ValType& val)
{
    A.push_back( val );
}

template<class ValType>
void RMQnlogn<ValType>::preprocess()
{
    size_t n = A.size();
    size_t log2_n = log2(n);
    M.assign(log2_n + 1, A);
    for(size_t i = 1; i <= log2_n; ++i)
        for(size_t j = 0; j < n; ++j)
        {
            M[i][j] = M[i-1][std::min(n-1, j + (1 << (i-1)))];
            if(M[i-1][j] < M[i][j])
                M[i][j] = M[i-1][j];
        }
}

template<class ValType>
ValType RMQnlogn<ValType>::query(size_t i, size_t j)
{
    if(i == j)  return M[0][i];
    if(i > j)   std::swap(i, j);
    size_t k = log2(j - i);
    return std::min(M[k][i], M[k][j+1-(1<<k)]);
}

}// namespace loon

#endif
