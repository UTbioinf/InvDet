/*
General RMQ algorithm
*/

#ifndef __RMQGENERAL_H
#define __RMQGENERAL_H

#include "RMQLinear.hpp"
#include "RMQnlogn.hpp"

namespace loon
{

template<class ValType>
class RMQ: public RMQnlogn<ValType>
{
private:
    RMQLinear<ValType> rmq_linear;
    size_t N_nlogn, N_linear;
    short cur_choice;
    
    ValType naive_query(size_t i, size_t j) const;
public:
    RMQ(size_t nlogn_n = 10, size_t linear_n = 10000);
    void clear();
    void push_back(const ValType& val);
    void preprocess();
    ValType query(size_t i, size_t j);
};

template<class ValType>
ValType RMQ<ValType>::naive_query(size_t i, size_t j) const
{
    if(i > j)   std::swap(i, j);
    ValType ret = RMQnlogn<ValType>::A[i];
    for(++i; i <= j; ++i)
        if(ret > RMQnlogn<ValType>::A[i])  ret = RMQnlogn<ValType>::A[i];
    return ret;
}

template<class ValType>
RMQ<ValType>::RMQ(size_t nlogn_n, size_t linear_n):
    N_nlogn(nlogn_n), N_linear(linear_n), cur_choice(0)
{}

template<class ValType>
void RMQ<ValType>::clear()
{
    cur_choice = 0;
    RMQnlogn<ValType>::clear();
    rmq_linear.clear();
}

template<class ValType>
void RMQ<ValType>::push_back(const ValType& val)
{
    if(cur_choice < 2)
    {
        size_t nn = RMQnlogn<ValType>::A.size();
        if(nn >= N_linear )
        {
            cur_choice = 2;
            rmq_linear.reserve( nn );
            for(size_t i = 0; i < nn; ++i)
                rmq_linear.push_back( RMQnlogn<ValType>::A[i] );
            rmq_linear.push_back( val );
            RMQnlogn<ValType>::A.clear();
        }
        else
        {
            RMQnlogn<ValType>::push_back( val );
        }
    }
    else
    {
        rmq_linear.push_back( val );
    }
}

template<class ValType>
void RMQ<ValType>::preprocess()
{
    if(cur_choice < 2)
    {
        if(RMQnlogn<ValType>::A.size() >= N_nlogn)
        {
            cur_choice = 1;
            RMQnlogn<ValType>::preprocess();
        }
    }
    else
    {
        rmq_linear.preprocess();
    }
}

template<class ValType>
ValType RMQ<ValType>::query(size_t i, size_t j)
{
    if(cur_choice == 0) return naive_query(i, j);
    else if(cur_choice == 1)    return RMQnlogn<ValType>::query(i, j);
    else    return rmq_linear.query(i, j);
}

}// namespace loon

#endif
