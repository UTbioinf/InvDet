namespace loon
{
/* relabel string */
template<class IndexType>
void RelabelString<IndexType>::clear()
{
    raw_to_new.clear();
    new_to_raw.clear();

    new_to_raw.reserve(1024);
}

template<class IndexType>
IndexType RelabelString<IndexType>::add_raw_id(const std::string& id)
{
    IndexType exist_id = raw_to_new.template exactMatchSearch<IndexType>( id.c_str(), id.length() );
    if(exist_id == -1)
    {
        exist_id = raw_to_new.update( id.c_str(), id.length(), new_to_raw.size() );
        new_to_raw.push_back( id );
    }
    return exist_id;
}

template<class IndexType>
const std::string& RelabelString<IndexType>::get_raw_id(IndexType id) const
{
    return new_to_raw.at( id );
}

template<class IndexType>
std::string& RelabelString<IndexType>::get_raw_id(IndexType id)
{
    return new_to_raw.at( id );
}

template<class IndexType>
IndexType RelabelString<IndexType>::get_new_id(const std::string& id) const
{
    IndexType ret_id = raw_to_new.template exactMatchSearch<IndexType>(id.c_str(), id.length());
    return ret_id;
}

template<class IndexType>
size_t RelabelString<IndexType>::size() const
{
    return new_to_raw.size();
}

/* relabel small positive integer */

template<class IntType, class IndexType>
void RelabelSmallPosInt<IntType, IndexType>::clear()
{
    raw_to_new.clear();
    new_to_raw.clear();

    new_to_raw.reserve(1024);
}

template<class IntType, class IndexType>
IndexType RelabelSmallPosInt<IntType, IndexType>::add_raw_id(IntType id)
{
    if(id >= raw_to_new.size())
        raw_to_new.resize( (id + 1) > (raw_to_new.size() << 1) ? (id+1) : (raw_to_new.size() << 1), static_cast<IndexType>(-1) );
    if(raw_to_new[ id ] == static_cast<IndexType>(-1))
    {
        raw_to_new[ id ] = new_to_raw.size();
        new_to_raw.push_back( id );
    }
    return raw_to_new[ id ];
}

template<class IntType, class IndexType>
IntType RelabelSmallPosInt<IntType, IndexType>::get_raw_id(IndexType id) const
{
    return new_to_raw.at( id );
}

template<class IntType, class IndexType>
IndexType RelabelSmallPosInt<IntType, IndexType>::get_new_id(IntType id) const
{
    return raw_to_new.at( id );
}

template<class IntType, class IndexType>
size_t RelabelSmallPosInt<IntType, IndexType>::size() const
{
    return new_to_raw.size();
}

/* relabel generic type */

template<class Label_T, class IndexType>
void Relabel<Label_T, IndexType>::clear()
{
    raw_to_new.clear();
    new_to_raw.clear();

    new_to_raw.reserve( 1024 );
}

template<class Label_T, class IndexType>
IndexType Relabel<Label_T, IndexType>::add_raw_id(const Label_T& id)
{
    typename std::map<Label_T, IndexType>::iterator it = raw_to_new.find( id );
    if(it == raw_to_new.end())
    {
        raw_to_new[ id ] = new_to_raw.size();
        new_to_raw.push_back( id );
        return (new_to_raw.size() - 1);
    }
    return (it->second);
}

template<class Label_T, class IndexType>
const Label_T& Relabel<Label_T, IndexType>::get_raw_id(IndexType id) const
{
    return new_to_raw.at( id );
}

template<class Label_T, class IndexType>
Label_T& Relabel<Label_T, IndexType>::get_raw_id(IndexType id)
{
    return new_to_raw.at( id );
}

template<class Label_T, class IndexType>
IndexType Relabel<Label_T, IndexType>::get_new_id(const Label_T& id) const
{
    typename std::map<Label_T, IndexType>::const_iterator it = raw_to_new.find( id );
    if(it == raw_to_new.end())
        return static_cast<IndexType>(-1);
    return (it->second);
}

template<class Label_T, class IndexType>
size_t Relabel<Label_T, IndexType>::size() const
{
    return new_to_raw.size();
}

}// namespace loon
