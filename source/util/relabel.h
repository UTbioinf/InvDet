#ifndef __UTIL_RELABEL_H
#define __UTIL_RELABEL_H

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "cedar.h"

namespace loon
{
/*! \ingroup Class_util
 * @{
 */

/*! \brief %Relabel strings to contiguous integers
 *
 * RelabelString is implemented by using a trie structure (c.f. cedar.h).
 *
 * In this implementation, the `IndexType` should always be at most 4 bytes.
 * That said, (unsigned) char/short/int are preferred. Using a type that is
 * larger than 4 bytes can lead to significant and unpredictable errors.
 * 
 * The constrains on at 4 bytes may not be a disadvantage, because if the
 * index is too large, the data itself will usually be much larger. But still,
 * we need to be careful, because there genome sequencing data can be extremely
 * large.
 */
template<class IndexType>
class RelabelString
{
private:
    cedar::da<IndexType> raw_to_new;
    std::vector<std::string> new_to_raw;
public:
    void clear();//!< clear the object

    /*!\brief Add a string to the relabeling object
     * 
     * Add a (new) string to the relabeling object. If the string already exists,
     * then nothing will done, otherwise, a new integer ID will be assigned to
     * the string
     *
     * \param [in] id The string that is to be relabelled.
     */
    IndexType add_raw_id(const std::string& id);

    /*!\brief Get the raw string
     *
     * Get the original string by it's relabelling ID.
     * 
     * This member function is used by a constant object. It also has a non-constant
     * version with the same kind of usage (c.f. std::string& get_raw_id(IndexType id)).
     * 
     * \param [in] id The relabelling ID of the string
     * \return The reference to the string with the specified relabelling ID.
     */
    const std::string& get_raw_id(IndexType id) const;
    std::string& get_raw_id(IndexType id);//!< Get the raw string by it's relabelling ID. \sa const std::string& get_raw_id(IndexType id) const

    /*!\brief Get relabelling ID of a string
     * 
     * \param [in] id The string that is to be queried for its relabelling ID
     * \return The relabelling ID of the string in query. If the string was not added
     * to the object, then -1, which is cast to type `IndexType`, will be returned.
     */
    IndexType get_new_id(const std::string& id) const;
    size_t size() const;//!< Number of different strings that have been relabelled.
};
/*! @} */

/*! \ingroup Class_util
 * @{
 */
/*!\brief %Relabel small positive integers to contiguous integers
 * 
 * RelabelSmallPosInt is use full when the largest integer is not very large.
 * In this case, a hash table is more efficient than a tree-based data 
 * structure.
 *
 * Notice:
 *  * `IntType` should always be non-negative
 *  * In this implementation, `IndexType` can be any type that is a positive 
 *    integer.
 * 
 * \todo Maybe a support for negative integer
 */
template<class IntType, class IndexType>
class RelabelSmallPosInt
{
private:
    std::vector<IndexType> raw_to_new;
    std::vector<IntType> new_to_raw;
public:
    void clear();//!< clear the object

    /*!\brief Add a small non-negative integer to the relabeling object
     * 
     * Add a (new) small integer to the relabeling object. If the small non-negative
     * integer already exists, then nothing will done, otherwise, a 
     * new integer ID will be assigned to the smaller non-negative integer
     *
     * \param [in] id The small integer that is to be relabelled.
     */
    IndexType add_raw_id(IntType id);

    /*!\brief Get the original small non-negative integer
     *
     * Get the original small non-negative integer by it's relabelling ID.
     * 
     * \param [in] id The relabelling ID of the original small integer
     * \return The corresponding small non-negative integer.
     */
    IntType get_raw_id(IndexType id) const;

    /*!\brief Get relabelling ID of a small non-negative integer 
     * 
     * \param [in] id The small non-negative integer that is to be queried for its relabelling ID
     * \return The relabelling ID of the small integer in query. If the small integer was not added
     * to the object, then -1, which is cast to type `IndexType`, will be returned.
     */
    IndexType get_new_id(IntType id) const;
    size_t size() const; //!< Number of different strings that have been relabelled.
};
/*! @} */

/*! \ingroup Class_util
 * @{
 */

/*! \brief %Relabel general object to contiguous integers
 * 
 * Relabel is the most generic implementation of the relabeling utility. It supports
 * relabelling all types of objects.
 * 
 * The `Label_T` can be any data type that has "<" implemented. It uses a map<>
 * STL, so it's slower than RelabelString. But the `IndexType` is not constrained
 * to be at most 4 bytes.
 */
template<class Label_T, class IndexType>
class Relabel
{
private:
    std::map<Label_T, IndexType> raw_to_new;
    std::vector<Label_T> new_to_raw;
public:
    void clear();//!< clear the object


    /*!\brief Add a generic obj to the relabeling object
     * 
     * Add a (new) generic obj to the relabeling object. If the generic obj already exists,
     * then nothing will done, otherwise, a new integer ID will be assigned to generic obj.
     *
     * \param [in] id The generic obj that is to be relabelled.
     */
    IndexType add_raw_id(const Label_T& id);

    /*!\brief Get the raw generic obj 
     *
     * Get the original generic obj by it's relabelling ID.
     * 
     * This member function is used by a constant object. It also has a non-constant
     * version with the same kind of usage (c.f. Label_T& get_raw_id(IndexType id)).
     * 
     * \param [in] id The relabelling ID of the generic obj
     * \return The reference to the generic obj with the specified relabelling ID.
     */
    const Label_T& get_raw_id(IndexType id) const;
    Label_T& get_raw_id(IndexType id);//!< Get the raw string by it's relabelling ID. \sa const Label_T& get_raw_id(IndexType id) const;


    /*!\brief Get relabelling ID of a generic obj
     * 
     * \param [in] id The generic obj that is to be queried for its relabelling ID
     * \return The relabelling ID of the generic in query. If the generic was not added
     * to the object, then -1, which is cast to type `IndexType`, will be returned.
     */
    IndexType get_new_id(const Label_T& id) const;
    size_t size() const;//!< Number of different generic objs that have been relabelled.
};

/*! @} */
}// namespace loon

#include "relabelImpl.h"

#endif
