#ifndef STLCONTAINERTRAITS_H
#define STLCONTAINERTRAITS_H

#include <internal/AllocatorTraits.h>
#include <vector>
#include <map>

namespace Internal {

template <typename T>
struct STLVectorTraits
{
    typedef std::vector<T, typename AllocatorTraits<T>::Alloc> Vector;
};

template <typename Key, typename T, typename Compare = std::less<Key>>
struct STLMapTraits
{
    typedef std::map<Key, T, Compare, typename AllocatorTraits<T>::Alloc> Map;
};

}
#endif
