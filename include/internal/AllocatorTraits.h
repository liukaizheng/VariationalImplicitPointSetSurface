#ifndef ALLOCATORTRAITS_H
#define ALLOCATORTRAITS_H

#include <Eigen/Dense>
#include <iosfwd>
#include <utility>

namespace Internal {

template <typename T>
struct AllocatorTraits
{
    typedef std::allocator<T> Alloc;
};

template <>
struct AllocatorTraits<Eigen::Vector2d>
{
    typedef Eigen::aligned_allocator<Eigen::Vector2d> Alloc;
};

template <>
struct AllocatorTraits<Eigen::Vector4d>
{
    typedef Eigen::aligned_allocator<Eigen::Vector4d> Alloc;
};

template <>
struct AllocatorTraits<Eigen::Vector2f>
{
    typedef Eigen::aligned_allocator<Eigen::Vector2f> Alloc;
};

template <>
struct AllocatorTraits<Eigen::Vector4f>
{
    typedef Eigen::aligned_allocator<Eigen::Vector4f> Alloc;
};

template <>
struct AllocatorTraits<Eigen::Matrix2d>
{
    typedef Eigen::aligned_allocator<Eigen::Matrix2d> Alloc;
};

template <>
struct AllocatorTraits<Eigen::Matrix4d>
{
    typedef Eigen::aligned_allocator<Eigen::Matrix4d> Alloc;
};

template <>
struct AllocatorTraits<Eigen::Matrix2f>
{
    typedef Eigen::aligned_allocator<Eigen::Matrix2f> Alloc;
};

template <>
struct AllocatorTraits<Eigen::Matrix4f>
{
    typedef Eigen::aligned_allocator<Eigen::Matrix4f> Alloc;
};

template<typename Key, typename T>
struct AllocatorTraits<std::pair<const Key, T>>
{
    typedef std::allocator<std::pair<const Key, T> > Alloc;
};
}

#endif
