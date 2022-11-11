//------------------------------------------------------------------------------------------------------------
//! \file   utils/bitmask_operators.hpp
//! \brief  Header file for enabling boolean operators on enum class types so that they may be used as
//!         bitfields.
//!
//! This code is based on the header distributed by Anthony Williams under the Boost Software License
//! (version 1.0, 17 August 2003) available <a href="https://www.justsoftwaresolutions.co.uk/cplusplus/using-enum-classes-as-bitfields.html">here</a>.
//! The code was forked on October 29, 2017 by Michael M. Crockatt and subsequently refactored into its
//! present form.
//!
//! \author Anthony Williams, Michael M. Crockatt
//! \date   October 2017
//
// (C) Copyright 2015 Just Software Solutions Ltd
//
// Distributed under the Boost Software License, Version 1.0.
//
// Boost Software License - Version 1.0 - August 17th, 2003
//
// Permission is hereby granted, free of charge, to any person or
// organization obtaining a copy of the software and accompanying
// documentation covered by this license (the "Software") to use,
// reproduce, display, distribute, execute, and transmit the
// Software, and to prepare derivative works of the Software, and
// to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:
//
// The copyright notices in the Software and this entire
// statement, including the above license grant, this restriction
// and the following disclaimer, must be included in all copies
// of the Software, in whole or in part, and all derivative works
// of the Software, unless such copies or derivative works are
// solely in the form of machine-executable object code generated
// by a source language processor.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
// KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
// COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE
// LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN
// CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
//------------------------------------------------------------------------------------------------------------

# ifndef __BITMASK_OPERATORS_HPP__
# define __BITMASK_OPERATORS_HPP__

# include <type_traits>


//------------------------------------------------------------------------------------------------------------
//! \brief  Used to constrain the template operators defined in this file only to scoped enumerations that are
//!         intended to act as bitmasks.
//!
//! To enable the use of these bitmask operators, the macro ENABLE_BITMASK_OPERATORS is provided to specialize
//! this template such that <code> enable = true </code>.
//!
//! \see    ENABLE_BITMASK_OPERATORS
//------------------------------------------------------------------------------------------------------------
template<typename E>
struct enable_bitmask_operators {  static const bool enable = false;  };


//------------------------------------------------------------------------------------------------------------
//! \brief  Macro used to specialize enable_bitmask_operators for the given enumeration type to enable bitmask
//!         operators.
//------------------------------------------------------------------------------------------------------------
# define ENABLE_BITMASK_OPERATORS(E) \
template<> \
struct enable_bitmask_operators<E> {  static const bool enable = true;  };


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns true if \pp{target} contains all of the bits set by \pp{desired} and false otherwise.
//!
//! \param[in]  target      The bitmask to determine the contents of.
//! \param[in]  desired     Bitmask containing bits to search for.
//------------------------------------------------------------------------------------------------------------
template<typename E>
typename std::enable_if< enable_bitmask_operators<E>::enable, bool >::type
BitmaskHasAll (

    const E target,
    const E desired
) {
    return ( target & desired ) == desired;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns true if \pp{target} contains any of the bits set by \pp{desired} and false otherwise.
//!
//! \param[in]  target      The bitmask to determine the contents of.
//! \param[in]  desired     Bitmask containing bits to search for.
//------------------------------------------------------------------------------------------------------------
template<typename E>
typename std::enable_if< enable_bitmask_operators<E>::enable, bool >::type
BitmaskHasAny (

    const E target,
    const E desired
) {
    typedef typename std::underlying_type<E>::type underlying;

    return static_cast<underlying>( target & desired ) != 0;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Bitwise OR between two bitmask values.
//------------------------------------------------------------------------------------------------------------
template<typename E>
typename std::enable_if< enable_bitmask_operators<E>::enable, E >::type
operator| (

    const E lhs,
    const E rhs
) {
    typedef typename std::underlying_type<E>::type underlying;

    return static_cast<E>(   static_cast<underlying>(lhs)
                           | static_cast<underlying>(rhs) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Bitwise AND between two bitmask values.
//------------------------------------------------------------------------------------------------------------
template<typename E>
typename std::enable_if< enable_bitmask_operators<E>::enable, E >::type
operator& (

    const E lhs,
    const E rhs
) {
    typedef typename std::underlying_type<E>::type underlying;

    return static_cast<E>(   static_cast<underlying>(lhs)
                           & static_cast<underlying>(rhs) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Bitwise XOR between two bitmask values.
//------------------------------------------------------------------------------------------------------------
template<typename E>
typename std::enable_if< enable_bitmask_operators<E>::enable, E >::type
operator^ (

    const E lhs,
    const E rhs
) {
    typedef typename std::underlying_type<E>::type underlying;

    return static_cast<E>(   static_cast<underlying>(lhs)
                           ^ static_cast<underlying>(rhs) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Bitwise NOT of bitmask value.
//------------------------------------------------------------------------------------------------------------
template<typename E>
typename std::enable_if< enable_bitmask_operators<E>::enable, E >::type
operator~ (

    const E lhs
) {
    typedef typename std::underlying_type<E>::type underlying;

    return static_cast<E>( ~static_cast<underlying>(lhs) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Bitwise OR assignment of bitmask values.
//------------------------------------------------------------------------------------------------------------
template<typename E>
typename std::enable_if< enable_bitmask_operators<E>::enable, E & >::type
operator|= (

    E & lhs,
    const E rhs
) {
    typedef typename std::underlying_type<E>::type underlying;

    lhs = static_cast<E>(   static_cast<underlying>(lhs)
                          | static_cast<underlying>(rhs) );
    return lhs;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Bitwise AND assignment of bitmask values.
//------------------------------------------------------------------------------------------------------------
template<typename E>
typename std::enable_if< enable_bitmask_operators<E>::enable, E & >::type
operator&= (

    E & lhs,
    const E rhs
) {
    typedef typename std::underlying_type<E>::type underlying;

    lhs = static_cast<E>(   static_cast<underlying>(lhs)
                          & static_cast<underlying>(rhs) );
    return lhs;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Bitwise XOR assignment of bitmask values.
//------------------------------------------------------------------------------------------------------------
template<typename E>
typename std::enable_if< enable_bitmask_operators<E>::enable, E & >::type
operator^= (

    E & lhs,
    const E rhs
) {
    typedef typename std::underlying_type<E>::type underlying;

    lhs = static_cast<E>(   static_cast<underlying>(lhs)
                          ^ static_cast<underlying>(rhs) );
    return lhs;
}


# endif // ifndef __BITMASK_OPERATORS_HPP__
