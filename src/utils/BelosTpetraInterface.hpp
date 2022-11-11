//------------------------------------------------------------------------------------------------------------
//! \file   utils/BelosTpetraInterface.hpp
//! \brief  Header for common definitions used by solver implementations using Belos and Tpetra packages.
//!
//! \author Michael Crockatt
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------

# ifndef __TRILINOS_INTERFACE_HPP__
# define __TRILINOS_INTERFACE_HPP__

# if defined (ENABLE_BELOS_TPETRA)


# include <BelosSolverManager.hpp>
# include <Tpetra_Operator.hpp>


//------------------------------------------------------------------------------------------------------------
// TPETRA TYPEDEFS
//------------------------------------------------------------------------------------------------------------

//!
//! \brief  Typedef of Tpetra::Operator.
//!
typedef Tpetra::Operator<> TpetraOperator;

//!
//! \brief  Scalar type for elements of vectors.
//!
typedef TpetraOperator::scalar_type TpetraScalar;

//!
//! \brief  Integral type for local indices.
//!
typedef TpetraOperator::local_ordinal_type TpetraLocalIndex;

//!
//! \brief  Integral type for global indices.
//!
typedef TpetraOperator::global_ordinal_type TpetraGlobalIndex;

//!
//! \brief  Type of Kokkos Node.
//!
typedef TpetraOperator::node_type TpetraNode;

//!
//! \brief  Typedef for vectors.
//!
typedef Tpetra::MultiVector< TpetraScalar, TpetraLocalIndex, TpetraGlobalIndex, TpetraNode > TpetraVec;

//!
//! \brief  Typedef for Tpetra::Map.
//!
typedef Tpetra::Map< TpetraLocalIndex, TpetraGlobalIndex, TpetraNode > TpetraMap;


//------------------------------------------------------------------------------------------------------------
// BELOS TYPEDEFS
//------------------------------------------------------------------------------------------------------------

//!
//! \brief  Typedef of Belos::LinearProblem.
//!
typedef Belos::LinearProblem< TpetraScalar, TpetraVec, TpetraOperator > BelosProblem;

//!
//! \brief  Typedef of Belos::SolverManager.
//!
typedef Belos::SolverManager< TpetraScalar, TpetraVec, TpetraOperator > BelosSolver;


# endif // if defined (ENABLE_BELOS_TPETRA)

# endif // ifndef __TRILINOS_INTERFACE_HPP__
