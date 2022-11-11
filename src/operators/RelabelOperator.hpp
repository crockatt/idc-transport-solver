//------------------------------------------------------------------------------------------------------------
//! \file   RelabelOperator.hpp
//! \brief  Header for RelabelOperator class.
//!
//! \author Michael M. Crockatt
//! \date   December 2017
//------------------------------------------------------------------------------------------------------------

# ifndef __RELABEL_OPERATOR_HPP__
# define __RELABEL_OPERATOR_HPP__


# include <cstdint>
# include <map>
# include <string>

# include "utils/Quadrule/Quadrule.hpp"
# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Class for performing relabeling operations mapping from one discrete ordinates quadrature set to
//!         another.
//------------------------------------------------------------------------------------------------------------
class RelabelOperator {

public:

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to specify the mathematical approach used to define the relabel operator.
    //--------------------------------------------------------------------------------------------------------
    enum class RelabelType {

        //!
        //! String descriptor: \e "none"
        //!
        //! Default null value.
        //!
        //! Does not yield a valid relabel operator. An error is thrown if this RelabelType is used to
        //! construct a RelabelOperator object.
        //!
        None

    # if SPACE_DIMS == 1

        //!
        //! String descriptor: \e "poly-interp"
        //!
        //! Defines a relabel operator on \f$ [-1,1] \f$ through polynomial interpolation.
        //!
    ,   PolyInterp

        //!
        //! String descriptor: \e "dbl-poly-interp"
        //!
        //! Defines a relabel operator on \f$ [-1,1] \f$ through polynomial interpolation over \f$ [-1,0] \f$
        //! and \f$ [0,1] \f$ separately.
        //!
    ,   DoublePolyInterp

    # elif SPACE_DIMS >= 2

        //!
        //! String descriptor: \e "sph-hyperinterp"
        //!
        //! Defines a relabel operator on \f$ \mathbb{S}^2 \f$ through hyperinterpolation by spherical
        //! harmonics.
        //!
    ,   SphHyperinterp

        //!
        //! String descriptor: \e "sph-ls"
        //!
        //! Defines a relabel operator on \f$ \mathbb{S}^2 \f$ through a least-squares reconstruction of the
        //! angular flux as an expansion of spherical harmonics.
        //!
    ,   SphLeastSquares

        //!
        //! String descriptor: \e "pwc-hierarchical"
        //!
        //! Applies a piecewise-constant hierarchical mapping between two compatible \f$ T_N \f$ quadrature
        //! sets.
        //!
    ,   PWC_Hierarchical

    # endif // if SPACE_DIMS == ?
    };


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to specify the BLAS algorithm used to apply the matrices defining the relabel operator.
    //--------------------------------------------------------------------------------------------------------
    enum class RelabelBLASOp {

        //!
        //! String descriptor: \e "none"
        //!
        //! Default null value.
        //!
        //! BLAS library routines are not used when applying the relabel operator.
        //!
        None

        //!
        //! String descriptor: \e "one-mat-gemv"
        //!
        //! Applies the relabel operator across each spatial coefficient separately using one matrix applied
        //! with the BLAS routine ?gemv.
        //!
    ,   OneMat_GEMV

        //!
        //! String descriptor: \e "one-mat-gemm"
        //!
        //! Applies the relabel operator across all spatial coefficients in a spatial cell simultaneously
        //! using one matrix applied with the BLAS routine ?gemm.
        //!
    ,   OneMat_GEMM

    # if SPACE_DIMS >= 2

        //!
        //! String descriptor: \e "two-mat-gemv"
        //!
        //! Applies the relabel operator across each spatial coefficient separately using two matrices applied
        //! sequentially with the BLAS routine ?gemv.
        //!
    ,   TwoMat_GEMV

        //!
        //! String descriptor: \e "two-mat-gemm"
        //!
        //! Applies the relabel operator across all spatial coefficients in a spatial cell simultaneously
        //! using two matrices applied sequentially with the BLAS routine ?gemm.
        //!
    ,   TwoMat_GEMM

    # endif // if SPACE_DIMS >= 2
    };


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert RelabelOperator::RelabelType values to their corresponding string descriptors.
    //--------------------------------------------------------------------------------------------------------
    static std::map< RelabelOperator::RelabelType, std::string > RelabelType_to_String;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert string descriptors to RelabelOperator::RelabelType values.
    //--------------------------------------------------------------------------------------------------------
    static std::map< std::string, RelabelOperator::RelabelType > String_to_RelabelType;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert RelabelOperator::RelabelBLASOp values to their corresponding string descriptors.
    //--------------------------------------------------------------------------------------------------------
    static std::map< RelabelOperator::RelabelBLASOp, std::string > RelabelBLASOp_to_String;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to convert string descriptors to RelabelOperator::RelabelBLASOp values.
    //--------------------------------------------------------------------------------------------------------
    static std::map< std::string, RelabelOperator::RelabelBLASOp > String_to_RelabelBLASOp;


    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete default constructor.
    //--------------------------------------------------------------------------------------------------------
    RelabelOperator( void ) = delete;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy constructor.
    //--------------------------------------------------------------------------------------------------------
    RelabelOperator( const RelabelOperator & ) = delete;


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Construct a RelabelOperator object using the given parameters.
    //--------------------------------------------------------------------------------------------------------
    RelabelOperator( const Quadrule::OrdinateSet::OrdinateType, const int64_t, const bool,
                     const Quadrule::OrdinateSet::OrdinateType, const int64_t, const bool );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Construct a RelabelOperator object using the given parameters.
    //--------------------------------------------------------------------------------------------------------
    RelabelOperator( const Quadrule::OrdinateSet &, const Quadrule::OrdinateSet & );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Destructor for class RelabelOperator.
    //--------------------------------------------------------------------------------------------------------
    ~RelabelOperator( void );


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Delete copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    RelabelOperator & operator=( const RelabelOperator & ) = delete;


    //========================================================================================================
    //=== RELABEL ROUTINES ===================================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes \f$ \vec{y} \gets \vec{y} + \alpha \mathcal{R} \vec{x} \f$.
    //--------------------------------------------------------------------------------------------------------
    void Relabel( double alpha, const RKDG::OrdinateFlux & x, RKDG::OrdinateFlux & y );


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes \f$ \vec{y} \gets \vec{y} + \alpha \mathcal{R} \vec{x} \f$.
    //--------------------------------------------------------------------------------------------------------
    void Relabel( double alpha, const STDG::OrdinateFlux & x, RKDG::OrdinateFlux & y );


    //========================================================================================================
    //=== ADDITIONAL MEMBER FUNCTIONS ========================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints a summary of the operator configuration to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( const std::string & = "  " ) const;


private:

    //========================================================================================================
    //=== PRIVATE MEMBER FUNCTIONS ===========================================================================
    //========================================================================================================


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Used to index into the array stored at RelabelOperator::R_mat_ptr.
    //!
    //! \param[in]  dst_q   Index of ordinate in destination quadrature set.
    //! \param[in]  src_q   Index of ordinate in source quadrature set.
    //!
    //! \return     Returns a reference to the specified element of the relabeling matrix.
    //!
    //! \see    RelabelOperator::SrcMat()
    //! \see    RelabelOperator::DstMat()
    //--------------------------------------------------------------------------------------------------------
    inline double & RMat( const int64_t dst_q, const int64_t src_q ) {

    # if defined (STRICT_CHECK)

        if ( R_mat_ptr == nullptr ) {

            std::string error_message = "Pointer R_mat_ptr is NULL in '" + std::string(__func__) + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::runtime_error( error_message );
        }

        if ( dst_q < 0 || dst_q >= dst.nq() ) {

            std::string error_message
                = "Value " + std::to_string(dst_q) + " for index dst_q in " + std::string(__func__)
                  + " outsize permissible range [ 0, " + std::to_string( dst.nq() ) + " ).\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

        if ( src_q < 0 || src_q >= src.nq() ) {

            std::string error_message
                = "Value " + std::to_string(src_q) + " for index src_q in " + std::string(__func__)
                  + " outsize permissible range [ 0, " + std::to_string( src.nq() ) + " ).\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

    # endif // if defined (STRICT_CHECK)

        return R_mat_ptr[ dst_q + dst.nq() * src_q ];
    }


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Indexes into the discrete-to-moment operator for the source ordinate set.
    //!
    //! \param[in]      k       Index of spherical harmonic.
    //! \param[in]      q       Index of ordinate.
    //!
    //! \return     Returns a reference to the specified element of the discrete-to-moment matrix.
    //!
    //! \see    RelabelOperator::RMat()
    //! \see    RelabelOperator::DstMat()
    //--------------------------------------------------------------------------------------------------------
    inline double & SrcMat( const int64_t k, const int64_t q ) {

    # if defined (STRICT_CHECK)

        if ( src_basis_ptr == nullptr ) {

            std::string error_message = "Pointer src_basis_ptr is NULL in '" + std::string(__func__) + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::runtime_error( error_message );
        }

        if ( q < 0 || q >= src.nq() ) {

            std::string error_message
                = "Value " + std::to_string(q) + " for index q in " + std::string(__func__)
                  + " outsize permissible range [ 0, " + std::to_string( src.nq() ) + " ).\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

    # endif // if defined (STRICT_CHECK)

        return src_basis_ptr[ k + stride * q ];
    }


    //--------------------------------------------------------------------------------------------------------
    //! \brief  Indexes into the moment-to-discrete operator for the destination ordinate set.
    //!
    //! \param[in]      k       Index of spherical harmonic.
    //! \param[in]      q       Index of ordinate.
    //!
    //! \return     Returns a reference to the specified element of the moment-to-discrete matrix.
    //!
    //! \see    RelabelOperator::RMat()
    //! \see    RelabelOperator::SrcMat()
    //--------------------------------------------------------------------------------------------------------
    inline double & DstMat( const int64_t k, const int64_t q ) {

    # if defined (STRICT_CHECK)

        if ( dst_basis_ptr == nullptr ) {

            std::string error_message = "Pointer dst_basis_ptr is NULL in '" + std::string(__func__) + "'.\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::runtime_error( error_message );
        }

        if ( q < 0 || q >= dst.nq() ) {

            std::string error_message
                = "Value " + std::to_string(q) + " for index q in " + std::string(__func__)
                  + " outsize permissible range [ 0, " + std::to_string( dst.nq() ) + " ).\n";

            PRINT_ERROR( error_message.c_str() )
            throw std::out_of_range( error_message );
        }

    # endif // if defined (STRICT_CHECK)

        return dst_basis_ptr[ k + stride * q ];
    }


    //========================================================================================================
    //=== PRIVATE MEMBER VARIABLES ===========================================================================
    //========================================================================================================


    RelabelType relabel_type;           //!< Specifies the approach used to define the relabel operator.
    RelabelBLASOp blas_op;              //!< Specifies the algorithm used to apply the relabel operator.

    Quadrule::OrdinateSet src;          //!< Ordinate set to map angular flux from.
    double * src_basis_ptr;             //!< Pointer to matrix mapping ordinates to moments.

    Quadrule::OrdinateSet dst;          //!< Ordinate set to map angular flux to.
    double * dst_basis_ptr;             //!< Pointer to matrix mapping moments to ordinates.

    int64_t interpolant_degree;         //!< Maximum degree polynomial used by relabeling operations employing (hyper)interpolation.
    int stride;                         //!< Row/column size in discrete-to-moment/moment-to-discrete matrices.

    double * R_mat_ptr;                 //!< Pointer to matrix mapping directly between ordinates.
};


# endif // ifndef __RELABEL_OPERATOR_HPP__
