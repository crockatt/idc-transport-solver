//------------------------------------------------------------------------------------------------------------
//! \file   linear_solvers/Abstract/SweepPattern/SweepPattern.cpp
//! \brief  Contains implementations and instantiations of methods from SweepPattern class template.
//!
//! \author Michael M. Crockatt
//! \date   January 2018
//------------------------------------------------------------------------------------------------------------

# include <typeinfo>

# include "linear_solvers/Abstract/SweepPattern/SweepPattern.hpp"


//============================================================================================================
//=== MEMBER DEFINITIONS FOR SweepPattern CLASS ==============================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Initializes a SweepPattern object.
//!
//! \param[in]  enclosing   Reference to instance of enclosing SweepOperator class.
// \param[in]  input_list  Contains parameters for initialization of object.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
SweepOperator<OrdinateFlux>::SweepPattern::SweepPattern (

    const SweepOperator<OrdinateFlux> & enclosing,
    const ParameterList & // input_list
) :
    sw_op{ enclosing }
{
    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPattern::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Virtual destructor.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
SweepOperator<OrdinateFlux>::SweepPattern::~SweepPattern ( void ) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPattern::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints a summary of the object's configuration to the logging interface.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::SweepPattern::Print (

    const std::string & // = "  "

) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPattern::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    /* empty */
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Sets parameters given in a ParameterList.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::SweepPattern::SetParameters (

    const ParameterList &
) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPattern::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    /* empty */
}


//--------------------------------------------------------------------------------------------------------
//! \brief  Construct lists for SIMD blocking across angular ordinates.
//!
//! \return     Returns a list of blocks (1D) or a list of blocks for each quadrant (2D) or octant (3D)
//!             for vector-blocking across the ordinates.
//--------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
# if SPACE_DIMS == 1
    std::vector< SIMD_BlkIdx<0> >
# else
    std::vector< std::vector< SIMD_BlkIdx<0> > >
# endif
SweepOperator<OrdinateFlux>::SweepPattern::ConstructSIMDAngleBlocks ( void ) const {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPattern::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

# if SPACE_DIMS == 1

    std::vector< SIMD_BlkIdx<0> > angle_blocks;

    for ( int64_t sign = 0; sign <= 1; ++sign ) {

        const int64_t q_min = this->sw_op.nq() * sign / 2;
        const int64_t q_max = this->sw_op.nq() * (sign + 1) / 2;

        for ( int64_t q = q_min; q < q_max; q += SWEEP_SIMD_LEN ) {

            // Full blocks.
            if ( q + SWEEP_SIMD_LEN <= q_max ) {

                angle_blocks.push_back(
                    SIMD_BlkIdx<0>{
                        /* .q   = */ q,
                        /* .i   = */ -1,
                        /* .len = */ SWEEP_SIMD_LEN
                    }
                );

            // Partial blocks.
            } else {

                // Determine number of angles in partial block.
                const int64_t q_remainder = q_max - q;

                // Peel remainder into decreasing powers of two.
                for ( int64_t simd_len = ( ((int64_t)(SWEEP_SIMD_LEN)) >> 1 );
                              simd_len > 0;
                              simd_len >>= 1
                ) {

                    if ( simd_len & q_remainder ) {

                        // Construct block of smaller vector length.
                        angle_blocks.push_back(
                            SIMD_BlkIdx<0>{
                                /* .q   = */ q,
                                /* .i   = */ -1,
                                /* .len = */ simd_len
                            }
                        );

                        // Increment q past block constructed above.
                        q += simd_len;
                    }

                } // end for simd_len.
            } // end else.

        } // end for q.
    } // end for sign.

# elif SPACE_DIMS == 2

    std::vector< std::vector< SIMD_BlkIdx<0> > > angle_blocks { std::vector< SIMD_BlkIdx<0> >{},
                                                                std::vector< SIMD_BlkIdx<0> >{},
                                                                std::vector< SIMD_BlkIdx<0> >{},
                                                                std::vector< SIMD_BlkIdx<0> >{} };

    for ( int64_t quad = 0; quad < 4; ++quad ) {

        const int64_t q_min = this->sw_op.Quadrants(quad);
        const int64_t q_max = this->sw_op.Quadrants(quad + 1);

        for ( int64_t q = q_min; q < q_max; q += SWEEP_SIMD_LEN ) {

            // Full blocks.
            if ( q + SWEEP_SIMD_LEN <= q_max ) {

                angle_blocks.at(quad).push_back(
                    SIMD_BlkIdx<0>{
                        /* .q        = */ q,
                        /* .i        = */ -1,
                        /* .j        = */ -1,
                        /* .len      = */ SWEEP_SIMD_LEN,
                        /* .quadrant = */ quad
                    }
                );

            // Partial blocks.
            } else {

                // Determine number of angles in partial block.
                const int64_t q_remainder = q_max - q;

                // Peel remainder into decreasing powers of two.
                for ( int64_t simd_len = ( ((int64_t)(SWEEP_SIMD_LEN)) >> 1 );
                              simd_len > 0;
                              simd_len >>= 1
                ) {

                    if ( simd_len & q_remainder ) {

                        // Construct block of smaller vector length.
                        angle_blocks.at(quad).push_back(
                            SIMD_BlkIdx<0>{
                                /* .q        = */ q,
                                /* .i        = */ -1,
                                /* .j        = */ -1,
                                /* .len      = */ simd_len,
                                /* .quadrant = */ quad
                            }
                        );

                        // Increment q back block constructed above.
                        q += simd_len;
                    }

                } // end for simd_len.
            } // end else.

        } // end for q.
    } // end for quad.

# elif SPACE_DIMS == 3

# endif // if SPACE_DIMS == ?

    return angle_blocks;
}



//============================================================================================================
//=== MEMBER DEFINITIONS FOR SweepPatternFactory CLASS =======================================================
//============================================================================================================


//!
//! \brief  Definition of SweepPatternFactory::factory_map_ptr.
//!
template< class OrdinateFlux >
std::map< std::string, typename SweepOperator<OrdinateFlux>::SweepPattern::SweepPatternFactory * > *
SweepOperator<OrdinateFlux>::SweepPattern::SweepPatternFactory::factory_map_ptr{ nullptr };


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a reference to the map stored at SweepPatternFactory::factory_map_ptr.
//!
//! Responsible for allocating the map if SweepPatternFactory::factory_map_ptr is null.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
std::map< std::string, typename SweepOperator<OrdinateFlux>::SweepPattern::SweepPatternFactory * > &
SweepOperator<OrdinateFlux>::SweepPattern::SweepPatternFactory::GetFactoryMap ( void ) {

    if ( factory_map_ptr == nullptr ) {

        factory_map_ptr = new std::map<
                std::string,
                typename SweepOperator<OrdinateFlux>::SweepPattern::SweepPatternFactory *
            >();
    }

    return *factory_map_ptr;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Adds a reference to a factory class to the map SweepPatternFactory::factory_map_ptr.
//!
//! Factory classes deriving from SweepPatternFactory should call this routine to add a mapping to their
//! singleton instance from the string descriptor of the SweepPattern subclass which they construct.
//!
//! \note   Because this function uses the \c insert routine of std::map each factory may only be added once.
//!
//! \param[in]  descriptor  String descriptor of the SweepPattern associated with the factory \pp{factory}.
//! \param[in]  factory     Factory object for creating objects of types deriving from SweepPattern.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
void SweepOperator<OrdinateFlux>::SweepPattern::SweepPatternFactory::AddFactory (

    const std::string & descriptor,
    SweepPatternFactory & factory
) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPattern::SweepPatternFactory::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    auto & factory_map = this->GetFactoryMap();

    if ( factory_map.count( descriptor ) ) {

        std::string error_message =   "Factory for SweepPattern '"
                                    + descriptor
                                    + "' already present in SweepOperator<"
                                    + std::string( typeid(OrdinateFlux).name() )
                                    + ">::SweepPattern::SweepPatternFactory::factory_map.\n";

        PRINT_ERROR( error_message.c_str() )
        throw std::runtime_error( error_message );
    }

    factory_map.insert( std::pair< std::string, SweepPatternFactory * >( descriptor, &factory ) );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Calls creator for desired SweepPattern implementation.
//!
//! Uses the map SweepPatternFactory::factory_map to determine which factory to call to create object.
//!
//! \param[in]  enclosing   Reference to instance of enclosing SweepOperator class.
//! \param[in]  input_list  Contains parameters for initialization of object.
//!
//! \return     Returns a pointer to the created object.
//------------------------------------------------------------------------------------------------------------
template< class OrdinateFlux >
typename SweepOperator<OrdinateFlux>::SweepPattern *
SweepOperator<OrdinateFlux>::SweepPattern::SweepPatternFactory::CreatePattern (

    const SweepOperator<OrdinateFlux> & enclosing,
    const ParameterList & input_list
) {

    PRINT_STATUS( "Executing SweepOperator<%s>::SweepPattern::SweepPatternFactory::%s.\n",
                  typeid(OrdinateFlux).name(), __func__ )

    std::string sweep_pattern_str;

    // Get string descriptor from input list.
    try {

        sweep_pattern_str = input_list.GetValue<std::string>( "sweep_pattern" );

    } catch (...) {

        PRINT_WARNING( "Using default sweep pattern.\n" )
        sweep_pattern_str = "default";
    }

    SweepPatternFactory * sweep_pattern_factory = nullptr;

    // Retrieve factory pointer for given string descriptor.
    try {

        auto & factory_map = GetFactoryMap();
        sweep_pattern_factory = factory_map.at( sweep_pattern_str );

    } catch (...) {

        std::string error_message =   "Failed to access factory for SweepPattern '"
                                    + sweep_pattern_str
                                    + "' in SweepOperator<"
                                    + std::string( typeid(OrdinateFlux).name() )
                                    + ">::SweepPattern:SweepPatternFactory::"
                                    + std::string(__func__)
                                    + ".\n";

        throw std::out_of_range( error_message );
    }

    // Create desired object using factory.
    return sweep_pattern_factory->Create( enclosing, input_list );
}


//============================================================================================================
//=== EXPLICIT INSTANTIATIONS ================================================================================
//============================================================================================================

# include "objects/RKDG/OrdinateFlux.hpp"
# include "objects/STDG/OrdinateFlux.hpp"


//!
//! \brief  Instantiation of SweepOperator<RKDG::OrdinateFlux>::SweepPattern class.
//!
template class SweepOperator<RKDG::OrdinateFlux>::SweepPattern;


//!
//! \brief  Instantiation of SweepOperator<STDG::OrdinateFlux>::SweepPattern class.
//!
template class SweepOperator<STDG::OrdinateFlux>::SweepPattern;
