//------------------------------------------------------------------------------------------------------------
//! \file   utils/DirectedGraph.hpp
//! \brief  Header for DirectedGraph class.
//!
//! \author Michael M. Crockatt
//! \date   December 2017
//------------------------------------------------------------------------------------------------------------

# ifndef __DIRECTED_GRAPH_HPP__
# define __DIRECTED_GRAPH_HPP__


# include <algorithm>
# include <list>
# include <map>
# include <set>
# include <sstream>
# include <vector>

# include "utils/CLog.hpp"


//------------------------------------------------------------------------------------------------------------
//! \brief  Implements a directed graph data structure with some basic graph algorithms.
//------------------------------------------------------------------------------------------------------------
template<typename NodeType>
class DirectedGraph {

public:

    //========================================================================================================
    //=== CONSTRUCTORS AND DESTRUCTOR ========================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Default constructor.
    //--------------------------------------------------------------------------------------------------------
    DirectedGraph( void ) = default;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Default copy constructor.
    //--------------------------------------------------------------------------------------------------------
    DirectedGraph( const DirectedGraph & ) = default;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Default destructor.
    //--------------------------------------------------------------------------------------------------------
    ~DirectedGraph( void ) = default;


    //========================================================================================================
    //=== OPERATOR OVERLOADS =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Default copy assignment operator.
    //--------------------------------------------------------------------------------------------------------
    DirectedGraph & operator=( const DirectedGraph & ) = default;


    //========================================================================================================
    //=== PUBLIC MEMBER FUNCTIONS ============================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Adds a directed edge to the graph between the two specified nodes. Nodes are added to the
    //!         graph if they are not already in the graph.
    //--------------------------------------------------------------------------------------------------------
    void AddEdge( const NodeType from_node, const NodeType to_node );

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Prints the adjacency list of the graph to the logging interface.
    //--------------------------------------------------------------------------------------------------------
    void Print( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Computes the in-degree of each node in the graph.
    //--------------------------------------------------------------------------------------------------------
    std::map<NodeType,int64_t> ComputeInDegrees( void ) const;

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Returns a std::list containing a topological sort of the nodes of the graph.
    //--------------------------------------------------------------------------------------------------------
    std::list<NodeType> TopologicalSort( void ) const;


private:

    //========================================================================================================
    //=== INTERNAL VARIABLES =================================================================================
    //========================================================================================================

    //--------------------------------------------------------------------------------------------------------
    //! \brief  Stores the adjacency list for the graph.
    //--------------------------------------------------------------------------------------------------------
    std::map< NodeType, std::set<NodeType> > adjacency_list;

};


//============================================================================================================
//=== PUBLIC MEMBER FUNCTIONS ================================================================================
//============================================================================================================


//------------------------------------------------------------------------------------------------------------
//! \brief  Adds a directed edge to the graph between the two specified nodes. Nodes are added to the graph if
//!         they are not already in the graph.
//------------------------------------------------------------------------------------------------------------
template<typename NodeType>
void DirectedGraph<NodeType>::AddEdge (

    const NodeType from_node,
    const NodeType to_node
) {

    // First check whether nodes are in the graph. If not, add them.
    if ( ! adjacency_list.count( from_node ) )
        adjacency_list.insert( std::make_pair( from_node, std::set<NodeType>() ) );

    if ( ! adjacency_list.count( to_node ) )
        adjacency_list.insert( std::make_pair( to_node, std::set<NodeType>() ) );

    // Then add edge to adjacency list.
    adjacency_list.at( from_node ).insert( to_node );
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Prints the adjacency list of the graph to the logging interface.
//------------------------------------------------------------------------------------------------------------
template<typename NodeType>
void DirectedGraph<NodeType>::Print( void ) const {

    PRINT_LOG( "Number of nodes: %" PRId64 "\n", (int64_t) adjacency_list.size() )
    PRINT_LOG( "\n" )

    for ( auto & node : adjacency_list ) {

        std::stringstream line;
        line << node.first << ": \t";

        for ( auto & edge_to : node.second )
            line << edge_to << "\t";

        PRINT_LOG( "%s\n", line.str().c_str() )
    }
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Computes the in-degree of each node in the graph.
//!
//! \return Returns an std::map object mapping each vertex of the graph to its in-degree.
//------------------------------------------------------------------------------------------------------------
template<typename NodeType>
std::map<NodeType,int64_t> DirectedGraph<NodeType>::ComputeInDegrees( void ) const {

    std::map<NodeType,int64_t> in_degrees;

    for ( auto & node : adjacency_list )
        in_degrees.insert( std::make_pair( node.first, 0 ) );

    for ( auto & node : adjacency_list ) {
    for ( auto & edge_to : node.second ) {

        in_degrees.at( edge_to )++;
    }}

    return in_degrees;
}


//------------------------------------------------------------------------------------------------------------
//! \brief  Returns a std::list containing a topological sort of the nodes of the graph.
//!
//! \attention  An exception is thrown if the graph contains a cycle.
//!
//! \return Returns a std::list containing a topological sort of the nodes of the graph.
//------------------------------------------------------------------------------------------------------------
template<typename NodeType>
std::list<NodeType> DirectedGraph<NodeType>::TopologicalSort( void ) const {

    std::list<NodeType> result;         // List containing result of topological sort.
    std::vector<NodeType> todo_stack;   // Stack of nodes to process.
    std::set<NodeType> processed;       // Contains processed nodes.

    auto in_degrees = ComputeInDegrees();

    // Push each node with in-degree zero onto the stack.
    for ( auto & node : in_degrees ) {

        if ( node.second == 0 )
            todo_stack.push_back( node.first );
    }

    // If no nodes with in-degree zero found, then topological sort is not possible.
    if ( todo_stack.size() == 0 ) {

        std::string error_message = "Topological sort of directed graph failed: "
                                    " graph contains a cycle.\n";

        throw std::runtime_error( error_message );
    }

    // Process nodes until none left.
    while ( todo_stack.size() > 0 ) {

        // If element at the top of the stack has not been processed.
        if ( ! processed.count( todo_stack.back() ) ) {

            processed.insert( todo_stack.back() );

            for ( auto & edge_to : adjacency_list.at( todo_stack.back() ) ) {

                if ( std::find( todo_stack.begin(), todo_stack.end(), edge_to ) != todo_stack.end() ) {

                    std::string error_message = "Topological sort of directed graph failed: "
                                                " graph contains a cycle.\n";

                    throw std::runtime_error( error_message );
                }

                if ( ! processed.count( edge_to ) )
                    todo_stack.push_back( edge_to );
            }

        // otherwise prepend value to topological sort.
        } else {

            result.push_front( todo_stack.back() );
            todo_stack.pop_back();
        }
    }

    return result;
}


# endif // ifndef __DIRECTED_GRAPH_HPP__
