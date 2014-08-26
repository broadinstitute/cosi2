//
//  Header: order_statistics.hpp
//
//  Adapted from order-statistics tree implementation by Szymon Wojciechowski
//  ( https://sourceforge.net/projects/orderstatistics/ )
//
//  Original code (C) Copyright Szymon Wojciechowski 2012-2013.
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef COSI_INCLUDE_ORDER_STATISTICS_HPP
#define COSI_INCLUDE_ORDER_STATISTICS_HPP

#include <cassert>
#include <cstdlib>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <stdexcept>
#include <iterator>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/utility/value_init.hpp>
#include <boost/call_traits.hpp>
#include <boost/concept_check.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/has_minus.hpp>
#include <boost/functional.hpp>
// #include <boost/pool/pool.hpp>
// #include <boost/pool/pool_alloc.hpp>
#include <boost/utility/enable_if.hpp>
#include <cosi/utils.h>

namespace cosi {
namespace util {

// Struct: IdentityFunctor
// A unary functor that returns its argument by value.
template <typename T>
  struct IdentityFunctor {
	 typedef T argument_type;
	 typedef T result_type;
	 typedef T value_type;
	 result_type operator()( const T& x) { return x; }
};

// Struct: IdentityReferenceFunctor
// A unary functor that returns its argument by mutable reference.
template <typename T>
  struct IdentityReferenceFunctor {
	 typedef T value_type;
	 typedef T& argument_type;
	 typedef T& result_type;
	 result_type operator()( T& x) { return x; }
};

namespace detail {

const char g_iteratorsOutOufBoundsException[] = "Iterator's operation exceeded bounds";

// Struct: DummyWeight
// Aux class for representing an _unweighted_ order statistics tree.
struct DummyWeight {
	 DummyWeight& operator+=( const DummyWeight& ) { return *this; }
	 DummyWeight& operator-=( const DummyWeight& ) { return *this; }
	 friend DummyWeight operator+( const DummyWeight&, const DummyWeight& ) { return DummyWeight(); }
	 friend DummyWeight operator-( const DummyWeight&, const DummyWeight& ) { return DummyWeight(); }
	 friend bool operator==( const DummyWeight&, const DummyWeight& ) { return true; }
   friend DummyWeight operator*( size_t, const DummyWeight& ) { return DummyWeight(); }
   friend DummyWeight operator*( const DummyWeight&, size_t ) { return DummyWeight(); }
};

// Struct: DummyWeightExtractor
// Aux class for representing an _unweighted_ order statistics tree.
template <typename ValueType>
struct DummyWeightExtractor {
	 typedef ValueType argument_type;
	 typedef DummyWeight result_type;
	 typedef DummyWeight value_type;
	 result_type operator()( const ValueType& ) const { return DummyWeight(); }
};

// Function: SystemSucceed2
// Execute given system command; 
inline void SystemSucceed2( std::string cmd ) {
	int exitCode = std::system( cmd.c_str() );
	std::cerr << "ran " << cmd << " exit code " << exitCode << "\n";
	if ( exitCode ) throw std::runtime_error( "error running command " + cmd );
}

namespace has_insertion_operator_impl {
typedef char no;
typedef char yes[2];

struct any_t {
	 template<typename T> any_t( T const& );
};

no operator<<( std::ostream const&, any_t const& );

yes& test( std::ostream& );
no test( no );

template<typename T>
struct has_insertion_operator {
	 static std::ostream &s;
	 static T const &t;
	 static bool const value = sizeof( test(s << t) ) == sizeof( yes );
};
}

template<typename T>
struct has_insertion_operator :
		 has_insertion_operator_impl::has_insertion_operator<T> {
};

}  // namespace detail


//
// Class: order_statistics_tree
//
// Stores an ordered set of objects.  Each object has a key (which may be the object itself),
// and an optional weight.  In addition to the usual map operations,
// supports computing partial sum of objects' weights for all objects with keys in a specified range,
// and finding the first object for which the sum of preceding keys is at least a given value.
//
// Template params:
//
//    ValueType - type of value stored in the container.
//    Comparator - functor specifying the less-than comparison operator between keys
//
// todo:
//
//   - option to use as intrusive container
//   - find distance betw pair of iterators faster than just subtracting positions
//   - option to cache the position when it is determined, and the timestamp for which it is valid.
//     then when determining position can go up and
//   - option to provide insert hints, emplace constructors, and difference operator.
//
//   - option to not store actual value, just store partial sum?
//     (can always get value as a diff)  -- like fenwick tree.
//     can get value in log time then.
//
  
template
<typename ValueType, typename Comparator = std::less<ValueType>,
 typename KeyExtractor = IdentityFunctor<ValueType>,
 typename WeightExtractor = detail::DummyWeightExtractor<ValueType> >
	class order_statistics_tree
	{
		 
	public:
		 typedef ValueType value_type;
		 typedef Comparator comparator_type;
		 typedef KeyExtractor key_extractor_type;
		 typedef WeightExtractor weight_extractor_type;
		 
		 typedef typename boost::unary_traits<key_extractor_type>::result_type key_type;
		 typedef typename boost::unary_traits<weight_extractor_type>::result_type weight_type;

		 typedef weight_type weight_difference_type;

		 typedef typename boost::binary_traits<comparator_type>::first_argument_type cmp_arg1_type;
		 typedef typename boost::binary_traits<comparator_type>::second_argument_type cmp_arg2_type;
		 typedef typename boost::binary_traits<comparator_type>::result_type cmp_result_type;
		 
		 BOOST_CONCEPT_ASSERT(( boost::Const_BinaryPredicate<comparator_type,key_type,key_type> ));
		 BOOST_CONCEPT_ASSERT(( boost::UnaryFunction<key_extractor_type,key_type,value_type> ));
		 
		 static void mpl_chk() {
			 
			 BOOST_MPL_ASSERT(( boost::is_same< typename boost::unary_traits<WeightExtractor>::argument_type, ValueType > ));
			 BOOST_MPL_ASSERT(( boost::is_same< typename boost::unary_traits<key_extractor_type>::argument_type, ValueType > ));
			 BOOST_MPL_ASSERT(( boost::is_same< cmp_arg1_type, cmp_arg2_type > ));
			 BOOST_MPL_ASSERT(( boost::is_convertible< key_type, cmp_arg1_type >));
			 BOOST_MPL_ASSERT(( boost::is_convertible< key_type, cmp_arg2_type >));

		 }

		 typedef ptrdiff_t difference_type;
		 typedef size_t size_type;

		 class value_compare
		 {   // in C++98, it is required to inherit binary_function<value_type,value_type,bool>
		 protected:
				comparator_type comp;
				value_compare (comparator_type c) : comp(c) {}  // constructed with map's comparison object
		 public:
				typedef bool result_type;
				typedef value_type first_argument_type;
				typedef value_type second_argument_type;
				bool operator() (const value_type& x, const value_type& y) const
					 {
						 return comp(x.first, y.first);
					 }
		 };

		 typedef comparator_type key_compare;
		 typedef value_type& reference;
		 typedef const value_type& const_reference;
		 typedef value_type *pointer;
		 typedef const value_type *const_pointer;

	private:

		 static WeightExtractor weightExtractor;
		 static key_extractor_type keyExtractor;

		// Class: tree_node
		// Represents one node of the tree.  Represents either a real node, or the special
		// sentinel node storing no value; each tree has one such node, stored in its <m_end> field.
		class tree_node
		{
			typedef tree_node* self_pointer;
			typedef const tree_node* self_const_pointer;

		public:

			tree_node(const ValueType& p_value) : 
				m_left(0),
				m_right(0),
				m_parent(0),
				m_isRed(false),
				m_subtreeSize(1),
				m_subtreeWeight( weightExtractor( p_value ) ),
				m_subtreeWeightDelta( weight_difference_type() ),
				m_value(p_value)
			{}

			 void init( const ValueType& p_value ) {
				 m_left = NULL;
				 m_right = NULL;
				 m_parent = NULL;
				 m_isRed = false;
				 m_subtreeSize = 1;
				 m_subtreeWeight = weightExtractor( p_value );
				 m_subtreeWeightDelta = weight_difference_type();
				 m_value = p_value;
			 }

			// Field: m_left
			// Left child of this node. 
			self_pointer		 m_left;

			// Field: m_right
			// Right child of this node. 
			self_pointer		 m_right;

			// Field: m_parent
			// Parent of this node. 
			self_pointer		 m_parent;

			// Field: m_isRed
			// Whether this node of the red-black tree is red (true), or black (false).
			bool			     m_isRed;

			// Field: m_subtreeSize
			// Number of nodes in this node's subtree (including this node). 
			size_t				 m_subtreeSize;

			// Field: m_subtreeWeight
			// Sum of weights of nodes in the subtree rooted at this node.
			// Includes the effect of any <m_delta> fields within this subtree,
			// but not above it. 
			weight_type m_subtreeWeight;

			// Field: m_subtreeWeightDelta
			// A weight value implicitly added to the weight of _every_ node in this node's
			// subtree.  The true weight of a node is obtained by adding to weightExtractor( m_value )
			// the m_subtreeWeightDelta of every node on the path to the root.
			weight_difference_type m_subtreeWeightDelta;

			// Field: m_value
			// The user value stored in the node.  May include the key, extractable
			// by the key extractor, and/or the weight, extractable by the weight extractor.
			ValueType			 m_value;

			size_t leftSize() const
			{
				if(m_left == 0) return 0;
				else return m_left->m_subtreeSize;
			}

			size_t rightSize() const
			{
				if(m_right == 0) return 0;
				else return m_right->m_subtreeSize;
			}

			weight_type leftWeight() const
			{
				assert( this );
				if(m_left == 0) return weight_type();
				else return m_left->m_subtreeWeight;
			}

			weight_type rightWeight() const
			{
				if(m_right == 0) return weight_type();
				else return m_right->m_subtreeWeight;
			}

			weight_type leftWeight( weight_difference_type p_parentDelta ) const
			{
				if(m_left == 0) return weight_type();
				else return m_left->m_subtreeWeight + m_left->m_subtreeSize * p_parentDelta;
			}

			weight_type rightWeight( weight_difference_type p_parentDelta ) const
			{
				if(m_right == 0) return weight_type();
				else return m_right->m_subtreeWeight + m_right->m_subtreeSize * p_parentDelta;
			}

			 
			// Method: subtreeWeightActual
			// Return the sum of actual weights of this subtree's nodes.
			 weight_type subtreeWeightActual( weight_difference_type p_deltasAbove = weight_difference_type() ) const {
				assert( !isEnd() );
				return m_subtreeWeight + ( p_deltasAbove * m_subtreeSize );
			}

			 weight_type weightWithLocalDelta() const {
				 return weightExtractor( m_value ) + m_subtreeWeightDelta;
			 }

			// Method: weightActual
			// Return the actual weight of this node, including any deltas above.  
  		weight_type weightActual() const {
				assert( !isEnd() );
				weight_type w = weightExtractor( m_value ) + m_subtreeWeightDelta;
				const tree_node *n = this;
				while ( n->m_parent && !n->m_parent->isEnd() ) {
					n = n->m_parent;
					w += n->m_subtreeWeightDelta;
				}
				return w;
			}

			// Method: ensureZeroSubtreeWeightDelta
			// Without changing any actual weights, ensure that this node has zero subtreeWeightDelta.
			void ensureZeroSubtreeWeightDelta() {
				assert( !isEnd() );
				if ( m_left ) {
					m_left->m_subtreeWeightDelta += m_subtreeWeightDelta;
					m_left->m_subtreeWeight += m_subtreeWeightDelta * m_left->m_subtreeSize;
				}
				if ( m_right ) {
					m_right->m_subtreeWeightDelta += m_subtreeWeightDelta;
					m_right->m_subtreeWeight += m_subtreeWeightDelta * m_right->m_subtreeSize;
				}
				weightExtractor( m_value ) += m_subtreeWeightDelta;
				m_subtreeWeightDelta = weight_type();
			}
			 

			self_pointer min() const
			{
				self_const_pointer l_result = this;
				while(l_result->m_left != 0)
					l_result = l_result->m_left;
				return const_cast<self_pointer>(l_result);
			}

			self_pointer max() const
			{
				self_const_pointer l_result = this;
				while(l_result->m_right != 0)
					l_result = l_result->m_right;
				return const_cast<self_pointer>(l_result);
			}

			self_pointer previous() const
			{
				self_const_pointer l_currentNode = this;
				if(l_currentNode->m_left != 0) return l_currentNode->m_left->max();

				self_const_pointer l_currentParent = l_currentNode ->m_parent;
				while(l_currentNode == l_currentParent->m_left)
				{
					l_currentNode = l_currentParent;
					l_currentParent = l_currentParent->m_parent;
				}

				if( !l_currentNode->isEnd() )
				{
					l_currentNode = l_currentParent;
				}

				return const_cast<self_pointer>(l_currentNode);
			}

			self_pointer next() const
			{
				self_const_pointer l_currentNode = this;
				if(l_currentNode->m_right != 0) return l_currentNode->m_right->min();

				self_const_pointer l_currentParent = l_currentNode ->m_parent;
				while(l_currentNode == l_currentParent->m_right)
				{
					l_currentNode = l_currentParent;
					l_currentParent = l_currentParent->m_parent;
				}

				if( !l_currentNode->isEnd() )
				{
					l_currentNode = l_currentParent;
				}

				return const_cast<self_pointer>(l_currentNode);
			}

			size_t position() const
			{
				if(m_parent == 0) return 0;
				if( isEnd() ) return m_parent->m_subtreeSize;

				size_t l_lesserValuesCounter = leftSize();
				self_const_pointer l_currentNode = this;

				while( !l_currentNode->isEnd() )
				{
					if(l_currentNode == l_currentNode->m_parent->m_right)
					{
						l_lesserValuesCounter += l_currentNode->m_parent->leftSize() + 1;
					}
					l_currentNode = l_currentNode->m_parent;
				}

				return l_lesserValuesCounter - l_currentNode->m_parent->m_subtreeSize - 1; // roots value was added unwillingly
			}

			size_t reversePosition() const
			{
				if(m_parent == 0) return 0;
				if( isEnd() ) return m_parent->m_subtreeSize;

				size_t l_lesserValuesCounter = rightSize();
				self_const_pointer l_currentNode = this;

				while( !l_currentNode->isEnd() )
				{
					if(l_currentNode == l_currentNode->m_parent->m_left)
					{
						l_lesserValuesCounter += l_currentNode->m_parent->rightSize() + 1;
					}
					l_currentNode = l_currentNode->m_parent;
				}

				return l_lesserValuesCounter - l_currentNode->m_parent->m_subtreeSize - 1; // roots value was added unwillingly
			}

			 template <typename WeightFieldExtractor>
			 typename WeightFieldExtractor::value_type cumulWeight( WeightFieldExtractor weightFieldExtractor =
																															IdentityFunctor<weight_type>() ) const
			{
				typedef typename WeightFieldExtractor::value_type weightField_t;
				
				if(m_parent == 0) return weightField_t();
				if( isEnd() ) return weightFieldExtractor( m_parent->m_subtreeWeight );

				weightField_t l_lesserValuesWeight = weightFieldExtractor( leftWeight() );
				size_t l_lesserValuesCounter = leftSize();
				self_const_pointer l_currentNode = this;

				while( !l_currentNode->isEnd() )
				{
					l_lesserValuesWeight +=
						 l_lesserValuesCounter * weightFieldExtractor( l_currentNode->m_subtreeWeightDelta );
					
					if(l_currentNode == l_currentNode->m_parent->m_right)
					{
						l_lesserValuesWeight += weightFieldExtractor( l_currentNode->m_parent->leftWeight() )
							 + weightFieldExtractor( weightExtractor( l_currentNode->m_parent->m_value ) );

						l_lesserValuesCounter += l_currentNode->m_parent->leftSize() + 1;
					}
					l_currentNode = l_currentNode->m_parent;
				}

				return l_lesserValuesWeight - weightFieldExtractor( l_currentNode->m_parent->m_subtreeWeight ); // roots value was added unwillingly
			}

			 weight_type cumulWeight() const {
				 return this->template cumulWeight< IdentityFunctor<weight_type> >();
			 }

			 // Method: addWeight
			 //
			 // Add a given value to the weight of this node; propagate subtree sums up the tree.
			 //
			 // Template params:
			 //
			 //    WeightFieldRefExtractor - if given, specifies which part of the weight to increment.
			 template <typename WeightFieldRefExtractor>
			 void addWeight( typename WeightFieldRefExtractor::value_type weightDelta,
											 WeightFieldRefExtractor refExtractor = IdentityReferenceFunctor<weight_type>() ) {
				 assert( !isEnd() );
				 //std::cerr << "addWeight: key=" << keyExtractor( m_value ) << " curWeight=" << weightExtractor( m_value ) << " weightDelta =" << weightDelta << std::endl;
				 refExtractor( weightExtractor( m_value ) ) += weightDelta;
				 self_pointer l_currentNode = this;
				 while ( !l_currentNode->isEnd() ) {
					 refExtractor( l_currentNode->m_subtreeWeight ) += weightDelta;
					 l_currentNode = l_currentNode->m_parent;
				 }
			 }

			 void addWeight( weight_type weightDelta ) {
				 this->addWeight( weightDelta, IdentityReferenceFunctor<weight_type>() );
			 }


			 void addWeightFrom( weight_difference_type weightDelta ) {
				 self_pointer l_currentNode = this;

				 bool arrivedFromLeft = true;
				 weight_type weightAdded = weight_type();
				 while ( !l_currentNode->isEnd() ) {
					 if ( arrivedFromLeft ) {
						 if ( l_currentNode->m_right ) {
							 l_currentNode->m_right->m_subtreeWeightDelta += weightDelta;
							 weight_type rightSubtreeWeightAdded = l_currentNode->m_right->m_subtreeSize * weightDelta;
							 l_currentNode->m_right->m_subtreeWeight += rightSubtreeWeightAdded;
							 weightAdded += rightSubtreeWeightAdded;
						 }
						 weightExtractor( l_currentNode->m_value ) += weightDelta;
						 weightAdded += weightDelta;
					 }
					 l_currentNode->m_subtreeWeight += weightAdded;
					 
					 arrivedFromLeft = ( l_currentNode->m_parent->m_left == l_currentNode );
					 l_currentNode = l_currentNode->m_parent;
				 }
			 }

			inline bool isEnd() const
			{
				return m_parent == m_right;
			}
		};

		typedef tree_node node;
		typedef tree_node* node_pointer;
		typedef const tree_node* const_node_pointer;
		 typedef order_statistics_tree<ValueType, comparator_type, key_extractor_type, WeightExtractor>& tree_reference;
		 typedef const order_statistics_tree<ValueType, comparator_type, key_extractor_type, WeightExtractor>& const_tree_reference;

	public:
	
		 class citerator: public boost::iterator_facade< citerator,
																										 /* value type = */ const ValueType,
																										 /* traversal category = */ boost::bidirectional_traversal_tag >
		{
			typedef citerator				self;
			typedef tree_node				node;
			 friend class order_statistics_tree<ValueType, comparator_type, key_extractor_type, WeightExtractor>;
		public:
			explicit citerator(const node* p_pointingNode = 0) : m_pointingNode(p_pointingNode) {}
			 
			size_t position() const
			{
				assert( m_pointingNode );
				return m_pointingNode->position();
			}
			 
			 template <typename WeightFieldExtractor>
			 typename WeightFieldExtractor::value_type cumulWeight( WeightFieldExtractor weightFieldExtractor =
																															IdentityFunctor<weight_type>() ) const
			{
				return m_pointingNode->cumulWeight( weightFieldExtractor );
			}
			 weight_type cumulWeight() const {
				 return this->template cumulWeight< IdentityFunctor<weight_type> >();
			 }

			 weight_type weightActual() const { assert( m_pointingNode ); return m_pointingNode->weightActual(); }
			 

		private:
			const node* m_pointingNode;

			 friend class boost::iterator_core_access;

			 void increment() {
				 if(m_pointingNode->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
				 m_pointingNode = m_pointingNode->next();
			 }
			 void decrement() {
				const node* l_previousNeighbor = m_pointingNode->previous();
				if(l_previousNeighbor->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
				m_pointingNode = l_previousNeighbor;
			 }

			 bool equal( citerator const& other ) const { return m_pointingNode == other.m_pointingNode; }

			 const ValueType& dereference() const {
				 assert( m_pointingNode && !m_pointingNode->isEnd() );
				 return m_pointingNode->m_value;
			 }

			 
		}; // const_iterator

		 class niterator: public boost::iterator_facade< niterator,
																										 ValueType,
																										 boost::bidirectional_traversal_tag >
																										 
		{
			typedef niterator				self;
			typedef tree_node				node;
			 friend class order_statistics_tree<ValueType, comparator_type, key_extractor_type, WeightExtractor>;
		public:
			// typedef std::bidirectional_iterator_tag		iterator_category;
			// typedef ptrdiff_t							difference_type;
			// typedef ValueType							value_type;
			// typedef ValueType*						pointer;
			// typedef ValueType&						reference;

			explicit niterator(node* p_pointingNode = 0) : m_pointingNode(p_pointingNode) {}

			// reference operator*() const
			// {
			// 	assert( m_pointingNode && !m_pointingNode->isEnd() );
			// 	return m_pointingNode->m_value;
			// }

			// pointer operator->() const
			// {
			// 	assert( m_pointingNode && !m_pointingNode->isEnd() );
			// 	return &(m_pointingNode->m_value);
			// }

			// inline bool operator==(const self& p_rhs) const
			// {
			// 	return m_pointingNode == p_rhs.m_pointingNode;
			// }

			// inline bool operator!=(const self& p_rhs) const
			// {
			// 	return m_pointingNode != p_rhs.m_pointingNode;
			// }

			// self& operator++()
			// {
			// 	if(m_pointingNode->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
			// 	m_pointingNode = m_pointingNode->next();
			// 	return *this;
			// }

			// self operator++(int)
			// {
			// 	if(m_pointingNode->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
			// 	self l_temporaryIterator = *this;
			// 	m_pointingNode = m_pointingNode->next();
			// 	return l_temporaryIterator;
			// }

			// self& operator--()
			// {
			// 	node* l_previousNeighbor = m_pointingNode->previous();
			// 	if(l_previousNeighbor->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
			// 	m_pointingNode = l_previousNeighbor;
			// 	return *this;
			// }

			// self operator--(int)
			// {
			// 	node* l_previousNeighbor = m_pointingNode->previous();
			// 	if(l_previousNeighbor->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
			// 	self l_temporaryIterator = *this;
			// 	m_pointingNode = l_previousNeighbor;
			// 	return l_temporaryIterator;
			// }

			size_t position() const
			{
				return m_pointingNode->position();
			}

			 weight_type weightActual() const { assert( m_pointingNode ); return m_pointingNode->weightActual(); }
			 

			// Method: cumulWeight
			// Cumulative weight of all lesser nodes
			 template <typename WeightFieldExtractor>
			 typename WeightFieldExtractor::value_type cumulWeight( WeightFieldExtractor weightFieldExtractor =
																														  IdentityFunctor<weight_type>()) const
			{
				return m_pointingNode->cumulWeight( weightFieldExtractor );
			}
			 weight_type cumulWeight() const {
				 return this->template cumulWeight< IdentityFunctor<weight_type> >();
			 }
			 

			operator citerator() const
			{
				return citerator(m_pointingNode);
			}

			 template <typename WeightFieldRefExtractor /*= */ >
			 void addWeight( typename WeightFieldRefExtractor::value_type weightDelta,
											 WeightFieldRefExtractor weightFieldRefExtractor = IdentityReferenceFunctor<weight_type>()
				 ) {
				 m_pointingNode->template addWeight< WeightFieldRefExtractor >( weightDelta, weightFieldRefExtractor );
			 }
			 
			 void addWeight( weight_difference_type weightDelta ) {
				 this->addWeight( weightDelta, IdentityReferenceFunctor<weight_type>() );
			 }

			 void addWeightFrom( weight_difference_type weightDelta ) {
				 m_pointingNode->addWeightFrom( weightDelta );
			 }


		private:
			node* m_pointingNode;

			 friend class boost::iterator_core_access;

			 void increment() {
				 if(m_pointingNode->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
				 m_pointingNode = m_pointingNode->next();
			 }
			 
			 void decrement() {
				node* l_previousNeighbor = m_pointingNode->previous();
				if(l_previousNeighbor->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
				m_pointingNode = l_previousNeighbor;
			 }

			 bool equal( niterator const& other ) const { return m_pointingNode == other.m_pointingNode; }

			 ValueType& dereference() const {
				 assert( m_pointingNode && !m_pointingNode->isEnd() );
				 return m_pointingNode->m_value;
			 }

			 
		};
		/// iterator

		class criterator
		{
			typedef criterator				self;
			typedef tree_node				node;
			 friend class order_statistics_tree<ValueType, comparator_type, key_extractor_type, WeightExtractor>;
		public:
			typedef std::bidirectional_iterator_tag		iterator_category;
			typedef ptrdiff_t							difference_type;
			typedef ValueType							value_type;
			typedef const ValueType*					pointer;
			typedef const ValueType&					reference;

			explicit criterator(const node* p_pointingNode = 0) : m_pointingNode(p_pointingNode) {}

			reference operator*() const
			{
				assert( m_pointingNode && !m_pointingNode->isEnd() );
				return m_pointingNode->m_value;
			}

			pointer operator->() const
			{
				assert( m_pointingNode && !m_pointingNode->isEnd() );
				return &(m_pointingNode->m_value);
			}

			inline bool operator==(const self& p_rhs) const
			{
				return m_pointingNode == p_rhs.m_pointingNode;
			}

			inline bool operator!=(const self& p_rhs) const
			{
				return m_pointingNode != p_rhs.m_pointingNode;
			}

			self& operator++()
			{
				if(m_pointingNode->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
				m_pointingNode = m_pointingNode->previous();
				return *this;
			}

			self operator++(int)
			{
				if(m_pointingNode->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
				self l_temporaryIterator = *this;
				m_pointingNode = m_pointingNode->previous();
				return l_temporaryIterator;
			}

			self& operator--()
			{
				const node* l_previousNeighbor = m_pointingNode->next();
				if(l_previousNeighbor->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
				m_pointingNode = l_previousNeighbor;
				return *this;
			}

			self operator--(int)
			{
				const node* l_previousNeighbor = m_pointingNode->next();
				if(l_previousNeighbor->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
				self l_temporaryIterator = *this;
				m_pointingNode = l_previousNeighbor;
				return l_temporaryIterator;
			}

			size_t position() const
			{
				return m_pointingNode->reversePosition();
			}

		private:
			const node* m_pointingNode;
		};
		/// const_reverse_iterator

		class riterator
		{
			typedef riterator				self;
			typedef tree_node				node;
			 friend class order_statistics_tree<ValueType, comparator_type, key_extractor_type, WeightExtractor>;
		public:
			typedef std::bidirectional_iterator_tag		iterator_category;
			typedef ptrdiff_t							difference_type;
			typedef ValueType							value_type;
			typedef const ValueType*					pointer;
			typedef const ValueType&					reference;

			explicit riterator(node* p_pointingNode = 0) : m_pointingNode(p_pointingNode) {}

			reference operator*() const
			{
				assert( m_pointingNode && !m_pointingNode->isEnd() );
				return m_pointingNode->m_value;
			}

			pointer operator->() const
			{
				assert( m_pointingNode && !m_pointingNode->isEnd() );
				return &(m_pointingNode->m_value);
			}

			inline bool operator==(const self& p_rhs) const
			{
				return m_pointingNode == p_rhs.m_pointingNode;
			}

			inline bool operator!=(const self& p_rhs) const
			{
				return m_pointingNode != p_rhs.m_pointingNode;
			}

			self& operator++()
			{
				if(m_pointingNode->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
				m_pointingNode = m_pointingNode->previous();
				return *this;
			}

			self operator++(int)
			{
				if(m_pointingNode->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
				self l_temporaryIterator = *this;
				m_pointingNode = m_pointingNode->previous();
				return l_temporaryIterator;
			}

			self& operator--()
			{
				node* l_previousNeighbor = m_pointingNode->next();
				if(l_previousNeighbor->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
				m_pointingNode = l_previousNeighbor;
				return *this;
			}

			self operator--(int)
			{
				node* l_previousNeighbor = m_pointingNode->next();
				if(l_previousNeighbor->isEnd()) throw std::runtime_error(detail::g_iteratorsOutOufBoundsException);
				self l_temporaryIterator = *this;
				m_pointingNode = l_previousNeighbor;
				return l_temporaryIterator;
			}

			size_t position() const
			{
				return m_pointingNode->reversePosition();
			}

			operator criterator() const
			{
				return criterator(m_pointingNode);
			}

			operator niterator() const
			{
				return niterator(m_pointingNode);
			}

		private:
			node* m_pointingNode;
		};
		/// reverse_iterator

		typedef niterator iterator;
		typedef citerator const_iterator;
		typedef riterator reverse_iterator;
		typedef criterator const_reverse_iterator;

		order_statistics_tree(const comparator_type& p_comparator = comparator_type()) : m_end(ValueType()), m_comparator(p_comparator), m_size(0)
		{
			setRoot(0);
		} 

		order_statistics_tree(const_tree_reference p_copyPattern) : m_end(ValueType()), m_size(0)
		{
			setRoot(0);
			this->operator=(p_copyPattern);
		}

		template <typename InputIterator>
		order_statistics_tree(InputIterator p_firstElement, const InputIterator & p_lastElement, const comparator_type& p_comparator = comparator_type()) : m_end(ValueType()), m_comparator(p_comparator), m_size(0)
		{
			setRoot(0);
			while(p_firstElement != p_lastElement)
			{
				insert(*p_firstElement);
				++p_firstElement;
			}
		}
	
		~order_statistics_tree()
		{
			clear();
			COSI_IF_DEBUG( freeListReporter.countDestr() );
		} 

		iterator begin()
		{
			return iterator(m_root != 0 ? m_root->min() : &m_end );
		}

		 

		iterator end()
		{
			return iterator( &m_end );
		}

		const_iterator cbegin() const
		{
			return const_iterator(m_root != 0 ? m_root->min() : &m_end );
		}

		 const_iterator begin() const { return cbegin(); }
		 const_iterator end() const { return cend(); }
		 const_iterator rbegin() const { return crbegin(); }
		 const_iterator rend() const { return crend(); }

		const_iterator cend() const
		{
			return const_iterator( &m_end );
		}

		reverse_iterator rbegin()
		{
			return reverse_iterator(m_root != 0 ? m_root->max() : &m_end );
		}

		reverse_iterator rend()
		{
			return reverse_iterator( &m_end );
		}

		const_reverse_iterator crbegin() const
		{
			return const_reverse_iterator(m_root != 0 ? m_root->max() : &m_end );
		}

		const_reverse_iterator crend() const
		{
			return const_reverse_iterator( &m_end );
		}
		 // typedef boost::fast_pool_allocator<node, boost::default_user_allocator_new_delete,
		 // 																		boost::details::pool::null_mutex > ost_alloc_t;
																				

		 static node_pointer freeList;
		 static size_t fromFreeList, fromMem, numDestroyed;

		 struct FreeListReporter {
				FreeListReporter() {  }
				~FreeListReporter() {
#if !defined(NDEBUG) && defined(COSI_DEV_PRINT)					
					std::cerr << " freelist info for " <<
						 typeid(*this).name() <<
						 " from free list: " << fromFreeList << " fromMem: "
																				<< fromMem << " numDestroyed=" << numDestroyed << std::endl;
#endif					
					while ( freeList ) {
						node_pointer p = freeList->m_left;
						delete freeList;
						freeList = p;
					}
				}

				void countDestr() { numDestroyed++;
					// std::cerr << " freelist dump for " <<
					// 	 typeid(*this).name() << std::endl;
				}
				
		 };
		 
		 static FreeListReporter freeListReporter;
		 
		 static node_pointer allocNode( const ValueType& p_insertingValue ) {
			//  node_pointer l_insertingNode = ost_alloc_t::allocate();
			// return new (l_insertingNode ) node(p_insertingValue);
			 if ( freeList ) {
				 node_pointer p = freeList;
				 freeList = freeList->m_left;
				 p->init( p_insertingValue );
				 COSI_IF_DEBUG( fromFreeList++ );
				 return p;
			 } else {
				 COSI_IF_DEBUG( fromMem++ );
				 return new node( p_insertingValue );
			 }
		 }
		 static void freeNode( node_pointer n ) {
			 n->m_left = freeList;
			 freeList = n;
			 // delete n; /*ost_alloc_t::deallocate( n );*/
		 }

		std::pair<iterator, bool> insert(const ValueType& p_insertingValue)
		{
			node_pointer l_insertingNode = allocNode( p_insertingValue );
			insertAndRebalance(l_insertingNode);

			return std::make_pair(iterator(l_insertingNode), true);
		}

		size_t erase(const key_type& p_deletingValue)
		{
			std::pair<iterator, iterator> l_equalRange = equal_range(p_deletingValue);
			size_t l_sizeBeforeDeleting = m_size;
			erase(l_equalRange.first, l_equalRange.second);
			return l_sizeBeforeDeleting - m_size;
		}

		iterator erase(const_iterator p_begin, const_iterator p_end)
		{
			while(p_begin != p_end)
			{
				node_pointer l_deletingNode = const_cast<node_pointer>(p_begin.m_pointingNode);
				++p_begin;
				deleteNode(l_deletingNode);
			}
			return iterator(const_cast<node_pointer>(p_end.m_pointingNode));
		}

		void erase(const_iterator l_deletingElement)
		{
			assert( l_deletingElement != end() );
			node_pointer l_deletingNode = const_cast<node_pointer>(l_deletingElement.m_pointingNode);
			deleteNode(l_deletingNode);
		}
		 

		size_t count (const key_type& p_searchingValue) const
		{
			std::pair<const_iterator, const_iterator> l_equalRange = equal_range(p_searchingValue);
			return l_equalRange.second.position() - l_equalRange.first.position();
		}

		const_iterator find(const key_type& p_searchingValue) const
		{
			const_iterator l_foundValue = lower_bound(p_searchingValue);
			return ( l_foundValue == end() || (*const_cast<comparator_type*>(&m_comparator))(p_searchingValue, keyExtractor(*l_foundValue)) ) ? end() : l_foundValue;
		}

		iterator find(const key_type& p_searchingValue)
		{
			iterator l_foundValue = lower_bound(p_searchingValue);
			return ( l_foundValue == end() || (*const_cast<comparator_type*>(&m_comparator))(p_searchingValue, keyExtractor(*l_foundValue)) ) ? end() : l_foundValue;
		}

		 
		iterator lower_bound(const key_type& p_searchingValue) const
		{
			return lower_bound_impl(p_searchingValue, m_root, &m_end);
		}

		iterator upper_bound(const key_type& p_searchingValue) const
		{
			return upper_bound_impl(p_searchingValue, m_root, &m_end);
		}

		std::pair<iterator,iterator> equal_range(const key_type& p_searchingValue)
		{
			const_node_pointer l_currentNode = m_root;
			const_node_pointer l_foundValue = &m_end;
			comparator_type* l_comparator = const_cast<comparator_type*>(&m_comparator);

			while(l_currentNode != 0)
			{
				if((*l_comparator)(keyExtractor(l_currentNode->m_value), p_searchingValue))
				{
					l_currentNode = l_currentNode->m_right;
				}
				else if((*l_comparator)(p_searchingValue, keyExtractor(l_currentNode->m_value)))
				{
					l_foundValue = l_currentNode;
					l_currentNode = l_currentNode->m_left;
				}
				else
				{
					const_node_pointer l_rootForUpperBound = l_currentNode;
					const_node_pointer l_searchingEnd = l_foundValue;
					l_foundValue = l_currentNode;
					l_currentNode = l_currentNode->m_left;
					l_rootForUpperBound = l_rootForUpperBound->m_right;

					return std::make_pair(lower_bound_impl(p_searchingValue, l_currentNode, l_foundValue),
									 upper_bound_impl(p_searchingValue, l_rootForUpperBound, l_searchingEnd));
				}
			}

			return std::make_pair(iterator(const_cast<node_pointer>(l_foundValue)), iterator(const_cast<node_pointer>(l_foundValue)));
		}


		std::pair<const_iterator,const_iterator> equal_range(const key_type& p_searchingValue) const
		{
			const_node_pointer l_currentNode = m_root;
			const_node_pointer l_foundValue = &m_end;
			comparator_type* l_comparator = const_cast<comparator_type*>(&m_comparator);

			while(l_currentNode != 0)
			{
				if((*l_comparator)(keyExtractor(l_currentNode->m_value), p_searchingValue))
				{
					l_currentNode = l_currentNode->m_right;
				}
				else if((*l_comparator)(p_searchingValue, keyExtractor(l_currentNode->m_value)))
				{
					l_foundValue = l_currentNode;
					l_currentNode = l_currentNode->m_left;
				}
				else
				{
					const_node_pointer l_rootForUpperBound = l_currentNode;
					const_node_pointer l_searchingEnd = l_foundValue;
					l_foundValue = l_currentNode;
					l_currentNode = l_currentNode->m_left;
					l_rootForUpperBound = l_rootForUpperBound->m_right;

					return std::make_pair(lower_bound_impl(p_searchingValue, l_currentNode, l_foundValue),
									 upper_bound_impl(p_searchingValue, l_rootForUpperBound, l_searchingEnd));
				}
			}

			return std::make_pair(const_iterator(const_cast<node_pointer>(l_foundValue)), const_iterator(const_cast<node_pointer>(l_foundValue)));
		}
		 

		iterator operator[](size_t p_order) const
		{
			return iterator(selectK(m_root, ++p_order));
		}

		 // Method: psum_inclusive_upper_bound
		 // Finds the first node for which the sum of weights of all nodes through and including n  is strictly greater than 'p_cumulWeight',
		 // or end() if no such node.
		 // Returns the node, and sets *residue to p_cumulWeight-(sum_of_weights_through_n).
		 const_iterator upper_bound_for_inclusive_partial_sum( weight_type p_cumulWeight, weight_type *residue) const
		{
			const_node_pointer p_root = m_root;
			
			if(p_root == 0 || p_root->m_subtreeWeight <= p_cumulWeight) {
				*residue = weight_type();
				return end();
			}

			weight_type l_deltaAboveRoot = weight_type();

			while(true)
			{
				assert( p_root );
				assert( !p_root->isEnd() );

				assert( static_cast< weight_type >( p_root->m_subtreeWeight +
																						l_deltaAboveRoot * p_root->m_subtreeSize )  > p_cumulWeight );
				
				weight_type l_leftWeight = p_root->leftWeight( l_deltaAboveRoot + p_root->m_subtreeWeightDelta );
				if ( p_cumulWeight < l_leftWeight ) {
					l_deltaAboveRoot += p_root->m_subtreeWeightDelta;
					p_root = p_root->m_left;
				} else {
					p_cumulWeight -= l_leftWeight;
					weight_type l_rootWeight = weightExtractor( p_root->m_value ) + l_deltaAboveRoot
						 + p_root->m_subtreeWeightDelta;

					if ( p_cumulWeight < l_rootWeight ) {
						*residue = l_rootWeight - p_cumulWeight;
						return const_iterator( p_root );
					} else {
						p_cumulWeight -= l_rootWeight;
						l_deltaAboveRoot += p_root->m_subtreeWeightDelta;
						p_root = p_root->m_right;
					}
				}
			}
		}
		 

		inline size_t size() const
		{
			return m_size;
		}

		 weight_type totalWeight() const { return !m_root ? weight_type() : m_root->m_subtreeWeight; }

		inline bool empty() const
		{
			return !m_size;
		}

		void clear()
		{
			deleteTree(m_root);
			setRoot(0);
			m_size = 0;
		}

		void swap (tree_reference p_swappingObject)
		{
			std::swap(p_swappingObject.m_size, m_size);
			std::swap(p_swappingObject.m_root, m_root);
			std::swap(p_swappingObject.m_end.m_left, m_end.m_left);
			std::swap(p_swappingObject.m_end.m_right, m_end.m_right);
			std::swap(p_swappingObject.m_end.m_parent, m_end.m_parent);
			std::swap(p_swappingObject.m_comparator, m_comparator);

			if(m_root) 
			{
				m_root->m_parent = &m_end;
			}

			if(p_swappingObject.m_root) 
			{
				p_swappingObject.m_root->m_parent = &p_swappingObject.m_end;
			}
		}

		tree_reference operator=(const order_statistics_tree<ValueType, comparator_type, key_extractor_type, WeightExtractor>& p_copyPattern)
		{
			if(this != &p_copyPattern)
			{
				clear();
				copy(p_copyPattern.m_root, m_root);
				setRoot(m_root);
				m_size = p_copyPattern.m_size;
				m_comparator = p_copyPattern.m_comparator;
			}

			return *this;
		}

		comparator_type key_comp() const
		{
			return m_comparator;
		}

		 bool is_less( const_iterator i1, const_iterator i2 ) const {
			 return ( i1 != this->end() ) &&
					( i2 == this->end()   || this->m_comparator( this->keyExtractor( *i1 ),
																											 this->keyExtractor( *i2 ) ) );
		 }

		 // Method: addWeight
		 // Add the specified delta to the weight of every item with key in the given range.
		 template <typename WeightFieldRefExtractor>
		 void addWeight( iterator beg, iterator end, typename WeightFieldRefExtractor::value_type delta,
										 WeightFieldRefExtractor weightFieldRefExtractor ) {
			 for ( iterator i = beg; i != end; i++ )
					i.addWeight( delta, weightFieldRefExtractor );
		 }
		 
		 void addWeight( iterator beg, iterator end, weight_type delta ) {
			 this->addWeight( beg, end, delta, IdentityReferenceFunctor<weight_type>() );
		 }


		 void addWeightFast( iterator beg, iterator end, weight_difference_type delta ) {
			 beg.addWeightFrom( delta );
			 end.addWeightFrom( -delta );
		 }

	private:

		void leftRotation(node_pointer p_rotationNode)
		{
			assert( p_rotationNode );
			node_pointer l_newSubTreeRoot = p_rotationNode -> m_right;

			p_rotationNode->ensureZeroSubtreeWeightDelta();
			if ( l_newSubTreeRoot ) l_newSubTreeRoot->ensureZeroSubtreeWeightDelta();
			
			p_rotationNode -> m_right = l_newSubTreeRoot->m_left;
			if(l_newSubTreeRoot->m_left != 0)
			{
				l_newSubTreeRoot->m_left ->m_parent = p_rotationNode;
			}
			l_newSubTreeRoot->m_parent = p_rotationNode->m_parent;
			if(p_rotationNode->m_parent == &m_end)
			{
				setRoot(l_newSubTreeRoot);
			}
			else
			{
				if(p_rotationNode == p_rotationNode->m_parent->m_left)
				{
					p_rotationNode->m_parent->m_left = l_newSubTreeRoot;
				}
				else
				{
					p_rotationNode->m_parent->m_right = l_newSubTreeRoot;
				}
			}
			l_newSubTreeRoot->m_left = p_rotationNode;
			p_rotationNode->m_parent = l_newSubTreeRoot;

			l_newSubTreeRoot->m_subtreeSize = p_rotationNode->m_subtreeSize;
			l_newSubTreeRoot->m_subtreeWeight = p_rotationNode->m_subtreeWeight;
			p_rotationNode->m_subtreeSize = p_rotationNode->leftSize() + p_rotationNode->rightSize() + 1;
			p_rotationNode->m_subtreeWeight = p_rotationNode->leftWeight() + p_rotationNode->rightWeight() +
				 weightExtractor( p_rotationNode->m_value ) ;
		}

		void rightRotation(node_pointer p_rotationNode)
		{
			assert( p_rotationNode );
			node_pointer l_newSubTreeRoot = p_rotationNode -> m_left;


			p_rotationNode->ensureZeroSubtreeWeightDelta();
			if ( l_newSubTreeRoot ) l_newSubTreeRoot->ensureZeroSubtreeWeightDelta();
			
			p_rotationNode -> m_left = l_newSubTreeRoot->m_right;
			if(l_newSubTreeRoot->m_right != 0)
			{
				l_newSubTreeRoot->m_right ->m_parent = p_rotationNode;
			}
			l_newSubTreeRoot->m_parent = p_rotationNode->m_parent;
			if(p_rotationNode->m_parent == &m_end)
			{
				setRoot(l_newSubTreeRoot);
			}
			else
			{
				if(p_rotationNode == p_rotationNode->m_parent->m_right)
				{
					p_rotationNode->m_parent->m_right = l_newSubTreeRoot;
				}
				else
				{
					p_rotationNode->m_parent->m_left = l_newSubTreeRoot;
				}
			}
			l_newSubTreeRoot->m_right = p_rotationNode;
			p_rotationNode->m_parent = l_newSubTreeRoot;

			l_newSubTreeRoot->m_subtreeSize = p_rotationNode->m_subtreeSize;
			l_newSubTreeRoot->m_subtreeWeight = p_rotationNode->m_subtreeWeight;
			p_rotationNode->m_subtreeSize = p_rotationNode->leftSize() + p_rotationNode->rightSize() + 1;
			p_rotationNode->m_subtreeWeight = p_rotationNode->leftWeight() + p_rotationNode->rightWeight()
				 + weightExtractor( p_rotationNode->m_value );
		}

		inline void setRoot(node_pointer p_newRoot)
		{
			m_root = p_newRoot;
			m_end.m_left = p_newRoot;
			m_end.m_right = p_newRoot;
			m_end.m_parent = p_newRoot;
			if(p_newRoot != 0)
			{
				p_newRoot->m_parent = &m_end;
			}
		}

		void insertAndRebalance(node_pointer p_checkingNode)
		{
			appendNode(p_checkingNode);
			p_checkingNode->m_isRed = true;
			node_pointer l_uncle;
			while(p_checkingNode != m_root && p_checkingNode->m_parent->m_isRed)
			{
				if(p_checkingNode->m_parent->m_parent != &m_end && p_checkingNode->m_parent == p_checkingNode->m_parent->m_parent->m_left)
				{
					l_uncle = p_checkingNode->m_parent->m_parent->m_right;
					if(l_uncle != 0 && l_uncle->m_isRed)
					{
						p_checkingNode->m_parent->m_isRed = false;
						l_uncle->m_isRed = false;
						p_checkingNode->m_parent->m_parent->m_isRed = true;
						p_checkingNode = p_checkingNode->m_parent->m_parent;
					}
					else
					{
						if(p_checkingNode == p_checkingNode->m_parent->m_right)
						{
							p_checkingNode = p_checkingNode->m_parent;
							leftRotation(p_checkingNode);
						}
						p_checkingNode->m_parent->m_isRed = false;
						if(p_checkingNode->m_parent->m_parent != &m_end)
						{
							p_checkingNode->m_parent->m_parent->m_isRed = true;
							rightRotation(p_checkingNode->m_parent->m_parent);
						}
					}
				}
				else
				{
					if(p_checkingNode->m_parent->m_parent != &m_end && (l_uncle = p_checkingNode->m_parent->m_parent->m_left) != 0 && l_uncle->m_isRed)
					{
						p_checkingNode->m_parent->m_isRed = false;
						l_uncle->m_isRed = false;
						p_checkingNode->m_parent->m_parent->m_isRed = true;
						p_checkingNode = p_checkingNode->m_parent->m_parent;
					}
					else
					{
						if(p_checkingNode == p_checkingNode->m_parent->m_left)
						{
							p_checkingNode = p_checkingNode->m_parent;
							rightRotation(p_checkingNode);
						}
						p_checkingNode->m_parent->m_isRed = false;
						if(p_checkingNode->m_parent->m_parent != &m_end)
						{
							p_checkingNode->m_parent->m_parent->m_isRed = true;
							leftRotation(p_checkingNode->m_parent->m_parent);
						}
					}
				}
			}
			m_root->m_isRed = false;
		}

		// Method: appendNode
		// Append the given node to the tree, maintaining the binary search tree invariant (but possibly
		// breaking red-black tree invariants).  If there are node(s) with keys equal to the key of the
		// inserted node, the inserted node is inserted AFTER all such nodes.
		void appendNode(const node_pointer p_inserting)
		{
			node_pointer l_insertionParent = &m_end;
			node_pointer l_insertionPlace = m_root;

			while(l_insertionPlace != 0)
			{
				l_insertionPlace->ensureZeroSubtreeWeightDelta();
				
				l_insertionParent = l_insertionPlace;
				++(l_insertionPlace->m_subtreeSize);
				l_insertionPlace->m_subtreeWeight += weightExtractor( p_inserting->m_value );
				if(m_comparator( keyExtractor( p_inserting->m_value ),
												 keyExtractor( l_insertionPlace->m_value ) ) )
				{
					l_insertionPlace = l_insertionPlace->m_left;
				}
				else
				{
					l_insertionPlace = l_insertionPlace->m_right;
				}
			}

			p_inserting->m_parent = l_insertionParent;
			if(l_insertionParent == &m_end)
			{
				setRoot(p_inserting);
			}
			else
			{
				if(m_comparator( keyExtractor( p_inserting->m_value ),
												 keyExtractor( l_insertionParent->m_value) ))
				{
					l_insertionParent->m_left = p_inserting;
				}
				else
				{
					l_insertionParent->m_right = p_inserting;
				}
			}
			++m_size;
		}

		node_pointer selectK(node_pointer p_root, size_t p_index) const
		{
			if(p_root == 0 || p_root->m_subtreeSize < p_index) return const_cast<node_pointer>(&m_end);

			while(true)
			{
				size_t l_lesserValuesCounter = p_root->leftSize() + 1;
				if(p_index == l_lesserValuesCounter) return p_root;
				else if(p_index < l_lesserValuesCounter)
				{
					p_root = p_root->m_left;
				}
				else
				{
					p_root = p_root->m_right;
					p_index -= l_lesserValuesCounter;
				}
			}
		}


		void deleteTree(node_pointer p_deletingNode)
		{
			if(p_deletingNode != 0)
			{
				deleteTree(p_deletingNode->m_left);
				deleteTree(p_deletingNode->m_right);
				freeNode( p_deletingNode );
				--m_size;
			}
		}

		iterator lower_bound_impl(const key_type& p_searchingValue, const_node_pointer p_currentNode, const_node_pointer p_foundValue) const
		{
			comparator_type* l_comparator = const_cast<comparator_type*>(&m_comparator);
			while(p_currentNode != 0)
			{
				if( (*l_comparator)(keyExtractor(p_currentNode->m_value), p_searchingValue))
				{
					p_currentNode = p_currentNode->m_right;
				}
				else
				{
					p_foundValue = p_currentNode;
					p_currentNode = p_currentNode->m_left;
				}
			}

			return iterator(const_cast< node_pointer >(p_foundValue));
		}

		iterator upper_bound_impl(const key_type& p_searchingValue, const_node_pointer p_currentNode, const_node_pointer p_foundValue) const
		{
			comparator_type* l_comparator = const_cast<comparator_type*>(&m_comparator);
			while(p_currentNode != 0)
			{
				if( (*l_comparator)(p_searchingValue, keyExtractor(p_currentNode->m_value)))
				{
					p_foundValue = p_currentNode;
					p_currentNode = p_currentNode->m_left;
				}
				else
				{
					p_currentNode = p_currentNode->m_right;
				}
			}

			return iterator(const_cast< node_pointer >(p_foundValue));
		}

		void copy(const tree_node* const p_copySource, tree_node* & p_copyTarget)
		{
			if(p_copySource != 0)
			{
				p_copyTarget = allocNode(p_copySource->m_value);

				p_copyTarget->m_isRed = p_copySource->m_isRed;
				p_copyTarget->m_subtreeSize = p_copySource->m_subtreeSize;
				p_copyTarget->m_subtreeWeight = p_copySource->m_subtreeWeight;

				copy(p_copySource->m_left, p_copyTarget->m_left);
				if(p_copyTarget->m_left != 0) p_copyTarget->m_left->m_parent = p_copyTarget;

				copy(p_copySource->m_right, p_copyTarget->m_right);
				if(p_copyTarget->m_right != 0) p_copyTarget->m_right->m_parent = p_copyTarget;
			}
		}

		std::vector<tree_node *> pathToRoot;

		void deleteNode(const node_pointer p_deletingNode)
		{
			assert( p_deletingNode );
			boost::value_initialized<ValueType> v;
			node null(v);
			null.m_subtreeSize = 0;
			null.m_subtreeWeight = weight_type();

			node_pointer y = 0;
			node_pointer x = 0;

			if(p_deletingNode->m_left == 0 || p_deletingNode->m_right == 0)
			{
				y = p_deletingNode;
			}
			else
			{
				y = p_deletingNode->next();
			}

			if(y->m_left != 0)
			{
				x = y->m_left;
			}
			else
			{
				x = y->m_right;
			}

			if(x == 0)
			{
				x = &null;
			}
			
			{
				pathToRoot.clear();
				for( tree_node *n = y; !n->isEnd(); n = n->m_parent )
					 pathToRoot.push_back( n );
				while ( !pathToRoot.empty() ) {
					pathToRoot.back()->ensureZeroSubtreeWeightDelta();
					pathToRoot.pop_back();
				}
			}

			// p_deletingNode->ensureZeroSubtreeWeightDelta();
			// if ( y != p_deletingNode ) y->ensureZeroSubtreeWeightDelta();

			//  std::map< std::string, node_pointer> names;
			//  names.clear();
			//  names.insert( std::make_pair( "x", x ) );
			//  names.insert( std::make_pair( "y", y ) );
			//  names.insert( std::make_pair( "del", p_deletingNode ) );
			//  names.insert( std::make_pair( "m_end", &m_end ) );
			//  names.insert( std::make_pair( "null", &null ) );
			// printSVG( "/home/unix/ilya/public_html/nmid", &names );

			//weight_type weightToSubtract = weightExtractor( y->m_value );
			weight_type y_weight = weightExtractor( y->m_value );

			node_pointer l_nodeToResizeItsSubtree = y;
			while(l_nodeToResizeItsSubtree != &m_end )
			{
				--(l_nodeToResizeItsSubtree->m_subtreeSize);
				// if ( l_nodeToResizeItsSubtree == p_deletingNode )
				// 	 weightToSubtract = weightExtractor( p_deletingNode->m_value );
				// cerr << "subtracting " << weightToSubtract << " from node with value " << l_nodeToResizeItsSubtree->m_value << endl;
//				l_nodeToResizeItsSubtree->m_subtreeWeight -= weightToSubtract;
				// cerr << "subtracted " << weightToSubtract << " from node with value " << l_nodeToResizeItsSubtree->m_value << " subtreeWeight here now " << l_nodeToResizeItsSubtree->m_subtreeWeight <<  endl;

				// weight_type rnodeWeight = weightExtractor( l_nodeToResizeItsSubtree->m_value );
				// if (  ) rnodeWeight = weight_type();
				// else if ( l_nodeToResizeItsSubtree == p_deletingNode ) rnodeWeight = y_weight;

				l_nodeToResizeItsSubtree->m_subtreeWeight = ( l_nodeToResizeItsSubtree == y ? weight_type() :
																											( l_nodeToResizeItsSubtree == p_deletingNode ? y_weight : weightExtractor( l_nodeToResizeItsSubtree->m_value ) ) )
					 + l_nodeToResizeItsSubtree->leftWeight() + l_nodeToResizeItsSubtree->rightWeight()
					 + l_nodeToResizeItsSubtree->m_subtreeWeightDelta * l_nodeToResizeItsSubtree->m_subtreeSize;

				// TODO: (option to) if weight didn't change here, and we're past p_deletingNode, can stop updating;
				//       will optimize e.g. for weight + operation being max or min.
				
				l_nodeToResizeItsSubtree = l_nodeToResizeItsSubtree->m_parent;
			}
			//			printSVG( "/home/unix/ilya/public_html/naft", &names );
		
			// names.clear();
			// names.insert( std::make_pair( "x", x ) );
			// names.insert( std::make_pair( "y", y ) );
			// names.insert( std::make_pair( "del", p_deletingNode ) );
			// printSVG( "/home/unix/ilya/public_html/mid", &names );

			if(p_deletingNode != y)
			{
				p_deletingNode->m_left->m_parent = y;
				y->m_left = p_deletingNode->m_left;

				if(y != p_deletingNode->m_right)
				{
					x->m_parent = y->m_parent;
					y->m_parent->m_left = x;
					y->m_right = p_deletingNode->m_right;
					p_deletingNode->m_right->m_parent = y;
				}
				else
				{
					x->m_parent = y;
				}

				if(p_deletingNode == m_root)
				{
					setRoot(y);
				}
				else if(p_deletingNode->m_parent->m_left == p_deletingNode)
				{
					p_deletingNode->m_parent->m_left = y;
				}
				else
				{
					p_deletingNode->m_parent->m_right = y;
				}

				y->m_parent = p_deletingNode->m_parent;
				std::swap(y->m_isRed, p_deletingNode->m_isRed);
				y->m_subtreeSize = p_deletingNode->m_subtreeSize;
				y->m_subtreeWeight = p_deletingNode->m_subtreeWeight;
				y = p_deletingNode;
			}
			else
			{
				x->m_parent = y->m_parent;
				if(m_root == p_deletingNode)
				{
					setRoot(x);
				}
				else if(p_deletingNode->m_parent->m_left == p_deletingNode)
				{
					p_deletingNode->m_parent->m_left = x;
				}
				else
				{
					p_deletingNode->m_parent->m_right = x;
				}
			}
			
			// names.clear();
			// names.insert( std::make_pair( "x", x ) );
			// names.insert( std::make_pair( "y", y ) );
			// names.insert( std::make_pair( "del", p_deletingNode ) );
			// printSVG( "/home/unix/ilya/public_html/bfix", &names );

			// cerr << "checking bfix" << endl;
			// check2();
			// cerr << "checked bfix" << endl;
			

			if(!y->m_isRed)
			{
				RBdeleteFix(x);
			}
			
			// names.clear();
			// names.insert( std::make_pair( "x", x ) );
			// names.insert( std::make_pair( "y", y ) );
			// names.insert( std::make_pair( "del", p_deletingNode ) );
			// printSVG( "/home/unix/ilya/public_html/afix", &names );

			if(x == &null)
			{
				if(m_root == &null)
				{
					setRoot(0);
				}
				else if(x->m_parent->m_left == x)
				{
					x->m_parent->m_left = 0;
				}
				else
				{
					x->m_parent->m_right = 0;
				}
			}

			// names.clear();
			// names.insert( std::make_pair( "x", x ) );
			// names.insert( std::make_pair( "y", y ) );
			// names.insert( std::make_pair( "del", p_deletingNode ) );
			// printSVG( "/home/unix/ilya/public_html/befdel", &names );
			

			--m_size;
			freeNode( y );
		}

		void RBdeleteFix(node_pointer x)
		{
			node_pointer w = 0;
			while(x != m_root && !x->m_isRed)
			{
				if(x == x->m_parent->m_left)
				{
					w = x->m_parent->m_right;
					if(w->m_isRed)
					{
						w->m_isRed = false;
						x->m_parent->m_isRed = true;
						leftRotation(x->m_parent);
						w = x->m_parent->m_right;
					}
					if((w->m_left == 0 || !w->m_left->m_isRed) && (w->m_right == 0 || !w->m_right->m_isRed))
					{
						w->m_isRed = true;
						x = x->m_parent;
					}
					else
					{
						if(w->m_right == 0 || !w->m_right->m_isRed)
						{
							w->m_left->m_isRed = false;
							w->m_isRed = true;
							rightRotation(w);
							w = x->m_parent->m_right;
						}
						w->m_isRed = x->m_parent->m_isRed;
						x->m_parent->m_isRed = false;
						if(w->m_right != 0)
							w->m_right->m_isRed = false;
						leftRotation(x->m_parent);
						x = m_root;
					}
				}
				else
				{
					w = x->m_parent->m_left;
					if(w->m_isRed)
					{
						w->m_isRed = false;
						x->m_parent->m_isRed = true;
						rightRotation(x->m_parent);
						w = x->m_parent->m_left;
					}
					if((w->m_right == 0 || !w->m_right->m_isRed) && (w->m_left == 0 || !w->m_left->m_isRed))
					{
						w->m_isRed = true;
						x = x->m_parent;
					}
					else
					{
						if(w->m_left == 0 || !w->m_left->m_isRed)
						{
							w->m_right->m_isRed = false;
							w->m_isRed = true;
							leftRotation(w);
							w = x->m_parent->m_left;
						}
						w->m_isRed = x->m_parent->m_isRed;
						x->m_parent->m_isRed = false;
						if(w->m_left != 0)
							w->m_left->m_isRed = false;
						rightRotation(x->m_parent);
						x = m_root;
					}
				}
			}
			x->m_isRed = false;
		}

		node			m_end;
		comparator_type		m_comparator;
		size_t			m_size;
		node_pointer	m_root;
//		 boost::pool::fast_pool_allocator<node> m_allocator;

	public:

		 size_type max_size() const { return std::map<key_type,value_type>().max_size(); }

		 void printSubtree( std::ostream& s, node_pointer p_node, std::string extraAttrs, const std::map< std::string, node_pointer> *nodeNames = NULL ) {
			 if ( p_node ) {
				 if ( p_node != &m_end ) {
					 printSubtree( s, p_node->m_left, extraAttrs, nodeNames );
					 printSubtree( s, p_node->m_right, extraAttrs, nodeNames );
				 }
				 std::string rootLbl;
				 if ( p_node == m_root )
						rootLbl = ",color=green";
				 s << "n" << (void *)p_node << " [ label=\"";
				 s << p_node->m_value;
					s << "; " <<
						p_node->m_subtreeSize << "; ";
				 s << p_node->m_subtreeWeight;
				 s << "; dl=" << p_node->m_subtreeWeightDelta;
				 s << "; " << ( p_node->m_isRed ? "R" : "B" ) << "\" " << extraAttrs <<
						rootLbl << " ] ;\n";
				 if ( p_node->m_left ) s << "n" << p_node << " -> n" << p_node->m_left << " [ label = \"L\" ];\n";
				 if ( p_node->m_right ) s << "n" << p_node << " -> n" << p_node->m_right << " [ label = \"R\" ];\n";
				 //if ( p_node->m_parent ) s << "n" << p_node << " -> n" << p_node->m_parent << " [ label = \"P\", style=dashed ];\n";
			 }
		 }
		 
		 void printDot( std::string fname, const std::map< std::string, node_pointer> *nodeNames = NULL ) {
			 std::ofstream f( fname.c_str() );
			 f << "digraph ost { ordering=out;\n";
			 printSubtree( f, m_root, "", nodeNames );
			 printSubtree( f, &m_end, ", color=red", nodeNames );
			 if ( nodeNames ) 
					for ( typename std::map< std::string, node_pointer >::const_iterator it = nodeNames->begin(); it != nodeNames->end(); it++ ) {
						f << it->first << " [ label=\"" << it->first << "\", shape=box, color=blue ];\n";
						f << it->first << " -> " << "n" << (void *)it->second << " [color=blue];\n";
					}
			 f << "}\n";
		 }

		 void printSVG( std::string fname, const std::map< std::string, node_pointer> *nodeNames = NULL ) {
			 std::string dotFile( fname + ".dot" );
			 printDot( dotFile, nodeNames );
			 detail::SystemSucceed2( std::string( "dot -Tsvg -o " ) + fname + ".svg " + dotFile );
		 }

		 void checkSubtree( node_pointer p_node, const key_type *mustBeBelow, const key_type *mustBeAbove,
												size_t *subtreeSize, weight_type *subtreeWeight, size_t *subtreeBlackHeight,
												weight_type weightDeltasAbove ) const {
			 assert( p_node );
			 assert( p_node->m_parent );
			 assert( !mustBeBelow || !( m_comparator( *mustBeBelow, keyExtractor( p_node->m_value ) ) ) );
			 assert( !mustBeAbove || !( m_comparator( keyExtractor( p_node->m_value ), *mustBeAbove ) ) );
			 assert( p_node == p_node->m_parent->m_left || p_node == p_node->m_parent->m_right );
			 assert( !p_node->m_isRed || ( ( !p_node->m_left || !p_node->m_left->m_isRed ) &&
																		 ( !p_node->m_right || !p_node->m_right->m_isRed ) ) );

			 size_t L_subtreeSize = 0;
			 weight_type L_subtreeWeight = weight_type();
			 size_t L_blackHeight = 0;
			 key_type keyHere = keyExtractor( p_node->m_value );
			 if ( p_node->m_left ) checkSubtree( p_node->m_left, &keyHere, mustBeAbove, &L_subtreeSize, &L_subtreeWeight, &L_blackHeight, weightDeltasAbove + p_node->m_subtreeWeightDelta );

			 size_t R_subtreeSize = 0;
			 weight_type R_subtreeWeight = weight_type();
			 size_t R_blackHeight = 0;
			 
			 if ( p_node->m_right ) checkSubtree( p_node->m_right, mustBeBelow, &keyHere, &R_subtreeSize, &R_subtreeWeight, &R_blackHeight, weightDeltasAbove + p_node->m_subtreeWeightDelta );

			 *subtreeSize = L_subtreeSize + 1 + R_subtreeSize;
			 assert( *subtreeSize == p_node->m_subtreeSize );
			 *subtreeWeight = L_subtreeWeight + weightExtractor( p_node->m_value ) + R_subtreeWeight
					+ p_node->m_subtreeWeightDelta * p_node->m_subtreeSize;
			 assert( *subtreeWeight == p_node->m_subtreeWeight );
			 assert( L_blackHeight == R_blackHeight );
			 *subtreeBlackHeight = L_blackHeight + ( p_node->m_isRed ? 0 : 1 );
		 }

		 // Method: check
		 // Check the internal consistency of the tree; cause assertion failure in case of error.
		 // If assert is disabled by NDEBUG conditional define, this method does nothing.
		 void check() const {
#ifndef NDEBUG			 
			 if ( !m_root ) {
				 assert( m_size == 0 );
				 assert( !m_end.m_left && !m_end.m_right && !m_end.m_parent && m_end.m_subtreeSize == 1 && m_end.m_subtreeWeight == weightExtractor( ValueType() ) );
			 } else {
				 assert( m_size == m_root->m_subtreeSize );
				 assert( m_root->m_parent == &m_end && m_end.m_left == m_root && m_end.m_right == m_root );
				 assert( !m_root->m_isRed );
				 size_t subtreeSize = 0;
				 weight_type subtreeWeight = weight_type();
				 size_t subtreeBlackHeight = 0;
				 checkSubtree( m_root, NULL, NULL, &subtreeSize, &subtreeWeight, &subtreeBlackHeight,
											 weight_type() );
			 }
#endif // #ifdef NDEBUG
		 }

	};

template<typename ValueType, typename Comparator,
				 typename KeyExtractor,
				 typename WeightExtractor >
typename order_statistics_tree<ValueType,Comparator,KeyExtractor,WeightExtractor>::node_pointer
order_statistics_tree<ValueType,Comparator,KeyExtractor,WeightExtractor>::freeList = NULL;

template<typename ValueType, typename Comparator,
				 typename KeyExtractor,
				 typename WeightExtractor >
size_t
order_statistics_tree<ValueType,Comparator,KeyExtractor,WeightExtractor>::fromFreeList = 0;

template<typename ValueType, typename Comparator,
				 typename KeyExtractor,
				 typename WeightExtractor >
size_t
order_statistics_tree<ValueType,Comparator,KeyExtractor,WeightExtractor>::fromMem = 0;


template<typename ValueType, typename Comparator,
				 typename KeyExtractor,
				 typename WeightExtractor >
size_t
order_statistics_tree<ValueType,Comparator,KeyExtractor,WeightExtractor>::numDestroyed = 0;

template<typename ValueType, typename Comparator,
				 typename KeyExtractor,
				 typename WeightExtractor >
typename order_statistics_tree<ValueType,Comparator,KeyExtractor,WeightExtractor>::FreeListReporter
order_statistics_tree<ValueType,Comparator,KeyExtractor,WeightExtractor>::freeListReporter;

template<typename ValueType, typename Comparator,
				 typename KeyExtractor,
				 typename WeightExtractor >
typename order_statistics_tree<ValueType,Comparator,KeyExtractor,WeightExtractor>::weight_extractor_type
order_statistics_tree<ValueType,Comparator,KeyExtractor,WeightExtractor>::weightExtractor;

template<typename ValueType, typename Comparator,
				 typename KeyExtractor,
				 typename WeightExtractor >
typename order_statistics_tree<ValueType,Comparator,KeyExtractor,WeightExtractor>::key_extractor_type
order_statistics_tree<ValueType,Comparator,KeyExtractor,WeightExtractor>::keyExtractor;

} // namespace util
} // namespace cosi

#endif  // #ifndef COSI_INCLUDE_ORDER_STATISTICS_HPP

