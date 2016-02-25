#ifndef INCLUDE_COSI_VECMAP_H
#define INCLUDE_COSI_VECMAP_H

#include <utility>
#include <vector>
#include <boost/optional.hpp>
#include <boost/container/stable_vector.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/typeof/typeof.hpp>

#include BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()


namespace cosi {
namespace util {

template <typename TKey, typename TVal,
					typename IdxMap>
class VecMap {
public:
	 typedef TKey key_type;
	 typedef std::pair< TKey, TVal > value_type;
	 typedef boost::container::stable_vector<value_type> vec_type;
	 typedef typename vec_type::const_iterator const_iterator;
	 typedef typename vec_type::iterator iterator;
	 typedef typename vec_type::const_pointer const_pointer;
	 typedef typename vec_type::const_reference const_reference;

	 const_iterator begin() const { return m_data.begin(); }
	 const_iterator end() const { return m_data.end(); }
	 std::size_t size() const { return m_data.size(); }

	 struct cmp_keys {
			bool operator()( const value_type& v1, const value_type& v2 ) const { return v1.first < v2.first; }
	 };

	 TVal& operator[]( TKey k ) {
		 std::size_t i = static_cast<std::size_t>( idxMap( k ) );
		 if ( i >= m_idx.size() )
				m_idx.resize( i+1 );
		 // std::cerr << "vecmap: k=" << k << " i=" << i << "\n"; 
		 // std::cerr << " m_idx size is " << m_idx.size() << "\n";
		 // std::cerr << " m_data size is " << m_data.size() << "\n";
		 assert( i < m_idx.size() );
				
		 if ( !m_idx[ i ] ) {
			 iterator it = m_data.begin();
			 while ( it != m_data.end() && it->first < k ) ++it;
			 m_idx[ i ] = m_data.emplace( it, k, TVal() );
				//boost::sort( m_data, cmp_keys() );
		 }
		 assert( m_idx[ i ] && get(m_idx[ i ])->first == k );
		 return get(m_idx[ i ])->second;
	 }

	 const_iterator find( const key_type& k ) const {
		 std::size_t i = static_cast<std::size_t>( idxMap( k ) );
		 
		 if ( i < m_idx.size() && m_idx[ i ] ) {
			 assert( get(m_idx[ i ])->first == k );
			 return get(m_idx[ i ]);
		 } else
				return m_data.end();
	 }

	 VecMap() { }

	 VecMap( VecMap const& m ):
		 m_data( m.m_data ) { reidx( m ); }

	 void reidx( VecMap const& m ) {
		 m_idx.resize( m.m_idx.size() );
		 for ( std::size_t i = 0; i < m_idx.size(); ++i )
				if ( m.m_idx[ i ] )
					 m_idx[ i ] = m_data.begin() + ( get(m.m_idx[ i ]) - m.m_data.begin() ); 
	 }

	 VecMap& operator=( VecMap const& m ) {
		 m_data = m.m_data;
		 reidx( m );
		 return *this;
	 }
	 
private:
	 vec_type m_data;
	 std::vector< boost::optional< iterator > > m_idx;
	 IdxMap idxMap;
};



template <typename K, typename V, typename I>
V const& at( const VecMap< K, V, I>& m, const K k ) {
	typename VecMap< K, V, I>::const_iterator it = m.find( k );
	cosi_chk( it != m.end(), "map lookup failed" );
	return it->second;
}

template <typename K1, typename K2, typename V, typename I1, typename I2> inline
V const& at( const VecMap< K1, VecMap< K2, V, I2>, I1 >& m, const K1 k1, const K2 k2 ) {
	return at( at( m, k1 ), k2 );
}


}  // namespace util
}  // namespace cosi

BOOST_TYPEOF_REGISTER_TEMPLATE(cosi::util::VecMap,3)


#endif // #ifndef INCLUDE_COSI_VECMAP_H
