#pragma once

#include <dvector.h>
#include <g3Debug.h>

// this is for debugging purposes only...
#include <string>
#include <strstream>


namespace g3 {


//   rcvector provides a vector of reference-counted objects,
//   with a built-in free list for re-using unreferenced indices.
//   For example in a mesh we might want to have an array of Vertex
//   elements, but if we delete a vertex we will need to leave a "hole"
//   in the index space. rcvector enables that kind of thing.
//
//   To use rcvector, your Type must contain the following functions:
#if 0
       int get_refcount() const;
       void set_refcount(int i);
       (some_time) printable() const;
#endif
//   Where the printable() function returns something safe to pass
//   to an ostrema via <<   (this can just be a dummy value if you want)
//
//   The reason for doing it this way is that it allows you to manage the 
//   actual refcount member. For example if a single used/not-used bit is
//   sufficient, you could pack that into your Type anywhere.
//
//   item_iterator and index_iterator are provided for iterating
//   over valid objects/indices. The standard begin()/end() functions 
//   provide item_iterators.
//
//   The internal vector type is templated. By default is a dvector,
//   this is a good choice. sparse_dvector can also be used if for
//   example there are guaranteed to be large gaps in the index space.
//
//   std::vector cannot be used directly because the templated VectorType
//   we do not have an an allocator, and the compilers do not allow the
//   nested template to inherit the default allocator. But you can easily
//   subclass std::vector<T, std::allocator<T>> to create a vector class
//   that will work.
//
//   TODO
//     - currently checking > 0 && != invalid in various funcs, but invalid = -1...

template<typename Type, template<typename> class VectorType = dvector >
class rcvector
{
public:
	static constexpr int invalid = -1;
		
	virtual ~rcvector();

	rcvector() : m_vFreeIndices(128) 
	{ 
		clear();
		Type t; t.get_refcount(); t.set_refcount(1);
	}
	rcvector( const rcvector<Type, VectorType> & copy ) {
		*this = copy;
	}
	rcvector( rcvector<Type, VectorType> && moved ) {
		*this = std::move(moved);
	}

	const rcvector<Type, VectorType> & 
	operator=(const rcvector<Type, VectorType> & copy ) {
		m_vData = copy.m_vData;
		m_nUsedCount = copy.m_nUsedCount;
		m_vFreeIndices = copy.m_vFreeIndices;
		return *this;
	}

	const rcvector<Type, VectorType> & 
	operator=(const rcvector<Type, VectorType> && moved ) {
		m_vData = std::move(moved.m_vData);
		m_nUsedCount = moved.m_nUsedCount;
		m_vFreeIndices = std::move(moved.m_vFreeIndices);
		return *this;
	}

	// check if element at nIndex has non-zero refcount
	inline bool isValid( unsigned int nIndex ) const {
		return ( nIndex < (int)m_vData.size() && m_vData[nIndex].get_refcount() > 0 );
	}

	// check if element at nIndex has non-zero refcount (no range check)
	inline bool isValidUnsafe( unsigned int nIndex ) const {
		return m_vData[nIndex].get_refcount() > 0;
	}

	// return refcount at nIndex
	inline int refCount( unsigned int nIndex ) const {
		gDevAssert( isValid(nIndex)  );
		return ( m_vData[nIndex].get_refcount() > 0 && m_vData[nIndex].get_refcount() != invalid )
			? m_vData[nIndex].get_refcount() : 0;
	}

	// increment refcount at nIndex. returns new refcount.
	inline int increment( unsigned int nIndex, int nIncrement = 1 ) {
		gDevAssert( isValid(nIndex)  );
		auto & e = m_vData[nIndex];
		e.set_refcount( e.get_refcount() + nIncrement );
		return e.get_refcount();
	}

	// decrement refcount at nIndex, move into free list if count hits 0
	inline void decrement( unsigned int nIndex, int nDecrement = 1 ) {
		gDevAssert( isValid(nIndex) );
		auto & e = m_vData[nIndex];
		e.set_refcount( e.get_refcount() - nDecrement );
		gDevAssert( e.get_refcount() >= 0 );
		if (e.get_refcount() == 0) {			// add to empty list
			m_vFreeIndices.push_back( nIndex );
			e.set_refcount(invalid);
			m_nUsedCount--;
		}
	}

	// insert value into vector, either using free list or appending to end.
	// returns insertion index
	inline unsigned int insert( const Type & t ) {
		m_nUsedCount++;
		if ( m_vFreeIndices.empty() ) {
			m_vData.push_back(t);
			m_vData.back().set_refcount(1);
			return (int)(m_vData.size() - 1);
		} else {
			unsigned int iFree = m_vFreeIndices.back();
			m_vFreeIndices.pop_back();
			m_vData[iFree] = t;
			m_vData[iFree].set_refcount(1);
			return iFree;
		}
	}

	// insert nIndex into free list)
	inline void remove( int nIndex ) {		// force remove
		gDevAssert( (unsigned int)nIndex < m_vData.size() );
		if ( m_vData[nIndex].get_refcount() < 0 || m_vData[nIndex].get_refcount() == invalid )
			return;			// already removed
		m_vFreeIndices.push_back(nIndex);
		m_vData[nIndex].set_refcount(invalid);
		m_nUsedCount--;
	}
		

	inline void clear( bool bFreeMem = false ) {
		if ( bFreeMem )
			m_vData.clear();
		else
			m_vData.resize(0);
		m_vFreeIndices.clear();
		m_nUsedCount = 0;
	}


	// return true if no elements are in use
	inline bool empty() const { return m_nUsedCount == 0; }

	// return count of elements with non-zero refcount
	inline unsigned int count() const { return m_nUsedCount; }

	// return highest index that has been used (which may currently be on free list)
	inline unsigned int max_index() const { return (unsigned int)m_vData.size(); }

	// return whether or not we have any free items
	inline bool is_dense() const { return m_vFreeIndices.empty(); }

	inline Type & operator[]( int nIndex ) {
		gDevAssert( nIndex < (int)m_vData.size() && m_vData[nIndex].get_refcount() > 0 && m_vData[nIndex].get_refcount() != invalid  );
		return m_vData[nIndex];
	}
	inline const Type & operator[]( int nIndex ) const {
		gDevAssert( nIndex < (int)m_vData.size() && m_vData[nIndex].get_refcount() > 0 && m_vData[nIndex].get_refcount() != invalid  );
		return m_vData[nIndex];
	}


	std::string toString(bool bPrintValues = true) const {
		std::ostrstream out;
		out << "[ Size: " <<  (int)m_vData.size() << " Used: " << m_nUsedCount << " Free: " << m_vFreeIndices.size() << " ]" << std::endl;

		for ( unsigned int i = 0; i < m_vData.size(); ++i ) {
			if ( m_vData[i].get_refcount() == invalid)
				out << "  [" << i << ",X]=";
			else
				out << "  [" << i << "," <<  m_vData[i].get_refcount() << "]=";
			if ( bPrintValues )
				out << m_vData[i].printable();
		}
		out << std::endl;
		for (unsigned int i = 0; i < m_vFreeIndices.size(); ++i) {
			out << "  [" << i << "," << m_vFreeIndices[i] << "]";
		}
		out << std::endl;
		out << '\0';
		return out.str();
	}



protected:
	VectorType< Type > m_vData;

	unsigned int m_nUsedCount;
	dvector<unsigned int> m_vFreeIndices;


public:


	// iterator for elements with non-zero refcount
	class base_iterator
	{
	public:
		inline base_iterator() { p = nullptr; m_nIndex = 0; m_nLast = 0;}

		inline bool operator==( const base_iterator & r2 ) const {
			return m_nIndex == r2.m_nIndex;
		}
		inline bool operator!=( const base_iterator & r2 ) const {
			return m_nIndex != r2.m_nIndex;
		}

	protected:
		inline void goto_next() {
			if ( m_nIndex != m_nLast ) 
				m_nIndex++;
			while ( m_nIndex != m_nLast && p->isValidUnsafe(m_nIndex) == false ) 
				m_nIndex++;
		}

		inline base_iterator( rcvector * pVector, int nIndex, int nLast )
		{
			p = pVector;
			m_nIndex = nIndex;
			m_nLast = nLast;
			if ( m_nIndex != m_nLast && p->isValidUnsafe(m_nIndex) == false )
				goto_next();		// initialize
		}
		rcvector * p;
		int m_nIndex;
		int m_nLast;
		friend class rcvector;
	};


	// iterator over item values at valid indices (ie non-zero refcount)
	// also provides index() function to retrieve current index
	class item_iterator : public base_iterator
	{
	public:
        inline int index() const {
            return m_nIndex;
        }
		inline Type & operator*() { 
			return (*p)[m_nIndex];
		}

		inline item_iterator & operator++() {		// prefix
			goto_next();
			return *this;
		}
		inline item_iterator operator++(int) {		// postfix
			item_iterator copy(*this);
			goto_next();
			return copy;
		}

	protected:
		inline item_iterator( rcvector * pVector, int nIndex, int nLast ) : base_iterator(pVector, nIndex, nLast)
		{}
		friend class rcvector;
	};

	inline item_iterator begin_items() {
		return item_iterator( this, (int)0, (int)m_vData.size() );
	}
	inline item_iterator end_items() {
		return item_iterator( this, (int)m_vData.size(), (int)m_vData.size() );
	}

	inline item_iterator begin() { return begin_items(); }
	inline item_iterator end() { return end_items(); }




	// iterator over valid indices (ie non-zero refcount)
	class index_iterator : public base_iterator
	{
	public:
		inline int operator*() const { 
			return m_nIndex;
		}

		inline index_iterator & operator++() {		// prefix
			goto_next();
			return *this;
		}
		inline index_iterator operator++(int) {		// postfix
			index_iterator copy(*this);
			goto_next();
			return copy;
		}

	protected:
		inline index_iterator( rcvector * pVector, int nIndex, int nLast ) : base_iterator(pVector, nIndex, nLast)
		{}
		friend class rcvector;
	};

	inline index_iterator begin_indices() {
		return index_iterator( this, (int)0, (int)m_vData.size() );
	}
	inline index_iterator end_indices() {
		return index_iterator( this, (int)m_vData.size(), (int)m_vData.size() );
	}


	class index_wrapper
	{
	public:
		typename rcvector * pVector;
		index_wrapper( typename rcvector * p ) { pVector = p; }
		typename rcvector::index_iterator begin() { return pVector->begin_indices(); }
		typename rcvector::index_iterator end() { return pVector->end_indices(); }
	};

};




template<typename Type, template<typename> class VectorType >
rcvector<Type,VectorType>::~rcvector() 
{
}




}  // end namespace g3
