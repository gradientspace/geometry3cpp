#pragma once

#include <vector>
#include <g3Debug.h>

namespace g3
{


template<class Type>
class sparse_dvector_segment {
public:
	Type * pData = nullptr;
	size_t nSize;
	size_t nCur;
	sparse_dvector_segment() = default;
	~sparse_dvector_segment() { }
};


template<class Type>
class sparse_dvector
{
public:
	sparse_dvector(const Type & defaultValue, unsigned int nSegmentSize = 0);
	sparse_dvector(const sparse_dvector & copy);
	sparse_dvector(sparse_dvector && moved);
	virtual ~sparse_dvector();

	const sparse_dvector & operator=( const sparse_dvector & copy );
	const sparse_dvector & operator=( sparse_dvector && moved );

	inline void clear( bool bFreeSegments = false );
	inline void resize( size_t nCount );
	inline void resize( size_t nCount, const Type & init_value );

	inline bool empty() const;
	inline size_t size() const;
	inline size_t allocated() const;

	inline void push_back( const Type & data );
	inline Type * push_back();
	inline void push_back( const sparse_dvector<Type> & data );
	inline void pop_back();

	inline Type & front();
	inline const Type & front() const;
	inline Type & back();
	inline const Type & back() const;

	inline Type & operator[]( unsigned int nIndex );
	inline const Type & operator[]( unsigned int nIndex ) const;

	// apply f() to each member sequentially
	template<typename Func>
	void apply(const Func & f);

	class iterator {
	public:
		inline const Type & operator*() const;
		inline Type & operator*();
		inline int index();
		inline void next();
		inline iterator & operator++();			// prefix
		inline iterator operator++(int);	// postfix
		inline bool operator==(const iterator & i2);
		inline bool operator!=(const iterator & i2);
	protected:
		sparse_dvector<Type> * pVector;
		int i;
		inline iterator(sparse_dvector<Type> * p, int iCur);
		friend class sparse_dvector<Type>;
	};

	iterator begin();
	iterator end();

protected:
	unsigned int m_nSegmentSize;
	unsigned int m_nCurSeg;
	size_t m_nAllocated;

	Type m_defaultValue;

	std::vector< sparse_dvector_segment<Type> > m_vSegments;

	Type * allocate_element();
	sparse_dvector_segment<Type> & get_or_allocate(unsigned int nSegIndex);
	bool is_allocated(unsigned int nSegIndex);

	friend class iterator;

	// parallel friend functions defined in sparse_dvector_util.h
	template<typename TypeX, typename Func>
	friend void parallel_apply( sparse_dvector<TypeX> & v, const Func & f);
};



template<class Type>
std::ostream& operator<<( std::ostream& os, sparse_dvector<Type> & dv )
{
	auto cur(dv.begin()), end(dv.end());
	while ( cur != end ) {
		os << cur.index() << "=" << *cur << " ";
		cur++;
	}
	return os;
}




template <class Type>
sparse_dvector<Type>::sparse_dvector(const Type & defaultValue, unsigned int nSegmentSize)
	: m_defaultValue(defaultValue)
{
	if ( nSegmentSize == 0 )
		m_nSegmentSize = 256;
	else
		m_nSegmentSize = nSegmentSize;

	m_vSegments.resize(1);
	m_vSegments[0].pData = nullptr;
	m_vSegments[0].nSize = m_nSegmentSize;
	m_vSegments[0].nCur = 0;

	m_nCurSeg = 0;
	m_nAllocated = 0;
}

template <class Type>
sparse_dvector<Type>::sparse_dvector(const sparse_dvector<Type> & copy) : sparse_dvector(0)
{
	*this = copy;
}

template <class Type>
sparse_dvector<Type>::sparse_dvector( sparse_dvector && moved )
{
	*this = std::move(moved);
}


template <class Type>
sparse_dvector<Type>::~sparse_dvector()
{
	size_t nCount = m_vSegments.size();
	for (unsigned int i = 0; i < nCount; ++i) {
		if (m_vSegments[i].pData != nullptr) {
			delete [] m_vSegments[i].pData;
			m_vSegments[i].pData = nullptr;
			m_nAllocated -= m_nSegmentSize;
		}
	}
}

template <class Type>
const sparse_dvector<Type> & sparse_dvector<Type>::operator=( const sparse_dvector & copy )
{
	m_defaultValue = copy.m_defaultValue;

	// if segments are the same size, we don't need to re-allocate any existing ones  (woot!)
	bool bSegmentsChanged = (m_nSegmentSize != copy.m_nSegmentSize);
	m_nSegmentSize = copy.m_nSegmentSize;
	if ( bSegmentsChanged )
		clear(true);

	m_nCurSeg = copy.m_nCurSeg;

	// allocate memory (or discard exisiting memory) for segments
	resize( copy.size() );

	// copy segment contents
	size_t nSegs = copy.m_vSegments.size();
	for ( unsigned int k = 0; k < nSegs; ++k ) {
		if (copy.m_vSegments[k].pData != nullptr) {
			get_or_allocate(k);
#ifdef WIN32	
			memcpy_s( m_vSegments[k].pData, m_nSegmentSize*sizeof( Type ), copy.m_vSegments[k].pData, m_nSegmentSize*sizeof( Type ) );
#else
			memcpy( m_vSegments[k].pData, copy.m_vSegments[k].pData, m_nSegmentSize*sizeof( Type ) );
#endif
			m_vSegments[k].nSize = m_nSegmentSize;
			m_vSegments[k].nCur = copy.m_vSegments[k].nCur;

		} else if ( is_allocated( k ) ) {
			delete [] m_vSegments[k].pData;
			m_vSegments[k].pData = nullptr;
		}
	}

	return *this;
}



template <class Type>
const sparse_dvector<Type> & sparse_dvector<Type>::operator=( sparse_dvector && moved)
{
	clear(this);
	this->m_nSegmentSize = moved.m_nSegmentSize;
	this->m_nCurSeg = moved.m_nCurSeg;
	m_vSegments = std::move(moved.m_vSegments);
	moved.m_nCurSeg = 0;
	return *this;
}


template <class Type>
void sparse_dvector<Type>::clear( bool bFreeSegments )
{
	size_t nCount = m_vSegments.size();
	for (unsigned int i = 0; i < nCount; ++i) {
		m_vSegments[i].nCur = 0;
		if (m_vSegments[i].pData != nullptr)
			for (unsigned int k = 0; k < m_nSegmentSize; ++k)
				m_vSegments[i].pData[k] = m_defaultValue;
	}

	if (bFreeSegments) {
		for (unsigned int i = 0; i < nCount; ++i) {
			if (m_vSegments[i].pData != nullptr) {
				delete [] m_vSegments[i].pData;
				m_vSegments[i].pData = nullptr;
				m_nAllocated -= m_nSegmentSize;
			}
		}

		m_vSegments.resize(1);
		m_vSegments[0].pData = nullptr;
		m_vSegments[0].nSize = m_nSegmentSize;
		m_vSegments[0].nCur = 0;
	}

	m_nCurSeg = 0;
}


template <class Type>
void sparse_dvector<Type>::resize( size_t nCount )
{
	// figure out how many segments we need
	unsigned int nNumSegs = 1 + (unsigned int)nCount / m_nSegmentSize;

	// figure out how many are currently allocated...
	size_t nCurCount = m_vSegments.size();

	// erase extra segments memory
	for (unsigned int i = nNumSegs; i < nCurCount; ++i) {
		if (m_vSegments[i].pData != nullptr) {
			delete [] m_vSegments[i].pData;
			m_vSegments[i].pData = nullptr;
			m_nAllocated -= m_nSegmentSize;
		}
	}

	// resize to right number of segments
	m_vSegments.resize(nNumSegs);

	// allocate new segments
	for (unsigned int i = (unsigned int)nCurCount; i < nNumSegs; ++i) {
		m_vSegments[i].pData = nullptr;
		m_vSegments[i].nSize = m_nSegmentSize;
		m_vSegments[i].nCur = 0;
	}

	// mark full segments as used
	for (unsigned int i = 0; i < nNumSegs-1; ++i)
		m_vSegments[i].nCur = m_nSegmentSize;

	// mark last segment
	m_vSegments[nNumSegs-1].nCur = nCount - (nNumSegs-1)*m_nSegmentSize;

	m_nCurSeg = nNumSegs-1;
}

template <class Type>
sparse_dvector_segment<Type> & sparse_dvector<Type>::get_or_allocate( unsigned int nSegIndex )
{
	gDevAssert( nSegIndex <= m_nCurSeg );
	sparse_dvector_segment<Type> & seg = m_vSegments[nSegIndex];
	if (seg.pData == nullptr) {
		seg.pData = new Type[m_nSegmentSize];
		for ( unsigned int k = 0; k < m_nSegmentSize; ++k )
			seg.pData[k] = m_defaultValue;
		m_nAllocated += m_nSegmentSize;
	}
	return seg;
}

template <class Type>
bool sparse_dvector<Type>::is_allocated( unsigned int nSegIndex )
{
	if ( nSegIndex <= m_nCurSeg && m_vSegments[nSegIndex].pData != nullptr )
		return true;
	return false;
}


template <class Type>
void sparse_dvector<Type>::resize( size_t nCount, const Type & init_value )
{
	size_t nCurSize = size();
	resize(nCount);

	// [TODO] could be more efficient!
	if (init_value != m_defaultValue) {
		for (size_t nIndex = nCurSize; nIndex < nCount; ++nIndex)
			get_or_allocate(nIndex / m_nSegmentSize).pData[nIndex % m_nSegmentSize] = init_value;
	}
}


template <class Type>
bool sparse_dvector<Type>::empty() const
{
	return ! ( m_nCurSeg > 0 || m_vSegments[0].nCur > 0 );
}

template <class Type>
size_t  sparse_dvector<Type>::size() const
{
	return m_nCurSeg*m_nSegmentSize + m_vSegments[m_nCurSeg].nCur;
}

template <class Type>
size_t  sparse_dvector<Type>::allocated() const
{
	return m_nAllocated;
}


template <class Type>
Type * sparse_dvector<Type>::allocate_element()
{
	sparse_dvector_segment<Type> & seg = m_vSegments[m_nCurSeg];
	if ( seg.nCur == seg.nSize ) {
		if ( m_nCurSeg == m_vSegments.size() - 1 ) {
			m_vSegments.resize( m_vSegments.size() + 1 );
			sparse_dvector_segment<Type> & newSeg = m_vSegments.back();
			newSeg.pData = nullptr;
			newSeg.nSize = m_nSegmentSize;
			newSeg.nCur = 0;
		}
		m_nCurSeg++;
	}

	sparse_dvector_segment<Type> & returnSeg = get_or_allocate(m_nCurSeg);
	return & returnSeg.pData[ returnSeg.nCur++ ];
}

template <class Type>
void sparse_dvector<Type>::push_back( const Type & data )
{
	Type * pNewElem = allocate_element();
	*pNewElem = data;
}

template <class Type>
Type * sparse_dvector<Type>::push_back()
{
	return allocate_element();
}

template <class Type>
void sparse_dvector<Type>::push_back( const sparse_dvector<Type> & data )
{
	// [RMS TODO] it would be a lot more efficient to use memcopies here...
	size_t nSize = data.size();
	for ( unsigned int k = 0; k < nSize; ++k )
		push_back( data[k] );
}


template <class Type>
void sparse_dvector<Type>::pop_back()
{
	if (m_vSegments[m_nCurSeg].nCur > 0) {
		m_vSegments[m_nCurSeg].nCur--;
	} else if (m_nCurSeg > 0) {
		m_nCurSeg--;
		m_vSegments[m_nCurSeg].nCur--;
	}
	// else we are in seg 0 and nCur is 0, so we have nothing to pop!
}


template <class Type>
Type & sparse_dvector<Type>::front()
{
	return get_or_allocate(0).pData[0];
}
template <class Type>
const Type & sparse_dvector<Type>::front() const
{
	return is_allocated(0) ? m_vSegments[0].pData[0] : m_defaultValue;
}

template <class Type>
Type & sparse_dvector<Type>::back()
{
	auto & seg = get_or_allocate(m_nCurSeg);
	return seg.pData[ seg.nCur-1 ];
}
template <class Type>
const Type & sparse_dvector<Type>::back() const
{
	if (is_allocated( m_nCurSeg )) {
		auto & seg = m_vSegments[m_nCurSeg];
		return seg.pData[seg.nCur-1];
	} else 
		return m_defaultValue;
}



template <class Type>
Type & sparse_dvector<Type>::operator[]( unsigned int nIndex )
{
	return get_or_allocate( nIndex / m_nSegmentSize ).pData[ nIndex % m_nSegmentSize ];
}

template <class Type>
const Type & sparse_dvector<Type>::operator[]( unsigned int nIndex ) const
{
	auto nSeg = nIndex / m_nSegmentSize;
	return is_allocated(nSeg) ? m_vSegments[nSeg].pData[ nIndex % m_nSegmentSize ] : m_defaultValue;
}



template<typename Type>
template<typename Func>
void sparse_dvector<Type>::apply( const Func & f )
{
	for (auto & seg : m_vSegments) {
		if (seg.pData != nullptr) {
			for (unsigned int i = 0; i < seg.nCur; ++i)
				f( seg.pData[i] );
		}
	}
}






template<typename Type>
const Type & sparse_dvector<Type>::iterator::operator*() const 
{
	return (*pVector)[i];
}

template<typename Type>
Type & sparse_dvector<Type>::iterator::operator*()
{
	return (*pVector)[i];
}

template<typename Type>
typename sparse_dvector<Type>::iterator & sparse_dvector<Type>::iterator::operator++()
{
	next();
	return *this;
}

template<typename Type>
typename sparse_dvector<Type>::iterator sparse_dvector<Type>::iterator::operator++(int i)
{
	iterator copy(*this);
	next();
	return copy;
}

template<typename Type>
int sparse_dvector<Type>::iterator::index()
{
	return i;
}

template<typename Type>
void sparse_dvector<Type>::iterator::next()
{
	int n = (int)pVector->size();
	if (i >= n) {
		i = n;
		return;
	}

	i++;
	if ( i == n )
		return;		// done!

	while ((*pVector)[i] == pVector->m_defaultValue) {
		unsigned int nSegment = i / pVector->m_nSegmentSize;
		if ( ! pVector->is_allocated( nSegment ) )
			i = (nSegment+1) * pVector->m_nSegmentSize;
		else
			i++;
		if (i >= n) {
			i = n;
			return;
		}
	}
}

template<typename Type>
bool sparse_dvector<Type>::iterator::operator==(const iterator & itr)
{
	return pVector == itr.pVector && i == itr.i;
}
template<typename Type>
bool sparse_dvector<Type>::iterator::operator!=(const iterator & itr)
{
	return pVector != itr.pVector || i != itr.i;
}


template<typename Type>
sparse_dvector<Type>::iterator::iterator( sparse_dvector<Type> * p, int iCur )
{
	pVector = p;
	i = iCur;
}

template<typename Type>
typename sparse_dvector<Type>::iterator sparse_dvector<Type>::begin()
{
	if ( empty() )
		return end();
	iterator itr(this,0);
	if ( *itr == m_defaultValue )
		itr.next();
	return itr;
}
template<typename Type>
typename sparse_dvector<Type>::iterator sparse_dvector<Type>::end()
{
	return iterator(this, (int)size());
}




} // end namespace g3


