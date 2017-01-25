#pragma once

#include <vector>

namespace g3
{


template<class Type>
class dvector_segment {
public:
	Type * pData;
	size_t nSize;
	size_t nCur;
	dvector_segment() { pData = NULL; }
	~dvector_segment() { }
};


template<class Type>
class dvector
{
public:
	dvector(unsigned int nSegmentSize = 0);
	dvector(const dvector & copy);
	dvector(dvector && moved);
	virtual ~dvector();

	const dvector & operator=( const dvector & copy );
	const dvector & operator=( dvector && moved );

	inline void clear( bool bFreeSegments = false );
	inline void resize( size_t nCount );
	inline void resize( size_t nCount, const Type & init_value );

	inline bool empty() const;
	inline size_t size() const;

	inline void push_back( const Type & data );
	inline Type * push_back();
	inline void push_back( const dvector<Type> & data );
	inline void pop_back();

	inline void insert( unsigned int nIndex, const Type & data );

	inline Type & front();
	inline const Type & front() const;
	inline Type & back();
	inline const Type & back() const;

	inline Type & operator[]( unsigned int nIndex );
	inline const Type & operator[]( unsigned int nIndex ) const;

	// apply f() to each member sequentially
	template<typename Func>
	void apply(const Func & f);

	// TODO insert function that handles resizing?

	class iterator {
	public:
		inline const Type & operator*() const;
		inline Type & operator*();
		inline iterator & operator++();			// prefix
		inline iterator operator++(int);	// postfix
		inline bool operator==(const iterator & i2);
		inline bool operator!=(const iterator & i2);
	protected:
		dvector<Type> * pVector;
		int i;
		inline iterator(dvector<Type> * p, int iCur);
		friend class dvector<Type>;
	};

	iterator begin();
	iterator end();

protected:
	unsigned int m_nSegmentSize;
	unsigned int m_nCurSeg;

	std::vector< dvector_segment<Type> > m_vSegments;

	Type * allocate_element();

	friend class iterator;

	// parallel friend functions defined in dvector_util.h
	template<typename TypeX, typename Func>
	friend void parallel_apply( dvector<TypeX> & v, const Func & f);
};


// stream-print operator
template<class Type>
std::ostream& operator<<( std::ostream& os, const dvector<Type> & dv )
{
	for ( unsigned int i = 0; i < dv.size(); ++i)
		os << i << "=" << dv[i] << " ";
	return os;
}



template <class Type>
dvector<Type>::dvector(unsigned int nSegmentSize)
{
	if ( nSegmentSize == 0 )
		m_nSegmentSize = (1 << 16) / sizeof(Type);		// 64k
	else
		m_nSegmentSize = nSegmentSize;

	m_vSegments.resize(1);
	m_vSegments[0].pData = new Type[ m_nSegmentSize ];
	m_vSegments[0].nSize = m_nSegmentSize;
	m_vSegments[0].nCur = 0;

	m_nCurSeg = 0;
}

template <class Type>
dvector<Type>::dvector(const dvector<Type> & copy) : dvector(0)
{
	*this = copy;
}

template <class Type>
dvector<Type>::dvector( dvector && moved )
{
	*this = std::move(moved);
}


template <class Type>
dvector<Type>::~dvector()
{
	size_t nCount = m_vSegments.size();
	for (unsigned int i = 0; i < nCount; ++i)
		delete [] m_vSegments[i].pData;
}

template <class Type>
const dvector<Type> & dvector<Type>::operator=( const dvector & copy )
{
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
#ifdef WIN32	
		memcpy_s( m_vSegments[k].pData, m_nSegmentSize*sizeof(Type), copy.m_vSegments[k].pData, m_nSegmentSize*sizeof(Type) );
#else
		memcpy( m_vSegments[k].pData, copy.m_vSegments[k].pData, m_nSegmentSize*sizeof(Type) );
#endif
		m_vSegments[k].nSize = m_nSegmentSize;
		m_vSegments[k].nCur = copy.m_vSegments[k].nCur;
	}
	return *this;
}



template <class Type>
const dvector<Type> & dvector<Type>::operator=( dvector && moved)
{
	clear(this);
	this->m_nSegmentSize = moved.m_nSegmentSize;
	this->m_nCurSeg = moved.m_nCurSeg;
	m_vSegments = std::move(moved.m_vSegments);
	moved.m_nCurSeg = 0;
	return *this;
}


template <class Type>
void dvector<Type>::clear( bool bFreeSegments )
{
	size_t nCount = m_vSegments.size();
	for (unsigned int i = 0; i < nCount; ++i) 
		m_vSegments[i].nCur = 0;

	if (bFreeSegments) {
		for (unsigned int i = 0; i < nCount; ++i)
			delete [] m_vSegments[i].pData;

		m_vSegments.resize(1);
		m_vSegments[0].pData = new Type[ m_nSegmentSize ];
		m_vSegments[0].nSize = m_nSegmentSize;
		m_vSegments[0].nCur = 0;
	}

	m_nCurSeg = 0;
}


template <class Type>
void dvector<Type>::resize( size_t nCount )
{
	// figure out how many segments we need
	unsigned int nNumSegs = 1 + (unsigned int)nCount / m_nSegmentSize;

	// figure out how many are currently allocated...
	size_t nCurCount = m_vSegments.size();

	// erase extra segments memory
	for ( unsigned int i = nNumSegs; i < nCurCount; ++i )
		delete [] m_vSegments[i].pData;

	// resize to right number of segments
	m_vSegments.resize(nNumSegs);

	// allocate new segments
	for (unsigned int i = (unsigned int)nCurCount; i < nNumSegs; ++i) {
		m_vSegments[i].pData = new Type[ m_nSegmentSize ];
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
void dvector<Type>::resize( size_t nCount, const Type & init_value )
{
	size_t nCurSize = size();
	resize(nCount);
	for ( size_t nIndex = nCurSize; nIndex < nCount; ++nIndex ) 
		m_vSegments[ nIndex / m_nSegmentSize ].pData[ nIndex % m_nSegmentSize ] = init_value;
}

template <class Type>
bool dvector<Type>::empty() const
{
	return ! ( m_nCurSeg > 0 || m_vSegments[0].nCur > 0 );
}

template <class Type>
size_t dvector<Type>::size() const
{
	return m_nCurSeg*m_nSegmentSize + m_vSegments[m_nCurSeg].nCur;
}


template <class Type>
Type * dvector<Type>::allocate_element()
{
	dvector_segment<Type> & seg = m_vSegments[m_nCurSeg];
	if ( seg.nCur == seg.nSize ) {
		if ( m_nCurSeg == m_vSegments.size() - 1 ) {
			m_vSegments.resize( m_vSegments.size() + 1 );
			dvector_segment<Type> & newSeg = m_vSegments.back();
			newSeg.pData = new Type[ m_nSegmentSize ];
			newSeg.nSize = m_nSegmentSize;
			newSeg.nCur = 0;
		}
		m_nCurSeg++;
	}
	dvector_segment<Type> & returnSeg = m_vSegments[m_nCurSeg];
	return & returnSeg.pData[ returnSeg.nCur++ ];
}

template <class Type>
void dvector<Type>::push_back( const Type & data )
{
	Type * pNewElem = allocate_element();
	*pNewElem = data;
}

template <class Type>
Type * dvector<Type>::push_back()
{
	return allocate_element();
}

template <class Type>
void dvector<Type>::push_back( const dvector<Type> & data )
{
	// [RMS TODO] it would be a lot more efficient to use memcopies here...
	size_t nSize = data.size();
	for ( unsigned int k = 0; k < nSize; ++k )
		push_back( data[k] );
}

template <class Type>
void dvector<Type>::pop_back()
{
	if (m_vSegments[m_nCurSeg].nCur > 0)
		m_vSegments[m_nCurSeg].nCur--;
	if (m_vSegments[m_nCurSeg].nCur == 0 && m_nCurSeg > 0 ) 
		m_nCurSeg--;
	// else we are in seg 0 and nCur is 0, so we have nothing to pop!
}


template <class Type>
void dvector<Type>::insert( unsigned int nIndex, const Type & data )
{
	//TODO this does not make sense if nIndex < size !!!
	abort();

	size_t s = size();
	if (nIndex == s) {
		push_back( data );
	} else {
		resize( nIndex+1 );
		(*this)[nIndex] = data;
	}
}


template <class Type>
Type & dvector<Type>::front()
{
	return m_vSegments[0].pData[0];
}
template <class Type>
const Type & dvector<Type>::front() const
{
	return m_vSegments[0].pData[0];
}

template <class Type>
Type & dvector<Type>::back()
{
	auto & seg = m_vSegments[m_nCurSeg];
	return seg.pData[ seg.nCur-1 ];
}
template <class Type>
const Type & dvector<Type>::back() const
{
	auto & seg = m_vSegments[m_nCurSeg];
	return seg.pData[ seg.nCur-1 ];
}



template <class Type>
Type & dvector<Type>::operator[]( unsigned int nIndex )
{
	return m_vSegments[ nIndex / m_nSegmentSize ].pData[ nIndex % m_nSegmentSize ];
}

template <class Type>
const Type & dvector<Type>::operator[]( unsigned int nIndex ) const
{
	return m_vSegments[ nIndex / m_nSegmentSize ].pData[ nIndex % m_nSegmentSize ];
}



template<typename Type>
template<typename Func>
void dvector<Type>::apply( const Func & f )
{
	for (auto & seg : m_vSegments) {
		for (unsigned int i = 0; i < seg.nCur; ++i) 
			f( seg.pData[i] );
	}
}







template<typename Type>
const Type & dvector<Type>::iterator::operator*() const 
{
	return (*pVector)[i];
}

template<typename Type>
Type & dvector<Type>::iterator::operator*()
{
	return (*pVector)[i];
}

template<typename Type>
typename dvector<Type>::iterator & dvector<Type>::iterator::operator++()
{
	i++;
	return *this;
}

template<typename Type>
typename dvector<Type>::iterator dvector<Type>::iterator::operator++(int i)
{
	iterator copy(*this);
	i++;
	return copy;
}

template<typename Type>
bool dvector<Type>::iterator::operator==(const iterator & itr)
{
	return pVector == itr.pVector && i == itr.i;
}
template<typename Type>
bool dvector<Type>::iterator::operator!=(const iterator & itr)
{
	return pVector != itr.pVector || i != itr.i;
}


template<typename Type>
dvector<Type>::iterator::iterator( dvector<Type> * p, int iCur )
{
	pVector = p;
	i = iCur;
}

template<typename Type>
typename dvector<Type>::iterator dvector<Type>::begin()
{
	return empty() ? end() : iterator(this,0);
}
template<typename Type>
typename dvector<Type>::iterator dvector<Type>::end()
{
	return iterator(this, (int)size());
}




} // end namespace g3


