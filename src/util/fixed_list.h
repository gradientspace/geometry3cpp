#pragma once

#include <object_pool.h>

namespace g3
{

extern void g3_testAssert(bool);
    
// fixed_index_list 

template<unsigned int N = 8>
class fixed_index_list
{
public:
	struct node {
		int index;
		node * pNext = nullptr;
	};
	typedef object_allocator<node> node_allocator;

	fixed_index_list();
	~fixed_index_list();

	// cannot do copy constructors as we do not allocate our own nodes
	fixed_index_list(const fixed_index_list<N> & copy ) = delete;
	const fixed_index_list<N> & operator=(const fixed_index_list<N> & copy ) = delete;

	// can do move semantics because pointers will not be invalidated
	fixed_index_list(const fixed_index_list<N> && moved );
	const fixed_index_list<N> & operator=(const fixed_index_list<N> && moved );

	void copy( const fixed_index_list & other, node_allocator * alloc );

	inline int num_fixed() const;
    inline size_t size() const;

    inline bool empty() const;
	inline bool is_full() const;

	inline void add(int index, node_allocator * alloc);
    inline bool remove(int index, node_allocator * alloc);
	inline void clear(node_allocator * alloc);

	// TODO SORT / SORTED INSERT ?

    inline bool contains(int index) const;
    inline int count(int index) const;
    
	template<typename Func>
	inline int find( const Func & f, int fail = -1 ) const;

	inline void get( std::vector<int> & v ) const;
    inline int front() const;

	//
	// index iterator
	//
	class iterator {
	public:
		inline int operator*() const;
		inline int & operator*();
		inline void next();
		inline bool is_index() const;
		inline int index() const;
		inline iterator & operator++();			// prefix
		inline iterator operator++(int);	// postfix
		inline bool operator==(const iterator & i2);
		inline bool operator!=(const iterator & i2);
	protected:
		fixed_index_list * pList = nullptr;
		int i = 0;
		node * pCur = nullptr;
		inline iterator(fixed_index_list * p, int iCur);
		inline iterator(fixed_index_list * p, node * pCur);
		friend class fixed_index_list<N>;
	};
    typedef const iterator const_iterator;

	iterator begin();
	iterator end();
    const_iterator begin() const;
    const_iterator end() const;
    
protected:
	int m_data[N];
	node * pFirst = nullptr;

	static constexpr int invalid_index = -9999;
	static node * invalid_node();

	friend class iterator;
    
public:
    void assert_valid(bool bCheckDuplicates = true);
};


// stream-print operator
template<unsigned int N>
std::ostream& operator<<( std::ostream& os, fixed_index_list<N> & l )
{
	int pi = N;
	for (auto cur = l.begin(); cur != l.end(); cur++) {
		if (cur.is_index())
			os << cur.index() << "=" << *cur << " ";
		else
			os << "p" << (pi++) << "=" << *cur << " ";
	}
	return os;
}



template<unsigned int N>
typename fixed_index_list<N>::node * fixed_index_list<N>::invalid_node()
{
	return (node *)1;
}



template<unsigned int N>
fixed_index_list<N>::fixed_index_list()
{
	for ( int k = 0; k < N; ++k )
		m_data[k] = invalid_index;
	pFirst = nullptr;
}

template<unsigned int N>
fixed_index_list<N>::~fixed_index_list()
{
}

template<unsigned int N>
fixed_index_list<N>::fixed_index_list( const fixed_index_list<N> && moved )
{
	*this = std::move(moved);
}

template<unsigned int N>
const fixed_index_list<N> & fixed_index_list<N>::operator=( const fixed_index_list<N> && moved )
{
	m_data = std::move(moved.m_data);
	pFirst = moved.pFirst;
	return *this;
}

template<unsigned int N>
void fixed_index_list<N>::copy( const fixed_index_list<N> & other, node_allocator * alloc )
{
	for ( int k = 0; k < N; ++k )
		m_data[k] = other.m_data[k];
	node * pOtherCur = other.pFirst;
	node * pCur = nullptr;
	while (pOtherCur != nullptr) {
		node * pNew = alloc->allocate();
		pNew->index = pOtherCur->index;
		if (pCur == nullptr) {
			pFirst = pCur = pNew;
		} else {
			pCur->pNext = pNew;
			pCur = pNew;
		}
		pOtherCur = pOtherCur->pNext;
	}
}


template<unsigned int N>
int fixed_index_list<N>::num_fixed() const
{
	return N;
}
    
template<unsigned int N>
size_t fixed_index_list<N>::size() const
{
    if ( pFirst != nullptr ) {
        size_t s = N;
        auto pCur = pFirst;
        while ( pCur != nullptr ) {
            s++;
            pCur = pCur->pNext;
        }
        return s;
    } else {
        for ( int k = N-1; k >= 0; k-- ) {
            if ( m_data[k] != invalid_index )
                return k+1;
        }
        return 0;
    }
}

template<unsigned int N>
bool fixed_index_list<N>::is_full() const
{
    return m_data[N-1] != invalid_index;
}

template<unsigned int N>
bool fixed_index_list<N>::empty() const
{
    return m_data[0] == invalid_index;
}


template<unsigned int N>
void fixed_index_list<N>::add( int index, node_allocator * alloc )
{
	if (pFirst == nullptr) {
		// [TODO] would be smarter to do binary search here...
		for (int k = 0; k < N; ++k) {
			if (m_data[k] == invalid_index) {
				m_data[k] = index;
				return;
			}
		}
	}
	node * pNode = alloc->allocate();
	pNode->index = index;
	if (pFirst == nullptr) {
		pFirst = pNode;
	} else {
		pNode->pNext = pFirst;
		pFirst = pNode;
	}
}

template<unsigned int N>
bool fixed_index_list<N>::remove( int index, node_allocator * alloc )
{
    // search index list
    for ( int k = 0; k < N; ++k ) {
        if ( m_data[k] == invalid_index )   // hit end of list
            return false;
        if ( m_data[k] == index ) {
            // shift fixed list to the left
            while ( k < N-1 && m_data[k+1] != invalid_index ) {
                m_data[k] = m_data[k+1];
                k++;
            }
            // shift overflow list if we have it
            if ( k == N-1 ) {
                if ( pFirst != nullptr ) {
                    m_data[N-1] = pFirst->index;
                    node * pTmp = pFirst;
                    pFirst = pFirst->pNext;
                    alloc->release(pTmp);
                } else {
                    m_data[N-1] = invalid_index;
                }
            } else
                m_data[k] = invalid_index;
            return true;
        }
    }
    
    // search overflow list
    node * pPrev = nullptr;
    node * pCur = pFirst;
    while ( pCur != nullptr ) {
        if ( pCur->index == index ) {
            node * pTmp = pCur;
            if ( pPrev == nullptr )
                pFirst = pCur->pNext;
            else
                pPrev->pNext = pCur->pNext;
            alloc->release(pTmp);
            return true;
        }
        pPrev = pCur;
        pCur = pCur->pNext;
    }
    
    // did not find
    return false;
}
    

template<unsigned int N>
void fixed_index_list<N>::clear( node_allocator * alloc )
{
	node * pCur = pFirst;
	while ( pCur != nullptr ) {
		node * pTmp = pCur;
		pCur = pCur->pNext;
		alloc->release(pTmp);
	}
	pFirst = nullptr;
	for ( int k = 0; k < N; ++k ) 
		m_data[k] = invalid_index;
}

template<unsigned int N>
int fixed_index_list<N>::iterator::operator*() const 
{
	return (pCur == nullptr) ? pList->m_data[i] : pCur->index;
}

template<unsigned int N>
int & fixed_index_list<N>::iterator::operator*()
{
	return (pCur == nullptr) ? pList->m_data[i] : pCur->index;
}

template<unsigned int N>
bool fixed_index_list<N>::iterator::is_index() const 
{
	return (pCur == nullptr);
}

template<unsigned int N>
int fixed_index_list<N>::iterator::index() const 
{
	return i;
}


template<unsigned int N>
typename fixed_index_list<N>::iterator & fixed_index_list<N>::iterator::operator++()
{
	next();
	return *this;
}

template<unsigned int N>
typename fixed_index_list<N>::iterator fixed_index_list<N>::iterator::operator++(int i)
{
	iterator copy(*this);
	next();
	return copy;
}

template<unsigned int N>
bool fixed_index_list<N>::iterator::operator==(const iterator & itr)
{
	return pList == itr.pList && i == itr.i && pCur == itr.pCur;
}
template<unsigned int N>
bool fixed_index_list<N>::iterator::operator!=(const iterator & itr)
{
	return pList != itr.pList || i != itr.i || pCur != itr.pCur;
}

template<unsigned int N>
void fixed_index_list<N>::iterator::next()
{
	if (pCur != nullptr) {
		if (pCur->pNext == nullptr)
			pCur = pList->invalid_node();
		else if ( pCur != pList->invalid_node() )
			pCur = pCur->pNext;
		return;
	} 

	if (i < N) 
		i++;
	if (i < N) {
		if (pList->m_data[i] == invalid_index) {
			pCur = pList->invalid_node();
			i = N;
			return;
		}
	} else {
		if (pList->pFirst == nullptr)
			pCur = pList->invalid_node();
		else
			pCur = pList->pFirst;
	}
}


template<unsigned int N>
fixed_index_list<N>::iterator::iterator( fixed_index_list * p, int iCur )
{
	pList = p;
	i = iCur;
}
template<unsigned int N>
fixed_index_list<N>::iterator::iterator( fixed_index_list * p, node * pCur )
{
	pList = p;
	this->pCur = pCur;
	i = N;
}

template<unsigned int N>
typename fixed_index_list<N>::iterator fixed_index_list<N>::begin()
{
	return (m_data[0] == invalid_index) ? end() : iterator(this, 0);
}
template<unsigned int N>
typename fixed_index_list<N>::iterator fixed_index_list<N>::end()
{
	return iterator(this, invalid_node());
}
    
template<unsigned int N>
const typename fixed_index_list<N>::iterator fixed_index_list<N>::begin() const
{
    return (m_data[0] == invalid_index) ? end()
    : iterator( const_cast<fixed_index_list<N> *>(this), 0);
}
template<unsigned int N>
const typename fixed_index_list<N>::iterator fixed_index_list<N>::end() const
{
    return iterator(const_cast<fixed_index_list<N> *>(this), invalid_node());
}

template<unsigned int N>
bool fixed_index_list<N>::contains(int index) const
{
    for (int k = 0; k < N; ++k) {
        if ( m_data[k] == invalid_index )
            return false;
        if ( m_data[k] == index )
            return true;
    }
    const node * pCur = pFirst;
    while (pCur != nullptr) {
        if (pCur->index == index)
            return true;
        pCur = pCur->pNext;
    }
    return false;
}

    
template<unsigned int N>
int fixed_index_list<N>::count(int index) const
{
    int n = 0;
    for (int k = 0; k < N; ++k) {
        if ( m_data[k] == invalid_index )
            return n;
        if ( m_data[k] == index )
            n++;
    }
    const node * pCur = pFirst;
    while (pCur != nullptr) {
        if (pCur->index == index)
            n++;
        pCur = pCur->pNext;
    }
    return n;
}
    

template<unsigned int N>
template<typename Func>
int fixed_index_list<N>::find( const Func & f, int fail ) const
{
	for (int k = 0; k < N; ++k) {
		if ( m_data[k] == invalid_index )
			return fail;
		if ( f(m_data[k]) )
			return m_data[k];
	}
	const node * pCur = pFirst;
	while (pCur != nullptr) {
		if ( f(pCur->index) )
			return pCur->index;
        pCur = pCur->pNext;
	}
	return fail;
}

template<unsigned int N>
int fixed_index_list<N>::front() const
{
    return m_data[0];
}

template<unsigned int N>
void fixed_index_list<N>::get( std::vector<int> & v ) const
{
	for (int k = 0; k < N; ++k) {
		if ( m_data[k] == invalid_index )
			return;
		v.push_back(m_data[k]);
	}
	const node * pCur = pFirst;
	while (pCur != nullptr) {
		v.push_back(pCur->index);
		pCur = pCur->pNext;
	}
}



template<unsigned int N>
void fixed_index_list<N>::assert_valid(bool bCheckDuplicates)
{
    for ( int k = 0; k < N; ++k ) {
        if ( m_data[k] == invalid_index ) {
            g3_testAssert( pFirst == nullptr );
            for ( int j = k+1; j < N; ++j )
                g3_testAssert( m_data[j] == invalid_index );
        }
    }
    if ( bCheckDuplicates ) {
        for ( int index : *this )
            g3_testAssert( count(index) == 1 );
    }
}
    
    
} // end namespace g3

