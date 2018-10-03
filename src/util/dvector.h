#pragma once

#include <vector>
#include <array>

namespace g3
{


template<class Type>
class dvector
{
public:
	dvector();
	dvector(const dvector & copy);
	dvector(dvector && moved);
	virtual ~dvector();

	const dvector & operator=( const dvector & copy );
	const dvector & operator=( dvector && moved );

	inline void clear();
	inline void fill(const Type & value);
	inline void resize( size_t nCount );
	inline void resize( size_t nCount, const Type & init_value );

	inline bool empty() const;
	inline size_t size() const;
	inline size_t length() const;
	inline int block_size() const;
	inline size_t byte_count() const;

	inline void add(const Type & data);
	inline void push_back( const Type & data );
	inline void push_back( const dvector<Type> & data );
	inline void pop_back();

	inline void insertAt(const Type & data, unsigned int nIndex);

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
	// [RMS] nBlockSize must be a power-of-two, so we can use bit-shifts in operator[]
	static constexpr int nBlockSize = 2048;   // (1 << 11)
	static constexpr int nShiftBits = 11;
	static constexpr int nBlockIndexBitmask = 2047;   // low 11 bits

	unsigned int iCurBlock;
	unsigned int iCurBlockUsed;
	
	using BlockType = std::array<Type, nBlockSize>;
	std::vector<BlockType> Blocks;

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
dvector<Type>::dvector()
{
	iCurBlock = 0;
	iCurBlockUsed = 0;
	Blocks.push_back(BlockType());
}

template <class Type>
dvector<Type>::dvector(const dvector<Type> & copy) : dvector()
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
}

template <class Type>
const dvector<Type> & dvector<Type>::operator=( const dvector & copy )
{
	Blocks = copy.Blocks;
	iCurBlock = copy.iCurBlock;
	iCurBlockUsed = copy.iCurBlockUsed;
	return *this;
}



template <class Type>
const dvector<Type> & dvector<Type>::operator=( dvector && moved)
{
	Blocks = std::move(moved.Blocks);
	iCurBlock = moved.iCurBlock;
	iCurBlockUsed = moved.iCurBlockUsed;
	return *this;
}


template <class Type>
void dvector<Type>::clear()
{
	Blocks.clear();
	iCurBlock = 0;
	iCurBlockUsed = 0;
	Blocks.push_back(BlockType());
}


template <class Type>
void dvector<Type>::fill(const Type & value)
{
	size_t nCount = Blocks.size();
	for (unsigned int i = 0; i < nCount; ++i)
		Blocks[i].fill(value);
}




template <class Type>
void dvector<Type>::resize( size_t nCount )
{
	if (length() == nCount)
		return;

	// figure out how many segments we need
	int nNumSegs = 1 + (int)nCount/nBlockSize;

	// figure out how many are currently allocated...
	size_t nCurCount = Blocks.size();

	// erase extra segments memory
	// [RMS] not necessary right? std::array will deallocate?
	//for (int i = nNumSegs; i < nCurCount; ++i)
	//	Blocks[i] = null;

	// resize to right number of segments
	if (nNumSegs >= Blocks.size()) {
		// allocate new segments
		for (int i = (int)nCurCount; i < nNumSegs; ++i) {
			Blocks.push_back(BlockType());
		}
	} else {
		//Blocks.RemoveRange(nNumSegs, Blocks.Count - nNumSegs);
		Blocks.resize(nNumSegs);
	}

	// mark last segment
	iCurBlockUsed = (unsigned int)(nCount - (nNumSegs-1)*nBlockSize);
	iCurBlock = nNumSegs-1;
}


template <class Type>
void dvector<Type>::resize( size_t nCount, const Type & init_value )
{
	size_t nCurSize = size();
	resize(nCount);
	for ( size_t nIndex = nCurSize; nIndex < nCount; ++nIndex ) 
		Blocks[ nIndex / m_nSegmentSize ].pData[ nIndex % m_nSegmentSize ] = init_value;
}



template <class Type>
bool dvector<Type>::empty() const
{
	return iCurBlock == 0 && iCurBlockUsed == 0;
}

template <class Type>
size_t dvector<Type>::size() const
{
	return iCurBlock * nBlockSize + iCurBlockUsed;
}
template <class Type>
size_t dvector<Type>::length() const
{
	return iCurBlock * nBlockSize + iCurBlockUsed;
}

template <class Type>
int dvector<Type>::block_size() const
{
	return nBlockSize;
}

template <class Type>
size_t dvector<Type>::byte_count() const
{
	int nb = (int)Blocks.size();
	return (nb == 0) ? 0 : nb * nBlockSize * sizeof(Type);
}



template <class Type>
void dvector<Type>::add(const Type & value)
{
	if (iCurBlockUsed == nBlockSize) {
		if (iCurBlock == Blocks.size() - 1)
			Blocks.push_back(BlockType());
		iCurBlock++;
		iCurBlockUsed = 0;
	}
	Blocks[iCurBlock][iCurBlockUsed] = value;
	iCurBlockUsed++;
}


template <class Type>
void dvector<Type>::push_back( const Type & data )
{
	add(data);
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
	if (iCurBlockUsed > 0)
		iCurBlockUsed--;
	if (iCurBlockUsed == 0 && iCurBlock > 0) {
		iCurBlock--;
		iCurBlockUsed = nBlockSize;
		// remove block ??
	}
}


template <class Type>
void dvector<Type>::insertAt(const Type & data, unsigned int nIndex)
{
	size_t s = size();
	if (nIndex == s) {
		push_back(data);
	}
	else if (nIndex > s) {
		resize(nIndex);
		push_back(data);
	}
	else {
		(*this)[nIndex] = data;
	}
}




template <class Type>
Type & dvector<Type>::front()
{
	return Blocks[0][0];
}
template <class Type>
const Type & dvector<Type>::front() const
{
	return Blocks[0][0];
}

template <class Type>
Type & dvector<Type>::back()
{
	return Blocks[iCurBlock][iCurBlockUsed - 1];
}
template <class Type>
const Type & dvector<Type>::back() const
{
	return Blocks[iCurBlock][iCurBlockUsed - 1];
}



template <class Type>
Type & dvector<Type>::operator[]( unsigned int i )
{
	return Blocks[i >> nShiftBits][i & nBlockIndexBitmask];
}

template <class Type>
const Type & dvector<Type>::operator[]( unsigned int i ) const
{
	return Blocks[i >> nShiftBits][i & nBlockIndexBitmask];
}



template<typename Type>
template<typename Func>
void dvector<Type>::apply( const Func & f )
{
	for (int bi = 0; bi < iCurBlock; ++bi) {
		auto block = Blocks[bi];
		for (int k = 0; k < nBlockSize; ++k)
			applyF(block[k], k);
	}
	auto lastblock = Blocks[iCurBlock];
	for (int k = 0; k < iCurBlockUsed; ++k)
		applyF(lastblock[k], k);
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


