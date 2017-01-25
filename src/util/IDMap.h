#pragma once

#include <g3types.h>
#include <g3Debug.h>

namespace g3 {

enum IDMapType 
{
	Dense = 0,
	Sparse = 1,
	Shift = 2
};


template <typename T>
class IDMap
{
public:

	IDMap() { m_eType = Sparse; }

	IDMap( IDMapType eType, int nInitial = 0 ) { 
		m_eType = eType; 
		if ( eType == Shift )
			m_nShift = nInitial;
		else if ( eType == Dense )
			Resize(nInitial, nInitial);
	}

	IDMap( unsigned int nOldSize, unsigned int nNewSize ) { 
		m_eType = Dense;
		Resize(nOldSize, nNewSize); 
	}

	IDMapType Type() const {
		return m_eType;
	}


	void Clear() {
		if ( m_eType == Shift )
			gDevAssert(false);
		else if ( m_eType == Sparse ) {
			vToNewMap.clear();
			vToOldMap.clear();
		} else {
			vToNew.resize(0);
			vToOld.resize(0);
		}
	}

	// for dense map, this sizes internal arrays appropriately
	inline void Resize( unsigned int nOldSize, unsigned int nNewSize ) {
		if ( m_eType == Shift || m_eType == Sparse ) {
			gDevAssert(false);
		} else {
			vToNew.resize(0);
			vToNew.resize(nOldSize, (T)InvalidID);
			vToOld.resize(0);
			vToOld.resize(nNewSize, (T)InvalidID);
		}
	}

	inline size_t OldSize() const {
		if ( m_eType == Shift || m_eType == Sparse ) {
			gDevAssert(false);
			return 0;
		} else
			return vToNewMap.size(); 
	}

	inline size_t NewSize() const {
		if ( m_eType == Shift || m_eType == Sparse ) {
			gDevAssert(false);
			return 0;
		} else
			return vToOldMap.size(); 
	}


	void SetShift( int nShift ) {
		if ( m_eType != Shift )
			gDevAssert(false);
		m_nShift = nShift;
	}


	inline void SetMap( T vOld, T vNew ) {
		if (m_eType == Shift) {
			gDevAssert( false );
		} else if (m_eType == Sparse) {
			vToNewMap[vOld] = vNew;
			if ( vNew != (T)InvalidID )
				vToOldMap[vNew] = vOld;
		} else {
			vToNew[vOld] = vNew;
			if ( vNew != (T)InvalidID )
				vToOld[vNew] = vOld;
		}
	}

	inline T GetNew( T vOld ) const {
		if (m_eType == Shift) {
			return vOld + m_nShift;
		} else if (m_eType == Sparse) {
			auto found = vToNewMap.find(vOld);
			return found != vToNewMap.end() ? found->second : (T)InvalidID;
		} else {
			return vToNew[vOld];
		}
	}
	inline T GetOld( T vNew ) const {
		if (m_eType == Shift) {
			return vNew - m_nShift;
		} else if (m_eType == Sparse) {
			auto found = vToOldMap.find(vNew);
			return found != vToOldMap.end() ? found->second : (T)InvalidID;
		} else {
			return vToOld[vNew];
		}
	}

	enum Direction {
		OldToNew = 0,
		NewToOld = 1
	};

	inline void Apply( std::vector<T> & v, Direction eDirection ) const {
		Apply( &v[0], v.size(), eDirection );
	}

	void Apply( T * v, size_t nCount, Direction eDirection ) const {
		if (m_eType == Shift) {
			if ( m_nShift == 0 )
				return;
			int nShift = (eDirection == OldToNew) ? m_nShift : -m_nShift;
			for ( unsigned int k = 0; k < nCount; ++k )
				v[k] += nShift;
		} else if ( m_eType == Sparse ) {
			const std::map<T,T> & m = (eDirection == OldToNew) ? vToNewMap : vToOldMap;
			for (unsigned int k = 0; k < nCount; ++k) {
				auto found = m.find( v[k] );
				v[k] = (found != m.end() ? found->second : (T)InvalidID);
			}
		} else {
			const std::vector<T> & m = (eDirection == OldToNew) ? vToNew : vToOld;
			for (unsigned int k = 0; k < nCount; ++k) 
				v[k] = m[v[k]];
		}
	}

	// [TODO] parallel apply?

private:
	IDMapType m_eType;
	std::vector<T> vToNew;
	std::vector<T> vToOld;
	std::map<T,T> vToNewMap;
	std::map<T,T> vToOldMap;
	int m_nShift = 0;
};


typedef IDMap<VertexID> VertexMap;
typedef IDMap<VertexID> TriangleMap;


} // namespace g3


