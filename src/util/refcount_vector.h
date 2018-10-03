#pragma once

#include <string>

#include <dvector.h>
#include <g3Debug.h>
#include <iterator_util.h>

namespace g3 {

/// <summary>
/// RefCountedvector is used to keep track of which indices in a linear index list are in use/referenced.
/// A free list is tracked so that unreferenced indices can be re-used.
///
/// The enumerator iterates over valid indices (ie where refcount > 0)
/// 
/// **refcounts are shorts** so the maximum count is 65536. 
/// No overflow checking is done in release builds.
/// 
/// </summary>
class refcount_vector
{
public:
    static constexpr short invalid = -1;

    dvector<short> ref_counts;
    dvector<int> free_indices;
    int used_count;

    refcount_vector()
    {
        ref_counts = dvector<short>();
        free_indices = dvector<int>();
        used_count = 0;
    }

    refcount_vector(const refcount_vector & copy)
    {
        ref_counts = dvector<short>(copy.ref_counts);
        free_indices = dvector<int>(copy.free_indices);
        used_count = copy.used_count;
    }

    //refcount_vector(short * raw_ref_counts, bool build_free_list = false)
    //{
    //    ref_counts = new dvector<short>(raw_ref_counts);
    //    free_indices = new dvector<int>();
    //    used_count = 0;
    //    if (build_free_list)
    //        rebuild_free_list();
    //}

    dvector<short> RawRefCounts() {
        return ref_counts; 
    }

    bool empty() const {
        return used_count == 0; 
    }
	size_t count() const {
        return used_count;
    }
    size_t max_index() const {
        return ref_counts.size(); 
    }
    bool is_dense() const {
        return free_indices.length() == 0;
    }


    bool isValid(int index) const {
        return ( index >= 0 && index < ref_counts.size() && ref_counts[index] > 0 );
    }
    bool isValidUnsafe(int index) const {
        return ref_counts[index] > 0;
    }


    int refCount(int index) const {
        int n = ref_counts[index];
        return (n == invalid) ? 0 : n;
    }
    int rawRefCount(int index) const {
        return ref_counts[index];
    }


    int allocate() {
        used_count++;
        if (free_indices.empty()) {
            // [RMS] do we need this branch anymore? 
            ref_counts.push_back(1);
            return (int)ref_counts.size() - 1;
        } else {
            int iFree = invalid;
            while (iFree == invalid && free_indices.empty() == false) {
                iFree = free_indices.back();
                free_indices.pop_back();
            }
            if (iFree != invalid) {
                ref_counts[iFree] = 1;
                return iFree;
            } else {
                ref_counts.push_back(1);
                return (int)ref_counts.size() - 1;
            }
        }
    }



    int increment(int index, short increment = 1) {
        gDevAssert( isValid(index)  );
        // debug check for overflow...
        gDevAssert(  (short)(ref_counts[index] + increment) > 0 );
        ref_counts[index] += increment;
        return ref_counts[index];       
    }

    void decrement(int index, short decrement = 1) {
        gDevAssert( isValid(index) );
        ref_counts[index] -= decrement;
        gDevAssert(ref_counts[index] >= 0);
        if (ref_counts[index] == 0) {
            free_indices.push_back(index);
            ref_counts[index] = invalid;
            used_count--;
        }
    }



    /// <summary>
    /// allocate at specific index, which must either be larger than current max index,
    /// or on the free list. If larger, all elements up to this one will be pushed onto
    /// free list. otherwise we have to do a linear search through free list.
    /// If you are doing many of these, it is likely faster to use 
    /// allocate_at_unsafe(), and then rebuild_free_list() after you are done.
    /// </summary>
    bool allocate_at(int index)
    {
        if (index >= ref_counts.size()) {
            int j = (int)ref_counts.size();
            while (j < index) {
                ref_counts.push_back(invalid);
                free_indices.push_back(j);
                ++j;
            }
            ref_counts.push_back(1);
            used_count++;
            return true;

        } else {
            if (ref_counts[index] > 0)
                return false;

            int N = (int)free_indices.size();
            for (int i = 0; i < N; ++i) {
                if ( free_indices[i] == index ) {
                    free_indices[i] = invalid;
                    ref_counts[index] = 1;
                    used_count++;
                    return true;
                }
            }
            return false;
        }
    }


    /// <summary>
    /// allocate at specific index, which must be free or larger than current max index.
    /// However, we do not update free list. So, you probably need to do 
    /// rebuild_free_list() after calling this.
    /// </summary>
    bool allocate_at_unsafe(int index)
    {
        if (index >= ref_counts.size()) {
            int j = (int)ref_counts.size();
            while (j < index) {
                ref_counts.push_back(invalid);
                ++j;
            }
            ref_counts.push_back(1);
            used_count++;
            return true;

        } else {
            if (ref_counts[index] > 0)
                return false;
            ref_counts[index] = 1;
            used_count++;
            return true;
        }
    }



    // [RMS] really should not use this!!
    void set_Unsafe(int index, short count)
    {
        ref_counts[index] = count;
    }

    // todo:
    //   remove
    //   clear


    void rebuild_free_list()
    {
        free_indices = dvector<int>();
        used_count = 0;

        int N = (int)ref_counts.length();
        for ( int i = 0; i < N; ++i ) {
            if (ref_counts[i] > 0)
                used_count++;
            else
                free_indices.add(i);
        }
    }


    void trim(int maxIndex)
    {
        free_indices = dvector<int>();
        ref_counts.resize(maxIndex);
        used_count = maxIndex;
    }



	/*
	 * base iterator for indices with valid refcount (skips zero-refcount indices)
	 */
	class base_iterator
	{
	public:
		inline base_iterator() { p = nullptr; m_nIndex = 0; m_nLast = 0; }

		inline bool operator==(const base_iterator & r2) const {
			return m_nIndex == r2.m_nIndex;
		}
		inline bool operator!=(const base_iterator & r2) const {
			return m_nIndex != r2.m_nIndex;
		}

	protected:
		inline void goto_next() {
			if (m_nIndex != m_nLast)
				m_nIndex++;
			while (m_nIndex != m_nLast && p->isValidUnsafe(m_nIndex) == false)
				m_nIndex++;
		}

		inline base_iterator(const refcount_vector * pVector, int nIndex, int nLast)
		{
			p = pVector;
			m_nIndex = nIndex;
			m_nLast = nLast;
			if (m_nIndex != m_nLast && p->isValidUnsafe(m_nIndex) == false)
				goto_next();		// initialize
		}
		const refcount_vector * p;
		int m_nIndex;
		int m_nLast;
		friend class refcount_vector;
	};


	/*
	 *  iterator over valid indices (ie non-zero refcount)
	 */
	class index_iterator : public base_iterator
	{
	public:
		inline index_iterator() : base_iterator() {}

		inline int operator*() const {
			return this->m_nIndex;
		}

		inline index_iterator & operator++() {		// prefix
			this->goto_next();
			return *this;
		}
		inline index_iterator operator++(int) {		// postfix
			index_iterator copy(*this);
			this->goto_next();
			return copy;
		}

	protected:
		inline index_iterator(const refcount_vector * pVector, int nIndex, int nLast) : base_iterator(pVector, nIndex, nLast)
		{}
		friend class refcount_vector;
	};


	inline index_iterator begin_indices() const {
		return index_iterator(this, (int)0, (int)ref_counts.size());
	}
	inline index_iterator end_indices() const {
		return index_iterator(this, (int)ref_counts.size(), (int)ref_counts.size());
	}



	/*
	 * enumerable object that provides begin()/end() semantics, so
	 * you can iterate over valid indices using range-based for loop
	 */
	class index_enumerable
	{
	public:
		const refcount_vector * pVector;
		index_enumerable() { pVector = nullptr; }
		index_enumerable(const refcount_vector * p) { pVector = p; }
		typename refcount_vector::index_iterator begin() { return pVector->begin_indices(); }
		typename refcount_vector::index_iterator end() { return pVector->end_indices(); }
	};

	/*
	 * returns iteration object over valid indices
	 * usage: for (int idx : indices()) { ... }
	 */
	inline index_enumerable indices() const {
		return index_enumerable(this);
	}


	/*
	 * enumerable object that maps indices output by index_iteration to a second type
	 */
	template<typename ToType>
	class mapped_enumerable
	{
	public:
		std::function<ToType(int)> map_func;
		index_enumerable enumerable;

		mapped_enumerable(const index_enumerable & enumerable, std::function<ToType(int)> map_func) {
			this->enumerable = enumerable;
			this->map_func = map_func;
		}

		typename mapped_iterator<int, ToType, index_iterator> begin() {
			return mapped_iterator<int, ToType, index_iterator>(enumerable.begin(), map_func);
		}

		typename mapped_iterator<int, ToType, index_iterator> end() {
			return mapped_iterator<int, ToType, index_iterator>(enumerable.end(), map_func);
		}
	};

	/*
	* returns iteration object over mapping applied to valid indices
	* eg usage: for (Vector3d v : mapped_indices(fn_that_looks_up_mesh_vtx_from_id)) { ... }
	*/
	template<typename ToType>
	inline mapped_enumerable<ToType> mapped_indices(std::function<ToType(int)> map_func) const {
		return mapped_enumerable<ToType>(indices(), map_func);
	}




	/*
	* iteration object that maps indices output by index_iteration to a second type
	*/
	class filtered_enumerable
	{
	public:
		std::function<bool(int)> filter_func;
		index_enumerable enumerable;

		filtered_enumerable(const index_enumerable & enumerable, std::function<bool(int)> filter_func) {
			this->enumerable = enumerable;
			this->filter_func = filter_func;
		}

		typename filtered_iterator<int, index_iterator> begin() {
			return filtered_iterator<int, index_iterator>(enumerable.begin(), enumerable.end(), filter_func);
		}

		typename filtered_iterator<int, index_iterator> end() {
			return filtered_iterator<int, index_iterator>(enumerable.end(), enumerable.end(), filter_func);
		}
	};

	inline filtered_enumerable filtered_indices(std::function<bool(int)> filter_func) const {
		return filtered_enumerable(indices(), filter_func);
	}





	// [RMS] how?
    //System.Collections.IEnumerator GetEnumerator()
    //{
    //    int nIndex = 0;
    //    int nLast = max_index;

    //    // skip leading empties
    //    while (nIndex != nLast && ref_counts[nIndex] <= 0)
    //        nIndex++;

    //    while (nIndex != nLast) {
    //        yield return nIndex;

    //        if (nIndex != nLast)
    //            nIndex++;
    //        while (nIndex != nLast && ref_counts[nIndex] <= 0)
    //            nIndex++;
    //    }
    //}


    std::string UsageStats() {
		std::ostringstream str;
		str << "RefCountSize " << ref_counts.size() 
			<< "  FreeSize" << free_indices.size() 
			<< " FreeMem " << (free_indices.byte_count() / 1024) << "kb";
		return str.str();
    }



    //std::string debug_print()
    //{
    //    string s = string.Format("size {0} used {1} free_size {2}\n", ref_counts.size, used_count, free_indices.size);
    //    for (int i = 0; i < ref_counts.size; ++i)
    //        s += string.Format("{0}:{1} ", i, ref_counts[i]);
    //    s += "\nfree:\n";
    //    for (int i = 0; i < free_indices.size; ++i)
    //        s += free_indices[i].ToString() + " ";
    //    return s;
    //}


};



}