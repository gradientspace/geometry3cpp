#pragma once

#include <string>

#include <dvector.h>
#include <g3Debug.h>


namespace g3
{
/// <summary>
/// small_list_set stores a set of short integer-valued variable-size lists.
/// The lists are encoded into a few large dvector buffers, with internal pooling,
/// so adding/removing lists usually does not involve any or delete ops.
/// 
/// The lists are stored in two parts. The first N elements are stored in a linear
/// subset of a dvector. If the list spills past these N elements, the extra elements
/// are stored in a linked list (which is also stored in a flat array).
/// 
/// Each list stores its count, so list-size operations are constant time.
/// All the internal "pointers" are 32-bit.
/// </summary>
class small_list_set
{
public:
	static constexpr int Null = -1;

	static constexpr int BLOCKSIZE = 8;
	static constexpr int BLOCK_LIST_OFFSET = BLOCKSIZE + 1;

    dvector<int> list_heads;        // each "list" is stored as index of first element in block-store (like a pointer)

    dvector<int> block_store;       // flat buffer used to store per-list initial block
                                    // blocks are BLOCKSIZE+2 long, elements are [CurrentCount, item0...itemN, LinkedListPtr]

    dvector<int> free_blocks;       // list of free blocks, indices into block_store
    int allocated_count = 0;

    dvector<int> linked_store;      // flat buffer used for linked-list elements,
                                    // each element is [value, next_ptr]

    int free_head_ptr;              // index of first free element in linked_store


    small_list_set()
    {
        list_heads = dvector<int>();
        linked_store = dvector<int>();
        free_head_ptr = Null;
        block_store = dvector<int>();
        free_blocks = dvector<int>();
    }


    small_list_set(const small_list_set & copy)
    {
        linked_store = dvector<int>(copy.linked_store);
        free_head_ptr = copy.free_head_ptr;
        list_heads = dvector<int>(copy.list_heads);
        block_store = dvector<int>(copy.block_store);
        free_blocks = dvector<int>(copy.free_blocks);
    }


    /// <summary>
    /// returns largest current list_index
    /// </summary>
    size_t Size() const {
        return list_heads.size();
    }

    /// <summary>
    /// resize the list-of-lists
    /// </summary>
    void Resize(int new_size)
    {
        int cur_size = (int)list_heads.size();
        if (new_size > cur_size) {
            list_heads.resize(new_size);
            for (int k = cur_size; k < new_size; ++k)
                list_heads[k] = Null;
        }
    }


    /// <summary>
    /// create a list at list_index
    /// </summary>
    void AllocateAt(int list_index)
    {
        if (list_index >= list_heads.size()) {
            int j = (int)list_heads.size();
            list_heads.insertAt(Null, list_index);
            // need to set intermediate values to null! 
            while (j < list_index) {
                list_heads[j] = Null;
                j++;
            }
        } else {
			gDevAssert(list_heads[list_index] == Null);
            //throw Exception("small_list_set: list at " + list_index + " is not empty!");
        }
    }


    /// <summary>
    /// insert val into list at list_index. 
    /// </summary>
    void Insert(int list_index, int val)
    {
        int block_ptr = list_heads[list_index];
        if ( block_ptr == Null ) {
            block_ptr = allocate_block();
            block_store[block_ptr] = 0;
            list_heads[list_index] = block_ptr;
        }

        int N = block_store[block_ptr];
        if (N < BLOCKSIZE) {
            block_store[block_ptr + N + 1] = val;
        } else {
            // spill to linked list
            int cur_head = block_store[block_ptr + BLOCK_LIST_OFFSET];

            if (free_head_ptr == Null) {
                // allocate linkedlist node
                int new_ptr = (int)linked_store.size();
                linked_store.add(val);
                linked_store.add(cur_head);
                block_store[block_ptr + BLOCK_LIST_OFFSET] = new_ptr;
            } else {
                // pull from free list
                int free_ptr = free_head_ptr;
                free_head_ptr = linked_store[free_ptr + 1];
                linked_store[free_ptr] = val;
                linked_store[free_ptr + 1] = cur_head;
                block_store[block_ptr + BLOCK_LIST_OFFSET] = free_ptr;
            }
        }

        // count element
        block_store[block_ptr] += 1;
    }



    /// <summary>
    /// remove val from the list at list_index. return false if val was not in list.
    /// </summary>
    bool Remove(int list_index, int val)
    {
        int block_ptr = list_heads[list_index];
        int N = block_store[block_ptr];


        int iEnd = block_ptr + std::min(N, BLOCKSIZE);
        for ( int i = block_ptr+1; i <= iEnd; ++i ) {

            if ( block_store[i] == val ) {
                for ( int j = i+1; j <= iEnd; ++j )     // shift left
                    block_store[j-1] = block_store[j];
                //block_store[iEnd] = -2;     // OPTIONAL

                if (N > BLOCKSIZE) {
                    int cur_ptr = block_store[block_ptr + BLOCK_LIST_OFFSET];
                    block_store[block_ptr + BLOCK_LIST_OFFSET] = linked_store[cur_ptr + 1];  // point to cur->next
                    block_store[iEnd] = linked_store[cur_ptr];
                    add_free_link(cur_ptr); 
                }

                block_store[block_ptr] -= 1;
                return true;
            }

        }

        // search list
        if ( N > BLOCKSIZE ) {
            if ( remove_from_linked_list(block_ptr, val) ) {
                block_store[block_ptr] -= 1;
                return true;
            }
        }

        return false;
    }



    /// <summary>
    /// move list at from_index to to_index
    /// </summary>
    void Move(int from_index, int to_index)
    {
		gDevAssert(list_heads[to_index] == Null);
		gDevAssert(list_heads[from_index] != Null);
        list_heads[to_index] = list_heads[from_index];
        list_heads[from_index] = Null;
    }






    /// <summary>
    /// remove all elements from list at list_index
    /// </summary>
    void Clear(int list_index)
    {
        int block_ptr = list_heads[list_index];
        if (block_ptr != Null) {
            int N = block_store[block_ptr];

            // if we have spilled to linked-list, free nodes
            if ( N > BLOCKSIZE ) {
                int cur_ptr = block_store[block_ptr + BLOCK_LIST_OFFSET];
                while (cur_ptr != Null) {
                    int free_ptr = cur_ptr;
                    cur_ptr = linked_store[cur_ptr + 1];
                    add_free_link(free_ptr);
                }
                block_store[block_ptr + BLOCK_LIST_OFFSET] = Null;
            }

            // free our block
            block_store[block_ptr] = 0;
            free_blocks.push_back(block_ptr);
            list_heads[list_index] = Null;
        }

    }


    /// <summary>
    /// return size of list at list_index
    /// </summary>
    int Count(int list_index) const
    {
        int block_ptr = list_heads[list_index];
        return (block_ptr == Null) ? 0 : block_store[block_ptr];
    }


    /// <summary>
    /// search for val in list at list_index
    /// </summary>
    bool Contains(int list_index, int val) const
    {
        int block_ptr = list_heads[list_index];
        if (block_ptr != Null) {
            int N = block_store[block_ptr];
            if (N < BLOCKSIZE) {
                int iEnd = block_ptr + N;
                for (int i = block_ptr + 1; i <= iEnd; ++i) {
                    if (block_store[i] == val)
                        return true;
                }
            } else {
                // we spilled to linked list, have to iterate through it as well
                int iEnd = block_ptr + BLOCKSIZE;
                for (int i = block_ptr + 1; i <= iEnd; ++i) {
                    if (block_store[i] == val)
                        return true;
                }
                int cur_ptr = block_store[block_ptr + BLOCK_LIST_OFFSET];
                while (cur_ptr != Null) {
                    if (linked_store[cur_ptr] == val)
                        return true;
                    cur_ptr = linked_store[cur_ptr + 1];
                }
            }
        }
        return false;
    }


    /// <summary>
    /// return the first item in the list at list_index (no zero-size-list checking)
    /// </summary>
    int First(int list_index) const
    {
        int block_ptr = list_heads[list_index];
        return block_store[block_ptr+1];
    }



	friend class value_iterator;

	// AAAHHH
	class value_iterator
	{
	public:
		inline value_iterator() { p = nullptr; mapF = nullptr; list_index = 0; }

		inline bool operator==(const value_iterator & r2) const {
			return p == r2.p && list_index == r2.list_index;
		}
		inline bool operator!=(const value_iterator & r2) const {
			return p != r2.p || list_index != r2.list_index || iCur != r2.iCur || cur_ptr != r2.cur_ptr;
		}

		inline int operator*() const {
			return (mapF == nullptr) ? cur_value : mapF(cur_value);
		}

		inline const value_iterator & operator++() {		// prefix
			this->goto_next();
			return *this;
		}
		//inline value_iterator operator++(int) {		// postfix
		//	index_iterator copy(*this);
		//	this->goto_next();
		//	return copy;
		//}


	protected:
		inline void goto_next() {
			if (N == 0)
				return;
			goto_next_overflow();
		}
		inline void goto_next_overflow() {
			if (iCur <= iEnd) {
				cur_value = p->block_store[iCur];
				iCur++;
			} else if (cur_ptr != Null) {
				cur_value = p->linked_store[cur_ptr];
				cur_ptr = p->linked_store[cur_ptr + 1];
			} else
				set_to_end();
		}

		inline value_iterator(const small_list_set * pVector, int list_index, bool is_end, 
			const std::function<int(int)> & value_mapper = nullptr )
		{
			p = pVector;
			mapF = value_mapper;
			this->list_index = list_index;
			if (is_end) {
				set_to_end();
			} else {
				block_ptr = p->list_heads[list_index];
				if (block_ptr != p->Null) {
					N = p->block_store[block_ptr];
					iEnd = (N < BLOCKSIZE) ? (block_ptr + N) : (block_ptr + BLOCKSIZE);
					iCur = block_ptr + 1;
					cur_ptr = (N < BLOCKSIZE) ? Null : p->block_store[block_ptr + BLOCK_LIST_OFFSET];
					goto_next();
				} else
					set_to_end();
			}
		}

		inline void set_to_end() {
			block_ptr = p->Null;
			N = 0;
			iCur = -1;
			cur_ptr = -1;
		}

		const small_list_set * p;
		std::function<int(int)> mapF;
		int list_index;
		int block_ptr;
		int N;
		int iEnd;
		int iCur;
		int cur_ptr;
		int cur_value;
		friend class small_list_set;
	};
	inline value_iterator begin_values(int list_index) const {
		return value_iterator(this, list_index, false);
	}
	inline value_iterator begin_values(int list_index, const std::function<int(int)> & value_mapper) const {
		return value_iterator(this, list_index, false, value_mapper);
	}
	inline value_iterator end_values(int list_index) const {
		return value_iterator(this, list_index, true);
	}


	class value_enumerable
	{
	public:
		const small_list_set * p;
		int list_index;
		std::function<int(int)> value_mapper;
		value_enumerable() {}
		value_enumerable(const small_list_set * p, int list_index, std::function<int(int)> value_mapper = nullptr) {
			this->p = p; 
			this->list_index = list_index; 
			this->value_mapper = value_mapper;
		}
		typename small_list_set::value_iterator begin() { return p->begin_values(list_index, value_mapper); }
		typename small_list_set::value_iterator end() { return p->end_values(list_index); }
	};
	inline value_enumerable values(int list_index) const {
		return value_enumerable(this, list_index);
	}
	inline value_enumerable values(int list_index, const std::function<int(int)> & value_mapper) const {
		return value_enumerable(this, list_index, value_mapper);
	}

/*
    /// <summary>
    /// iterate over the values of list at list_index
    /// </summary>
    IEnumerable<int> ValueItr(int list_index)
    {
        int block_ptr = list_heads[list_index];
        if (block_ptr != Null) {
            int N = block_store[block_ptr];
            if ( N < BLOCKSIZE ) {
                int iEnd = block_ptr + N;
                for (int i = block_ptr + 1; i <= iEnd; ++i)
                    yield return block_store[i];
            } else {
                // we spilled to linked list, have to iterate through it as well
                int iEnd = block_ptr + BLOCKSIZE;
                for (int i = block_ptr + 1; i <= iEnd; ++i)
                    yield return block_store[i];
                int cur_ptr = block_store[block_ptr + BLOCK_LIST_OFFSET];
                while (cur_ptr != Null) {
                    yield return linked_store[cur_ptr];
                    cur_ptr = linked_store[cur_ptr + 1];
                }
            }
        }
    }
*/

    /// <summary>
    /// search for findF(list_value) == true, of list at list_index, and return list_value
    /// </summary>
    int Find(int list_index, const std::function<bool(int)> & findF, int invalidValue = -1 ) const
    {
        int block_ptr = list_heads[list_index];
        if (block_ptr != Null) {
            int N = block_store[block_ptr];
            if (N < BLOCKSIZE) {
                int iEnd = block_ptr + N;
                for (int i = block_ptr + 1; i <= iEnd; ++i) {
                    int val = block_store[i];
                    if (findF(val))
                        return val;
                }
            } else {
                // we spilled to linked list, have to iterate through it as well
                int iEnd = block_ptr + BLOCKSIZE;
                for (int i = block_ptr + 1; i <= iEnd; ++i) {
                    int val = block_store[i];
                    if (findF(val))
                        return val;
                }
                int cur_ptr = block_store[block_ptr + BLOCK_LIST_OFFSET];
                while (cur_ptr != Null) {
                    int val = linked_store[cur_ptr];
                    if (findF(val))
                        return val;
                    cur_ptr = linked_store[cur_ptr + 1];
                }
            }
        }
        return invalidValue;
    }





    /// <summary>
    /// search for findF(list_value) == true, of list at list_index, and replace with new_value.
    /// returns false if not found
    /// </summary>
    bool Replace(int list_index, const std::function<bool(int)> & findF, int new_value)
    {
        int block_ptr = list_heads[list_index];
        if (block_ptr != Null) {
            int N = block_store[block_ptr];
            if (N < BLOCKSIZE) {
                int iEnd = block_ptr + N;
                for (int i = block_ptr + 1; i <= iEnd; ++i) {
                    int val = block_store[i];
                    if (findF(val)) {
                        block_store[i] = new_value;
                        return true;
                    }
                }
            } else {
                // we spilled to linked list, have to iterate through it as well
                int iEnd = block_ptr + BLOCKSIZE;
                for (int i = block_ptr + 1; i <= iEnd; ++i) {
                    int val = block_store[i];
                    if (findF(val)) {
                        block_store[i] = new_value;
                        return true;
                    }
                }
                int cur_ptr = block_store[block_ptr + BLOCK_LIST_OFFSET];
                while (cur_ptr != Null) {
                    int val = linked_store[cur_ptr];
                    if (findF(val)) {
                        linked_store[cur_ptr] = new_value;
                        return true;
                    }
                    cur_ptr = linked_store[cur_ptr + 1];
                }
            }
        }
        return false;
    }

protected:


    // grab a block from the free list, or allocate a one
    int allocate_block()
    {
        int nfree = (int)free_blocks.size();
        if ( nfree > 0 ) {
            int ptr = free_blocks[nfree - 1];
            free_blocks.pop_back();
            return ptr;
        }
        int nsize = (int)block_store.size();
        block_store.insertAt(Null, nsize + BLOCK_LIST_OFFSET);
        block_store[nsize] = 0;
        allocated_count++;
        return nsize;
    }


    // push a link-node onto the free list
    void add_free_link(int ptr)
    {
        linked_store[ptr + 1] = free_head_ptr;
        free_head_ptr = ptr;
    }


    // remove val from the linked-list attached to block_ptr
    bool remove_from_linked_list(int block_ptr, int val)
    {
        int cur_ptr = block_store[block_ptr + BLOCK_LIST_OFFSET];
        int prev_ptr = Null;
        while (cur_ptr != Null) {
            if (linked_store[cur_ptr] == val) {
                int next_ptr = linked_store[cur_ptr + 1];
                if (prev_ptr == Null) {
                    block_store[block_ptr + BLOCK_LIST_OFFSET] = next_ptr;
                } else {
                    linked_store[prev_ptr + 1] = next_ptr;
                }
                add_free_link(cur_ptr);
                return true;
            }
            prev_ptr = cur_ptr;
            cur_ptr = linked_store[cur_ptr + 1];
        }
        return false;
    }


public:
    std::string MemoryUsage()
    {
		std::ostringstream s;
		s << "ListSize " << list_heads.size() 
			<< "  Blocks Count " << allocated_count 
			<< " Free " << (free_blocks.size() * sizeof(int) / 1024) 
			<< " Mem " << block_store.size() 
			<< "kb  Linked Mem " << (linked_store.size() * sizeof(int) / 1024) << "kb";
		return s.str();
    }


};

}