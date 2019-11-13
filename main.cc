#include <assert.h>
#include <iostream>

#include <emu_c_utils/emu_c_utils.h>
#include <cassert>
#include <cstddef>
#include <memory>

namespace emu {

    class local_arena
    {
    protected:
	typedef unsigned char uchar;
	// Value originally returned from mw_mallocrepl
	void * buffer;
	// Pointer to next available byte on this nodelet
	uchar * next_chunk;
	// Size of each stripe, in elements
	size_t size;
    public:

	// @param stripe_size: Number of bytes to reserve on each nodelet
	local_arena(size_t stripe_size)
	{
            assert(stripe_size > 0);
            // Allocate a stripe on each nodelet
            void * ptr = mw_mallocrepl(stripe_size);
            assert(ptr); // FIXME
            // Save the pointer in all copies of the allocator
            // since we don't know which one will get freed
            mw_replicated_init((long*)&buffer, (long)ptr);
            mw_replicated_init((long*)&size, (long)stripe_size);

            // Let the nth copy of next_chunk point to the stripe on the nth nodelet
            for (long nlet = 0; nlet < NODELETS(); ++nlet) {
                uchar * nth_chunk = (uchar*)mw_get_nth(buffer, nlet);
                uchar ** nth_ptr = (uchar**)mw_get_nth(&next_chunk, nlet);
                *nth_ptr = nth_chunk;
            }
	}

	~local_arena()
	{
            mw_free(buffer);
	}

	void *
        allocate(size_t n, const void * hint)
	{
            if (hint) { MIGRATE((void*)hint); }
            // Serial version
            void * ptr = next_chunk;
            next_chunk += n;
            // Parallel-safe version
            // T * ptr = (T*)ATOMIC_ADDMS((long*)next_chunk, (long)(n * sizeof(T)));
            return ptr;
	}
    };

// Reserve 2GB on each nodelet for satisfying allocations
    extern replicated local_arena g_arena;

    template<typename T>
    class local_arena_allocator
    {
    protected:
	local_arena & arena;
    public:
	typedef T value_type;

        // Declare friendship with all templated versions of this class
        template<typename U>
        friend class local_arena_allocator;

        // Default constructor, uses the global replicated arena
        local_arena_allocator() : arena(g_arena) {}

        // Alternate constructor for specifying a different arena to use
	explicit
        local_arena_allocator(local_arena & arena) : arena(arena) {}

        // Templated copy-constructor, for use with rebind
        template<typename U>
        local_arena_allocator(const local_arena_allocator<U>& other) : arena(other.arena) {}

        template<typename U>
        bool operator== (const local_arena_allocator<U>& other)
        {
	    return arena == other.arena;
        }

	T *
        allocate(size_t n)
	{
            return static_cast<T*>(
                arena.allocate(n * sizeof(T), nullptr)
		);
	}

	T *
        allocate(size_t n, const void * hint)
	{
            return static_cast<T*>(
                arena.allocate(n * sizeof(T), hint)
		);
	}

	void
        deallocate(T *, size_t)
	{
            // Not implemented
	}
    };

} // end namespace emu

// Reserve 2GB on each nodelet for satisfying allocations
replicated emu::local_arena emu::g_arena((1UL<<31));

//#include <tuple>
#include <vector>

#include <cilk.h>
//#include <memoryweb.h>

typedef long Index_t;
//typedef long Scalar_t;
//typedef std::vector<Index_t> IndexArray_t;
//typedef std::tuple<Index_t, Scalar_t> Pair_t;
//#if 0
//typedef std::vector<Pair_t> Row_t;
//#else
////typedef std::vector<Pair_t, emu::local_arena_allocator<Pair_t>>
////Row_t;
typedef std::vector<Index_t, emu::local_arena_allocator<Index_t>> Row_t;
//#endif

typedef Row_t * pRow_t;
typedef pRow_t * ppRow_t;

static inline Index_t n_map(Index_t i) { return i % NODELETS(); }
static inline Index_t r_map(Index_t i) { return i / NODELETS(); }
// inverse mapping
//static inline Index_t nr_inv(Index_t n, Index_t r)
//{ return r * NODELETS() + n; }

/*
 * Overrides default new to always allocate replicated storage for instances
 * of this class. repl_new is intended to be used as a parent class for
 * distributed data structure types.
 */
class repl_new
{
public:
    // Overrides default new to always allocate replicated storage for
    // instances of this class
    static void *
    operator new(std::size_t sz)
    {
        return mw_mallocrepl(sz);
    }

    // Overrides default delete to safely free replicated storage
    static void
    operator delete(void * ptr)
    {
        mw_free(ptr);
    }
};

class rMatrix_t : public repl_new
{
public:
    static rMatrix_t * create(Index_t nrows)
    {
        return new rMatrix_t(nrows);
    }

    rMatrix_t() = delete;
    rMatrix_t(const rMatrix_t &) = delete;
    rMatrix_t & operator=(const rMatrix_t &) = delete;
    rMatrix_t(rMatrix_t &&) = delete;
    rMatrix_t & operator=(rMatrix_t &&) = delete;

//    void build(Index_t nl_id,
//               IndexArray_t::iterator i_it,
//               IndexArray_t::iterator j_it,
//               IndexArray_t::iterator v_it,
//               Index_t nedges)
//    {
////        for (Index_t ix = 0; ix < nedges; ++ix)
////        {
////            setElement(*i_it, *j_it, *v_it);
////            ++i_it; ++j_it; ++v_it; // increment iterators
////        }
////
////        // count max degree
////        max_degree_ = 0;
////        for (Index_t row_idx = 0; row_idx < nrows_per_nodelet_; ++row_idx)
////        {
////            Index_t irow = nr_inv(nl_id, row_idx); // absolute row index
////            pRow_t r = getrow(irow);
////            Index_t deg = r->size();
////            max_degree_ = (deg > max_degree_) ? deg : max_degree_;
////        }
//    }
//
//    void set_max_degree()
//    {
////        Index_t max = 0;
////        for (Index_t i = 0; i < NODELETS(); ++i)
////        {
////            Index_t lmax = *(Index_t*)mw_get_nth(&this->max_degree_, i);
////            max = (lmax > max) ? lmax : max;
////        }
////        mw_replicated_init(&this->max_degree_, max);
//    }

//    Index_t nrows() { return nrows_; }
//    Index_t nrows() const { return nrows_; }
//
//    Index_t nrows_nl() { return nrows_per_nodelet_; }
//    Index_t nrows_nl() const { return nrows_per_nodelet_; }

    pRow_t getrow(Index_t i) { return rows_[n_map(i)] + r_map(i); }
//    pRow_t getrow(Index_t i) const { return rows_[n_map(i)] + r_map(i); }
//
//    Index_t * row_addr(Index_t i)
//    {
//        return (Index_t *)(rows_ + n_map(i));
//    }
//
//    Index_t max_degree() { return max_degree_; }
//
private:
    rMatrix_t(Index_t nrows) : nrows_(nrows)
    {
        nrows_per_nodelet_ = r_map(nrows_);
        if (n_map(nrows_) != 0) nrows_per_nodelet_++; // add empty rows

        rows_ = (ppRow_t)mw_malloc2d(NODELETS(),
                                     nrows_per_nodelet_ * sizeof(Row_t));

        // replicate the class across nodelets
        for (Index_t i = 1; i < NODELETS(); ++i)
        {
            memcpy(mw_get_nth(this, i), mw_get_nth(this, 0), sizeof(*this));
        }

        for (Index_t i = 0; i < NODELETS(); ++i)
        {
            cilk_migrate_hint(rows_ + i);
            cilk_spawn allocateRows(i);
        }
        cilk_sync;
    }

    void allocateRows(Index_t i)
    {
        for (Index_t row_idx = 0; row_idx < nrows_per_nodelet_; ++row_idx)
        {
            new(rows_[i] + row_idx) Row_t();
        }
    }

//    void setElement(Index_t irow, Index_t icol, Index_t const &val)
//    {
//        pRow_t r = rows_[n_map(irow)] + r_map(irow);
//
//        if (r->empty()) // empty row
//        {
//            r->push_back(std::make_tuple(icol, val));
//        }
//        else // insert into row
//        {
//            Row_t::iterator it = r->begin();
//            while (it != r->end() and std::get<0>(*it) < icol)
//            {
//                ++it;
//            }
//            if (it == r->end())
//            {
//                r->push_back(std::make_tuple(icol, val));
//            }
//            else
//            {
//                it = r->insert(it, std::make_tuple(icol, val));
//            }
//        }
//    }

    Index_t nrows_;
    Index_t nrows_per_nodelet_;
    ppRow_t rows_;
    //Index_t max_degree_;
};
typedef rMatrix_t * prMatrix_t;

////static inline
//bool index_exists(pRow_t r, Index_t icol)
//{
//    bool result = false;
//    Row_t::iterator rit = r->begin();
//    while (rit != r->end())
//    {
//        if (icol == std::get<0>(*rit))
//        {
//            result = true;
//            break;
//        }
//        ++rit;
//    }
//    return result;
//}
//
////static inline
//bool dot(Scalar_t & ans, pRow_t a, pRow_t b) // no semiring
//{
//    bool result = false;
//    Row_t::iterator ait = a->begin();
//    Row_t::iterator bit = b->begin();
//
//    ans = 0;
//    while (ait != a->end() && bit != b->end())
//    {
//        Index_t a_idx = std::get<0>(*ait);
//        Index_t b_idx = std::get<0>(*bit);
//
//        if (a_idx == b_idx)
//        {
//            ans += std::get<1>(*ait) * std::get<1>(*bit);
//            result = true;
//            ++ait;
//            ++bit;
//        }
//        else if (a_idx < b_idx)
//        {
//            ++ait;
//        }
//        else
//        {
//            ++bit;
//        }
//    }
//    return result;
//}
//
////static inline
//void row_kernel(Index_t irow,
//                prMatrix_t C,
//                prMatrix_t const M,
//                prMatrix_t const A,
//                prMatrix_t const B)
//{
//    // return for empty row of A
//    if (!A->getrow(irow)) return;
//
//    std::tuple<Index_t, Scalar_t> tmp;
//    // loop over columns
//    for (Index_t icol = 0; icol < A->nrows(); ++icol)
//    {
//        // continue for empty column of B
//        if (!B->getrow(icol)) continue;
//        // apply mask
//        if (!index_exists(M->getrow(irow), icol)) continue;
//
//        // compute the dot
//        Scalar_t ans;
//        if (dot(ans, A->getrow(irow), B->getrow(icol)))
//        {
//            std::get<0>(tmp) = icol;
//            std::get<1>(tmp) = ans;
//            C->getrow(irow)->push_back(tmp);
//        }
//    }
//}
//
////static inline
//void multi_row_kernel(Index_t nl_id,
//                      Index_t t,
//                      Index_t nrow,
//                      prMatrix_t C,
//                      prMatrix_t const M,
//                      prMatrix_t const A,
//                      prMatrix_t const B)
//{
//    for (Index_t j = t*nrow; j < (t+1)*nrow; ++j)
//    {
//        // absolute row index
//        Index_t irow = nr_inv(nl_id, j);
//        row_kernel(irow, C, M, A, B);
//    }
//
//}

//typedef std::vector<Index_t, emu::local_arena_allocator<Index_t>> myVec_t;

//static inline
//void ABT_Mask_NoAccum_kernel(
//    Index_t nl_id,              // nodelet_id
//    prMatrix_t C)//,               // output matrix
//    prMatrix_t const M,         // mask matrix
//    // SemiringT,               // semiring
//    prMatrix_t const A,         // Input matrix 1
//    prMatrix_t const B,         // Input matrix 2
//    bool replace_flag = false)  // put the answer in place?
//{

//    std::tuple<Index_t, Scalar_t> tmp;
//    std::get<0>(tmp) = 2;
//    std::get<1>(tmp) = 3;
//    C->getrow(0)->push_back(tmp);
//    C->getrow(0)->push_back(1);
//    myVec_t a;
//    a.push_back(10);

//    // making use of the fact we know that B equals L^T
//
//    // compute rows per thread
//    Index_t threads_per_nodelet = THREADS_PER_NODELET;
//    Index_t nrows_per_thread = A->nrows_nl() / threads_per_nodelet;
//    // if nrows_nl < threads_per_nodelet, all rows are remainder rows
//    Index_t nremainder_rows = A->nrows_nl() % threads_per_nodelet;
//
//    // spawn threads_per_nodelet threads
//    if (nrows_per_thread)
//    {
//        for (Index_t t = 0; t < threads_per_nodelet; ++t)
//        {
//            cilk_spawn multi_row_kernel(nl_id, t, nrows_per_thread,
//                                        C, M, A, B);
//        }
//        cilk_sync;
//    }
//
//    // spawn nremainder_rows threads
//    if (nremainder_rows)
//    {
//        Index_t offset = nrows_per_thread * threads_per_nodelet;
//
//        for (Index_t t = 0; t < nremainder_rows; ++t)
//        {
//            // absolute row index
//            Index_t irow = nr_inv(nl_id, t + offset);
//            cilk_spawn row_kernel(irow, C, M, A, B);
//        }
//        cilk_sync;
//    }
//}

//Scalar_t reduce(prMatrix_t A)
//{
//    Scalar_t sum = 0;
////    for (Index_t irow = 0; irow < A->nrows(); ++irow)
////    {
////        pRow_t pArow = A->getrow(irow);
////        Row_t::iterator ait = pArow->begin();
////        while (ait != pArow->end())
////        {
////            sum += std::get<1>(*ait);
////            ++ait;
////        }
////    }
//    return sum;
//}

//void initialize(Index_t nl_id, std::string const & filename, prMatrix_t M,
//                Index_t const nnodes, Index_t const nedges)
//{
//    Index_t tmp;
//    FILE *infile = mw_fopen(filename.c_str(), "r", &tmp);
//    mw_fread(&tmp, sizeof(Index_t), 1, infile);
//    //assert(tmp == nnodes); // needed with 19.09
//    mw_fread(&tmp, sizeof(Index_t), 1, infile);
//    //assert(tmp == nedges); // needed with 19.09
//
//    // thread local storage to read into
//    IndexArray_t iL(nedges);
//    IndexArray_t jL(nedges);
//    mw_fread(reinterpret_cast<void *>(iL.data()),
//             //sizeof(Index_t), iL.size(), infile);
//             1, sizeof(Index_t)*iL.size(), infile); // bug work around
//    mw_fread(reinterpret_cast<void *>(jL.data()),
//             //sizeof(Index_t), jL.size(), infile);
//             1, sizeof(Index_t) * jL.size(), infile); // bug work around
//    mw_fclose(infile);
//
//    // remove edges where i is a row not owned by this nodelet.
//    IndexArray_t iL_nl;
//    IndexArray_t jL_nl;
//    Index_t nedges_nl = 0;
//
//    for (Index_t e = 0; e < iL.size(); ++e)
//    {
//        Index_t i = iL[e];
//        Index_t j = jL[e];
//        if (n_map(i) == nl_id)
//        {
//            iL_nl.push_back(i);
//            jL_nl.push_back(j);
//            ++nedges_nl;
//        }
//    }
//
//    // build matrix
//    IndexArray_t v_nl(iL_nl.size(), 1);
//    M->build(nl_id, iL_nl.begin(), jL_nl.begin(), v_nl.begin(), nedges_nl);
//}

//Row_t testfn(prMatrix_t C)
//{
//    C->getrow(0)->push_back(1);
//    return C->getrow(0)[0];
//}

void testfn(Row_t & r)
{
    r.push_back(1);
}

int main(int argc, char* argv[])
{
//    if (argc != 2)
//    {
//        std::cerr << "Requires binary edge list." << std::endl;
//        std::cerr << "Usage: ./llt input.bin" << std::endl;
//        exit(1);
//    }
//
//    Index_t nnodes, nedges;
//    std::string filename = std::string(argv[1]);
//
//    // open file to get number of nodes and edges, then close
//    FILE *infile = mw_fopen(filename.c_str(), "r", &nnodes);
//    if (!infile)
//    {
//        fprintf(stderr, "Unable to open file: %s\n", filename.c_str());
//        exit(1);
//    }
//    mw_fread(&nnodes, sizeof(Index_t), 1, infile);
//    mw_fread(&nedges, sizeof(Index_t), 1, infile);
//    mw_fclose(infile);
//
//    std::cerr << "nnodes: " << nnodes << std::endl;
//    std::cerr << "nedges: " << nedges << std::endl;

//    Index_t nnodes = 1000;
//    prMatrix_t L = rMatrix_t::create(nnodes);

//    // spawn threads on each nodelet to read and build
//    for (Index_t i = 0; i < NODELETS(); ++i)
//    {
//        cilk_migrate_hint(L->row_addr(i));
//        cilk_spawn initialize(i, filename, L, nnodes, nedges);
//    }
//    cilk_sync;
//    L->set_max_degree();
//
//    std::cerr << "Initialization complete." << std::endl;
//    std::cerr << "Max degree: " << L->max_degree() << std::endl;
//
//    // answer matrix
//    prMatrix_t C = rMatrix_t::create(nnodes);
//
//#ifdef __PROFILE__
//    hooks_region_begin("6.1_llt");
//#endif
//    // solve L * L^T using ABT kernel
//    for (Index_t i = 0; i < NODELETS(); ++i)
//    {
//        cilk_migrate_hint(L->row_addr(i));
//        //cilk_spawn ABT_Mask_NoAccum_kernel(i, C, L, L, L);
//        cilk_spawn ABT_Mask_NoAccum_kernel(i, C);
    //cilk_spawn testfn(L);


    //pRow_t p = new Row_t;
    Row_t r;
//    pRow_t p = &r;
//    p->push_back(1);
    testfn(r);

//    }
//    cilk_sync;
//
//    // reduce
//    std::cerr << "Start reduction." << std::endl;
//    Scalar_t nTri = reduce(C);
//    std::cerr << "nTri: " << nTri << std::endl;
//
//    // clean up matrices
//    delete L;
//    delete C;
//
//#ifdef __PROFILE__
//    hooks_region_end();
//#endif

    return 0;
}
