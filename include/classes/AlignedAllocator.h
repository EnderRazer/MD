#include <cstddef>
#include <new>
#include <malloc/_malloc.h>

template <typename T, size_t Alignment>
struct AlignedAllocator {
    using value_type = T;
    
    template <typename U>
    struct rebind {
        using other = AlignedAllocator<U, Alignment>;
    };
    T* allocate(size_t n) {
        void* ptr = nullptr;
        if (posix_memalign(&ptr, Alignment, n * sizeof(T)) != 0)
            throw std::bad_alloc();
        return reinterpret_cast<T*>(ptr);
    }
    void deallocate(T* p, size_t) noexcept { free(p); }
};