#include <cstddef>
#include <new>

#if defined(__APPLE__) && defined(__MACH__)
// macOS / iOS
#include <malloc/_malloc.h>
#elif defined(__linux__)
// Linux
#include <malloc.h>
#endif

template <typename T, size_t Alignment> struct AlignedAllocator {
  using value_type = T;

  template <typename U> struct rebind {
    using other = AlignedAllocator<U, Alignment>;
  };
  T *allocate(size_t n) {
    void *ptr = nullptr;
    if (posix_memalign(&ptr, Alignment, n * sizeof(T)) != 0)
      throw std::bad_alloc();
    return reinterpret_cast<T *>(ptr);
  }
  void deallocate(T *p, size_t) noexcept { free(p); }

  bool operator==(const AlignedAllocator &) const noexcept { return true; }
  bool operator!=(const AlignedAllocator &) const noexcept { return false; }

  using propagate_on_container_swap = std::true_type;
};