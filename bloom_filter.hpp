#ifndef _BLOOM_FILTER_HPP
#define _BLOOM_FILTER_HPP

#include <algorithm>
#include <array>
#include <cstring>
#include <sdsl/bit_vectors.hpp>

#include "xxhash.hpp"

class BF {

private:
  uint64_t _get_hash(const uint64_t kmer) const {
    return xxh::xxhash<64>(&kmer, sizeof(uint64_t), 0);
  }

public:
  BF() : _mode(false), _bf(0, 0) { _size = 0; };
  BF(const size_t size) : _mode(false), _bf(size, 0) { _size = size; }
  ~BF() {}

  void add_key(const uint64_t kmer) {
    uint64_t hash = _get_hash(kmer);
    _bf[hash % _size] = 1;
  }

  bool test_key(const uint64_t kmer) {
    uint64_t hash = _get_hash(kmer);
    return _bf[hash % _size];
  }

  void switch_mode() {
    _mode = true;
    _brank = sdsl::bit_vector::rank_1_type(&_bf);
    _counts = sdsl::int_vector<16>(_brank(_size), 0, 16);
  }

  bool set_count(const uint64_t kmer, const uint16_t counter) {
    if (!_mode)
      return false;
    uint64_t hash = _get_hash(kmer);
    size_t bf_idx = hash % _size;
    if (_bf[bf_idx]) {
      size_t cnts_idx = _brank(bf_idx);
      uint32_t new_value = _counts[cnts_idx] + counter;
      _counts[cnts_idx] = new_value;
    }
    return true;
  }

  uint16_t get_count(const uint64_t kmer) const {
    if (_mode) {
      uint64_t hash = _get_hash(kmer);
      size_t bf_idx = hash % _size;
      if (_bf[bf_idx])
        return _counts[_brank(bf_idx)];
    }
    return 0;
  }

private:
  bool _mode; // false = write, true = read
  size_t _size;
  sdsl::bit_vector _bf;
  sdsl::bit_vector::rank_1_type _brank;
  sdsl::int_vector<16> _counts;
};

#endif
