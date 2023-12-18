#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <zlib.h>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "bloom_filter.hpp"
#include "kmc_api/kmc_file.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

static const char *const USAGE_MESSAGE =
    "Usage: MFCNV <reference.fa> <control_kmc_db> <case_kmc_db> \n"
    "      -m <INT>   minimum weight for kmers (default: 0)\n"
    "      -M <INT>   maximum weight for kmers (default: 65535)\n"
    "      -w <INT>   bin size (default: 1000)\n"
    "      -b <INT>   Bloom filter size in GB (default: 1)\n"
    "      -a         use average instead of median for normalization\n"
    // "      -@ <INT>   set threads (default: 1)\n"
    // "      -v         verbose mode\n"
    "      -h         display this help and exit\n"
    "\n";

static const uint8_t to_int[128] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 0
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 10
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 20
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 30
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 40
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 50
                                    0, 0, 0, 0, 0, 1, 0, 2, 0, 0, // 60
                                    0, 3, 0, 0, 0, 0, 0, 0, 0, 0, // 70
                                    0, 0, 0, 0, 4, 0, 0, 0, 0, 0, // 80
                                    0, 0, 0, 0, 0, 0, 0, 1, 0, 2, // 90
                                    0, 0, 0, 3, 0, 0, 0, 0, 0, 0, // 100
                                    0, 0, 0, 0, 0, 0, 4, 0, 0, 0, // 110
                                    0, 0, 0, 0, 0, 0, 0, 0};      // 120

inline uint8_t reverse_char(const uint8_t c) { return ((~c) & 3); }

uint64_t revcompl(uint64_t kmer, const uint8_t k) {
  uint64_t rckmer = 0;
  kmer = ~kmer;
  for (uint8_t i = 0; i < k; ++i) {
    rckmer = (rckmer << 2) | (kmer & 3);
    kmer >>= 2;
  }
  return rckmer;
}

inline uint64_t lsappend(const uint64_t kmer, const uint64_t c,
                         const uint64_t k) { // left shift and append
  return ((kmer << 2) | c) & ((1UL << 2 * k) - 1);
}

inline uint64_t rsprepend(const uint64_t kmer, const uint64_t c,
                          const uint64_t k) { // right shift and prepend
  return (kmer >> 2) | (c << (2 * k - 2));
}

uint64_t kmer2d(char *kmer, uint8_t k) {
  uint64_t kmer_d = 0;
  for (uint8_t i = 0; i < k; ++i)
    kmer_d = (kmer_d << 2) | (to_int[kmer[i]] - 1);
  // std::cerr << kmer << " " << kmer_d << " " << counter << std::endl;
  return kmer_d;
}

std::map<uint64_t, uint16_t> extract_kmers_from_db(char *kmc_path,
                                                   uint32_t *klen,
                                                   uint16_t min_w,
                                                   uint16_t max_w) {
  std::map<uint64_t, uint16_t> KMERS;
  CKMCFile kmer_db;
  if (!kmer_db.OpenForListing(kmc_path)) {
    std::cerr << "ERROR: cannot open " << kmc_path << std::endl;
    exit(1);
  }
  uint32_t mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64_t tot_kmers, max_c;
  kmer_db.Info(*klen, mode, min_counter, pref_len, sign_len, min_c, max_c,
               tot_kmers);
  CKmerAPI kmer_obj(*klen);
  char kmer[*klen + 1];
  std::cerr << "Parsing KMC database.." << std::endl;
  while (kmer_db.ReadNextKmer(kmer_obj, counter)) {
    if (counter >= min_w && counter <= max_w) {
      kmer_obj.to_string(kmer);
      KMERS[kmer2d(kmer, *klen)] = counter;
    }
  }
  return KMERS;
}

float average(std::map<uint16_t, uint32_t> hist) {
  // TODO
  (void)(hist); // suppress unused parameter warning
  return 0.0;
}

float median(std::map<uint16_t, uint32_t> hist) {
  // CHECKME
  float med = 0.0;
  uint64_t total = 0;
  for (const auto &x : hist)
    total += x.first * x.second;
  uint64_t sum = 0;
  for (const auto &x : hist) {
    sum += x.first * x.second;
    if (sum > total / 2) {
      med = x.first;
      break;
    }
  }
  return med;
}

float pass1(char *kmc_path, BF *bf, uint32_t *klen, uint16_t min_w,
            uint16_t max_w, bool do_avg) {
  std::map<uint16_t, uint32_t> khist;
  CKMCFile kmer_db;
  if (!kmer_db.OpenForListing(kmc_path)) {
    std::cerr << "ERROR: cannot open " << kmc_path << std::endl;
    return -1;
  }

  uint32 mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64 tot_kmers, max_c;
  kmer_db.Info(*klen, mode, min_counter, pref_len, sign_len, min_c, max_c,
               tot_kmers);
  CKmerAPI kmer_obj(*klen);
  char kmer[*klen + 1];
  while (kmer_db.ReadNextKmer(kmer_obj, counter)) {
    if (counter >= min_w && counter <= max_w) {
      kmer_obj.to_string(kmer);
      bf->add_key(kmer2d(kmer, *klen));
      ++khist[counter];
    }
  }
  return do_avg ? average(khist) : median(khist);
}

void pass2(char *kmc_path, BF *bf, uint32_t *klen, uint16_t min_w,
           uint16_t max_w) {
  CKMCFile kmer_db;
  if (!kmer_db.OpenForListing(kmc_path)) {
    std::cerr << "ERROR: cannot open " << kmc_path << std::endl;
    exit(1);
  }

  uint32 mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64 tot_kmers, max_c;
  kmer_db.Info(*klen, mode, min_counter, pref_len, sign_len, min_c, max_c,
               tot_kmers);
  CKmerAPI kmer_obj(*klen);
  char kmer[*klen + 1];
  while (kmer_db.ReadNextKmer(kmer_obj, counter)) {
    if (counter >= min_w && counter <= max_w) {
      kmer_obj.to_string(kmer);
      bf->set_count(kmer2d(kmer, *klen), counter);
    }
  }
}

int main(int argc, char *argv[]) {
  spdlog::set_default_logger(spdlog::stderr_color_st("stderr"));

  uint16_t wsize = 1000; // window size
  uint16_t min_w = 0;
  uint16_t max_w = -1;
  uint64_t bf_size = ((uint64_t)0b1 << 33); // 1GB
  bool do_avg = false;
  int a;
  while ((a = getopt(argc, argv, "w:m:M:b:avh")) >= 0) {
    switch (a) {
    case 'w':
      wsize = atoi(optarg);
      continue;
    case 'm':
      min_w = atoi(optarg);
      continue;
    case 'M':
      max_w = atoi(optarg);
      continue;
    case 'a':
      do_avg = true;
      continue;
    case 'b':
      // Let's consider this as GB
      bf_size = atoi(optarg) * ((uint64_t)0b1 << 33);
      break;
    // case '@':
    //   threads = atoi(optarg);
    //   continue;
    case 'v':
      spdlog::set_level(spdlog::level::debug);
      continue;
    case 'h':
      std::cerr << USAGE_MESSAGE;
      return 0;
    default:
      std::cerr << USAGE_MESSAGE;
      return 1;
    }
  }

  if (argc - optind < 3) {
    std::cerr << USAGE_MESSAGE;
    return 1;
  }

  char *fa_path = argv[optind++];       // reference
  char *kcontrol_path = argv[optind++]; // control KMC database
  char *kcase_path = argv[optind];      // case KMC database

  spdlog::info("Allocating Bloom filters");
  BF control_bf(bf_size);
  BF case_bf(bf_size);

  // TODO: fail if k is different between databases
  uint32_t klen;

  spdlog::info("First pass on {}", kcontrol_path);
  float control_norm =
      pass1(kcontrol_path, &control_bf, &klen, min_w, max_w, do_avg);
  spdlog::info("First pass on {}", kcase_path);
  float case_norm = pass1(kcase_path, &case_bf, &klen, min_w, max_w, do_avg);

  spdlog::info("Switching mode on control BF");
  control_bf.switch_mode();
  spdlog::info("Fill ratio for control BF:\t{}", control_bf.get_ratio());
  spdlog::info("Switching mode on case BF");
  case_bf.switch_mode();
  spdlog::info("Fill ratio for case BF:\t{}", case_bf.get_ratio());

  spdlog::info("Second pass on {}", kcontrol_path);
  pass2(kcontrol_path, &control_bf, &klen, min_w, max_w);

  spdlog::info("Second pass on {}", kcase_path);
  pass2(kcase_path, &case_bf, &klen, min_w, max_w);

  char kmer[klen + 1];   // first kmer on sequence (plain)
  int l;                 // chromosome size
  uint64_t kmer_d = 0;   // kmer
  uint64_t rckmer_d = 0; // reverse and complemented kmer
  uint64_t ckmer_d = 0;  // canonical kmer
  uint8_t c;             // new character to append
  int p = 0;             // current position on chromosome
  int bin_p = 0;         // current position in bin

  std::vector<uint16_t> control_bin(wsize);
  std::vector<uint16_t> case_bin(wsize);
  std::vector<float> ratio_bin(wsize);

  gzFile fa = gzopen(fa_path, "r");
  kseq_t *seq = kseq_init(fa);
  while ((l = kseq_read(seq)) >= 0) {
    spdlog::info("Iterating over {}", seq->name.s);

    strncpy(kmer, seq->seq.s, klen);
    kmer_d = kmer2d(kmer, klen);
    rckmer_d = revcompl(kmer_d, klen);
    ckmer_d = std::min(kmer_d, rckmer_d);

    control_bin[bin_p] = control_bf.get_count(ckmer_d);
    case_bin[bin_p] = case_bf.get_count(ckmer_d);
    ratio_bin[bin_p] =
        (float)(case_bin[bin_p] == 0 ? 0.0001 : case_bin[bin_p] / case_norm) /
        (float)(control_bin[bin_p] == 0 ? 0.0001
                                        : control_bin[bin_p] / control_norm);

    // std::cerr << seq->name.s << "\t" << p + 1 << "\t" << control_bin[bin_p]
    // << "\t" << case_bin[bin_p] << "\t" << ratio_bin[bin_p] << std::endl;
    ++bin_p;

    for (p = klen; p < l; ++p) {
      c = to_int[seq->seq.s[p]] - 1; // A is 1 but it should be 0
      kmer_d = lsappend(kmer_d, c, klen);
      rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
      ckmer_d = std::min(kmer_d, rckmer_d);

      control_bin[bin_p] = control_bf.get_count(ckmer_d);
      case_bin[bin_p] = case_bf.get_count(ckmer_d);
      ratio_bin[bin_p] =
          (float)(case_bin[bin_p] == 0 ? 0.0001 : case_bin[bin_p] / case_norm) /
          (float)(control_bin[bin_p] == 0 ? 0.0001
                                          : control_bin[bin_p] / control_norm);

      // std::cerr << seq->name.s << "\t" << p - klen + 2 << "\t" <<
      // case_bin[bin_p] << "\t" << control_bin[bin_p] << "\t" <<
      // ratio_bin[bin_p] << std::endl;
      ++bin_p;
      if (bin_p >= wsize) {
        std::cout << seq->name.s << ":" << p - klen + 2 - wsize << "-"
                  << p - klen + 1 << "\t";
        std::cout << case_bin[0];
        for (int _ = 1; _ < wsize; ++_)
          std::cout << "," << case_bin[_] / case_norm;
        std::cout << "\t";
        std::cout << control_bin[0];
        for (int _ = 1; _ < wsize; ++_)
          std::cout << "," << control_bin[_] / control_norm;
        std::cout << "\t";

        std::cout << ratio_bin[0];
        for (int _ = 1; _ < wsize; ++_)
          std::cout << "," << ratio_bin[_];
        std::cout << "\t";
        size_t n = wsize / 2; // CHECKME: what about even wsize?
        nth_element(ratio_bin.begin(), ratio_bin.begin() + n, ratio_bin.end());
        std::cout << ratio_bin[n];
        bin_p = 0;
        std::cout << std::endl;
      }
    }
    // TODO: manage last short (<wsize) bin
    // while (p - klen < l - 1) {
    //   std::cout << seq->name.s << "\t" << p - klen + 2 << "\t" << 0
    //             << std::endl;
    //   ++p;
    // }
  }

  kseq_destroy(seq);
  gzclose(fa);

  return 0;
}
