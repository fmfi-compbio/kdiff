#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <zlib.h>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "kmc_api/kmc_file.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;

static const char *const USAGE_MESSAGE =
    "Usage: MFCNV <reference.fa> <control_kmc_db> <case_kmc_db> \n"
    "      -m <INT>   minimum weight for kmers (default: 0)\n"
    "      -M <INT>   maximum weight for kmers (default: 65535)\n"
    "      -w <INT>   bin size (default: 1000)\n"
    "      -a         use average instead of median for normalization\n"
    // "      -@ <INT>   set threads (default: 1)\n"
    "      -v         verbose mode\n"
    "      -h         display this help and exit\n"
    "\n";

float average(map<uint16_t, uint32_t> hist) {
  // TODO
  // (void)(hist); // suppress unused parameter warning
  float total = 0.0;
  for (const auto &x : hist)
    total += x.first * x.second;
  return total / 248956422.0; // FIXME
}

float median(map<uint16_t, uint32_t> hist) {
  // CHECKME
  float med = 0.0;
  float total = 0.0;
  for (const auto &x : hist)
    total += x.second;
  uint64_t sum = 0;
  for (const auto &x : hist) {
    sum += x.second;
    if (sum >= ceil(total / 2)) {
      med = x.first;
      break;
    }
  }
  return med;
}

float get_norm_factor(char *kmc_path, uint16_t min_w, uint16_t max_w,
                      bool do_avg) {
  map<uint16_t, uint32_t> khist;
  CKMCFile kmer_db;
  if (!kmer_db.OpenForListing(kmc_path))
    return -1;

  kmer_db.SetMinCount(min_w);
  kmer_db.SetMaxCount(max_w);

  uint32 klen, mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64 tot_kmers, max_c;
  kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c,
               tot_kmers);
  CKmerAPI kmer_obj(klen);
  while (kmer_db.ReadNextKmer(kmer_obj, counter))
    ++khist[counter];
  return do_avg ? average(khist) : median(khist);
}

int main(int argc, char *argv[]) {
  spdlog::set_default_logger(spdlog::stderr_color_st("stderr"));

  uint16_t wsize = 1000; // window size
  uint16_t min_w = 0;
  uint16_t max_w = -1;
  bool do_avg = false;
  bool verbose = false;
  int a;
  while ((a = getopt(argc, argv, "w:b:m:M:@:avh")) >= 0) {
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
    case 'v':
      verbose = true;
      continue;
    case 'h':
      cerr << USAGE_MESSAGE;
      return 0;
    default:
      cerr << USAGE_MESSAGE;
      return 1;
    }
  }

  if (argc - optind < 3) {
    cerr << USAGE_MESSAGE;
    return 1;
  }

  char *fa_path = argv[optind++];    // reference
  char *kctrl_path = argv[optind++]; // control KMC database
  char *kcase_path = argv[optind];   // case KMC database

  if (verbose) {
    spdlog::set_level(spdlog::level::debug);
    spdlog::debug("Setting threads to 1 for verbose mode");
  }

  spdlog::info("Computing normalization factor on control db");
  float ctrl_norm = get_norm_factor(kctrl_path, min_w, max_w, do_avg);
  if (ctrl_norm < 0) {
    spdlog::critical("Something went wrong while reading {}", kctrl_path);
    return 1;
  }
  spdlog::info("Computing normalization factor on case db");
  float case_norm = get_norm_factor(kcase_path, min_w, max_w, do_avg);
  if (case_norm < 0) {
    spdlog::critical("Something went wrong while reading {}", kcase_path);
    return 1;
  }

  CKMCFile ctrl_db;
  vector<uint32_t> ctrl_counts;
  ctrl_db.OpenForRA(kctrl_path);
  ctrl_db.SetMinCount(min_w);
  ctrl_db.SetMaxCount(max_w);

  CKMCFile case_db;
  vector<uint32_t> case_counts;
  case_db.OpenForRA(kcase_path);
  case_db.SetMinCount(min_w);
  case_db.SetMaxCount(max_w);
  if (ctrl_db.KmerLength() != case_db.KmerLength()) {
    spdlog::critical("Control and case databases have different k ({} vs {})",
                     ctrl_db.KmerLength(), case_db.KmerLength());
    return 1;
  }
  vector<float> ratios(wsize);
  gzFile fa = gzopen(fa_path, "r");
  kseq_t *seq = kseq_init(fa);
  int l; // chromosome size
  int p; // position on chromosome
  char *binseq = (char *)malloc(wsize);
  while ((l = kseq_read(seq)) >= 0) {
    spdlog::info("Processing {}", seq->name.s);
    p = 0;
    while (p < l - wsize) {
      strncpy(binseq, seq->seq.s + p, wsize);
      binseq[wsize] = '\0';
      ctrl_db.GetCountersForRead(binseq, ctrl_counts);
      case_db.GetCountersForRead(binseq, case_counts);
      for (size_t i = 0; i < ratios.size(); ++i) {
        ratios[i] =
            (float)(case_counts[i] == 0 ? 0.0001 : case_counts[i] / case_norm) /
            (float)(ctrl_counts[i] == 0 ? 0.0001 : ctrl_counts[i] / ctrl_norm);
      }
      // CHECKME: what about even wsize?
      nth_element(ratios.begin(), ratios.begin() + wsize / 2, ratios.end());
      cout << seq->name.s << ":" << p << "-" << p + wsize << "\t"
           << ratios[wsize / 2] << endl;
      p += wsize;
    }
    // last part
    strncpy(binseq, seq->seq.s + p, l - p);
    binseq[l - p] = '\0';
    ctrl_db.GetCountersForRead(binseq, ctrl_counts);
    case_db.GetCountersForRead(binseq, case_counts);
    for (size_t i = 0; i < ratios.size(); ++i) {
      ratios[i] =
          (float)(case_counts[i] == 0 ? 0.0001 : case_counts[i] / case_norm) /
          (float)(ctrl_counts[i] == 0 ? 0.0001 : ctrl_counts[i] / ctrl_norm);
    }
    // CHECKME: what about even wsize?
    nth_element(ratios.begin(), ratios.begin() + wsize / 2, ratios.end());
    cout << seq->name.s << ":" << p << "-" << l << "\t" << ratios[wsize / 2]
         << endl;
  }
  kseq_destroy(seq);
  gzclose(fa);

  return 0;
}

int main_test(int argc, char *argv[]) {
  (void)(argc); // suppress unused parameter warning
  spdlog::set_default_logger(spdlog::stderr_color_st("stderr"));

  char *kmc_path = argv[1];

  CKMCFile kmer_db_ra;
  kmer_db_ra.OpenForRA(kmc_path);

  CKMCFile kmer_db;
  kmer_db.OpenForListing(kmc_path);
  uint32 klen, mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64 tot_kmers, max_c;
  kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c,
               tot_kmers);
  CKmerAPI kmer_obj(klen);
  char kmer[klen + 1];
  bool found;
  uint32_t counter_ra;
  while (kmer_db.ReadNextKmer(kmer_obj, counter)) {
    kmer_obj.to_string(kmer);
    found = kmer_db_ra.CheckKmer(kmer_obj, counter_ra);
    assert(counter == counter_ra);
    spdlog::info("{}: {} - {} {}", string(kmer), found, counter, counter_ra);
  }

  kmer_obj.from_string("AAC");
  found = kmer_db_ra.CheckKmer(kmer_obj, counter_ra);
  spdlog::info("AAC: {}", counter_ra);

  kmer_obj.from_string("GTT");
  found = kmer_db_ra.CheckKmer(kmer_obj, counter_ra);
  spdlog::info("GTT: {}", counter_ra);

  vector<uint32_t> counts;
  kmer_db_ra.GetCountersForRead("AACNGTT", counts);
  for (const auto &c : counts)
    cout << c << " ";
  cout << endl;

  return 0;
}
