#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <zlib.h>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "kmc_api/kmc_file.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;

static const char *const USAGE_MESSAGE =
    "Usage: MFCNV <reference.fa> <reference_kmc_db> <control_kmc_db> "
    "<case_kmc_db> \n"
    "      -m <INT>     minimum weight for kmers (default: 0)\n"
    "      -M <INT>     maximum weight for kmers (default: 65535)\n"
    "      -w <INT>     bin size (default: 1000)\n"
    "      -z <FLOAT>   filter windows with less than this percentage of\n"
    "                   non-zeros counts (value in [0,1], default: 0.0,\n"
    "                   i.e., no filtering)\n"
    "      -r <BED>     regions file in BED format (default: None, i.e., full "
    "genome)\n"
    "      -a           use average instead of median for normalization\n"
    // "      -@ <INT>   set threads (default: 1)\n"
    "      -v           verbose mode\n"
    "      -h           display this help and exit\n"
    "\n";

float average(map<uint16_t, uint32_t> hist, uint64_t &tot_refkmers) {
  float total = 0.0;
  for (const auto &x : hist)
    total += x.first * x.second;
  return total / tot_refkmers;
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
                      bool do_avg, uint64_t &tot_refkmers) {
  map<uint16_t, uint32_t> khist;
  CKMCFile kmer_db;
  if (!kmer_db.OpenForListing(kmc_path))
    return -1;

  kmer_db.SetMinCount(min_w);
  kmer_db.SetMaxCount(max_w);

  uint32_t klen, mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64_t tot_kmers, max_c;
  kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c,
               tot_kmers);
  CKmerAPI kmer_obj(klen);
  while (kmer_db.ReadNextKmer(kmer_obj, counter))
    ++khist[counter];
  kmer_db.Close();
  return do_avg ? average(khist, tot_refkmers) : median(khist);
}

float get_score(char *binseq, int wsize, CKMCFile &ctrl_db,
                const float &ctrl_norm, float eps1, CKMCFile &case_db,
                const float &case_norm, float eps2, vector<float> &ratios,
                float nzt) {
  vector<uint32_t> ctrl_counts;
  vector<uint32_t> case_counts;
  ctrl_db.GetCountersForRead(binseq, ctrl_counts);
  // Check for 0s in the control
  int ctrl_zeros = 0;
  for (const auto &w : ctrl_counts)
    ctrl_zeros += w == 0;
  float wscore = -1.0;
  if (ctrl_zeros <= (1.0 - nzt) * wsize) {
    case_db.GetCountersForRead(binseq, case_counts);
    for (size_t i = 0;
         i < case_counts.size() /* wsize - ctrl_db.KmerLength() + 1*/; ++i) {
      ratios[i] =
          (float)(case_counts[i] == 0 ? eps2
                                      : (float)case_counts[i] / case_norm) /
          (float)(ctrl_counts[i] == 0 ? eps1
                                      : (float)ctrl_counts[i] / ctrl_norm);
    }

    nth_element(ratios.begin(), ratios.begin() + wsize / 2, ratios.end());
    if (wsize % 2 == 0) {
      const auto x = max_element(ratios.begin(), ratios.begin() + (wsize / 2));
      wscore = (ratios[wsize / 2] + *x) / 2;
    } else {
      wscore = ratios[wsize / 2];
    }
  }
  return wscore;
}

int main(int argc, char *argv[]) {
  spdlog::set_default_logger(spdlog::stderr_color_st("stderr"));

  uint16_t wsize = 1000;  // window size
  uint16_t min_w = 0;     // minimum weight for kmers
  uint16_t max_w = -1;    // maximum weight for kmers
  float nzt = 0.0;        // filtering ratio for 0s in a region
  string regions_fn = ""; // list of regions to analyze
  bool do_avg = false;
  bool verbose = false;
  int a;
  while ((a = getopt(argc, argv, "w:b:m:M:z:r:avh")) >= 0) {
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
    case 'z':
      nzt = atof(optarg);
      continue;
    case 'r':
      regions_fn = optarg;
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

  if (argc - optind < 4) {
    cerr << USAGE_MESSAGE;
    return 1;
  }

  char *fa_path = argv[optind++];    // reference
  char *kref_path = argv[optind++];  // reference KMC database
  char *kctrl_path = argv[optind++]; // control KMC database
  char *kcase_path = argv[optind];   // case KMC database

  if (verbose) {
    spdlog::set_level(spdlog::level::debug);
    // spdlog::debug("Setting threads to 1 for verbose mode");
  }

  spdlog::info("Opening reference db");
  CKMCFile kmer_db;
  if (!kmer_db.OpenForRA(kref_path))
    return -1;
  uint64_t tot_refkmers = kmer_db.KmerCount();
  kmer_db.Close();

  spdlog::info("Computing normalization factor on control db");
  float ctrl_norm =
      get_norm_factor(kctrl_path, min_w, max_w, do_avg, tot_refkmers);
  if (ctrl_norm < 0) {
    spdlog::critical("Something went wrong while reading {}", kctrl_path);
    return 1;
  }
  if (verbose)
    spdlog::debug("Normalization factor on control db: {}", ctrl_norm);

  spdlog::info("Computing normalization factor on case db");
  float case_norm =
      get_norm_factor(kcase_path, min_w, max_w, do_avg, tot_refkmers);
  if (case_norm < 0) {
    spdlog::critical("Something went wrong while reading {}", kcase_path);
    return 1;
  }
  if (verbose)
    spdlog::debug("Normalization factor on case db: {}", case_norm);

  // CHECKME
  float eps1 = 0.0001 * ctrl_norm / (ctrl_norm + case_norm);
  float eps2 = 0.0001 * case_norm / (ctrl_norm + case_norm);

  CKMCFile ctrl_db;
  ctrl_db.OpenForRA(kctrl_path);
  ctrl_db.SetMinCount(min_w);
  ctrl_db.SetMaxCount(max_w);

  CKMCFile case_db;
  case_db.OpenForRA(kcase_path);
  case_db.SetMinCount(min_w);
  case_db.SetMaxCount(max_w);
  if (ctrl_db.KmerLength() != case_db.KmerLength()) {
    spdlog::critical("Control and case databases have different k ({} vs {})",
                     ctrl_db.KmerLength(), case_db.KmerLength());
    return 1;
  }

  map<string, vector<pair<int, int>>> regions;
  if (regions_fn != "") {
    ifstream infile(regions_fn);
    string line;
    if (infile.is_open()) {
      while (getline(infile, line)) {
        istringstream iss(line);
        string name;
        int start, end;
        if (!(iss >> name >> start >> end)) {
          spdlog::critical("Error at line {}\t{}\t{}", name, start, end);
          return 1;
        }
        regions[name].push_back(make_pair(start, end));
      }
    } else {
      spdlog::critical("Can't open regions file");
      return 1;
    }
  }

  vector<float> ratios(wsize - ctrl_db.KmerLength() + 1);
  gzFile fa = gzopen(fa_path, "r");
  kseq_t *seq = kseq_init(fa);
  int l; // chromosome size
  int p; // position on chromosome
  char *binseq = (char *)malloc((wsize + 1) * sizeof(char));
  int binseq_c = wsize;
  float wscore;
  while ((l = kseq_read(seq)) >= 0) {
    if (regions.empty()) {
      spdlog::info("Processing {}", seq->name.s);
      p = 0;
      while (p < l - wsize) {
        strncpy(binseq, seq->seq.s + p, wsize);
        binseq[wsize] = '\0';
        wscore = get_score(binseq, wsize, ctrl_db, ctrl_norm, eps1, case_db,
                           case_norm, eps2, ratios, nzt);
        cout << seq->name.s << ":" << p << "-" << p + wsize << "\t" << wscore
             << endl;
        p += wsize;
      }

      // last window with length < wsize
      strncpy(binseq, seq->seq.s + p, l - p);
      binseq[l - p] = '\0';
      ratios.resize(l - p);
      wscore = get_score(binseq, l - p, ctrl_db, ctrl_norm, eps1, case_db,
                         case_norm, eps2, ratios, nzt);
      cout << seq->name.s << ":" << p << "-" << l << "\t" << wscore << endl;
    } else {
      // Parse regions
      if (regions.find(seq->name.s) != regions.end()) {
        spdlog::info("Processing {} ({} regions)", seq->name.s,
                     regions.at(seq->name.s).size());
        for (const pair<int, int> &region : regions.at(seq->name.s)) {
          wsize = region.second - region.first;
          if (wsize > binseq_c) {
            spdlog::debug("Reallocating bin sequence from {} to {}", binseq_c,
                          wsize);
            binseq = (char *)realloc(binseq, wsize * 2);
            binseq_c = wsize * 2;
          }
          if (wsize > ratios.size())
            ratios.reserve(wsize - ctrl_db.KmerLength() + 1);

          strncpy(binseq, seq->seq.s + region.first, wsize);
          binseq[wsize] = '\0';
          wscore = get_score(binseq, wsize, ctrl_db, ctrl_norm, eps1, case_db,
                             case_norm, eps2, ratios, nzt);
          cout << seq->name.s << ":" << region.first << "-" << region.second
               << "\t" << wscore << endl;
        }
      }
    }
  }
  kseq_destroy(seq);
  gzclose(fa);
  ctrl_db.Close();
  case_db.Close();
  free(binseq);

  return 0;
}
