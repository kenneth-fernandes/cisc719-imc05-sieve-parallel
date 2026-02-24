// ============================================================
// sieve_openmp.cpp — OpenMP Parallel Segmented Sieve (Beginner-Friendly)
// Shared-memory parallel version
// PCAM: Range partitioned across threads; base primes shared read-only
// Output: N=<N> threads=<T> count=<count> time_sec=<time>
// ============================================================

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <omp.h>

using namespace std;
using namespace chrono;

// Segment size (numbers, not just odds)
// You can tune this later for performance experiments
static const long long SEG_SIZE = 1 << 20; // ~1M numbers per segment

// Toggle verbose prints here (true = learning/debugging, false = benchmarking)
static const bool VERBOSE = true;

// ------------------------------------------------------------
// Step 1: Sequential simple sieve up to sqrt(N)
// Returns list of base primes
// ------------------------------------------------------------
vector<int> simple_sieve(int limit) {
    if (VERBOSE) {
        cout << "\nStarting simple_sieve(limit = " << limit << ")" << endl;
    }

    vector<bool> is_prime(limit + 1, true);

    if (limit >= 0) is_prime[0] = false;
    if (limit >= 1) is_prime[1] = false;

    for (int i = 2; (long long)i * i <= limit; i++) {
        if (is_prime[i]) {
            if (VERBOSE && i <= 20) {
                cout << "Found base prime " << i
                     << ". Marking multiples starting at " << (1LL * i * i) << "." << endl;
            }

            for (long long j = 1LL * i * i; j <= limit; j += i) {
                is_prime[(int)j] = false;
            }
        }
    }

    vector<int> primes;
    for (int i = 2; i <= limit; i++) {
        if (is_prime[i]) primes.push_back(i);
    }

    if (VERBOSE) {
        cout << "simple_sieve complete. Number of base primes = " << primes.size() << endl;
        cout << "First few base primes: ";
        for (size_t i = 0; i < primes.size() && i < 15; i++) {
            cout << primes[i] << (i + 1 < primes.size() && i < 14 ? ", " : "");
        }
        cout << endl;
    }

    return primes;
}

// ------------------------------------------------------------
// Step 2: OpenMP parallel segmented sieve over [2..N]
// - Base primes computed once (sequential)
// - Segments distributed across threads
// - Each thread uses its own local segment buffer
// ------------------------------------------------------------
long long sieve_openmp(long long N, int num_threads) {
    if (N < 2) return 0;

    omp_set_num_threads(num_threads);

    if (VERBOSE) {
        cout << "\nStarting sieve_openmp(N = " << N
             << ", num_threads = " << num_threads << ")" << endl;
    }

    // Count prime 2 separately (odd-only segment representation later)
    long long total_count = 1; // prime = 2

    // Base primes up to sqrt(N), computed once
    int limit = (int)floor(sqrt((long double)N));
    vector<int> base_primes = simple_sieve(limit);

    // Number of segments for [3..N]
    long long first_value = 3;
    if (first_value > N) return total_count;

    long long total_numbers = N - first_value + 1;
    long long num_segments = (total_numbers + SEG_SIZE - 1) / SEG_SIZE;

    if (VERBOSE) {
        cout << "Computed parameters:" << endl;
        cout << "  floor(sqrt(N)) = " << limit << endl;
        cout << "  first_value = " << first_value << endl;
        cout << "  total_numbers = " << total_numbers << endl;
        cout << "  SEG_SIZE = " << SEG_SIZE << endl;
        cout << "  num_segments = " << num_segments << endl;
        cout << "Starting OpenMP parallel loop over segments..." << endl;
    }

    // Parallelize over segments
    // dynamic scheduling helps load balancing if segments vary
    #pragma omp parallel for schedule(dynamic) reduction(+:total_count)
    for (long long seg_id = 0; seg_id < num_segments; seg_id++) {
        int tid = omp_get_thread_num();

        long long low  = first_value + seg_id * SEG_SIZE;
        long long high = min(low + SEG_SIZE - 1, N);

        // Make low odd (we store only odd numbers)
        if (low % 2 == 0) low++;

        // If segment became invalid after adjusting, skip
        if (low > high) {
            if (VERBOSE) {
                #pragma omp critical
                {
                    cout << "Thread " << tid
                         << " skipped segment " << seg_id
                         << " because low > high after odd adjustment." << endl;
                }
            }
            continue;
        }

        // Number of odd values in [low..high]
        long long odd_count = ((high - low) / 2) + 1;

        if (VERBOSE) {
            #pragma omp critical
            {
                cout << "Thread " << tid
                     << " processing segment " << seg_id
                     << " with range [" << low << ", " << high << "]"
                     << " (odd_count = " << odd_count << ")" << endl;
            }
        }

        // Local segment buffer (thread-local because created inside loop)
        vector<bool> segment((size_t)odd_count, true);

        // Mark composites in this segment using base primes
        int printed_mark_actions = 0; // limit prints per segment

        for (int p : base_primes) {
            if (p == 2) continue; // segment stores only odd numbers

            long long p64 = (long long)p;
            long long p2  = p64 * p64;

            // If p^2 > high, no need to continue for this segment
            if (p2 > high) break;

            // First multiple of p within [low..high]
            long long start = max(p2, ((low + p64 - 1) / p64) * p64);

            // Ensure start is odd (segment stores odd numbers only)
            if ((start % 2) == 0) start += p64;

            if (VERBOSE && seg_id < 2 && p <= 19) {
                #pragma omp critical
                {
                    cout << "Thread " << tid
                         << " | Segment " << seg_id
                         << " | Using prime p = " << p
                         << " | p^2 = " << p2
                         << " | first odd multiple = " << start
                         << " | step = " << (2 * p64) << endl;
                }
            }

            // Mark odd multiples only (step = 2p)
            for (long long x = start; x <= high; x += 2 * p64) {
                long long idx = (x - low) / 2;
                segment[(size_t)idx] = false;

                // Print only a few actual mark actions (to avoid huge output)
                if (VERBOSE && seg_id < 1 && printed_mark_actions < 12) {
                    #pragma omp critical
                    {
                        cout << "Thread " << tid
                             << " marked composite number " << x
                             << " (segment index " << idx << ") using prime " << p << endl;
                    }
                    printed_mark_actions++;
                }
            }
        }

        // Count remaining primes in this segment
        long long local_count = 0;
        for (long long i = 0; i < odd_count; i++) {
            if (segment[(size_t)i]) {
                long long num = low + 2 * i;
                if (num <= N) local_count++;

                // Print first few primes found in first segment
                if (VERBOSE && seg_id == 0 && local_count <= 10) {
                    #pragma omp critical
                    {
                        cout << "Thread " << tid
                             << " found surviving prime candidate: " << num << endl;
                    }
                }
            }
        }

        if (VERBOSE) {
            #pragma omp critical
            {
                cout << "Thread " << tid
                     << " finished segment " << seg_id
                     << " with local prime count = " << local_count << endl;
            }
        }

        total_count += local_count;
    }

    if (VERBOSE) {
        cout << "OpenMP loop finished. Total prime count (including 2) = "
             << total_count << endl;
    }

    return total_count;
}

int main(int argc, char* argv[]) {
    long long N = (argc > 1) ? atoll(argv[1]) : 100000000LL;
    int threads = (argc > 2) ? atoi(argv[2])  : 4;

    if (VERBOSE) {
        cout << "Program started." << endl;
        cout << "Input N = " << N << ", threads = " << threads << endl;
    }

    auto t0 = high_resolution_clock::now();
    long long count = sieve_openmp(N, threads);
    auto t1 = high_resolution_clock::now();

    double elapsed = duration<double>(t1 - t0).count();

    if (VERBOSE) {
        cout << "Execution finished in " << elapsed << " seconds." << endl;
    }

    printf("N=%lld threads=%d count=%lld time_sec=%.6f\n", N, threads, count, elapsed);

    return 0;
}
