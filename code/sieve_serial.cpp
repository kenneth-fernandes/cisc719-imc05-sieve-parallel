// ============================================================
// sieve_serial.cpp — Beginner-Friendly Serial Segmented Sieve (Verbose Version)
// - Uses Sieve of Eratosthenes
// - Uses segmentation (processes in chunks)
// - Uses odd-only optimization (skips even numbers except 2)
// - Includes print statements (toggle VERBOSE=true/false)
// Output: N=<N> count=<count> time_sec=<time>
// ============================================================

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <algorithm>

using namespace std;
using namespace chrono;

// Toggle verbose prints here
static const bool VERBOSE = true;

// Segment size = how many numbers we process at once
static const long long SEG_SIZE = 1 << 20; // about 1 million numbers per segment

// ------------------------------------------------------------
// Step 1: Build base primes up to sqrt(N) using simple sieve
// ------------------------------------------------------------
vector<int> simple_sieve(int limit) {
    if (VERBOSE) {
        cout << "\nEntering simple_sieve(limit = " << limit << ")" << endl;
    }

    // is_prime[x] = true means x is currently assumed prime
    vector<bool> is_prime(limit + 1, true);

    // 0 and 1 are not prime
    if (limit >= 0) is_prime[0] = false;
    if (limit >= 1) is_prime[1] = false;

    // Standard Sieve of Eratosthenes
    for (int i = 2; (long long)i * i <= limit; i++) {
        if (is_prime[i]) {
            if (VERBOSE && i <= 20) {
                cout << "Found base prime " << i
                     << ". Marking multiples starting from " << (1LL * i * i)
                     << " with step " << i << "." << endl;
            }

            // Start from i*i (optimization)
            for (long long j = 1LL * i * i; j <= limit; j += i) {
                is_prime[(int)j] = false;
            }
        }
    }

    // Collect primes into a list
    vector<int> primes;
    for (int i = 2; i <= limit; i++) {
        if (is_prime[i]) {
            primes.push_back(i);
        }
    }

    if (VERBOSE) {
        cout << "Total base primes found up to sqrt(N): " << primes.size() << endl;
        cout << "First few base primes: ";
        for (size_t i = 0; i < primes.size() && i < 15; i++) {
            cout << primes[i];
            if (i + 1 < primes.size() && i < 14) cout << ", ";
        }
        cout << endl;
    }

    return primes;
}

// ------------------------------------------------------------
// Step 2: Segmented sieve over [2..N]
// We process the full range in smaller chunks (segments)
// ------------------------------------------------------------
long long sieve_serial(long long N) {
    if (VERBOSE) {
        cout << "\nStarting sieve_serial(N = " << N << ")" << endl;
        cout << "Segment size = " << SEG_SIZE << endl;
    }

    if (N < 2) {
        if (VERBOSE) {
            cout << "N is less than 2, so there are no primes." << endl;
        }
        return 0;
    }

    // Count prime 2 separately (because we skip evens later)
    long long prime_count = 1; // prime number 2

    // We only need base primes up to sqrt(N)
    int limit = (int)floor(sqrt((long double)N));
    if (VERBOSE) {
        cout << "floor(sqrt(N)) = " << limit << endl;
    }

    vector<int> base_primes = simple_sieve(limit);

    if (VERBOSE) {
        cout << "\nStarting segmented sieve over odd numbers in range [3.." << N << "]" << endl;
    }

    long long segment_number = 0;

    // Process [3..N] in segments
    // We will only store/check odd numbers in each segment
    for (long long low = 3; low <= N; low += SEG_SIZE) {
        long long high = min(low + SEG_SIZE - 1, N);
        segment_number++;

        // Make sure segment starts on an odd number
        if (low % 2 == 0) low++;
        if (high < low) continue;

        // Number of odd values in this segment
        long long odd_count = ((high - low) / 2) + 1;

        if (VERBOSE) {
            cout << "\nSegment " << segment_number
                 << ": range [" << low << ", " << high << "]"
                 << " | odd_count = " << odd_count << endl;
        }

        // segment[i] corresponds to number = low + 2*i
        // true = prime candidate, false = composite
        vector<bool> segment((size_t)odd_count, true);

        int printed_prime_steps = 0;  // limit prints per segment

        // Mark composites using base primes
        for (int p : base_primes) {
            if (p == 2) continue; // segment stores odd numbers only

            long long p64 = (long long)p;
            long long p2 = p64 * p64;

            // If p^2 > high, no need to mark multiples for this segment
            if (p2 > high) {
                if (VERBOSE) {
                    cout << "Stopping marking for this segment at p = " << p
                         << " because p^2 = " << p2
                         << " is greater than segment high = " << high << endl;
                }
                break;
            }

            // Find first multiple of p inside [low..high]
            long long start = max(p2, ((low + p64 - 1) / p64) * p64);

            // Make start odd, because segment stores only odd numbers
            if (start % 2 == 0) {
                start += p64;
            }

            if (VERBOSE && printed_prime_steps < 8) {
                cout << "Using base prime p = " << p
                     << " | p^2 = " << p2
                     << " | first odd multiple in segment = " << start
                     << " | step = " << (2 * p64) << endl;
                printed_prime_steps++;
            }

            // Mark odd multiples of p in this segment
            long long marks_for_p = 0;
            for (long long x = start; x <= high; x += 2 * p64) {
                long long idx = (x - low) / 2;
                if (idx >= 0 && idx < odd_count) {
                    if (segment[(size_t)idx]) { // count only first-time changes
                        segment[(size_t)idx] = false;
                        marks_for_p++;
                    }
                }
            }

            if (VERBOSE && p <= 19) {
                cout << "Prime " << p << " marked " << marks_for_p
                     << " odd composite numbers in this segment." << endl;
            }
        }

        // Count remaining true values in this segment (these are primes)
        long long segment_prime_count = 0;
        if (VERBOSE) {
            cout << "Counting surviving prime candidates in segment " << segment_number << "..." << endl;
        }

        for (long long i = 0; i < odd_count; i++) {
            if (segment[(size_t)i]) {
                long long num = low + 2 * i;
                if (num <= N) {
                    prime_count++;
                    segment_prime_count++;

                    // Print only first few primes in segment
                    if (VERBOSE && segment_prime_count <= 10) {
                        cout << "  Prime found in this segment: " << num << endl;
                    }
                }
            }
        }

        if (VERBOSE) {
            cout << "Segment " << segment_number
                 << " complete. Primes in this segment = " << segment_prime_count
                 << " | Running total = " << prime_count << endl;
        }
    }

    if (VERBOSE) {
        cout << "\nFinished sieve_serial. Final prime count = " << prime_count << endl;
    }

    return prime_count;
}

int main(int argc, char* argv[]) {
    long long N = (argc > 1) ? atoll(argv[1]) : 1000000LL;

    if (VERBOSE) {
        cout << "Program started." << endl;
        cout << "Input N = " << N << endl;
    }

    auto t0 = high_resolution_clock::now();
    long long count = sieve_serial(N);
    auto t1 = high_resolution_clock::now();

    double elapsed = duration<double>(t1 - t0).count();

    if (VERBOSE) {
        cout << "Execution time = " << elapsed << " seconds" << endl;
    }

    // Machine-readable output for benchmark parser
    printf("N=%lld count=%lld time_sec=%.6f\n", N, count, elapsed);

    return 0;
}
