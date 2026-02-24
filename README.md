# IMC05 — Parallel Sieve of Eratosthenes

**Course:** CISC 719 — Contemporary Computing Systems Modeling Algorithms (CCSM)  
**Program:** Ph.D. in Computational Sciences  
**University:** Harrisburg University of Science and Technology  
**Instructor:** Professor Majid Shaalan, Ph.D.  
**Student:** Kenneth Peter Fernandes  
**Semester:** Spring 2026  

---

## About

This assignment implements and benchmarks the Sieve of Eratosthenes using three parallel strategies, applying the PCAM (Partitioning, Communication, Agglomeration, Mapping) model to each design. The goal is to compare performance across a serial baseline, a shared-memory parallel implementation, and a GPU implementation across problem sizes up to N = 10⁸.

The assignment also includes a real-world extension applying sieve-generated prime data to cryptographic primality analysis — specifically identifying Sophie Germain primes and safe primes used in Diffie-Hellman and DSA protocols.

**Implementations:**
- Serial C++ — segmented, odd-only, bit-packed baseline
- OpenMP C++ — shared-memory parallel, `schedule(dynamic)`, 1/2/4 threads
- Numba CUDA (Python) — GPU parallel using Numba JIT kernels

---

## Project Structure

```
cisc719-imc05-sieve-parallel/
├── notebooks/
│   └── imc05_sieve_parallel.ipynb   # Main notebook — all implementations, benchmarks, crypto extension
├── code/
│   ├── sieve_serial.cpp             # Serial C++ sieve (standalone copy)
│   └── sieve_openmp.cpp             # OpenMP C++ sieve (standalone copy)
├── docs/
│   ├── src/
│   │   ├── report.tex               # Full technical report (LaTeX)
│   │   ├── memo.tex                 # Reflective memo (LaTeX)
│   │   └── references.bib          # Bibliography (BibTeX)
│   ├── pdf/
│   │   ├── report.pdf               # Compiled report
│   │   └── memo.pdf                 # Compiled memo
│   └── results/
│       ├── results.csv              # Benchmark timing data (all runs)
│       ├── table7_sg_safe_counts.csv # Sophie Germain / safe prime counts
│       ├── Benchmarks.png           # Combined 6-panel benchmark figure
│       └── performance_analysis.png # Additional performance plot
├── ai-log/
│   └── ai_usage_log.csv             # AI tool usage log (required by assignment)
└── README.md
```

---

## Running the Notebook on Google Colab

The entire implementation runs in a single Jupyter notebook on Google Colab (T4 GPU runtime).

**Steps:**

1. Open [Google Colab](https://colab.research.google.com)
2. Upload `notebooks/imc05_sieve_parallel.ipynb` via **File → Upload notebook**
3. Enable GPU: **Runtime → Change runtime type → T4 GPU**
4. Run all cells in order: **Runtime → Run all**

**Runtime environment used:**
- CPU: Intel Xeon, ~2 physical cores
- GPU: NVIDIA T4, 16 GB GDDR6, 2,560 CUDA cores
- RAM: 12 GB
- Python: 3.12
- Numba: 0.60.0
- CUDA: 13.0

**Notebook sections:**

| Section | Description |
|---------|-------------|
| 0 | Title and student info |
| 1 | Environment setup (GPU/CPU info) |
| 2 | Problem overview and PCAM design table |
| 3 | Benchmark configuration |
| 4 | Shared utilities |
| 5 | Serial C++ sieve (compile + run) |
| 6 | OpenMP C++ sieve (compile + run) |
| 7 | Numba GPU kernel |
| 8 | Correctness validation |
| 9 | Benchmark runner → exports `results.csv` |
| 10 | Performance plots → exports `Benchmarks.png` |
| 11 | Crypto extension — Sophie Germain and safe primes |
| 12 | Reflective memo draft |
| 13 | Export artifacts to Google Drive |

---

## Generating the PDF Documents

### Requirements

Requires a TeX Live installation with `pdflatex` and `biber`.

**macOS (Homebrew):**
```bash
brew install --cask mactex
```

**Ubuntu / Debian:**
```bash
sudo apt install texlive-full biber
```

**Verify installation:**
```bash
pdflatex --version
biber --version
```

### Build Command

Run from the project root:

```bash
cd docs/src && \
pdflatex -interaction=nonstopmode report.tex && \
biber report && \
pdflatex -interaction=nonstopmode report.tex && \
pdflatex -interaction=nonstopmode report.tex && \
mv report.pdf ../pdf/report.pdf && \
pdflatex -interaction=nonstopmode memo.tex && \
mv memo.pdf ../pdf/memo.pdf
```

Output PDFs are written to `docs/pdf/`. All auxiliary files (`.aux`, `.bbl`, `.log`, `.toc`, etc.) remain in `docs/src/`.

**Why three pdflatex passes for the report?**
The first pass generates the `.aux` file. `biber` reads it to resolve citations and writes the `.bbl` file. The second pass incorporates the bibliography. The third pass resolves any remaining cross-references (table of contents, figure labels).

The memo has no bibliography and requires only one pass.
