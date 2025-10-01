# ProtSpaM (Sequential Version in C++)

This is a **sequential** adaptation of [ProtSpaM](https://github.com/jschellh/ProtSpaM), a tool for estimating phylogenetic distances between proteins based on **spaced-word matches**.  

The original version includes parallelism with **OpenMP**.  
This version removes all parallelization to run in a **strictly sequential** manner, while keeping the same core functionality.  

---

## ✨ Changes Compared to the Original Version
- A new file `main_sequential.cpp` was created instead of the original `main.cpp`.  
- Removed all **OpenMP** dependencies (`#include <omp.h>`, `omp_get_wtime`, `omp_set_num_threads`, `#pragma omp parallel for`).  
- The timing system now uses `std::chrono` instead of OpenMP functions.  
- The `-t` parameter (number of threads) was removed since it is not needed in sequential code.  

---

## 🔧 Compilation
Make sure you have a compiler compatible with **C++11 or later** (e.g., `g++` or `clang++`).  

```bash
g++ -std=c++11 -O3 -o protspam_sequential main_sequential.cpp *.cpp
```

---

## 📂 Usage

The program accepts the same parameters as the original version.

Example with multiple input files:

```bash
./protspam_sequential -w 6 -d 40 -m 5 -l input_files.txt -o DMat
```

---
## ⚙️ Main Options

```bash
-w <int> : Pattern weight (default: 6)

-d <int> : Number of “don’t-care” positions (default: 40)

-s <int> : Minimum score to consider a spaced-word match as homologous (default: 0)

-m <int> : Number of patterns (default: 5)

-o <file> : Name of the output file containing the distance matrix (default: DMat)

-l <file> : List of input files in multifasta format

-p <file> : Load a predefined pattern set

-z : Save generated patterns into patterns.txt

-r : Generate individual scores for each sequence pair (spamograms)
```

---
## 📊 Output
The program produces the following results:

- Distance matrix → output file (DMat by default).
- Generated patterns → saved in patterns.txt if the -z option is used.
- Spamograms (scores) → saved in the scores/ directory if the -r option is used.


### Contact:
ana.izaguirre@udc.es
