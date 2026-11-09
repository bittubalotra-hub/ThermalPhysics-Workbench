# C++ reference implementation

The C++17 sample is a **standalone original teaching example** for a uniform 1D slab with equal fixed-temperature wall conditions. It uses a conservative cell-centered finite-volume discretization and guarded explicit Euler integration.

No libraries beyond the C++17 standard library are needed.

```bash
# From repository root (Linux/macOS/Windows with g++)
g++ -std=c++17 -O2 -Wall -Wextra -pedantic cpp/heat1d_explicit.cpp -o heat1d
./heat1d 80 8000 0.0005 cpp_profile.csv
python scripts/compare_cpp_python.py cpp_profile.csv
```

On Windows PowerShell, run the executable as `.\heat1d.exe 80 8000 0.0005 cpp_profile.csv`.

The cross-language verification compares every spatial cell's final temperature against the Python solver, using identical time stepping, boundary treatment and material parameters. The fixture values are *generic illustrative parameters*; this is not a model of any spacecraft or proprietary material.
