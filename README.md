This is the Fortran source code in the Julia package [FuzzifiED](https://github.com/FuzzifiED/FuzzifiED.jl) intended for exact diagonalisation (ED) calculations on fuzzy sphere. 

The code can be compiled by `make`. Three options are required : 
- `LIB_DIR` : the directory of the product libraries, `.` by default ; 
- `DLEXT` : the extension of the shared library on the target system, `so` for Linux (by default), `dll` for Windows, and `dylib` for macOS ;
- `NBITS` : the number of bits of the target architecture, `64` (by default) or `32`. 

This package is developed by Zheng Zhou (周正) and contributors. If you have any questions, please contact [fuzzified@zhengzhou.page](mailto:fuzzified@zhengzhou.page).
