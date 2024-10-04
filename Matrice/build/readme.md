# Compiled binaries and dependencies

## Usage
>   - Unzip **matrice_x64[d].zip** to get the *lib* files, where the postflag *d* indicates debug version.
>   - Put them together in one folder, such as *lib*.

## 3rdparty lib
Matrice is speeded up by the Intel MKL library when the ```MATRICE_MATH_KERNEL == MATRICE_USE_MKL``` is enabled (default).
To make your code run correctly, copy the files in the folder *bin* to the folder where the exe resides. 
