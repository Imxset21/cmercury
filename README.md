## Synopsis

A simple C99 rewrite of the original Mercury software initially
created by John E. Chambers. It is the NBody integrator based on
Bulirsh-Stoer, Everhart and other methods. 

## Building

Building requires `gcc (version 4.8+)` and `make`. Once downloaded,
open a shell in the directory and type:
```
make all
```

A libmercury.so file will be produced. Warning: this is non-portable
(specific to your platform) by default!

## Tests

Testing is simple. Simply run:
```
make test
```

and a small test executable (default: `main.c` -> `main.bin`) will execute.

## Original Contributors

* John E. Chambers
* Hal Levison 
* Martin Duncan

## License

CMercury is published under the LGPLv3. See COPYING for details.

## Citation Guidelines

```
If you publish the results of calculations using (C)MERCURY, please
reference the package using J.E.Chambers (1999) "A Hybrid
Symplectic Integrator that Permits Close Encounters between
Massive Bodies". Monthly Notices of the Royal Astronomical
Society, vol 304, pp793-799. 
```
