Lattice-Liouville
=================

The code generates a discrete Gaussian free field on a periodic square lattice. From this field a discrete Liouville measure can be constructed, using a regularization method if desired. Geodesic distances can be computed w.r.t. the Liouville measure using either a simple Dijkstra algorithm or a Fast Marching solution to the Eikonal equation. 

For more information check:
* J. Ambjorn, T. Budd, *Geodesic distances in quantum Liouville gravity*, [arXiv:1405.3424](http://arxiv.org/abs/1405.3424)
* [Timothy Budd's academic webpage](http://www.nbi.dk/~budd/)

Pre-requisites
--------------

* BOOST C++ libraries (tested with versions 1.41 and 1.55)
* FFTW 3 (see http://www.fftw.org/)
* FFTW++ (tested with version 1.13, see http://fftwpp.sourceforge.net/)
* g++ or Visual C++ compiler 

To compile with g++ on linux:

g++ shortest_cycle.cpp liouville.cpp fftw++.cc -lfftw3 -o CYCLE

Usage
-----

When running ./CYCLE a several of inputs are requested. All inputs can also be given as command line arguments.
For example:

```
width = 512           // 512 x 512 lattice
gamma = 1.412         // scaling exponent 
averaging method (0=BOX_DIVISION,1=BOX_AVERAGING,2=DISK_AVERAGING,3=DISK_DIFFUSION) = 2
delta min = 0.0001    // perform measurements for values delta = 0.0001, 0.0002, ..., 0.0064
delta max = 0.0064
delta factor = 2.0
normalize volume (1=yes,0=no) = 0 
measurements per delta per output = 50   // determines how often the program writes its output
output path = output                     // puts the data in folder ./output/
```

License
-------

This source code may be reused without limitations according to the MIT license. If you use this code in your scientific publication, please cite the accompanying paper:

J. Ambjorn, T. Budd, *Geodesic distances in quantum Liouville gravity*, preprint, [arXiv:1405.3424](http://arxiv.org/abs/1405.3424)
