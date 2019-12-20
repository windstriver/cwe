# Benchmark ERCOFTAC UFR 2-02

[Source](https://openfoamwiki.net/index.php/Benchmark_ercoftac_ufr2-02)

## Flow past a square cylinder
### Description
For a full description of the case see:
http://qnet-ercoftac.cfms.org.uk/index.php?title=Flow_past_cylinder

or

[here](http://cfd.mace.manchester.ac.uk/cgi-bin/cfddb/prpage.cgi?43&EXP&database/cases/case43/Case_data&database/cases/case43&cas43_head.html&cas43_desc.html&cas43_meth.html&cas43_data.html&cas43_refs.html&cas43_rsol.html&1&0&0&0&0)

U_in = 0.54 m/s with 2% turbulence

Re_D = 22 000

Size of the square cylinder is 0.04 m wide and 0.392 m high.

Its situated in a tunnel with a width of 14 D = 0.56 m.

### Results

Experimental results

The Strouhal number is around 0.132

St = n*D/U

Drag coefficient is around 1.9-2.2

### Notes

In the case directory there is a perl script called generateMeshDict.pl.

This script will generate blockMeshDict file given input arguments.

### References

(http://www.tfd.chalmers.se/~lada/postscript_files/ahmad_corsica.pdf)

(http://www.tfd.chalmers.se/~lada/postscript_files/ahmad_ASME-JSME_paper.pdf)
