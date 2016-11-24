# SPOT
Sequential Phase Optimization Technique (SPOT) package


Contents of this README file:

1. Warning
2. License
3. Package contents
4. Getting started


Contact: m.taylor@sbs.uq.edu.au

============
1. Warning
============

This code is a work in progress and may contain bugs. I do not guarantee the 
accuracy of all results. The expansions used can be susceptible to numerical 
errors due to the addition of very many numbers, particularly for very large 
particle size. Please let me know of any bugs you find. 

Any version updates will be posted at
[web address]

Contact me if you require help getting this to work.



THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 


==========
2. License
==========

Copyright 2016 The University of Queensland.

This package and its component files are copyright 2016 by The University of 
Queensland. Not-for-profit re-distribution of the package is permitted. The 
package may be used free-of-charge for research, teaching, or personal use. 

This package includes code from the Optical Tweezers Toolbox version 1.3, 
Copyright 2007-2014 The University of Queensland, reused with permission. 
The Optical Tweezers Toolbox is hosted at: 
http://www.physics.uq.edu.au/people/nieminen/software.html
 
If results obtained using this optimization package are published, both the 
package and the Optical Tweezers Toolbox should be appropriately referenced.


Further, it also includes code from the Quadrant Detection Package, reused
with permission. The Quadrant Detection Package is hosted at:



References

*The Optical Tweezers Toolbox package:
T. A. Nieminen, V. L. Y. Loke, A. B. Stilgoe, G. Knoener,
A. M. Branczyk, N. R. Heckenberg, H. Rubinsztein-Dunlop,
"Optical tweezers computational toolbox",
Journal of Optics A 9, S196-S203 (2007)
http://www.physics.uq.edu.au/people/nieminen/software.html


*The Quadrant Detection Package:
M. A. Taylor and W. P. Bowen, “A computational tool to characterize particle 
tracking measurements in optical tweezers,” J. Opt. 15, 085701 (2013). 


*SPOT optimization:



===================
3. Package contents
===================
The code was written as an add-on to the Optical Tweezers Toolbox, and most 
of the files included in the package are copied unmodified from the Optical 
Tweezers Toolbox version 1.3. 

Aside from these, the SPOT package files consist of:
- SPOT_optimize_kx.m
- SPOT_optimize_Fx.m
- SPOT_optimize_neg_Fz.m
- SPOT_optimize_kz.m
- SPOT_optimize_SNR.m
- Test_phase_solution.m
- Test_phase_solution_Fz.m

- farfield_matrix2.m (from Quadrant Detection Package)
- Quadrant_measurement.m (from Quadrant Detection Package)

- spharm.m (modified from Optical Tweezers Toolbox)








File descriptions:

- SPOT_optimize_kx.m: SPOT optimization of lateral trap stiffness. It 
contains the core SPOT algorithm and everything associated with that. 
It is quite well commented, by my standards.

- SPOT_optimize_Fx.m: SPOT optimization of lateral force magnitude. 

- SPOT_optimize_neg_Fz.m: SPOT optimization of pulling force magnitude. 

- SPOT_optimize_kz.m: SPOT optimization of axial trap stiffness. 

- SPOT_optimize_SNR.m: SPOT optimization of lateral displacement 
sensitivity, assuming measurement with a quadrant detector (for details
see Quadrant Detection Package). It can also calculate shot-noise limit 
to displacement sensitivity.

- Test_phase_solution.m: This was written to verify the results. It takes a 
calculated phase profile and re-evaluates the trap forces using whatever 
angular grid you choose and any Nmax; this can be used to check for numerical 
problems. It can also be used to change particle parameters. 

- Test_phase_solution_Fz.m: Similar to Test_phase_solution.m, but calculates 
both axial force, and the lateral force taken at the axial equilibrium. For 
all other calculateions, lateral forces are calculated from z=0 rather than 
the equilibrium point.

- farfield_matrix2: Calculates matrices that convert between the E(theta,phi) 
basis and vector spherical harmonics.

- Quadrant_measurement.m: Calculates quadrant detection signal. This is taken
straight from Quadrant Detection Package, and may be helpful when optimizing
signal-to-noise ratio (SNR).

- spharm.m: modified to improve numerical stability



==================
4. Getting started
==================


(a) Firstly, read the paper [CITATION] as well as its supplementary information. 
This will outline the basic principle of the SPOT algorithm, and how it is 
implemented in this software package.

(b) Install the package. Unzip all of the files into a directory. Either work in
 that directory, or add it to your MATLAB path.

(c) Play around with the code. Start with Optimize_Matrix_phase.m. I suggest 
starting with moderately easy parameters, e.g. 2um silica particle, as this is 
big enough to allow some interesting features in the optimum while also running 
quickly. The computation time scales strongly with particle size.


(d) All length units are specified in the medium wavelength, which in water is 
free-space wavelength/1.33.
    
(e) All force calculations are in units of photon momentum per photon. To convert to
    SI units:
                force_SI = force_Q * n * P/c
    where n is the refractive index of the surrounding medium,
          P is the beam power in watts,
          c is the speed of light in free space,

Although it works in photon units, all calculated spring constants are converted 
into SI units of N/m/W prior to reporting. As such, renormalization of the spring 
constant should not be necessary.
