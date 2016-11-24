% Optical Tweezers Toolbox 1.3
% 
% Contents of this README file:
% 
% 1. A warning!
% 2. License
% 3. Installation
% 4. Getting started
% 5. Miscellany
% 6. References
% 
% This README file is supplied as a text file (README.txt)
% and as an m-file (README.m).
% 
% Contact: timo@physics.uq.edu.au
% 
% ============
% 1. A warning
% ============
% 
% This release incorporates many of the previous bug fixes. However, as
% always this is a work in progress and may contain bugs from the new code.
% Bugs you discover may prevent the correct results from being obtained. 
% Let us know what you find.
% 
% See Section 4 Getting started for some basics,
% and check out the example files. Contact us if you can't figure
% it out, and we'll provide what help we can.
% 
% Anyway, check
% http://www.physics.uq.edu.au/people/nieminen/software.html
% or
% http://www.physics.uq.edu.au/omg/Links.html
% to see if a new version is out.
% 
% We do not guarantee that results in all cases will be correct and disavow
% any liability for any damages this code may inflict. However, please keep
% us in the loop so that we can endeavour to fix any issues you experience.
% 
% ==========
% 2. License
% ==========
% 
% Copyright 2007-2013 The University of Queensland.
% 
% This package and its component files are copyright 2007-2013 by 
% The University of Queensland. Not-for-profit re-distribution of 
% the unmodified complete package.
% 
% The package may be used free-of-charge for research, teaching,
% or personal use. If results obtained using the package are
% published, the package should be appropriately referenced.
% 
% ===============
% 3. Installation
% ===============
% 
% The easy way: unzip all of the files into a directory, and work
% in there.
% 
% The better way: unzip all of the files into a directory, and add
% that directory to your MATLAB path. This way, you can keep your
% files separate from the package files.
% 
% ==================
% 4. Getting started
% ==================
% 
% (a) Read the paper (ott_preprint.pdf, available from website above).
%     Read the user Guide (ottug_1.2.pdf). It will outline the basic 
%     operation of the toolbox.
% 
% (b) Copy the examples to your working directory, and play with
%     them. Start with example_gaussian.m.
% 
% (c) It's best to use length units of the wavelength in the trapping
%     medium, usually free-space wavelength/1.33.
%     
% (d) The examples calculate the force and torque efficiencies. These are
%     the force and torque per photon, in photon units. To convert to
%     SI units:
%                 force_SI = force_Q * n * P/c
%                torque_SI = torque_Q * P/w
%     where n is the refractive index of the surrounding medium,
%           P is the beam power in watts,
%           c is the speed of light in free space,
%           w is the angular optical frequency, in radians/s.
% 
% =============
% 5. Miscellany
% =============
% 
% Plans for the future?
% (a) Make the package less user-unfriendly. A good start would be to
%     include automatic choice of Nmax.
% (b) A cool GUI interface would be nice, but (a) and (b) are much
%     higher priorities.
% (c) Add T-matrix routines for essentially arbitrary particles.
% 
% Want to contribute? Feel free to do so, either in the form of code
% or suggestions. If you contribute code, you'll need to assign
% copyright to The University of Queensland.
% 
% =============
% 6. References
% =============
% 
% The package:
% 
% T. A. Nieminen, V. L. Y. Loke, A. B. Stilgoe, Y. Hu, G. Knoener,
% A. M. Branczyk,
% "Optical tweezers toolbox 1.1",
% http://www.physics.uq.edu.au/people/nieminen/software.html
%  
% 
% Descriptions of the package:
% 
% T. A. Nieminen, V. L. Y. Loke, A. B. Stilgoe, G. Knoener,
% A. M. Branczyk, N. R. Heckenberg, H. Rubinsztein-Dunlop,
% "Optical tweezers computational toolbox",
% Journal of Optics A 9, S196-S203 (2007)
% 
% T. A. Nieminen, V. L. Y. Loke, G. Knoener, A. M. Branczyk,
% "Toolbox for calculation of optical forces and torques",
% PIERS Online 3(3), 338-342 (2007)
% 
% 
% More about computational modelling of optical tweezers:
% 
% T. A. Nieminen, N. R. Heckenberg, H. Rubinsztein-Dunlop,
% "Computational modelling of optical tweezers",
% Proc. SPIE 5514, 514-523 (2004)
% 
% 
% More about our beam multipole expansion algorithm:
% 
% T. A. Nieminen, H. Rubinsztein-Dunlop, N. R. Heckenberg,
% "Multipole expansion of strongly focussed laser beams",
% Journal of Quantitative Spectroscopy and Radiative Transfer 79-80,
% 1005-1017 (2003)
% 
% More about our T-matrix algorithm:
% 
% T. A. Nieminen, H. Rubinsztein-Dunlop, N. R. Heckenberg,
% "Calculation of the T-matrix: general considerations and
% application of the point-matching method",
% Journal of Quantitative Spectroscopy and Radiative Transfer 79-80,
% 1019-1029 (2003)
% 
% 
% The multipole rotation matrix algorithm we used:
% 
% C. H. Choi, J. Ivanic, M. S. Gordon, K. Ruedenberg,
% "Rapid and stable determination of rotation matrices between
% spherical harmonics by direct recursion"
% Journal of Chemical Physics 111, 8825-8831 (1999)
% 
% 
% The multipole translation algorithm we used:
% 
% G. Videen,
% "Light scattering from a sphere near a plane interface",
% pp 81-96 in:
% F. Moreno and F. Gonzalez (eds),
% Light Scattering from Microstructures, LNP 534,
% Springer-Verlag, Berlin, 2000
% 
% 
% More on optical trapping landscapes:
% 
% A. B. Stilgoe, T. A. Nieminen, G. Knoener, N. R. Heckenberg, H. 
% Rubinsztein-Dunlop, "The effect of Mie resonances on trapping in 
% optical tweezers", Optics Express, 15039-15051 (2008)
% 
% 
% Multi-layer sphere algorithm:
% 
% W. Yang, "Improved recursive algorithm for light scattering by a
%  multilayered sphere", Applied Optics 42(9), (2003)
