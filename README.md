# SIDM

SIDM_pde.py calculates the log of the density of DM in the presence of some baryon potential

it includes optional use of command line arguments:
--output: type = str, default = None,
    help = 'names (prefix) of output files to write')
--r_len: type = int, default=50,
    help = 'size of array along r direction')
--z_len: type = int, default=50,
    help = 'size of array along r direction')
--scale: type = float, default=-4,
    help = 'exponential term to describe innermost point r_inner = r1 * 10^scale')
--r1: type = float, default=5.,
    help = 'the boundary of the array in kpc')
--rho0: type = float, default=5.e6,
    help = 'scale density of the solution in M_sun/kpc')
--num_iter: type = int, default=1000000,
    help = 'number of iterations in the relaxation')
--no_baryons: action='store_true',
    help = 'type --no_baryons if you want an isothermal core')

example command line:
python SIDM_pde.py --r_len 100 --z_len 100 --scale -5. --r1 10. --rho0 1.e7 --num_iter 4000000 --no_baryons

also writes a .txt file with the parameters used

no longer includes no baryon case as a seperate file since that functionality is now in an if statement    
    
as of 2/15/17 all the plots with _iso suffix comes from the nobaryons solution