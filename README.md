# snlines.py
supernova line identification tool

#########################################################
 Code for spectral line identification
 
 Thanks to the Bloodhound Gang...
 
 Daniel Kasen (kasen@berkeley.edu)
 
 improvements by Josiah Schwab

 usage: snlines.py spectrum.dat
 
 where spectrum.dat is a spectrum file with two
 
 column format:
 
    wavelength(angstroms)   flux


Line data from Kurucz cd23, stored in file kurucz.pkl

 Press ? at the prompt for commands:

  === List of Commands ===

 n, new      add species
 
 c, cycle    cycle the active species
 
 k, kill     remove active species
 
 e, species  list all species

 d, incv     increase velocity (blueshift)

 f, decv     decrease velocity (redshift)
 
 v, setv     set the velocity
 

 a, more     add line(s) to active species
 
 r, less     remove line(s) from active species
 
 l, lines    list all current lines

 -, load     load & plot a datafile
 
 -, oload    load & overplot a datafile
 
 -, uload    unload the most recently loaded datafile

 z, setz     set cosomological redshift
 
 q, quit     exit SNLines
 
 !, shell    pass a command to the shell
 
#########################################################
