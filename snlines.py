#!/usr/bin/env python
import os
import cmd
import readline
import optparse
import numpy as np
import pylab as py
import cPickle as pickle
from collections import deque

#########################################################
# Code for spectral line identification
# Thanks to the Bloodhound Gang...
# Daniel Kasen (kasen@berkeley.edu)
# improvements by Josiah Schwab
#
# usage: snlines.py spectrum.dat
# where spectrum.dat is a spectrum file with two
# column format:
#    wavelength(angstroms)   flux
#########################################################

romannums = ('I','II','III', 'IV')
elements = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni','Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'I', 'Te', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U')

def species_name(specid):
    return "{} {}".format(elements[specid[0]-1], romannums[specid[1]])

### Define SNConsole ###

class SNConsole(cmd.Cmd):

    def __init__(self, state = None):
        cmd.Cmd.__init__(self)
        self.prompt = "SN> "
        self.intro = "Welcome to SNLines"

        # start with a state, or if not, create a blank one
        if state is None:
            self.state = SNLinesState()
        else:
            self.state = state

    ## Internal functions to aid with argument conversions ##
    def _try_int(self,args, default = 1):
        try:
            value = int(args)
        except (TypeError, ValueError):
            value = default
        return value

    def _try_float(self,args, default = 0):
        try:
            value = float(args)
        except (TypeError, ValueError):
            value = default
        return value

    ## Command definitions ##
    def do_exit(self, args):
        """Exit SNLines
        
        SN> exit
        """
        return -1

    def do_new(self,args):
        """Add a new species 

        SN> new <element> <ionstate>

        Either as 
        SN> new Ca II
        or 
        SN> new 20 1
        """
        try:
            e, i = args.split(' ')
        except:
            e = 0
            i = 0

        if e in elements and i in romannums:
            e = elements.index(e) + 1
            i = romannums.index(i)
        else:
            try:
                e = int(e)
                i = int(i)                
            except ValueError:
                e = 0
                i = 0

        specid = (e,i)
        return self.state.add_species(specid)

    def do_cycle(self, args):
        """Cycle through the species

        SN> cycle [<n> = 1]
        """
        n = self._try_int(args)
        return self.state.cycle_species(n)

    def do_kill(self,args):
        """Remove the active species

        SN> kill
        """
        return self.state.remove_species()

    def do_load(self,args):
        """Load and plot a new data file

        SN> load <filename>
        """
        filename = args
        return self.state.read_data(filename)

    def do_oload(self,args):
        """Load and overplot a new data file

        SN> oload <filename>
        """
        filename = args
        return self.state.read_data(filename, append = True)

    def do_uload(self,args):
        """Unload the most recent data file
        
        SN> uload
        """
        filename = args
        return self.state.remove_data()

    def do_incv(self,args):
        """Increase velocity of active speciess

        SN> incv [<deltav> = 1e3]
        """
        active = self.state.get_active_species()
        if active is not None:
            dv = self._try_float(args, default = 1e3)
            active.v += dv

    def do_decv(self,args):
        """Decrease the velocity of the active species

        SN> decv [<deltav> = 1e3]
        """
        active = self.state.get_active_species()
        if active is not None:
            dv = self._try_float(args, default = 1e3)
            active.v -= dv

    def do_setv(self,args):
        """Set the velocity of the active speciess

        SN> setv [<v> = 1e4]
        """
        active = self.state.get_active_species()
        if active is not None:
            v = self._try_float(args, default = 1e4)
            active.v = v

    def do_more(self,args):
        """Show more lines of the active species

        SN> more [<deltan> = 1]
        """
        active = self.state.get_active_species()
        if active is not None:
            dn = self._try_int(args)
            active.n += dn

    def do_less(self,args):
        """Show fewer lines of the active species

        SN> less [<deltan> = 1]
        """
        active = self.state.get_active_species()
        if active is not None:
            dn = self._try_int(args)
            active.n -= dn

    def do_setn(self,args):
        """Show (# = 1) lines"""
        active = self.state.get_active_species()
        if active is not None:
            n = self._try_int(args, default = active.n)
            active.n = n

    def do_setz(self,args):
        """Set the cosmological redshift

        SN> setz [<z> = 0]
        """
        z = self._try_float(args,default = 0)
        self.state. z = z

    def do_lines(self,args):
        """List displayed lines

        SN> lines
        """
        return self.state.describe_lines()

    def do_species(self,args):
        """List displayed species

        SN> species
        """
        return self.state.describe_species()


    ## Command aliases ##

    def do_n(self, args):
        """alias for new"""
        return self.do_new(args)

    def do_c(self, args):
        """alias for cycle"""
        return self.do_cycle(args)

    def do_k(self, args):
        """alias for kill"""
        return self.do_kill(args)

    def do_e(self, args):
        """alias for species"""
        return self.do_species(args)

    def do_f(self,args):
        """alias for decv"""
        return self.do_decv(args)

    def do_d(self,args):
        """alias for incv"""
        return self.do_incv(args)

    def do_v(self,args):
        """alias for setv"""
        return self.do_setv(args)

    def do_a(self,args):
        """alias for more"""
        return self.do_more(args)

    def do_r(self,args):
        """alias for incv"""
        return self.do_less(args)

    def do_l(self,args):
        """alias for lines"""
        return self.do_lines(args)

    def do_z(self,args):
        """alias for setz"""
        return self.do_setz(args)

    def do_q(self,args):
        """alias for exit"""
        return self.do_exit(args)

    def do_quit(self,args):
        """alias for exit"""
        return self.do_exit(args)

    ## Override cmd Commands ##
    def do_shell(self, args):
        """Pass command to a system shell when line begins with '!'"""
        os.system(args)

    def do_help(self, args):
        """Get help on commands
           'help <command>' or '? <command>' gives help on <command>
        """
        if args == "":
        # for no argument, list our commands

            print("for more detailed help, type ? <command>")
            print("this will display usage information of the form")
            print("  SN> command <required> [<optional> = default]")
            print("")

            print("=== List of Commands ===")
            print("")

            print("n, new      add species")
            print("c, cycle    cycle the active species")
            print("k, kill     remove active species")
            print("e, species  list all species")

            print("")

            print("d, incv     increase velocity (blueshift)")
            print("f, decv     decrease velocity (redshift)")
            print("v, setv     set the velocity")

            print("")

            print("a, more     add line(s) to active species")
            print("r, less     remove line(s) from active species")
            print("l, lines    list all current lines")
            
            print("")

            print("-, load     load & plot a datafile")
            print("-, oload    load & overplot a datafile")
            print("-, uload    unload the most recently loaded datafile")
            
            print("")

            print("z, setz     set cosomological redshift")
            print("q, quit     exit SNLines")

            print("!, shell    pass a command to the shell")

        else:
        # use built-in method to fetch the doc string
            cmd.Cmd.do_help(self, args)


    def postcmd(self, stop, line):
        """After command has been processed, update the screen"""
        self.state.show()
        return stop


class Species(object):

    def __init__(self, name = "active", velocity = 1e4, nshow = 1):

        self._name = name
        self._v = velocity
        self._nshow = nshow

    @property
    def name(self):
        """Species name"""
        return self._name

    @property
    def v(self):
        """Species velocity"""
        return self._v

    @v.setter
    def v(self, value):
        self._v = value

    @property
    def n(self):
        """Species number of lines"""
        return self._nshow

    @n.setter
    def n(self, value):
        self._nshow = value
        
    def describe(self):
        return "{} ; vel = {} km/s".format(self._name, self._v)


        
class SNLinesState:

    def __init__(self, pklfile = "/Users/kasen/Dropbox/codes/my_python/kurucz.pkl", defaultv = 1e4):
        
        # read in line data, pickled by pickedata.py
        with open(pklfile, 'rb') as f:
            self.linedict = pickle.load(f)

        # empty state
        self.species = {}
        self.specids = deque()
        self._z = 0 
        self.xdata = []
        self.ydata = []
        self.xrange = [3000, 9000]
        self.yrange = [0,1]
        self.defaultv = defaultv

    @property
    def z(self):
        """Cosmological redshift"""
        return self._z

    @z.setter
    def z(self, value):
        self._z = value

    def get_active_species(self):
        if len(self.specids) > 0:
            return self.species[self.specids[-1]]
        else:
            return None

    def read_data(self, filename, append = False, rescale = True):
        """read in the data file"""
        try:
            ddata  = np.loadtxt(filename) #, unpack = True)
            newxdata = ddata[:,0]
            newydata = ddata[:,1]
        except IOError:
            print("*** Unable to read file: {}".format(filename))
        else:
            if rescale:
                self.yrange = (0,1.1*max(newydata))
                self.xrange = (min(newxdata),max(newxdata))
            if append:
                self.xdata.append(newxdata)
                self.ydata.append(newydata)
            else:
                self.xdata = [newxdata]
                self.ydata = [newydata]

    def remove_data(self):
        if len(self.xdata) > 0:
            self.xdata.pop()
            self.ydata.pop()

    def add_species(self, specid):
        """add a species"""
        if specid in self.specids:
            print("Species already exists")
            return

        if specid in self.linedict:
            
            active = self.get_active_species()
            if active is not None:
                v = active.v
            else:
                v = self.defaultv
            
            self.species[specid] = Species(name = species_name(specid),
                                           velocity = v)
            self.specids.append(specid)
        else:
            print("Species not in datafile")

    def remove_species(self):
        """take species out of states and specids"""
        if len(self.specids) > 0:
            self.species.pop(self.specids.pop())
        
    def cycle_species(self, n):
        """shift the species deque"""
        self.specids.rotate(n)

    def describe_species(self):
        """print which species we're using"""
        for id in self.specids:
            print(self.species[id].describe())

    def describe_lines(self):
        """print which lines we're displaying"""
        for id in self.specids:
            s = self.species[id]
            ld = self.linedict[id]
            print(s.describe())
            print("   lambda        gf       E_low      tau")

            # lam_obs = ld['lam'] * (1 - s.v/3e5)
            # visible = (lam_obs > self.xrange[0]) & (lam_obs < self.xrange[1])
            # nshow = max(s.n, len(visible))
            
            nshow = min(s.n, len(ld['lam']))
            for i in range(nshow):                
                print(" %10.3e %10.3e %10.3e %10.3e" % (ld['lam'][i],
                                                        ld['gf'][i],
                                                        ld['El'][i],
                                                        ld['tau'][i]))


    def show(self):
        """Plots the spectral data and overlays the lines"""
        py.clf()
        
        for xpts,ypts in zip(self.xdata,self.ydata):
            py.plot(xpts /(1+self.z), ypts)

        py.xlabel('wavelength',size=15)
        py.ylabel('flux',size=15)

        py.ylim(self.yrange)
        py.xlim(self.xrange)

        # make sure that things are set ok for active region
        try:
            active = self.specids[-1]
        except IndexError:
            active = None

        # overplot lines
        for id in self.specids:

            ld = self.linedict[id]
            s = self.species[id]

            if id == active:
                color = 'red'
                py.title(s.describe(),size=15)
            else:
                color = 'black'

            lam_obs = ld['lam'] * (1 - s.v/3e5)
            # visible = (lam_obs > self.xrange[0]) & (lam_obs < self.xrange[1])

            nshow = min(s.n, len(lam_obs))
            for lam in lam_obs[:nshow]:
                py.axvline(lam, color=color)

        py.show()


if __name__ == "__main__":
    
    # get command line arguments
    parser = optparse.OptionParser()
    parser.add_option("--xr",dest="xrange")
    parser.add_option("--yr",dest="yrange")
    parser.add_option("-v", dest="vel", type="float", default=1e4)

    # set defaults, if not specified
    (opts, args) = parser.parse_args()

    state = SNLinesState(defaultv = opts.vel)
    for arg in args:
        state.read_data(arg, append = True)

    # override defaults with values from command line
    if opts.xrange is not None:
        xx = opts.xrange.split(',')
        state.xrange = [float(x) for x in xx]

    if opts.yrange is not None:
        yy = opts.yrange.split(',')
        state.yrange = [float(y) for y in yy]

    # turn interactive mode on & show first plot
    py.ion()
    state.show()

    # get a console object and enter the main loop
    console = SNConsole(state)
    console.cmdloop() 

    print("Exiting...")
