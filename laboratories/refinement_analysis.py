import numpy as np
import matplotlib.pyplot as plt

class Error(Exception):
    '''Class used to create errors'''
    pass

class gridsNotOrdered(Error):
    '''Raised when no the grids are not ordered with the finest grid first'''


class refinement_analysis:
    '''Grid refinement analysis tools

    A series of tools and a driver function for performing a grid
    refinement study using Python.  The tools are based on the
    RefinenementAnalysis.m Matlab programme written by D Ingram based
    on the NASA NPARC Fortran code published as part of the grid
    convergence assessment collaboration.

    The equations used are those given by Patric Roache in his 1998
    Book.

    David Ingram, School of Engingeering
    (c) 2021 The University of Edinburgh 
    License: CC-BY, Creative Commons by Attribution
    '''

    def __init__(self,x,f):
        ''' Create a new instance of the refinement analysis using
        data obtained from a mesh refinement study.  The two
        parameters given are arrays of the same length containing the
        grid spacing (x) and the function values (f) that have been
        evaluated on each grid.
        '''
        N = len(x) # number of grids used in the study
        # store the mesh spacing and the function values in the object
        self.f = f
        self.x = x
        
        # check the grid spacings are in the correct order and if not
        # throw an exception.
        if all(x[i] < x[i+1] for i in range(N-1)) is False:
            raise gridsNotOrdered

        # now compute the mesh spacings
        self.r = np.divide(x[1:],x[:-1])

        # finally estimate the the order of convergence using the
        # first three data pairs and assuming that the grid refinement
        # ratio is constant, r(1) = r(2). This is done using
        # Eqn. 5.10.6.1 of (Roache, 1998).  This needs three grids

        try:
            self.p = abs(np.log(abs(f[2]-f[1])/abs(f[1]-f[0]))/np.log(self.r[0]))
        except:
            print("Warning: Only 2 grids are provided, 2nd order is assumed")
            self.p = 2.0
            
    def richardson(self):
        '''calculate the Richardson extrapolation Estimate of the
        value of f on a mesh with delta_x = 0. '''

        f_exact = self.f[0] + (self.f[0]-self.f[1]) / (
            self.r[0]**self.p - 1.0)
        return f_exact

    def plot(self, Richardson = False, f_name = 'function'):
        '''produce a plot of f against dx.  If Richardson is True then
        the value of f with delta_x = 0 will be plotted.  f_name is
        used to determine the name of the function value to be used on
        the y-axis.'''

        # create a plot in Matplotlib
        fig, ax = plt.subplots(1)
        ax.plot(self.x,self.f, color = 'blue', linestyle=':',
            marker='o' )

        # add axis titles
        ax.set(xlabel = r'$\Delta x$',
            ylabel = f_name)
        
        # If Richardson is True add a symbol marker for Delta_x = 0 on
        # the y-axis
        if Richardson is True:
            ax.plot(0.0,self.richardson(),'rx')
        else:
            ax.set_xlim(left=0.0)

        # We don't want a box round the plot and we want to ensure the
        # x-axis sits exactly on x=0.0  The following commands turn
        # off the top and right "spines" and enforce the left edge at
        # x=0.0. 
        ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.yaxis.tick_left()

        ax.spines['top'].set_color('none')
        ax.xaxis.tick_bottom()

        # Display the plot.
        plt.show()
        
    def graph(self, Richardson = False, f_name = 'function', filename='refinement_analysis.pdf'):
        '''output a plot of f against dx.  If Richardson is True then
        the value of f with delta_x = 0 will be plotted.  f_name is
        used to determine the name of the function value to be used on
        the y-axis.'''

        # create a plot in Matplotlib
        fig, ax = plt.subplots(1)
        ax.plot(self.x,self.f, color = 'blue', linestyle=':',
            marker='o' )

        # add axis titles
        ax.set(xlabel = r'$\Delta x$',
            ylabel = f_name)
        
        # If Richardson is True add a symbol marker for Delta_x = 0 on
        # the y-axis
        if Richardson is True:
            ax.plot(0.0,self.richardson(),'rx')
        else:
            ax.set_xlim(left=0.0)

        # We don't want a box round the plot and we want to ensure the
        # x-axis sits exactly on x=0.0  The following commands turn
        # off the top and right "spines" and enforce the left edge at
        # x=0.0. 
        ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.yaxis.tick_left()

        ax.spines['top'].set_color('none')
        ax.xaxis.tick_bottom()

        # Display the plot.
        plt.savefig(filename)


    def GCI(self,A,B):
        ''' Compute Grid Convergence Index (GCI) between grids A and
        B using Eqn. 5.6.1 from Roache's book. Using a factor of
        safety as recommended on page 123. The parameters A and B
        are grid numbers, so grid 1 is the finest grid and grid N is
        the coarsest.''' 

        # How many grids have been used
        N = len(self.x)

        # Set the satfey factor, if there are more than 2 grids then
        # this is 1.25, if there are only 2 grids then the saftey
        # factor is 3.0
        if (N>2):
              f_safe = 1.25
        else:
              f_safe = 3.0

        # Now compute the relative error in f between the grids
        epsilon = abs(self.f[B-1] - self.f[A-1]) / abs(self.f[A-1])

        # and return the GCI
        return f_safe * epsilon / (self.r[0]**self.p - 1.0)

    def GCI_all(self):
        ''' Compute the GCI on the sequence of all grid pairs'''
        
        # How many grids have been used
        N = len(self.x)

        # Loop over the grids calculating the GCI(A,B) with A going
        # from 1 to N-1 and B from 2 to N.
        gci = []
        for i in range(N-1):
            gci.append(self.GCI(i+1,i+2))
            
        return gci

    def GCI_ratio(self,tol=1.0e-2):
        '''Examine if asymptotic range has been achieved by examining
        ratio of Eqn. 5.10.5.2 of Roache's book is one.  tol is used
        as a tollerance with np.isclose(ratio,1.0,tol) to see if the
        computed ratio is close to 1.0'''

        # Calculate the GCI for all the grids
        gci = self.GCI_all()

        # find out how many GCI's we have calculated
        N = len(gci)

        if (N>2):
            # Loop over the grids calculating the ratios
            ratio = []
            for i in range(N-1):
                ratio.append(self.r[0]**self.p * gci[i]/gci[i+1])
        else:
            print("Warning: Only 2 grids used, ratios cannot be calculated.")
            ratio = [1000.0]

        return [ratio, np.isclose(ratio, 1.0, tol)]

    def report(self,f_name = 'function', tol =1.0e-2):
        ''' print a report on the grid convergence study. f_name is
        the description of the function being measured on the grid and
        tol is the tolerance used to see if the GCI ratio is close to
        one. '''

        # How many grids
        N = len(self.x)
        print(f"Refinement Analysis on {N} grids.\n")

        # Tabulate the grids and the function values
        heading = f"Grid Delta x ratio  {f_name}"
        print(heading)
        print("-"*len(heading))
        for i in range(N-1):
            print("{:4d} {:7.4g} {:5.2g} {:9.6g}"
                      .format(i+1, self.x[i], self.r[i], self.f[i]))
        print("{:4d} {:7.4g}   --- {:9.6g}"
                  .format(N, self.x[-1], self.f[-1]))
        print("-"*len(heading),"\n")

        # report the order of accuracy
        print("order of convergence, p = {:6.2f}".format(self.p))

        # report the Grid Convergence Index
        gci = self.GCI_all()
        print("\nGrids       GCI\n---------------")
        for i in range(len(gci)):
            print("{:2d} {:2d} {:9.4g}".format(i+1,i+2, gci[i]))
        print("---------------\n")

        # report the GCI Ratios
        if (len(gci)>2):
            ratio, converged = self.GCI_ratio(tol)
            print(" Grid Step  GCI Ratio Converged")
            print("-------------------------------")
            for i in range(len(ratio)):
                print("{:2d}{:2d} {:2d}{:2d}     {:7.4g} {}"
                          .format(i+1,i+2,i+2,i+3,ratio[i],converged[i]))        
            print("-------------------------------\n")

    # These last two functions enable python to describe the
    # refinement analysis object (using repr) and to give a short
    # human readable summayry (using str).

    def __repr__(self):
        ''' describe the object for the command line'''
        return f'refinement analysis(x={self.x}, f={self.f})'

    def __str__(self):
        ''' describe the object for a print command'''

        N = len(self.x)

        return f'''Refinement Analysis on {N} grids.
        
  dx = {self.x}
  f  = {self.f}

  refinement ration, r = {self.r}
  order of convergence, p = {self.p}

  Grid Convergence Indicies: {self.GCI_all()}

  GCI ratios : {self.GCI_ratio()}
'''
    
        
