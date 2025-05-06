from __future__ import print_function
import numpy as np
from numpy import pi
from . import model
from scipy.integrate import solve_bvp


class SQGModel(model.Model):
    """Surface quasigeostrophic model."""

    def __init__(
        self,
        beta=0.,                    # gradient of coriolis parameter
        Nb = 1.,                    # Buoyancy frequency
        rd=0.,                      # deformation radius
        H = 1.,                     # depth of layer
        U=0.,                       # max vel. of base-state
        finite=-1,                  # 0 is infinite SQG, 1 is SQG with vanishing bottom streamfunction
        sigma2=None,
        slip=None,
        alpha_exp = 0,              # alpha turbulence exponent when finite = 3
        nvgrid = 200, # Number of grid points in vertical (only for custom stratification)
        load_mk = False, # load mk from file
        mk_file = None, # name of .npy file with mk array
        **kwargs
        ):
        """
        Parameters
        ----------

        beta : number
            Gradient of coriolis parameter. Units: meters :sup:`-1`
            seconds :sup:`-1`
        Nb : number
            Buoyancy frequency. Units: seconds :sup:`-1`.
        U : number
            Background zonal flow. Units: meters.
        """

        # physical
        self.beta = beta
        self.Nb = Nb
        #self.rek = rek
        self.rd = rd
        self.H = H
        self.Hi = np.array(H)[np.newaxis,...]
        self.U = U
        self.nvgrid = nvgrid
        
        self.finite = finite
        #self.filterfac = filterfac
        self.load_mk = load_mk
        self.mk_file = mk_file

        self.nz = 1        
       
        self.sigma2 = sigma2
        # set Nb from sigma2 profile
        if self.sigma2 and not self.load_mk: 
            self.Nb = np.sqrt(self.sigma2(0))
        self.slip = slip
        
        self.alpha_exp = alpha_exp

        super(SQGModel, self).__init__(**kwargs)

        # initial conditions: (PV anomalies)
        self.set_q(1e-3*np.random.rand(1,self.ny,self.nx))

    ### PRIVATE METHODS - not meant to be called by user ###

    def _initialize_background(self):
        """Set up background state (zonal flow and PV gradients)."""

        # the meridional PV gradients in each layer
        self.Qy = np.asarray(self.beta)[np.newaxis, ...]

        # background vel.
        self.set_U(self.U)

        # complex versions, multiplied by k, speeds up computations to pre-compute
        self.ikQy = self.Qy * 1j * self.k

        self.ilQx = 0.
    
    def solve_bvp_params(self,k2,nv=200):
        z1 = -self.H
        z2 = 0.
        Z = np.linspace(z1, z2, nv)
        Y = np.ones((2, Z.size))
        
        self.Z = Z
    
        def in_dYdz(Z,Y):
            sigma2Z = np.asarray([self.sigma2(z) for z in Z])
            return np.vstack([sigma2Z*Y[1,:], k2*Y[0,:]])
    
        if self.slip==0: # no slip 
            def in_Ybc(Ya,Yb):
                return np.array([Ya[0],Yb[0]-1])
        elif self.slip==1: # slip 
            def in_Ybc(Ya,Yb):
                return np.array([Ya[1],Yb[0]-1])
        else:
            print("WARNING: SLIP MUST BE 0 OR 1.")
            return None
    
        sol = solve_bvp(in_dYdz, in_Ybc, Z, Y)
        sigma2Z = np.asarray([self.sigma2(z) for z in Z])
        
        y = sol.sol(Z)[0]
        yp  = sol.sol(Z)[1]*sigma2Z
        return (y,yp)
    
    def find_vertical_wavenumber(self,nv=200):
        print("Calculating vertical wavenumbers...")
        self.mk = np.zeros_like(self.wv2)
        for j in range(len(self.ll)):
            if j%10==0:
                print("{0:.1f}% complete".format(100*j/len(self.ll)))
            for i in range(len(self.kk)):
                k2 = self.ll[j]**2+ self.kk[i]**2
                y,yp = self.solve_bvp_params(k2,nv=nv)
                self.mk[j,i] = yp[-1]/self.Nb**2
        print("Done.")

    def _initialize_inversion_matrix(self):
        """ the inversion """
        # The sqg model is diagonal. The inversion is simply qh = -kappa**2 ph
        
        print("NSQGv1.1\n") 
        
        if self.sigma2 == None:
            if self.finite ==0:
                self.a = -np.asarray(self.Nb*np.sqrt(self.wv2i))[np.newaxis, np.newaxis, :, :]
            elif self.finite==1: # SQG in finite depth
                print("Running finite SQG.")
                tanh = np.tanh(self.Nb*np.sqrt(self.wv2)*self.H)
                itanh = tanh !=0.
                coth = np.zeros_like(tanh)
                coth[itanh] = tanh[itanh]**-1
                self.a = -np.asarray(self.Nb*coth*np.sqrt(self.wv2i))[np.newaxis, np.newaxis, :, :]     
            elif self.finite==2: # SQG with vanishing streamfunction
                print("Running topographic SQG.")
                tanh = np.tanh(self.Nb*np.sqrt(self.wv2)*self.H)
                self.a = -np.asarray(self.Nb*tanh*np.sqrt(self.wv2i))[np.newaxis, np.newaxis, :, :]
            elif self.finite==3: # alpha turbulence; alpha must be specified
                 self.a = -np.asarray(self.Nb*np.sqrt(self.wv2i)**self.alpha_exp)[np.newaxis, np.newaxis, :, :]
            elif self.load_mk: # sigma2=None and m(k) is prescribed
                print("Loading m(k) with sigma2=None.")
                self.mk = np.load(self.mk_file)
                imk = self.mk !=0.
                mki = np.zeros_like(self.mk)
                mki[imk] = self.mk[imk]**(-1)
                self.a = -np.asarray(mki)[np.newaxis,np.newaxis,:,:]
            else:
                print("Invalid option.Choose finite=n for the following options.\n 0 - SQG\n 1 - fSQG\n 2 - tSQG\n 3")            
        else: # if sigma2 is given    
            print("Running custom stratification")
            if self.load_mk:
                self.mk = np.load(self.mk_file)
            else: # problem: find_vertical_wavenumber is not parallized.
                self.find_vertical_wavenumber(self.nvgrid)
            imk = self.mk !=0.
            mki = np.zeros_like(self.mk)
            mki[imk] = self.mk[imk]**(-1)
            self.a = -np.asarray(mki)[np.newaxis,np.newaxis,:,:]

    def _initialize_forcing(self):
        pass

    def set_U(self, U):
        """Set background zonal flow"""
        self.Ubg = np.asarray(U)[np.newaxis,...]

    def _calc_diagnostics(self):
        # here is where we calculate diagnostics
        if (self.t>=self.dt) and (self.tc%self.taveints==0):
            self._increment_diagnostics()

    ### All the diagnostic stuff follows. ###
    def _calc_cfl(self):
        return np.abs(
            np.hstack([self.u + self.Ubg, self.v])
        ).max()*self.dt/self.dx


    # calculate KE: this has units of m^2 s^{-2}
    def _calc_ke(self):
        ke = .5*self.spec_var(self.wv*self.ph)
        return ke.sum()

    # calculate eddy turn over time
    # (perhaps should change to fraction of year...)
    def _calc_eddy_time(self):
        """ estimate the eddy turn-over time in days """
        ens = .5*self.H * spec_var(self, self.wv2*self.ph)
        return 2.*pi*np.sqrt( self.H / ens ) / year
    
    def _calc_derived_fields(self):
        self.p = self.ifft(self.ph)