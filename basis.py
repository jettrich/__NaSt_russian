'''
HINTS:
q=pylab.quiver(a.xu,a.yu,a.u[-1],a.v[-1])
pylab.show(q)
'''


from numpy import pi, shape, ones, zeros
from numpy import meshgrid, linspace
# misc.derivative
# from scipy.misc import *
# from scipy.integrate import *
from scipy import ndimage as nd
import pylab


class NS: 
    '''
    DESCRIPTION:
    solve Navier-Stokes equation:
       d(ee)/dt = -u*d(ee)/dx - v*d(ee)/dy + Re*laplace(ee)
       ee = laplace(ps)
       d(ps)/dy=u,
       d(ps)/dx=-v
    for tube with barrier (rectangle) inside
    
    WARNING:
    for more accurate solution time steps (TN) shold be
    more 100. For my PC it spread 3 hours.
    Despite of  memory error while it solves you can view
    solution at last working step due to save history in
    instance object (for that reason we use NS class instead
    just branch of functions). Just call instance.plot() function.
    '''

    def __init__(self, Nx=60, Ny=60, TN=22, TN_pde1=1, c=1, x01=3, y01=2):
        '''
        INPUT:
        Nx,Ny- count of points between [0,x01] and [0,y01]
        TN-count of time steps
        others see below
        '''
        # barrier rotation (not tested yet)
        self.a = 0

        # length of barrier
        self.lenght = 43

        # size of tube
        self.Nx = Nx
        self.Ny = Ny

        # width of barier
        self.hh = Ny/2

        # count of steps
        self.TN = TN

        # steps for convergence
        self.TN_pde1 = TN_pde1*int(2*(max(Nx, Ny))**2/pi**2)

        self.x0 = [0, Nx]
        self.y0 = [0, Ny]
        self.x0[1] = x01
        self.y0[1] = y01

        self.ht0 = 0.001
        self.u0 = 0.1
        self.u00 = 1
        
        # Renolds
        self.Re = 1 # 1

        # arrays for barrier (internal and external)
        self.D1 = zeros((Nx, Ny))
        self.D2 = zeros((self.Nx, self.Ny))
        
        # f will be const for test
        self.f_test = lambda x, y: x-x+4
        self.psG_test = lambda x, y: x**2+y**2+2*x*y

        # create states
        self.ee = zeros((self.TN, self.Nx, self.Ny))
        self.ps = zeros((self.TN, self.Nx, self.Ny))
        self.psG = zeros((self.Nx, self.Ny))
        self.f = zeros((self.Nx, self.Ny))
        self.u = zeros((self.TN, self.Nx, self.Ny))
        self.v = zeros((self.TN, self.Nx, self.Ny))

        # width of tube for init conditions
        self.H = abs(self.y0[1]-self.y0[0])

        # integrating constant for init conditions
        self.c = c

        self.x = linspace(self.x0[0], self.x0[1], self.Nx)
        self.y = linspace(self.y0[0], self.y0[1], self.Ny)
        
        # for createing simmetry in init conditions restrictions
        self.constPs0 = (self.c*(self.y[int(Ny/2)]**2*self.H/2.
                                 - self.y[int(Ny/2)]**3/3.))
        self.constPs1 = self.u00*self.y[int(Ny/2)]

        # type of conditions
        self.cn_u = 0
        self.cn_ps = 0
        self.cn_ee = 0
        
        # t = linspace(0, T, TN)
        self.hx = abs(self.x[1]-self.x[0])
        self.hy = abs(self.y[1]-self.y[0])
        # ht = abs(t[1]-t[0])

        if self.hy < self.hx:  # for convergence
            self.ht = (self.hy**2)/4.
        else:
            self.ht = (self.hx**2)/4.

        self.yu, self.xu = meshgrid(self.y, self.x)
        self.set_S()

    def test_default(self, timeSteps=22):
        '''
        DESCRIPTION:
        Get solution phase field with default values.
        '''
        self.TN = timeSteps
        self.clear()
        self.pde2()
        self.plot()

    def plot(self):
        q = pylab.quiver(self.xu, self.yu, self.u[-1], self.v[-1])
        pylab.show(q)
        
    def clear(self):
        '''
        DESCRIPTION:
        Clear instanse.
        Probably better create new object for second time.
        '''
        self.ee = zeros((self.TN, self.Nx, self.Ny))
        self.ps = zeros((self.TN, self.Nx, self.Ny))
        self.psG = zeros((self.Nx, self.Ny))
        self.f = zeros((self.Nx, self.Ny))
        self.u = zeros((self.TN, self.Nx, self.Ny))
        self.v = zeros((self.TN, self.Nx, self.Ny))
        # self.D1 = zeros((self.Nx, self.Ny))
        # self.D2 = zeros((self.Nx, self.Ny))
        self.set_S()

    def set_S(self):
        '''
        DESCRIPTION:
        Create barrier like two rectange one inside other
        (two- because we need differences at any border)
        for binary_dilation see:

        https://en.wikipedia.org/wiki/Dilation_%28morphology%29
        http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.ndimage.morphology.binary_dilation.html
        '''

        Nx = self.Nx
        Ny = self.Ny
        x = self.x
        y = self.y
        hx = self.hx
        hy = self.hy
        D1 = self.D1
        D2 = self.D2
        # d = nd.generate_binary_structure()

        D1[:] = zeros((Nx, Ny))[:]
        D2[:] = zeros((Nx, Ny))[:]
        D1[int(Nx/2), int(self.hh)] = 1
        D1[:] = nd.binary_dilation(D1,
                                   ones((3, self.lenght))).astype(D1.dtype)[:]
        D2[:] = nd.binary_dilation(D1,
                                   ones((3, 3))).astype(D2.dtype)[:]

        if(self.a != 0):
            D1[:] = nd.rotate(D1, self.a, reshape=False)[:]
            D2[:] = nd.rotate(D2, self.a, reshape=False)[:]
            t1 = max([max(D1[i][:]) for i in range(shape(D1)[0])])
            t2 = max([max(D2[i][:]) for i in range(shape(D2)[0])])
            D1[:] = (D1 >= t1/10).astype(D1.dtype)[:]
            D2[:] = (D2 >= t2/10).astype(D2.dtype)[:]

    def set_borders_ps(self, m, cn_ps):
        '''
        DESCRIPTION:
          B3
        B4  B6
          B1
        
        B4 is steam input. It weak near edges (coasts) and strong near center
        so it use parabola equation for speed V=(u(y),0) where u(y)=C*y*(H-y).
        for ps=integrate(u(s), s, 0, y)
        B6 is same like B4
        B1 and B3 we use
        u=dps/dy,
        v=-dps/dx
        =>
        ps=v*x+c=c1*x+c=u0*x+c (u0 is const for init)
        
        INPUT:
        m is time step (for differ m==0 and others)
        cn_ps it choice condition type (use only 0)
        
        '''
        # self.f = ones((self.Nx, self.Ny))
        self.psG = zeros((self.Nx, self.Ny))
        if cn_ps == 0:

            # enter B4
            self.psG[0] = (self.c*(self.y**2*self.H/2.-self.y**3/3.)
                           - ones(self.Ny)*self.constPs0)

            self.psG[-1] = (self.c*(self.y**2*self.H/2.-self.y**3/3.)
                            + self.u0*self.x[-1]*ones(self.Ny)
                            - ones(self.Ny)*self.constPs0)
            # or self.c*(self.y**2*self.H/2.-self.y**3/3.)

            # low B1
            self.psG.T[0] = self.u0*self.x - ones(self.Nx)*self.constPs0

            # hight B3
            self.psG.T[-1] = self.u0*self.x + self.psG[0][-1]*ones(self.Nx)
            # or -ones(self.Nx)*self.constPs0
            # or self.psG.T[-1] = ones(self.Nx)*self.c*self.H**3/6.
        
        if m == 0:
            self.ps[m][:] = self.psG[:]  # for copy, not reference

    def set_borders_u(self, m, cn_u=0):
        '''
        DESCRIPTION:
        See DESCRIPTION for set_borders_ps.
        '''

        u = self.u[m]
        v = self.v[m]
        Nx = self.Nx
        Ny = self.Ny
        x = self.x
        y = self.y
        hx = self.hx
        hy = self.hy

        # low B1
        u.T[0] = self.u0*ones(Nx)
        
        if m == 0:

            # hight B3
            u.T[-1] = self.u0*ones(Nx)
            
        # v.T[0] = zeros(Nx)
        # u.T[Ny/2] = zeros(Nx)
        
        if cn_u == 0:
            
            # enter B4
            u[0][1:-1] = (self.c*y[1:-1]*(ones(Ny)[1:-1]*self.H-y[1:-1])
                          + self.u0*ones(Ny)[1:-1])
  
            if m == 0:
                u[-1][1:-1] = (self.c*y[1:-1]*(ones(Ny)[1:-1]*self.H-y[1:-1])
                               + self.u0*ones(Ny)[1:-1])
            
            # u.T[Ny/2] = zeros(Nx)
            if m != 0:

                # for bad point
                u[0][-1] = self.u[m-1][1][-1]

    def set_borders_ee(self, m, cn_ee=0):
        '''
        DESCRIPTION:
        in case 1 (cn_ee=0)
        for all borders Bi it
        used deffinition of function of vorticity
        ee = du/dy-dv/dx
        '''

        Nx = self.Nx
        Ny = self.Ny
        x = self.x
        y = self.y
        hx = self.hx
        hy = self.hy
        
        # vorticity
        if m != 0:
            # case 1
            if cn_ee == 0:

                # low B1
                self.ee[m].T[0][1:-1] = ((self.u[m-1].T[0][2:]
                                          - self.u[m-1].T[0][1:-1])/float(hy)
                                         - (self.v[m-1].T[1][1:-1]
                                            - self.v[m-1].T[0][1:-1])/float(hx))

                # height B3
                self.ee[m].T[-1][1:-1] = ((self.u[m-1].T[-1][2:]
                                           - self.u[m-1].T[-1][1:-1])/float(hy)
                                          - (self.v[m-1].T[-1][1:-1]
                                             - self.v[m-1].T[-2][1:-1])/float(hx))

                # self.ee[m][0][1:-1] = zeros(Ny)[1:-1]
                # self.ee[m][-1][1:-1] = zeros(Ny)[1:-1]

                # enter B4
                self.ee[m][0][1:-1] = ((self.u[m-1][0][2:]
                                        - self.u[m-1][0][1:-1])/float(hy)
                                       - (self.v[m-1][1][1:-1]
                                          - self.v[m-1][0][1:-1])/float(hx))
                
                # exit B4
                self.ee[m][-1][1:-1] = ((self.u[m-1][-1][2:]
                                         - self.u[m-1][-1][1:-1])/float(hy)
                                        - (self.v[m-1][-1][1:-1]
                                           - self.v[m-1][-2][1:-1])/float(hx))
                
                # changing only hight points where ee[] cannot be calculs
                self.ee[m].T[0][0] = self.ee[m].T[0][1]
                self.ee[m].T[-1][0] = self.ee[m].T[-1][1]
                self.ee[m].T[0][-1] = self.ee[m].T[0][-2]
                self.ee[m].T[-1][-1] = self.ee[m].T[-1][-2]
                self.ee[m][self.D2.astype(bool)] = 2*self.ps[m-1][self.D2.astype(bool)]/float(self.hy)
            
    def test_pde1(self):
        x = self.x
        y = self.y

        # create a grid for F:R^2->R after which
        # f(xu,yu) will be see like
        # f[0]=f(0,y)= 1line,
        # f[1]=f(1,y) = 2line
        # and so on;
        # f.T[0]=f(x,0)=1colum,
        # f.T[1]=f(x,1)=2colum
        # and so on
        yu, xu = meshgrid(y, x)
 
        # it will be transpose view and
        # don't need a F.T in next steps (see previous line)
        # u0 = phi(xu, yu)
        # phi[0] = phi(0, y);
        # phi[1]=phi(1,y);
        # ...
        psG = self.psG_test(xu, yu)
        self.psG_test_m = psG
        self.psG[0] = psG[0]
        self.psG[-1] = psG[-1]
        self.psG.T[0] = psG.T[0]
        self.psG.T[-1] = psG.T[-1]
        self.f = self.f_test(xu, yu)
        self.pde1(0)

    def pde2(self):
        '''
        DESCRIPTION:
        It is main function for all program.
        first it init borders and init conditions.
        loop n:
             step1: it solve ee (using ee,u,v from step n-1).
             step2: find ps from equation laplace(ps) = ee (ee from step1).
             step3: then find u,v (u=d(ps)/dy, v=-d(ps)/dx) use ps from step 2.
        '''
        Re = self.Re  # 50
        ee = self.ee
        ps = self.ps
        Nx = self.Nx
        Ny = self.Ny  
        x = self.x
        y = self.y
        hx = self.hx
        hy = self.hy
        # ht = self.ht
        ht = self.ht0
        u = self.u
        v = self.v

        self.set_borders_ps(0, self.cn_ps)
        self.set_borders_u(0, self.cn_u)
        # self.set_borders_ee(0)
        
        for n in range(self.TN)[1:]:
            
            self.n = n
            if n % 10 == 0:
                print("n = %d" % n)

            # solve step 1
            self.set_borders_ee(n, self.cn_ee)
            for i in range(Nx)[1:]:  # i
                for j in range(Ny)[1:]:  # j
                    if (i != (Nx-1) and j != (Ny-1)):

                        # barier conditions
                        if self.D1[i, j] != 0:
                            ee[n][i][j] = 0
                        else:
                            if self.D2[i, j] != 0:
                                pass
                                # ee[n][i][j] = 2*ps[n-1][i][j]/float(hy)
                            else:
                                dx = (ee[n-1][i][j]-ee[n-1][i-1][j])/float(hx)
                                dy = (ee[n-1][i][j]-ee[n-1][i][j-1])/float(hy)
                                Dxx = (ee[n-1][i-1][j]-2*ee[n-1][i][j]+ee[n-1][i+1][j])/float(hx**2)
                                Dyy = (ee[n-1][i][j-1]-2*ee[n-1][i][j]+ee[n-1][i][j+1])/float(hy**2)
                                ee[n][i][j] = (ee[n-1][i][j]
                                               + ht*(- u[n-1][i][j]*dx
                                                     - v[n-1][i][j]*dy
                                                     + Re*(Dxx+Dyy)))
            
            # solve step 2
            self.f = ee[n]
            self.set_borders_ps(n-1, self.cn_ps)
            
            self.pde1(n)

            # solve step 3
            self.set_borders_u(n, self.cn_u)

            for i in range(Nx)[1:]:  # i
                for j in range(Ny)[1:]:  # j
                    # if (i != (Nx-1) and j != (Ny-1)):
                    if self.D1[i, j] != 0:
                        v[n][i][j] = u[n][i][j] = 0
                    else:
                        v[n][i][j] = -(ps[n][i][j]-ps[n][i-1][j])/float(hx)
                        u[n][i][j] = (ps[n][i][j]-ps[n][i][j-1])/float(hy)
        # mesh(v[-1][-1]) or quiver(v[0][0], v[0][1], v[1][0], v[1][1])

    def pde1(self, m):
        '''
        DESCRIPTION:
        Solve Dxx+Dyy=f(x,y)
        TN >= 2*N**2/pi**2
        '''
        Nx = self.Nx
        Ny = self.Ny
        hx = self.hx
        hy = self.hy
        ht = self.ht
        TN = self.TN_pde1
        u = zeros((TN, Nx, Ny))
        uG = self.psG
                
        u[0] = uG
        for n in range(TN)[1:]:
            for i in range(Nx):
                if i == 0 or i == Nx-1:  # or i==1 or i==Nx-2#!!!!!
                    # contain borders with -1
                    u[n][i] = uG[i]  # -ones(Ny)
                else:
                    for j in range(Ny):
                        if j == 0 or j == Ny-1 or self.D1[i, j] != 0:
                            # contain borders with -1
                            u[n][i][j] = uG[i][j]  # -1
                        else:
                            Dxx = (u[n-1][i-1][j] - 2*u[n-1][i][j] + u[n-1][i+1][j])/float(hx**2)
                            Dyy =(u[n-1][i][j-1] - 2*u[n-1][i][j] + u[n-1][i][j+1])/float(hy**2) 
                            u[n][i][j] = (u[n-1][i][j]
                                          + ht*(Dxx+Dyy-self.f[i][j]))
        self.ps[m] = u[-1]
        # return ((xu,yu),(hx,hy,ht),u[TN-1],uG)
        # mesh(r[4])
