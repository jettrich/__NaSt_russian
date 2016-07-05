from basis import *

class NSExtend(basis.NS):
    '''DESCRIPTION:
    some extend cases for borders
    for each method see description in 
    parent class
    '''
    def set_borders_ps(self,m,cn_ps):
        #self.f =ones((self.Nx,self.Ny))
        self.psG = zeros((self.Nx,self.Ny))
        if cn_ps==0:
            self.psG[0] =self.c*(self.y**2*self.H/2.-self.y**3/3.)-ones(self.Ny)*self.constPs0#enter B4
            self.psG[-1]=self.c*(self.y**2*self.H/2.-self.y**3/3.)+self.u0*self.x[-1]*ones(self.Ny)-ones(self.Ny)*self.constPs0#self.c*(self.y**2*self.H/2.-self.y**3/3.)
            self.psG.T[0]=self.u0*self.x-ones(self.Nx)*self.constPs0#low B1
            self.psG.T[-1]=self.u0*self.x+self.psG[0][-1]*ones(self.Nx)#-ones(self.Nx)*self.constPs0#hight B3
            #self.psG.T[-1]=ones(self.Nx)*self.c*self.H**3/6.
           
        if cn_ps==1:
            self.psG[0][1:-1] =self.u00*self.y[1:-1]-self.constPs1*ones(self.Ny)[1:-1]
            
            self.psG.T[0]=self.u0*self.x-self.constPs1*ones(self.Nx)#low B1
            #self.psG.T[self.Ny/2]=zeros(self.Nx)
            self.psG.T[-1]=self.u0*self.x+self.u00*self.y[-1]*ones(self.Nx)-self.constPs1*ones(self.Nx)#hight B3
            self.psG[-1][1:-1] =self.u00*self.y[1:-1]+self.psG.T[0][-1]*ones(self.Ny)[1:-1]#-self.constPs1*ones(self.Ny)[1:-1]
            
            #self.psG.T[-1]=ones(self.Nx)*self.c*self.H**3/6.
        #self.psG[self.D2.astype(bool)]=0
        ##self.psG[1][1:-1]=self.psG[0][1:-1]#!!!
        ##self.psG[-2][1:-1]=self.psG[-1][1:-1]#!!!
        if m==0:               
            self.ps[m][:] = self.psG[:]#for copy, not reference

    
    def set_borders_u(self,m,cn_u=0):
        u = self.u[m]
        v = self.v[m]
        Nx=self.Nx
        Ny=self.Ny
        x = self.x
        y=self.y
        hx = self.hx
        hy = self.hy
        u.T[0]=self.u0*ones(Nx)# low B1
        if m==0:
            u.T[-1]=self.u0*ones(Nx)# hight B3
            
        #v.T[0]=zeros(Nx)
        #u.T[Ny/2]=zeros(Nx)
        if cn_u==0:
            u[0][1:-1]=self.c*y[1:-1]*(ones(Ny)[1:-1]*self.H-y[1:-1])+self.u0*ones(Ny)[1:-1]# enter B4  
            if m==0:
                u[-1][1:-1]=self.c*y[1:-1]*(ones(Ny)[1:-1]*self.H-y[1:-1])+self.u0*ones(Ny)[1:-1]
            
            #u.T[Ny/2]=zeros(Nx)
            if m!=0:
                u[0][-1]=self.u[m-1][1][-1]#for bad point
        if cn_u==1:
            u[0][1:-1]=self.u00*ones(Ny)[1:-1]# enter B4
            if m==0:
                u[-1][1:-1]=self.u00*ones(Ny)[1:-1]
            if m!=0:
                u[0][-1]=self.u[m-1][1][-1]#for bad point
            #u.T[Ny/2]=u.T[Ny/2-1]
            #u[-1][1:-1]=self.u00*ones(Ny)[1:-1]

    
    def set_borders_ee(self,m,cn_ee=0):
        Nx=self.Nx
        Ny=self.Ny
        x = self.x
        y=self.y
        hx = self.hx
        hy = self.hy
        #vorticity
        if m!=0:
            #case 1
            if cn_ee ==0:
                self.ee[m].T[0][1:-1]=(self.u[m-1].T[0][2:]-self.u[m-1].T[0][1:-1])/float(hy)-(self.v[m-1].T[1][1:-1]-self.v[m-1].T[0][1:-1])/float(hx)#low B1
                self.ee[m].T[-1][1:-1]=(self.u[m-1].T[-1][2:]-self.u[m-1].T[-1][1:-1])/float(hy)-(self.v[m-1].T[-1][1:-1]-self.v[m-1].T[-2][1:-1])/float(hx)#height B3
                ##self.ee[m][0][1:-1]=zeros(Ny)[1:-1]
                ##self.ee[m][-1][1:-1]=zeros(Ny)[1:-1]
                self.ee[m][0][1:-1]=(self.u[m-1][0][2:]-self.u[m-1][0][1:-1])/float(hy)-(self.v[m-1][1][1:-1]-self.v[m-1][0][1:-1])/float(hx)#enter B4
                self.ee[m][-1][1:-1]=(self.u[m-1][-1][2:]-self.u[m-1][-1][1:-1])/float(hy)-(self.v[m-1][-1][1:-1]-self.v[m-1][-2][1:-1])/float(hx)
                #changing only hight points where ee[] cannot be calculs
                self.ee[m].T[0][0]=self.ee[m].T[0][1]
                self.ee[m].T[-1][0]=self.ee[m].T[-1][1]
                self.ee[m].T[0][-1]=self.ee[m].T[0][-2]
                self.ee[m].T[-1][-1]=self.ee[m].T[-1][-2]
                self.ee[m][self.D2.astype(bool)]=2*self.ps[m-1][self.D2.astype(bool)]/float(self.hy)
            if cn_ee==1:
                #case 2        
                #self.ee[m][0]=(self.ps[m][0]+self.ps[m][2]-2*self.ps[m][1])/float(hx**2)+ones(Ny)*self.uu0#enter B4        
                #self.ee[m][-1]=self.ee[m-1][-2]#(2*ps[n][1])/float(hx**2)##exit B6
                self.ee[m][0]=self.c*(self.H-2*y)-(self.v[m-1][1]-self.v[m-1][0])/float(hx)#+(self.ps[m-1][0]+self.ps[m-1][2]-2*self.ps[m-1][1])/float(hx**2)#enter B4
                self.ee[m][-1]=self.c*(self.H-2*y)-(self.v[m-1][-2]-self.v[m-1][-1])/float(hx)
                self.ee[m].T[0][1:-1]=2*(self.ps[m-1].T[1][1:-1]-self.ps[m-1].T[0][1:-1]-self.u0*hy*ones(Nx)[1:-1])/float(hy**2)-(self.v[m-1].T[1][1:-1]-self.v[m-1].T[0][1:-1])/float(hx)#low B1
                self.ee[m].T[-1][1:-1]=2*(self.ps[m-1].T[-2][1:-1]-self.ps[m-1].T[-1][1:-1]-self.u0*hy*ones(Nx)[1:-1])/float(hy**2)-(self.v[m-1].T[-1][1:-1]-self.v[m-1].T[-2][1:-1])/float(hx)#height B3  
                
                #self.ee[m][-1]=self.ee[m-1][-2]#(2*ps[n][1])/float(hx**2)##exit B6
            
        #case 2        
        #if m!=0:
            #self.ee[m][0]=(3*self.ps[m][1])/float(hx**2)-self.ee[m-1][1]/2.
            #self.ee[m][-1]=(3*self.ps[m][-2])/float(hx**2)-self.ee[m-1][-2]/2.        
            ##self.ee[m].T[0][1:]=(u[n-1].T[1][1:]-u[n-1].T[0][1:])/float(hy)
            ##self.ee[m].T[-1]=(3*self.ps[m].T[-2])/float(hy**2)-self.ee[m-1].T[-2]/2.
        #ee[n][-1]=ee[n-1][-2]#(-2*ps[n][-1])/float(hx**2)        
        #ee[n].T[0][:]=(2*ps[n].T[0][:])/float(hx**2)#ee[n-1].T[1][:]##
        #ee[n].T[0][:]=ee[n-1].T[1][:]
        #ee[n].T[-1][:]=(2*ps[n].T[-1][:])/float(hx**2)#ee[n-1].T[-2][:]#
        #ee[n].T[-1][:]=ee[n-1].T[-2][:]
