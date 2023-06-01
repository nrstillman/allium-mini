import numpy as np
class Parameters(object):
    def __init__(self, p):
        for key, values in p.items():
            setattr(self, key, values)

# global param #required for pickling
class param:
    def __init__(self, framerate, Lx, Ly, R ):
        self.Lx = Lx
        self.Ly = Ly
        self.R = R
        self.framerate = framerate #given in hours
        self.dt = 0.001
        self.output_time = int(self.framerate/self.dt)

class ExperimentData:

    def gettypes(self, readtypes, frames):
        # 0 is posn data (only)
        # 1 is full tracks
        # 2 is non-full tracks
        # 4 is array padding 
        boollabels =  np.isin(self.ptype[frames],readtypes)
        return np.nonzero(boollabels)[0]
        
    def truncateto(self,start, endtime):
        if not self.truncated:
            self.Nsnap = endtime - start
            self.flag =  self.flag[start:endtime]
            self.rval = self.rval[start:endtime]
            self.vval = self.vval[start:endtime]
            # self.radius = self.radius[start:endtime]
            self.ptype = self.ptype[start:endtime]
            self.truncated = True
            self.Nvals = self.Nvals[start:endtime]
        else:
            print("Already truncated. Skipping this step")

    # Subtracting avg posn and vel and each timepoint
    def takeDrift(self):
        if not self.drift_removed:

            for t in range(self.Nsnap):
                rtype = self.gettypes([0,1,2],t)
                rdrift = self.rval[t,rtype,:].mean(axis=0) 
                self.rval[t,rtype,:] -= rdrift

            for t in range(self.Nsnap-1):
                vtype = self.gettypes([1,2],t)
                vdrift = self.vval[t,vtype,:].mean(axis=0)
                self.vval[t,vtype,:] -= vdrift

            self.drift_removed = True
        else:
            print("Drift already removed. Skipping this step")

    def __init__(self,data,properties, takeDrift=True,framerate=0.166, Lx = 800, Ly =800, R = 10, umpp=0.8,minlen=10,debug=False):

        self.Nvariable = True
        self.truncated = False
        self.drift_removed = False
        self.param = param(framerate, Lx, Ly, R )
        self.umpp = umpp #microns per pixel
        self.sigma = R

        #shift all flags by one to standardise w python syntax
        data[:,0] -= 1

        self.Nvals = []
        for t in range(properties['t'].min(),properties['t'].max()):
            print(t,end='\r')
            self.Nvals.append(len(np.where(data[:,1] == t)[0]))
        
        flags = np.unique(data[:,0]).astype(int)  
        # count number of occurances of each flag
        self.flag_lifespan = np.bincount(data[:,0].astype(int))

        self.Nparticles = len(flags)
        self.Tinit = properties['t'].min()
        self.Tfinal = properties['t'].max()
        self.Nsnap = (self.Tfinal - self.Tinit) +1

        # number of flags per timestep
        # self.flags_per_timestep = np.asarray([len(flag) for flag in flag_new])
        # self.maxN = self.flags_per_timestep.max()

        self.flag = np.zeros((self.Nsnap, self.Nparticles)).astype(int)
        self.rval = np.zeros((self.Nsnap, self.Nparticles,2))
        self.vval = np.zeros((self.Nsnap-1, self.Nparticles,2))

        # tracking data (set all to 3 to account for padding)
        self.ptype = np.ones((self.Nsnap, self.Nparticles)).astype(int)*4
        self.Ntracers = 0
        self.Nvuse = 0
        for counter, f in enumerate(flags):
            
            # get all rvalues related to a flag and turn into a numpy array
            rval = data[:,2:][data[:,0] == f]*umpp 
            rval[:,0] -= Lx/2
            rval[:,1] -= Ly/2
            # get all time values related to a flag
            time = (data[:,1][data[:,0] == f] - self.Tinit).astype(int)
            if len(time)<1:
                print("Error: empty track detected!")
                break

            self.rval[time, counter] = rval
            self.flag[time, counter] = f
            # set everyone to the ptype for position existing
            self.ptype[time, counter] = 0
            
            if len(time) > 1:
                # time snap difference between ends
                dtime = (time[-1] - time[0])
                # if this is one shorter than the array, this is a complete track
                # i.e. diff will result in a valid velocity field
                vel = np.diff(rval, axis=0)/framerate
                self.vval[time[:-1], counter] = vel
                if (len(time)-1-dtime) == 0:
                    # locally useful velocity data on array one short
                    self.ptype[time[:-1] , counter] = 2
                    self.Nvuse+=1
                # complete track, will be used for time correlations
                # label 1, is a tracer
                if len(time)==self.Nsnap:
                    self.ptype[time[:-1], counter] = 1
                    self.Ntracers += 1
                        
        print(f"Total number of frames = {self.Nsnap}")
        print(f"Total number of tracers/particles = {self.Ntracers}/{self.Nparticles} ({self.Ntracers/self.Nparticles:.3})\n")
        print(f"Total number of velocity use/particles = {self.Nvuse}/{self.Nparticles} ({self.Nvuse/self.Nparticles:.3})\n")

class SimData:
    def ApplyPeriodic2d(self,dr):
        dr[:,0]-=self.param.Lx*np.round(dr[:,0]/self.param.Lx)
        dr[:,1]-=self.param.Ly*np.round(dr[:,1]/self.param.Ly)
        return dr

    def checkTypes(readtypes,data):
        #check which particles to load 
        if len(readtypes) > 0:
            usetypes = np.isin(data[:,-1],readtypes)
        else:
            usetypes = [True]*len(data)
        return usetypes

    # Data object for summary statistics
    def __init__(self,**kwargs):
        # check for debugging
        try:
            self.debug = kwargs['debug']
            if self.debug:
                print('kwargs: ', kwargs)
        except:
            self.debug = False
        # check for disp_field opt
        try:
            self.disp_field = kwargs['disp_field']
        except:
            self.disp_field = True
        # check for specific loadtimes
        try:    
            self.start = kwargs["loadtimes"][0]
            self.end = kwargs["loadtimes"][1]
            self.multiopt = True
        except:
            self.multiopt = False
        # check for specific types
        try:
            self.readtypes = kwargs["readtypes"]
        except:
            self.readtypes = []
        # load parameters
        try:    
            self.param = Parameters(kwargs['params'])
        except:
            print('Error! Parameters must be a dictionary')
            return 1
        # load multiple simulation snapshots
        if self.multiopt:

            data = kwargs['data']

            all_flags = []
            for dt in data:
                all_flags = np.append(all_flags,dt[:,0])
            all_flags = np.unique(all_flags) + 1

            self.Nvuse = 0
            self.Ntracers = 0
            self.Nsnap = len(data)
            self.Nparticles = int(max(all_flags) + 1)

            self.ptype = np.ones((self.Nsnap, self.Nparticles)).astype(int)*4
            self.flag = np.zeros((self.Nsnap, self.Nparticles))
            self.rval = np.zeros((self.Nsnap, self.Nparticles,2))
            self.vval = np.zeros((self.Nsnap, self.Nparticles,2))
            self.theta = np.zeros((self.Nsnap, self.Nparticles))
            self.radius = np.zeros((self.Nsnap, self.Nparticles))
            self.Z = np.zeros((self.Nsnap, self.Nparticles))
            self.sigma = self.param.R
            self.Nvals = [len(dt) for dt in data]
            
            for t, dt in enumerate(data):
                self.flag[t, dt[:,0].astype(int)] = dt[:,0] + 1
                self.rval[t, dt[:,0].astype(int), :] = dt[:,1:3]
                self.theta[t, dt[:,0].astype(int)] = dt[:,5].reshape(-1)
                self.radius[t, dt[:,0].astype(int)] = dt[:,6].reshape(-1)
                self.Z[t, dt[:,0].astype(int)] = dt[:,8].reshape(-1)
                if not self.disp_field:
                    self.vval[t, dt[:,0].astype(int), :] = dt[:,3:5]

            if self.disp_field:
                for counter, f in enumerate(all_flags):
                    # get all rvalues related to a flag and turn into a numpy array
                    ind = (self.flag == f).any(axis=1)

                    tmprval = self.rval[ind,int(f)-1,:]

                    # get all time values related to a flag
                    time = np.linspace(0,self.Nsnap -1,self.Nsnap)[ind].astype(int)

                    if len(time)<1:
                        print("Error: empty track detected!")
                        break

                    #set everyone to the ptype for position existing
                    self.ptype[time, counter] = 0
                    if len(time) > 1:
                        # time snap difference between ends
                        dtime = (time[-1] - time[0])
                        # if this is one shorter than the array, this is a complete track
                        # i.e. diff will result in a valid velocity field
                        dr = np.diff(tmprval, axis=0)
                        vel = self.ApplyPeriodic2d(dr)/self.param.framerate
                        self.vval[time[:-1], int(f)-1] = vel
                        if (len(time)-1-dtime) == 0:
                            # locally useful velocity data on array one short
                            self.ptype[time[:-1], int(f)-1] = 2
                            self.Nvuse+=1
                        # complete track, will be used for time correlations
                        # label 1, is a tracer
                        if len(time) == self.Nsnap:
                            self.ptype[time, int(f)-1] = 1
                            self.Ntracers += 1

        # or a single snapshot
        else:
            # only get particles we're interestsed in
            usetypes = SimData.checkTypes(self.readtypes, kwargs['data'])
            self.Ntrack = sum(usetypes)
            #check whether data is old or new style
            if kwargs['data'].shape[1] > 4:
                #new output includes v,theta,radius
                self.flag =  kwargs['data'][usetypes,0]
                self.rval = kwargs['data'][usetypes,1:3]
                self.vval = kwargs['data'][usetypes,3:5]
                self.theta = kwargs['data'][usetypes,5]
                self.nval = np.array([np.cos(self.theta), np.sin(self.theta)]).T
                self.radius = kwargs['data'][usetypes,6]
                self.ptype = kwargs['data'][usetypes,7]
                self.Z = kwargs['data'][usetypes,8]
            else:
                #old output only contains flag, r and type
                self.flag =  kwargs['data'][usetypes,0]
                self.rval = kwargs['data'][usetypes,1:3]
                self.ptype = kwargs['data'][usetypes, 3]
        
                # For defect tracking
                self.vnorm = np.sqrt(self.vval[:,0]**2 + self.vval[:,1]**2+self.vval[:,2]**2)
                self.vhat = self.vval / np.outer(vnorm,np.ones((3,)))
                
                self.N = len(radius)
                self.sigma = np.mean(radius)
                print("New sigma is " + str(self.sigma))

    def gettypes(self, readtypes, frames):
        return np.isin(self.ptype[frames],readtypes)
        
    def truncateto(self,start, endtime):
        self.Nsnap = endtime - start
        self.flag =  self.flag[start:endtime]
        self.rval = self.rval[start:endtime]
        self.vval = self.vval[start:endtime]
        self.theta = self.theta[start:endtime]
        self.radius = self.radius[start:endtime]
        self.ptype = self.ptype[start:endtime]
        self.Nvals = self.Nvals[start:endtime]
        
    def spatialcut(self,minL=-400, maxL=400, dim=0):
        for t in range(self.Nsnap):
            cut_indices = (self.rval[t][:,dim] < minL) | (self.rval[t][:,dim] > maxL)
            self.flag[t][cut_indices] = [0]
            self.vval[t][cut_indices,:] = [0,0]
            self.theta[t][cut_indices] = [0]
            self.nval[t][cut_indices,:] = [0,0]
            self.radius[t][cut_indices] = [0]
            self.ptype[t][cut_indices] = [0]
            self.rval[t][cut_indices,:] = [0,0]
            self.Z[t][cut_indices,:] = [0,0]