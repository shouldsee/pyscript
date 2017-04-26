### STOCHASTIC march forward
def march(yiter,ytnow,dt,params):
    ynow,tnow=ytnow;
    yout=yiter(ynow,tnow,params);
    tout=tnow+dt;
    return (yout,tout);


class state():
    def __init__(self,yiter,yt,dt):
        self.yiter=yiter
        self.dt=dt;
        self.yt=np.array(yt);
#         self.t=0;
        self.ys=[yt[0]];
        self.ts=[yt[1]];
#         self.ts=np.array();
    def get_params(self):
        pass
    def forward(self,dur):
        ts=np.arange(self.yt[1],self.yt[1]+dur,self.dt);
        for t in ts:
            self.yt=march(self.yiter,
                          self.yt,
                          self.dt,
                          self.get_params());
#             print(self.ys)
            self.ys=np.vstack((np.array(self.ys) ,
                                np.array(self.yt[0]) ))
        self.ts =self.ts + list(ts);
                    
        
    
def yiter(ynow,t,params):
    fmrna, nfmrna = ynow;
    react= int(random.random() < fmrna/(fmrna+nfmrna));
#     print fmrna/(fmrna+nfmrna);
    fmrna +=-react;
    nfmrna+=react;
    
    
    return  (fmrna,nfmrna)