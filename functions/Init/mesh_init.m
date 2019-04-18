%awesome constructor for my mesh
classdef mesh_init
    
    properties
    Lz,dz,z,zlength,fmax,df,fmin,f,flength,fmid,freal,t,dt,R,dr,r,rlength,wmid,wvl,dwvl,rmin,fbound,rmid
    end
    
   methods 
       function s=mesh_init(beam,dim)
      
        s.Lz=1000e-3;%[m]
        s.dz=10e-3;
        s.z=0:s.dz:s.Lz;
        s.zlength=length(s.z);

        s.fmax=beam.f0*40;%[1/s]
        s.df=2e11;
        s.fmin=s.df*2;
        s.f=-s.fmax:s.df:s.fmax;
        s.flength=length(s.f);
        s.fmid=s.f(round(s.flength/2));
        s.freal=s.f+beam.f0;
        s.fbound=find(s.f>0,1);
        
        s.wvl=const.c./s.f(s.fbound:end);
        s.dwvl=abs(s.wvl(2)-s.wvl(1));
        
        s.t=linspace(-1/(2*s.df),1/(2*s.df),s.flength);%[s]
        s.dt=abs(s.t(2)-s.t(1));
        
        switch dim
            case 1
                s.r=0;
                
            case 2
                s.R=400e-6;%800e-6;%[m]
                s.dr=2e-6;
                s.rmin=-400e-6;
                s.r=s.rmin:s.dr:s.R;%start at r0=3*dr to avoid singularity at r0=0! 
        end
                s.rlength=length(s.r); 
                s.rmid=round(s.rlength/2);
       end

   end
end