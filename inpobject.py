# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 11:31:42 2013

@author: jing
"""

# ------------- lissajous test generator --------#
from numpy import *
from pylab import *

def Tmore(T,p,q,omiga,bo,bar_l,bar_w,no):    
    t = linspace(0,T,T)

    a = 100
    b = 100 
    lx = a*sin(p*t + random()*2*pi)+150
    ly = b*sin(q*t + random()*2*pi)+150 
    xx = floor(lx)
    yy = floor(ly)
    x = ones([bo,bo])*range(0,bo)
    y = x.T
    alpha0 = random()*2*pi
    alpha =1e-20+omiga*t+alpha0  
    patterns = zeros((T,size(x)))

    for i in range(T):
        a =  alpha[i]
        
        dl = abs(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-bar_l-1e-20
        dw = abs((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-bar_w-1e-20
        
        dl = -(dl/abs(dl)-1)/2.
        dw = -(dw/abs(dw)-1)/2.
        hinp = abs(dl*dw)
    
        r = abs(cos(a)+1e-20)/(cos(a)+1e-20)
        vdl1 = r*(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-(bar_l) -1e-20
        vdl2 =  r*(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-(bar_l-2*bar_w) -1e-20
        vdw = abs((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-bar_l-1e-20
    
        vdl = vdl1*vdl2
        vdl = -(vdl/abs(vdl)-1)/2.
        vdw = -(vdw/abs(vdw)-1)/2.
        vinp = abs(vdl*vdw)
        
        inp = ((hinp + vinp)>0)*1
        patterns[i] = reshape(inp,size(inp))      
    return patterns
    
def Lmore(T,p,q,omiga,bo,bar_l,bar_w,no):   
    t = linspace(0,T,T)

    a = 100
    b = 100
    lx = a*sin(p*t + random()*2*pi)+150
    ly = b*sin(q*t + random()*2*pi)+150 

    xx = floor(lx)
    yy = floor(ly)
    x = ones([bo,bo])*range(0,bo)
    y = x.T
    alpha0 = random()*2*pi
    alpha =1e-20+omiga*t+alpha0 
    patterns = zeros((T,size(x)))

    for i in range(T):
        a =  alpha[i]
        
        r = abs(cos(a)+1e-20)/(cos(a)+1e-20)
        q = abs(sin(a)+1e-20)/(sin(a)+1e-20)
        hdl = abs(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-bar_l-1e-20
        hdw1 = q*((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-(bar_l)-1e-20
        hdw2 = q*((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-(bar_l-2*bar_w)-1e-20
    
    
        hdl = -(hdl/abs(hdl)-1)/2.
        hdw = hdw1*hdw2    
        hdw = -(hdw/abs(hdw)-1)/2.
        hinp = abs(hdl*hdw)
        
    
        vdl1 = r*(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-(bar_l) -1e-20
        vdl2 =  r*(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-(bar_l-2*bar_w) -1e-20
        vdw = abs((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-bar_l-1e-20
    
        vdl = vdl1*vdl2
        vdl = -(vdl/abs(vdl)-1)/2.
        vdw = -(vdw/abs(vdw)-1)/2.
        vinp = abs(vdl*vdw)
    
        inp = ((hinp + vinp)>0)*1
        patterns[i] = reshape(inp,size(inp))  
    return patterns
 

def Umore(T,p,q,omiga,bo,bar_l,bar_v,bar_w,no):   
    t = linspace(0,T,T)
    a = 100
    b = 100
    lx = a*sin(p*t + random()*2*pi)+150
    ly = b*sin(q*t + random()*2*pi)+150 

    xx = floor(lx)
    yy = floor(ly)
    x = ones([bo,bo])*range(0,bo)
    y = x.T
    alpha0 = random()*2*pi
    alpha =1e-20+omiga*t+alpha0     
    patterns = zeros((T,size(x)))

    for i in range(T):
        a =  alpha[i]
        
        q = abs(sin(a)+1e-20)/(sin(a)+1e-20)
        hdl = abs(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-bar_l-1e-20
        hdw1 = q*((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-(bar_l)-1e-20
        hdw2 = q*((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-(bar_l-2*bar_w)-1e-20
       
        hdl = -(hdl/abs(hdl)-1)/2.
        hdw = hdw1*hdw2    
        hdw = -(hdw/abs(hdw)-1)/2.
        hinp = abs(hdl*hdw)        
    
        vdl1 = abs(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-(bar_v) -1e-20
        vdl2 = abs(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-(bar_v-2*bar_w) -1e-20
        vdw = abs((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-bar_l-1e-20
    
        vdl = vdl1*vdl2
        vdl = -(vdl/abs(vdl)-1)/2.
        vdw = -(vdw/abs(vdw)-1)/2.
        vinp = abs(vdl*vdw)
    
        inp = hinp + vinp -hinp*vinp
        patterns[i] = reshape(inp,size(inp))  
    return patterns

def Arrmore(T,p,q,omiga,bo,bar_l,bar_w,no):   
    t = linspace(0,T,T)
    a = 100
    b = 100
    lx = a*sin(p*t + random()*2*pi)+150
    ly = b*sin(q*t + random()*2*pi)+150 

    xx = floor(lx)
    yy = floor(ly)
    x = ones([bo,bo])*range(0,bo)
    y = x.T
    alpha0 = random()*2*pi
    alpha =1e-20+omiga*t+alpha0 
 
    patterns = zeros((T,size(x)))

    for i in range(T):
        a =  alpha[i]
        
        q = abs(sin(a)+1e-20)/(sin(a)+1e-20)
        hdl = abs(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-bar_l-1e-20
        hdw1 = q*((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-(bar_l)-1e-20
        hdw2 = q*((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-(bar_l-2*bar_w)-1e-20
    
    
        hdl = -(hdl/abs(hdl)-1)/2.
        hdw = hdw1*hdw2    
        hdw = -(hdw/abs(hdw)-1)/2.
        hinp = abs(hdl*hdw)
        
        r = abs(cos(a)+1e-20)/(cos(a)+1e-20)
        vdl1 = r*(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-(bar_l) -1e-20
        vdl2 =  r*(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-(bar_l-2*bar_w) -1e-20
        vdw = abs((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-bar_l-1e-20
        
        vdl = vdl1*vdl2
        vdl = -(vdl/abs(vdl)-1)/2.
        vdw = -(vdw/abs(vdw)-1)/2.
        vinp = abs(vdl*vdw)        
        
        b = a - pi/4.    
        tdl = abs(tan(b)*x-y+ yy[i]-tan(b)*xx[i])/sqrt(tan(b)**2+1)-(bar_l)-1e-20
        tdw = abs((-1/tan(b))*x-y+yy[i]+xx[i]/tan(b))/sqrt(1/tan(b)**2+1)-bar_w-1e-20
    
        tdl = -(tdl/abs(tdl)-1)/2.
        tdw = -(tdw/abs(tdw)-1)/2.
        tinp = abs(tdl*tdw)
    
        inp = ((hinp + vinp +tinp)>0)*1
        patterns[i] = reshape(inp,size(inp))  
    return patterns    
def Kmore(T,p,q,omiga,bo,bar_l,bar_w,no):   
    t = linspace(0,T,T)
    a = 100
    b = 100
    lx = a*sin(p*t + random()*2*pi)+150
    ly = b*sin(q*t + random()*2*pi)+150 

    xx = floor(lx)
    yy = floor(ly)
    x = ones([bo,bo])*range(0,bo)
    y = x.T
    alpha0 = random()*2*pi
    alpha =1e-20+omiga*t+alpha0   
    patterns = zeros((T,size(x)))

    for i in range(T):
        a =  alpha[i]
        
        q = abs(sin(a)+1e-20)/(sin(a)+1e-20)
        hdl = abs(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-bar_l-1e-20
        hdw1 = q*((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-(bar_l)-1e-20
        hdw2 = q*((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-(bar_l-2*bar_w)-1e-20
    
    
        hdl = -(hdl/abs(hdl)-1)/2.
        hdw = hdw1*hdw2    
        hdw = -(hdw/abs(hdw)-1)/2.
        hinp = abs(hdl*hdw)
        
        r = abs(cos(a)+1e-20)/(cos(a)+1e-20)
        vdl1 = r*(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-(bar_l) -1e-20
        vdl2 =  r*(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-(bar_l-2*bar_w) -1e-20
        vdw = abs((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-bar_l-1e-20
        
        vdl = vdl1*vdl2
        vdl = -(vdl/abs(vdl)-1)/2.
        vdw = -(vdw/abs(vdw)-1)/2.
        vinp = abs(vdl*vdw)        
        
        b = a - pi/4.
        s = abs(cos(b)+1e-20)/(cos(b)+1e-20)

        tdl1 = s*(tan(b)*x-y+ yy[i]-tan(b)*xx[i])/sqrt(tan(b)**2+1)-(sqrt(2)*bar_l) -1e-20
        tdl2 =  s*(tan(b)*x-y+ yy[i]-tan(b)*xx[i])/sqrt(tan(b)**2+1)-(sqrt(2)*bar_l-2*bar_w) -1e-20
        tdw = abs((-1/tan(b))*x-y+yy[i]+xx[i]/tan(b))/sqrt(1/tan(b)**2+1)-sqrt(2)*bar_l-1e-20
        
        tdl = tdl1*tdl2
        tdl = -(tdl/abs(tdl)-1)/2.
        tdw = -(tdw/abs(tdw)-1)/2.
        tinp = abs(tdl*tdw)

    
        inp = ((hinp + vinp +tinp)>0)*1
        patterns[i] = reshape(inp,size(inp))  
    return patterns
def Trimore(T,p,q,omiga,bo,bar_l,bar_w,no):   
    t = linspace(0,T,T)
    a = 100
    b = 100
    lx = a*sin(p*t + random()*2*pi)+150
    ly = b*sin(q*t + random()*2*pi)+150 

    xx = floor(lx)
    yy = floor(ly)
    x = ones([bo,bo])*range(0,bo)
    y = x.T
    alpha0 = random()*2*pi
    alpha =1e-20+omiga*t+alpha0 
    patterns = zeros((T,size(x)))

    for i in range(T):
        a =  alpha[i]

        q = abs(sin(a)+1e-20)/(sin(a)+1e-20)
        hdl = abs(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-bar_l-1e-20
        hdw1 = q*((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-(bar_l)-1e-20
        hdw2 = q*((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-(bar_l-2*bar_w)-1e-20


        hdl = -(hdl/abs(hdl)-1)/2.
        hdw = hdw1*hdw2
        hdw = -(hdw/abs(hdw)-1)/2.
        hinp = abs(hdl*hdw)

        r = abs(cos(a)+1e-20)/(cos(a)+1e-20)
        vdl1 = r*(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-(bar_l) -1e-20
        vdl2 =  r*(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-(bar_l-2*bar_w) -1e-20
        vdw = abs((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-bar_l-1e-20

        vdl = vdl1*vdl2
        vdl = -(vdl/abs(vdl)-1)/2.
        vdw = -(vdw/abs(vdw)-1)/2.
        vinp = abs(vdl*vdw)

        b = a+pi/4.
        tdl = abs(tan(b)*x-y+ yy[i]-tan(b)*xx[i])/sqrt(tan(b)**2+1)-(bar_l)-1e-20
        tdw = abs((-1/tan(b))*x-y+yy[i]+xx[i]/tan(b))/sqrt(1/tan(b)**2+1)-bar_w-1e-20

        tdl = -(tdl/abs(tdl)-1)/2.
        tdw = -(tdw/abs(tdw)-1)/2.
        tinp = abs(tdl*tdw)

        inp = ((hinp + vinp +tinp)>0)*1
        patterns[i] = reshape(inp,size(inp))
    return patterns
def Xmore(T,p,q,omiga,bo,bar_l,bar_w,no):    
    t = linspace(0,T,T)
    a = 100
    b = 100
    lx = a*sin(p*t + random()*2*pi)+150
    ly = b*sin(q*t + random()*2*pi)+150

    xx = floor(lx)
    yy = floor(ly)
    x = ones([bo,bo])*range(0,bo)
    y = x.T
    alpha0 = random()*2*pi
    alpha =1e-20+omiga*t+alpha0   
    patterns = zeros((T,size(x)))

    for i in range(T):
        a =  alpha[i]
        a2 = a+5*pi/12.  
        dl = abs(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-bar_l-1e-20
        dw = abs((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-bar_w-1e-20
        
        dl = -(dl/abs(dl)-1)/2.
        dw = -(dw/abs(dw)-1)/2.
        hinp = abs(dl*dw)
    
        vdl =  abs(tan(a2)*x-y+ yy[i]-tan(a2)*xx[i])/sqrt(tan(a2)**2+1)-bar_l-1e-20   
        vdw =  abs((-1/tan(a2))*x-y+yy[i]+xx[i]/tan(a2))/sqrt(1/tan(a2)**2+1)-bar_w-1e-20
        vdl = -(vdl/abs(vdl)-1)/2.
        vdw = -(vdw/abs(vdw)-1)/2.
        vinp = abs(vdl*vdw)
        
        inp = hinp + vinp -hinp*vinp
        patterns[i] = reshape(inp,size(inp))      
    nmatrix = random((T,size(x)))<no
    patterns = abs(patterns - nmatrix)
    return patterns
    
def Hmore(T,p,q,omiga,bo,bar_l,bar_v,bar_w,no):
    t = linspace(0,T,T)
    a = 100
    b = 100
    lx = a*sin(p*t + random()*2*pi)+150
    ly = b*sin(q*t + random()*2*pi)+150 

    xx = floor(lx)
    yy = floor(ly)
    x = ones([bo,bo])*range(0,bo)
    y = x.T
    alpha0 = random()*2*pi
    alpha =1e-20+omiga*t+alpha0  
    patterns = zeros((T,size(x)))

    for i in range(T):
        a = alpha[i]
        dl = abs(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-bar_l-1e-20
        dw = abs((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-bar_w-1e-20

        dl = -(dl/abs(dl)-1)/2.
        dw = -(dw/abs(dw)-1)/2.
        hinp = abs(dl*dw)

        vdl1 = abs(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-(bar_l+bar_w) -1e-20
        vdl2 =  abs(tan(a)*x-y+ yy[i]-tan(a)*xx[i])/sqrt(tan(a)**2+1)-(bar_l-bar_w) -1e-20
        vdw = abs((-1/tan(a))*x-y+yy[i]+xx[i]/tan(a))/sqrt(1/tan(a)**2+1)-bar_v-1e-20

        vdl = vdl1*vdl2
        vdl = -(vdl/abs(vdl)-1)/2.
        vdw = -(vdw/abs(vdw)-1)/2.
        vinp = abs(vdl*vdw)
        inp = hinp + vinp -hinp*vinp
        patterns[i] = reshape(inp,size(inp))
    return patterns

