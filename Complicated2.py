# -*- coding: utf-8 -*-
"""
Created on Wed May 27 21:34:15 2020

@author: stnav
"""

def complicated2():
    import numpy as np
    from numpy import sin, cos
    import matplotlib.pyplot as plt
    from math import pi
    import matplotlib.animation as animation
    import matplotlib.patches as patches
    
    g = 9.81
    m1 = 275
    m2 = 20
    I1 = 28.34
    I2 = 0.94
    
    R = 0.5
    h = 3/8*R
    H = 0.9
    l = 0.73
    rs = 0.094
    
    Amp = 30/180*pi
    T = 2.18
    
    def f(r,t,b):
        beta = b[0]
        fbeta = b[1]
        ffbeta = b[2]
        theta = r[0]
        ftheta = r[1]
        nom1 = -m1*h*sin(theta)*(g+R*ftheta**2) + m2*(H*sin(theta)+l/2*sin(theta+beta))*(g+R*ftheta**2) + m2*l/2*sin(theta+beta)*ftheta*fbeta*R
        nom2 = -ffbeta*(m2*(l/2*R*cos(theta+beta)+l/2*H*cos(beta)+(l/2)**2)+I2)
        nom3 = -ftheta**2*(m1*(2*h*R*sin(theta))+m2*(-2*H*R*sin(theta)-2*l/2*R*sin(theta+beta)))
        nom4 = -fbeta**2*(m2*(-l/2*R*sin(theta+beta)-l/2*H*sin(beta)))
        nom5 = -ftheta*fbeta*(m2*(-2*l/2*R*sin(theta+beta)-2*H*l/2*sin(beta)-2*l/2*R*sin(theta+beta)))
        denom = m1*(R**2+h**2-2*h*R*cos(theta))+I1+m2*(R**2+H**2+(l/2)**2+2*H*R*cos(theta)+2*l/2*R*cos(theta+beta)+2*H*l/2*cos(beta))+I2
        fftheta = (nom1+nom2+nom3+nom4+nom5)/denom
        return np.array([ftheta,fftheta], float)

    tmin = 0
    tmax = 10
    N = 20*(tmax-tmin) 

    theta0 = 0/180*pi
    ftheta0 = 0

    step = (tmax-tmin)/N
    t_pts = np.arange(tmin,tmax,step)
    beta_pts = Amp*sin(2*pi/T*t_pts)
    fbeta_pts = Amp*(2*pi/T)*cos(2*pi/T*t_pts)
    ffbeta_pts = -Amp*(2*pi/T)**2*sin(2*pi/T*t_pts)
    theta_pts = np.array([],float)
    ftheta_pts = np.array([],float)
    fftheta_pts = np.array([],float)
    r = np.array([theta0, ftheta0], float)
    i = 0
    for t in t_pts:
        theta_pts = np.append(theta_pts, r[0])
        ftheta_pts = np.append(ftheta_pts, r[1])
        b = np.array([beta_pts[i], fbeta_pts[i], ffbeta_pts[i]], float)
        fftheta_pts = np.append(fftheta_pts, f(r,t,b)[1])
        i += 1
        k1 = step*f(r,t,b)
        k2 = step*f(r+0.5*k1,t+0.5*step,b)
        k3 = step*f(r+0.5*k2,t+0.5*step,b)
        k4 = step*f(r+k3,t+step,b)
        r += (k1+2*k2+2*k3+k4)/6  

#    x1 = np.array(theta_pts)*(-R) + h*sin(theta_pts)
#    y1 = R - h*cos(theta_pts)
#    x2 = np.array(theta_pts)*(-R) - H*sin(theta_pts) - l/2*sin(theta_pts+beta_pts)
#    y2 = R + H*cos(theta_pts) + l/2*cos(theta_pts+beta_pts)
    
    xo = np.array(theta_pts)*(-R)
    yo = np.ones(N)*R
    xa = xo - rs*cos(theta_pts) 
    ya = yo - rs*sin(theta_pts)
    xb = xo - H*sin(theta_pts) - rs*cos(theta_pts+beta_pts)
    yb = yo + H*cos(theta_pts) - rs*sin(theta_pts+beta_pts)
        
    Fx = m2*(-R*fftheta_pts-1/2*l*cos(theta_pts+beta_pts)*(fftheta_pts+ffbeta_pts)+H*sin(theta_pts)*ftheta_pts**2-H*cos(theta_pts)*fftheta_pts+1/2*l*sin(theta_pts+beta_pts)*(ftheta_pts+fbeta_pts)**2)
    Fy = m2*(g-1/2*l*sin(theta_pts+beta_pts)*(fftheta_pts+ffbeta_pts)-H*cos(theta_pts)*ftheta_pts**2-H*sin(theta_pts)*fftheta_pts-1/2*l*cos(theta_pts+beta_pts)*(ftheta_pts+fbeta_pts)**2)
    
    plt.figure(1, figsize=(6,4))
    plt.title("Angle vs. t")
    plt.ylabel(r"$\theta(t)$")
    plt.xlabel(r"$t$")
    plt.plot(t_pts,theta_pts,'*-')

    plt.figure(2, figsize=(6,4))
    plt.title("Relative angle vs. t")
    plt.ylabel(r"$\beta(t)$")
    plt.xlabel(r"$t$")
    plt.plot(t_pts,beta_pts,'*-')

    plt.figure(3, figsize=(6,4))
    plt.title("Horizontal force vs. t")
    plt.ylabel(r"$F_x(t)$")
    plt.xlabel(r"$t$")
    plt.plot(t_pts,Fx,'*-')

    plt.figure(4, figsize=(6,4))
    plt.title("Vertical force vs. t")
    plt.ylabel(r"$F_y(t)$")
    plt.xlabel(r"$t$")
    plt.plot(t_pts,Fy,'*-')
    
    fig = plt.figure(5)
    ax = fig.add_subplot(111, autoscale_on=False, ylim=[-0.5,R+H+l+0.5], xlim=[-(R+H+l),(R+H+l)], title="Animation of Roly Poly Toy")
    ax.set_aspect('equal')
    ax.grid()

    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    
    w1 = patches.Wedge((xo[0],yo[0]), R, 180+theta_pts[0]/pi*180, theta_pts[0]/pi*180,color='b')
    r1 = patches.Rectangle((xa[0],ya[0]), 2*rs, H,color='b')
    r2 = patches.Rectangle((xb[0],yb[0]), 2*rs, l,color='r')

    def init():
        time_text.set_text('')
        ax.plot([-10,10],[0,0],'-',color='black')
        ax.add_patch(w1)
        ax.add_patch(r1)
        ax.add_patch(r2)
        return w1,r1,r2,time_text


    def animate(i):
        w1.set_center((xo[i],yo[i]))
        w1.theta1 = 180 + theta_pts[i]/pi*180
        w1.theta2 = theta_pts[i]/pi*180
        r1.xy = xa[i], ya[i]
        r2.xy = xb[i], yb[i]
        r1.angle = theta_pts[i]/pi*180
        r2.angle = (theta_pts[i]+beta_pts[i])/pi*180
        time_text.set_text(time_template % (i*step))
        return w1,r1,r2,time_text
    
    ani = animation.FuncAnimation(fig, animate, range(1, N),
                              interval=step*1000, blit=True, init_func=init)
    ani.save('RolyPoly_complicated2.mp4')
    
complicated2()