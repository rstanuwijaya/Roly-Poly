# -*- coding: utf-8 -*-
"""
Created on Wed May 27 20:54:49 2020

@author: stnav
"""

##simplified model of Roly Poly Toy
def simplified():
    import numpy as np
    from numpy import sin, cos
    import matplotlib.pyplot as plt
    from math import pi
    import matplotlib.animation as animation
    import matplotlib.patches as patches


    R = 0.5
    m = 250
    h = 3/8*R
    Ic = 83/320*m*R**2
    g = 9.81

    def f(r,t):
        theta = r[0]
        ftheta = r[1]
#         fftheta = -(h*m*sin(theta)*(g+R*ftheta**2))/(Ic+m*(R**2+h**2-2*h*R*cos(theta)))
        fftheta = -(h*m*sin(theta)*(g+R*ftheta**2)+2*h*m*R*sin(theta))/(Ic+m*(R**2+h**2-2*h*R*cos(theta)))
        return np.array([ftheta,fftheta], float)

    tmin = 0.0
    tmax = 10.0
    N = 400 

    theta0 = 60/180*pi
    ftheta0 = 0

    step = (tmax-tmin)/N
    t_pts = np.arange(tmin,tmax,step)
    theta_pts = []
    ftheta_pts = []
    r = np.array([theta0, ftheta0], float)
    for t in t_pts:
        theta_pts.append(r[0])
        ftheta_pts.append(r[1])
        k1 = step*f(r,t)
        k2 = step*f(r+0.5*k1,t+0.5*step)
        k3 = step*f(r+0.5*k2,t+0.5*step)
        k4 = step*f(r+k3,t+step)
        r += (k1+2*k2+2*k3+k4)/6  

    hc_pts = R - h*cos(theta_pts)
    xo = np.array(theta_pts)*(-R)
    yo = np.ones(N)*R

    
    plt.figure(1, figsize=(6,4))
    plt.title("Angle vs. t")
    plt.ylabel(r"$\theta(t)$")
    plt.xlabel(r"$t$")
    plt.plot(t_pts,theta_pts,'*-')

    plt.figure(2, figsize=(6,4))
    plt.title("Angular velocity vs. t")
    plt.ylabel(r"$\dot{\theta}(t)$")
    plt.xlabel(r"$t$")
    plt.plot(t_pts,ftheta_pts,'*-')

    plt.figure(3, figsize=(6,4))
    plt.title("Height of cm vs. t")
    plt.ylabel(r"$h_c(t)$")
    plt.xlabel(r"$t$")
    plt.plot(t_pts,hc_pts,'*-')
    
    fig = plt.figure(4)
    ax = fig.add_subplot(111, autoscale_on=False, ylim=[-1,R+1], xlim=[-(R+1),(R+1)], title="Animation of Roly Poly Toy")
    ax.set_aspect('equal')
    ax.grid()

#    line, = ax.plot([], [], 'o-', lw=2)
    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    
    w1 = patches.Wedge((xo[0],yo[0]), R, 180+theta_pts[0]/pi*180, theta_pts[0]/pi*180,color='b')

    def init():
        time_text.set_text('')
        ax.plot([-10,10],[0,0],'-',color='black')
        ax.add_patch(w1)
        return w1,time_text


    def animate(i):
        w1.set_center((xo[i],yo[i]))
        w1.theta1 = 180 + theta_pts[i]/pi*180
        w1.theta2 = theta_pts[i]/pi*180
        time_text.set_text(time_template % (i*step))
        return w1,time_text
    
    ani = animation.FuncAnimation(fig, animate, range(1, N),
                              interval=step*1000, blit=True, init_func=init)
    ani.save('RolyPoly_simplified.mp4')
#    plt.show()
    
simplified()