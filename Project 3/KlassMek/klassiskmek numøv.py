import matplotlib.pyplot as plt
import numpy as np

from matplotlib.backends.backend_pdf import PdfPages

def calc_r(x_1,y_1,x_2,y_2):
    return np.sqrt((x_1 - x_2)**2 + (y_1 - y_2)**2)

def calc_F(x_1,y_1,x_2,y_2):
    beta = 2
    r = calc_r(x_1,y_1,x_2,y_2)
    GM = (4*(np.pi**2))
    F_x = (GM*(x_1 - x_2))/(r**(beta+1))
    F_y = (GM*(y_1 - y_2))/(r**(beta+1))
    return F_x,F_y

#M_1 is the mass of the object influenced by the force, in astronomical units (earth masses)
#M_2 is the mass of the object influencing M_1
def calc_F_general(x_1,y_1,x_2,y_2,M_1,M_2):
    M_s = 333480 #Mass of the sun in earth masses
    beta = 2
    r = calc_r(x_1,y_1,x_2,y_2)
    GM = (4*(np.pi**2))*(M_2/M_s)
    F_x = (GM*(x_1 - x_2))/(r**(beta+1))
    F_y = (GM*(y_1 - y_2))/(r**(beta+1))
    return F_x,F_y

def Euler_cromer_earth(x_1,y_1,x_2,y_2,vx,vy,dt):
    F_x,F_y = calc_F(x_1,y_1,x_2,y_2)
    vx_ny = vx - F_x*dt
    vy_ny = vy - F_y*dt
    
    x_ny = x_1 + vx_ny*dt
    y_ny = y_1 + vy_ny*dt
    
    return x_ny,y_ny,vx_ny,vy_ny
    
def Euler_cromer_general(x_1,y_1,x_2,y_2,M_1,M_2,vx,vy,dt):
    F_x_1,F_y_1 = calc_F_general(x_1,y_1,x_2,y_2,M_1,M_2)
    F_x_2,F_y_2 = calc_F_general(x_1,y_1,0,0,M_1,333480)
    
    vx_ny = vx - F_x_1*dt - F_x_2*dt 
    vy_ny = vy - F_y_1*dt - F_y_2*dt 
    
    x_ny = x_1 + vx_ny*dt
    y_ny = y_1 + vy_ny*dt
    
    return x_ny,y_ny,vx_ny,vy_ny

#For earth: eccentricity e = 0.0167
# a = semimajor axis of orbit (in astronomical units, for earth a = 1.000)
# e = eccentricity
# M_p = mass of planet (in earth masses)
# M_s = mass of the sun (in earth masses, = 333 480)
def get_initial_cond(a,e,M_p,M_s):
    GM = 4*(np.pi**2)
    v_p = np.sqrt(GM)*np.sqrt(((1+e)/(a*(1-e)))*(1 + M_p/M_s))#Max, at perihilion
    v_a = np.sqrt(GM)*np.sqrt(((1-e)/(a*(1+e)))*(1 + M_p/M_s))#Min,at aphelion
    
    R_p = a*(1-e) #Min, at perihilion
    R_a = a*(1+e) #Max, at aphelion
    
    return v_p,v_a,R_p,R_a

def plot_radius(r,t,label):
    ax2 = plt.xlim(0,max(t))
    ax2 = plt.ylim(min(r)*0.9,max(r)*1.1)
    ax2 = plt.xlabel('$t$ (years)',fontsize=15)
    ax2 = plt.ylabel('$r$ (AU)',fontsize=15)
    
    ax2 = plt.plot(t,r,'-r',Linewidth=2,label=label)
    ax2 = plt.legend()
    
    r_1 = r[0]
    r_2 = r[-1:]
    
    delta_r = np.sqrt((r_2 - r_1)**2)
    
    plt.show()
    print(min(r))

def plot_radius_tbp(r_1,r_2,t,label_1,label_2):
    ax2 = plt.xlim(0,max(t))
    ax2 = plt.ylim(min(r_1)*0.9,max(r_2)*1.25)
    ax2 = plt.xlabel('$t$ (years)',fontsize=15)
    ax2 = plt.ylabel('$r$ (AU)',fontsize=15)
    
    ax2 = plt.plot(t,r_1,'-b',Linewidth=2,label=label_1)
    ax2 = plt.plot(t,r_2,'-r',Linewidth=2,label=label_2)
    ax2 = plt.legend()
    
    plt.show()
    #print(min(r))

def plot_area(dA,t):
    ax3 = plt.xlim(0,max(t))
    ax3 = plt.ylim(min(dA)*0.9,max(dA)*1.1)
    ax3 = plt.xlabel('$t$ (years)',fontsize=14)
    ax3 = plt.ylabel('$dA$ (AU$^2$)',fontsize=14)
    
    ax3 = plt.plot(t,dA,'-k',label='Area swept per $\Delta t$')
    ax3 = plt.legend()
    
    plt.show()

def plot_conservationtest(E,L,t):
    ax4 = plt.xlim(0,max(t))
    ax4 = plt.ylim(min(E)*0.9,max(E)*1.1)
    ax4 = plt.xlabel('t (years)',fontsize=14)
    ax4 = plt.ylabel('Energy ()',fontsize=14)
    
    ax4 = plt.plot(t,E,'-g',label='Energy')
    ax4 = plt.legend()    
    
    plt.show()       
    
    ax5 = plt.xlim(0,max(t))
    ax5 = plt.ylim(min(L)*0.9,max(L)*1.1)
    ax5 = plt.xlabel('t (years)',fontsize=14)
    ax5 = plt.ylabel('Angular Momentum ()',fontsize=14) 
    
    ax5 = plt.plot(t,L,'-y',label='Angular Momentum')
    ax5 = plt.legend()
    
    plt.show()

def keplers_second_law(r,vx,vy,dt):
    v = np.sqrt((vx**2) + (vy**2))
    d_theta = (v/r)*dt
    dA = 0.5*(r**2)*d_theta
    return dA
    
def check_for_conservation(r,vx,vy):
    m = 1    
    
    GM = (4*(np.pi**2))    
    
    T = (0.5)*m*((vx**2)+(vy**2))
    V = (-GM)/(r)
    
    E = T + V
    
    L = m*r*np.sqrt((vx**2) + (vy**2))
    
    return E , L

def earth_orbit():
    #Time in years, 1 year gives one period
    dt = 0.0001 
    t = 0
    end = 1
    
    #Computes values for velocity and radius at perihelion and apohelion
    v_p,v_a,R_p,R_a = get_initial_cond(1,0.0167,1,333480)    
    
    #Setting initial conditions
    x = R_a
    y = 0
    vx = 0
    vy = -v_a
    
    ax = plt.subplot2grid((6, 6), (0, 0), colspan=6, rowspan=6)    
    
    limit = 1.3
    ax = plt.xlim(-limit,limit)
    ax = plt.ylim(-limit,limit)
    ax = plt.xlabel('$x$ (AU)',fontsize=15)
    ax = plt.ylabel('$y$ (AU)',fontsize=15)
    
    ax = plt.plot(0,0,'oy',MarkerSize=10)
    
    plt.text(1,1,'t=%d')
    
    x_list = [x]
    y_list = [y]
    
    r_list = [np.sqrt((x**2)+(y**2))] #List of values for radius r
    t_list = [t] #List of values for time t
    dA_list = [keplers_second_law(np.sqrt((x**2)+(y**2)),vx,vy,dt)] #List of values for every part area swept for check of Keples second law
    E,L = check_for_conservation(np.sqrt((x**2)+(y**2)),vx,vy) #Finds values for total mechanical energy and angular momentum
    E_list = [E]
    L_list = [L]
    
    while t < end:
        
        r_1 = calc_r(x,y,0,0)         #Hvorfor denne utregningen?
        
        x,y,vx,vy = Euler_cromer_earth(x,y,0,0,vx,vy,dt)
        
        r_2 = calc_r(x,y,0,0)        #Vet ikke ass
        
        t += dt
        x_list.append(x)
        y_list.append(y)
        
        r_list.append(np.sqrt((x**2) + (y**2)))
        t_list.append(t)
        dA_list.append(keplers_second_law(np.sqrt((x**2)+(y**2)),vx,vy,dt))
        E,L = check_for_conservation(np.sqrt((x**2)+(y**2)),vx,vy)
        E_list.append(E)
        L_list.append(L)
    with PdfPages('earth orbit.pdf') as pdf:
        plt.plot(x_list,y_list,'-b',linewidth=2,label='Earth orbit')
        plt.legend(bbox_to_anchor=(1.05,1))
        pdf.savefig()
    plt.show()
    
    with PdfPages('earth radius.pdf') as pdf:
        fig2 = plt.figure()
        plot_radius(r_list,t_list,'Earth radius') #Plot of radius as function of time for check of eccentricity
        pdf.savefig(fig2)
    plt.show()
    
    with PdfPages('kepler.pdf') as pdf:
        fig2 = plt.figure()
        plot_area(dA_list,t_list) #Plot of radius as function of time for check of eccentricity
        pdf.savefig(fig2)
    plt.show()
    #plot_area(dA_list,t_list) #Plot of area swept every dt for check of Keplers second law
    #plot_conservationtest(E_list,L_list,t_list)
    
def three_body_problem():
    
    dt = 0.0001 
    t = 0
    end = 20
    
    #Mass of jupiter
    M_j = 317.89    
    
    #Mass of earth and mars
    M_e = 1
    #M_e = M_j
    M_m = 0.1074
    #M_m = M_j
    
    #epihilion and perihilion velocities and radius of earth
    ve_p,ve_a,Re_p,Re_a = get_initial_cond(1,0.0167,M_e,333480)
    
    #epihilion and perihilion velocities and radius of mars
    vm_p,vm_a,Rm_p,Rm_a = get_initial_cond(1.5237,0.0934,M_m,333480)
    
    #initial conditions for earth
    x_e = Re_a
    y_e = 0
    vx_e = 0
    vy_e = -ve_a
    
    #initial conditions for mars
    x_m = Rm_a
    y_m = 0
    vx_m = 0
    vy_m = -vm_a
    
    ax = plt.subplot2grid((6, 7), (0, 0), colspan=6, rowspan=7)    
    
    limit = 2
    ax = plt.xlim(-limit,limit)
    ax = plt.ylim(-limit,limit)
    ax = plt.xlabel('$x$ (AU)',fontsize=14)
    ax = plt.ylabel('$y$ (AU)',fontsize=14)
    
    text = '$t$ = %d year' % end
    text += '\n' + r'$\Delta t$' + ' = %d' % dt
    #plt.text(1.47,-2,text, bbox={'facecolor':'white'})
    
    
    ax = plt.plot(0,0,'oy',MarkerSize=10)
    
    x_e_list = [x_e]
    y_e_list = [y_e]    
    
    x_m_list = [x_m]
    y_m_list = [y_m]     
    
    r_e_list = [np.sqrt((x_e**2)+(y_e**2))]
    r_m_list = [np.sqrt((x_m**2)+(y_m**2))]  
    
    t_list = [t]    
    
    while t < end:
        x_e_ny,y_e_ny,vx_e,vy_e = Euler_cromer_general(x_e,y_e,x_m,y_m,M_e,M_m,vx_e,vy_e,dt)
        x_m_ny,y_m_ny,vx_m,vy_m = Euler_cromer_general(x_m,y_m,x_e,y_e,M_m,M_e,vx_m,vy_m,dt)
        
        t += dt
        
        x_e = x_e_ny
        y_e = y_e_ny  
        x_m = x_m_ny
        y_m = y_m_ny
        
        x_e_list.append(x_e)
        y_e_list.append(y_e)
        x_m_list.append(x_m)
        y_m_list.append(y_m)
        
        r_e_list.append(np.sqrt((x_e**2)+(y_e**2)))
        r_m_list.append(np.sqrt((x_m**2)+(y_m**2)))
        t_list.append(t)
    
    with PdfPages('threebodyproblem.pdf') as pdf:  
        ax = plt.plot(x_e_list,y_e_list,'-b',linewidth=1,label='Earth orbit')
        ax = plt.plot(x_m_list,y_m_list,'-r',linewidth=1,label='Mars orbit')
    
        #ax = plt.title('Simulation of planetary orbits',fontsize=15)
        ax = plt.legend(fontsize=10,bbox_to_anchor=(1.01,1))
        
        pdf.savefig()
    plt.show()
    
    with PdfPages('tbpradius.pdf') as pdf:
        ax2 = plt.figure()
        plot_radius_tbp(r_e_list,r_m_list,t_list,'Earth radius','Mars radius')
        pdf.savefig(ax2)
        plt.show()
        
#with PdfPages('numøv.pdf') as pdf:
#    ax = plt.figure()
#    earth_orbit()
#    pdf.savefig(ax)
#    plt.show()

three_body_problem()