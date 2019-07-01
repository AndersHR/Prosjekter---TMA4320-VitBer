import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

from sys import path
from os.path import dirname

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('solar system.pdf') 
def calculate_r(x_1,y_1,x_0,y_0):
    return np.sqrt((x_1 - x_0)**2 + (y_1 - y_0)**2)

def get_initial_conditions(a,e,M_p,M_s):
    GM = 4*(np.pi**2)
    v_p = np.sqrt(GM)*np.sqrt(((1+e)/(a*(1-e)))*(1 + M_p/M_s))#Max, at perihilion
    v_a = np.sqrt(GM)*np.sqrt(((1-e)/(a*(1+e)))*(1 + M_p/M_s))#Min,at aphelion
    
    R_p = a*(1-e) #Min, at perihilion
    R_a = a*(1+e) #Max, at aphelion
    
    return v_p,v_a,R_p,R_a

def get_initial_lists(a_list,e_list,mass):
    M_sun = 333480
    positions_perihelion = []
    positions_aphelion = []
    velocities_perihelion = []
    velocities_aphelion = []
    even = True
    for n in range(len(mass)):
        v_p,v_a,R_p,R_a = get_initial_conditions(a_list[n],e_list[n],mass[n],M_sun)
        if even == True:
            a = 1
            even = False
        else:
            a = -1
            even = True
            
        positions_perihelion.append(a*R_p)
        positions_perihelion.append(0)
        positions_aphelion.append(a*R_a)
        positions_aphelion.append(0)
        velocities_perihelion.append(0)
        velocities_perihelion.append(-(a*v_p))
        velocities_aphelion.append(0)
        velocities_aphelion.append(-(a*v_a))
    
    return velocities_perihelion,velocities_aphelion,positions_perihelion,positions_aphelion
        

def calculate_centralforce(x_1,y_1,x_0,y_0,M_1,M_0):
    M_s = 333480 #Mass of the sun in earth masses
    beta = 2
    r = calculate_r(x_1,y_1,x_0,y_0)
    GM = (4*(np.pi**2))*(M_1/M_s)
    F_x = (GM*(x_0 - x_1))/(r**(beta+1))
    F_y = (GM*(y_0 - y_1))/(r**(beta+1))
    return F_x,F_y

def Euler_Cromer_method(i,positions,velocities,mass,dt):
    x_0 = positions[2*i]
    y_0 = positions[2*i + 1]    
    vx_0 = velocities[2*i]
    vy_0 = velocities[2*i + 1]
    
    M_0 = mass[i]
    M_sun = 333480
    
    F_x = 0
    F_y = 0    
    
    for k in range(len(mass)):
        if k != i:
            x_1 = positions[2*k]
            y_1 = positions[2*k + 1]
            M_1 = mass[k]
            dF_x, dF_y = calculate_centralforce(x_1,y_1,x_0,y_0,M_1,M_0)
            F_x += dF_x
            F_y += dF_y
        else:
            dF_x, dF_y = calculate_centralforce(0,0,x_0,y_0,M_sun,M_0)
            F_x += dF_x
            F_y += dF_y
    
    vx_new = vx_0 - F_x*dt
    vy_new = vy_0 - F_y*dt 
    
    x_new = x_0 + vx_new*dt
    y_new = y_0 + vy_new*dt
    
    return x_new,y_new,vx_new,vy_new
    
def solar_system():
    # 0:Mercury, 1:Venus, 2:Earth, 3:Eros(asteroid), 4:Mars, 5:Jupiter, 6:Saturn, 7:Uranus
    # 8: Neptune, 9:Pluto, 10:Halley's comet
    semimajor_axis = [0.3871,0.7233,1.000,1.4583,1.5237,5.2028,9.5388,19.191,30.061,39.529,18]
    eccentricity = [0.2056,0.0068,0.0167,0.2230,0.0934,0.0483,0.0560,0.046,0.0100,0.2484,0.967]
    mass = [0.0553,0.8150,1.000,2*(10**(-9)),0.1074,317.89,95.16,14.56,17.15,0.002,10**(-10)]
    
    velocities_perihelion,velocities_aphelion,positions_perihelion,positions_aphelion = get_initial_lists(semimajor_axis,
                                                                                                          eccentricity,mass)
    current_velocities = velocities_perihelion
    current_positions = positions_perihelion
    
    total_positions = [current_positions]    
    total_velocities = [current_velocities] 
    
    t = 0
    dt = 0.001
    end = 500
    
    time = [t]    
    
    while t < end:
        t += dt
        time.append(t)
        
        new_positions = []
        new_velocities = []
        
        for i in range(len(mass)):
            x_new,y_new,vx_new,vy_new = Euler_Cromer_method(i,current_positions,current_velocities,mass,dt)
            new_positions += [x_new,y_new]
            new_velocities += [vx_new,vy_new]
            
        total_positions.append(new_positions)
        total_velocities.append(new_velocities)
        
        current_positions = new_positions
        current_velocities = new_velocities
        
    all_x_coordinates, all_y_coordinates = process_positionlists(total_positions)
    #plot_solar_system(all_x_coordinates,all_y_coordinates,end,dt)
    save_plot(all_x_coordinates,all_y_coordinates,end,dt)
                                                                                              
def process_positionlists(total_positions):
    all_x_coordinates = []
    all_y_coordinates = []    
    
    planets = int(len(total_positions[0])/2)
    
    for k in range(planets):
        all_x_coordinates.append([])
        all_y_coordinates.append([])    

    for positions in total_positions:
        for i in range(0,planets):
            all_x_coordinates[i].append(positions[2*i])
            all_y_coordinates[i].append(positions[2*i + 1])
            
    return all_x_coordinates,all_y_coordinates
        


# 0:Mercury, 1:Venus, 2:Earth, 3:Eros(asteroid), 4:Mars, 5:Jupiter, 6:Saturn, 7:Uranus
# 8: Neptune, 9:Pluto, 10:Halley's comet
def plot_solar_system(all_x_coordinates,all_y_coordinates,end,dt):
    
    fig = plt.figure(dpi=300)

    ax = plt.subplot2grid((6, 7), (0, 0), colspan=6, rowspan=7)
    #ax = plt.title('Numerical Simulation of the Solar System')

    ax = plt.plot(0,0,'oy',MarkerSize=5)
    colors = ['red','yellow','cyan','black','red','orange','green','blue','purple','gray','indigo','.r','.r']
    labels = ['Mercury','Venus','Earth','Eros (asteroid)','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto','''Halley's comet''','test','test']

    time_box = dict(boxstyle='square,pad=0.3',fc='white')
    ax = plt.text(22,-45,('Time lapsed = %d' % end) +
    ' years',fontsize=6,bbox=time_box)

    limit = 50
    ax = plt.xlim(-limit,limit)
    ax = plt.ylim(-limit,limit)
    ax = plt.xlabel('$x$ (AU)',fontsize=8)
    ax = plt.ylabel('$y$ (AU)',fontsize=8)
    ax = plt.rc('xtick',labelsize=8)
    ax = plt.rc('ytick',labelsize=8)

    np.savez('Solar_System_coords',x=all_x_coordinates,y=all_y_coordinates)

    for planet in range(len(all_x_coordinates)):
        x = all_x_coordinates[planet]
        y = all_y_coordinates[planet]
        
        ax = plt.plot(x,y,linestyle='-',color=colors[planet],label=labels[planet])
        ax = plt.plot(x[-1],y[-1],linestyle='-',color=colors[planet],linewidth=1)
    ax = plt.legend(bbox_to_anchor=(1.26,1),fontsize=6)
    
    pp.savefig(ax = plt.gcf)

def save_plot(all_x_coordinates,all_y_coordinates,end,dt):
    with PdfPages('simulation.pdf') as pdf:
        plot_solar_system(all_x_coordinates,all_y_coordinates,end,dt)
        pdf.savefig()
    plt.show()
        
solar_system()
    