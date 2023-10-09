import math
import matplotlib.pyplot as plt
import numpy as np

#----------------ORBIT-DETAILS
G = 6.6732*10**(-11)
Sun = ["Sun", 1.989*10**30, 696340*10**3]
Earth = ["Earth", 5.98*10**24, 6370*10**3]
Mars = ["Mars", 6.417*10**23, 3389.5*10**3]
Celestial_bodies = [Sun,Earth,Mars]

print("Please chhoose your orbiting Celestial body")
for b in range(len(Celestial_bodies)):
    print(str(b)+"- " + Celestial_bodies[b][0])

answer = int(input("Answer: "))
M_e = Celestial_bodies[answer][1]
r_e = Celestial_bodies[answer][2]
miu = M_e * G

#----------------ORBIT-DETAILS

# M_e = 1.989*10**30
# miu = M_e * G
# r_e = 696340*10**3

#--- Earth Specs:
# M_e = 5.98*10**24
# miu = M_e * G
# r_e = 6370000

#--- Mars specs
# M_e = 6.417*10**23
# miu = M_e * G
# r_e = 3389.5*10**3

#------------- ANALYTICAL METHOD TO FIND VELOCITY FOR GIVEN ORBIT
print("This program has two modes: \n 1 – Gives velocities for a given apogee and perigee \n 2 - Gives orbital parameters and energy relations for a chosen velocity at the perigee")
answer = str(input("Please make your choice: "))

if answer == "1":
    Given_perigee = float(input("Insert the perigee in km: "))*10**3
    Target_apogee = float(input("Insert the target apogee in km: "))*10**3

    
    Right_extreme_distance= Given_perigee
    Right_extreme_velocity = 10**3
    error = Target_apogee

    while True:

        E = 1/2*Right_extreme_velocity**2 - miu/Right_extreme_distance
        H = Right_extreme_distance * Right_extreme_velocity
        C = 1/Right_extreme_distance *(1- miu/(Right_extreme_distance * Right_extreme_velocity**2))
        e = C*H**2/miu

        Left_extreme_distance = (-C+miu/H**2)**-1
        Left_extreme_velocity = (2*(E + miu/Left_extreme_distance))**0.5

        if e > 1:
            break

        if abs(Left_extreme_distance-Target_apogee) <= error:
            Best_apogee = Left_extreme_distance
            Best_velocity_perigee = Right_extreme_velocity
            Best_velocity_apogee = Left_extreme_velocity
    
        error = abs(Left_extreme_distance-Target_apogee)
        Right_extreme_velocity += 0.005*10**3

    print("Best value for apogee distance: ", round(Best_apogee/1000,1), "km")
    print("Best value for velocity at perigee: ", round(Best_velocity_perigee/1000,4), "km/s")
    print("Best value for velocity at apogee: ", round(Best_velocity_apogee/1000,4), "km/s")
    key = input("Press enter do leave...")
    quit()
#---------------------- METHOD TO FIND ORBIT DETAILS GIVEN PERIGEE VELOCITY
print("\n----------------------------\n")
print("Right extreme and left extreme terms are used just to identify  the sides of the orbit which might be the perigee or the apogee depending on the specifications you give. It is only necessary to keep it consistent for the correct use of this program. \n") 

Right_extreme_velocity= 1000 * float(input("Insert the right extreme velocity [km/s]: "))
Right_extreme_distance= 1000 * float(input("Insert the right extreme distance [km]: "))
Inclination = math.radians(float(input("Insert the orbital inclination in degrees: ")))


E = 1/2*Right_extreme_velocity**2 - miu/Right_extreme_distance
H = Right_extreme_distance * Right_extreme_velocity
C = 1/Right_extreme_distance *(1- miu/(Right_extreme_distance * Right_extreme_velocity**2))
e = C*H**2/miu

#--- PLOTTING INTERVALS AND FIRST ORBIT DETAILS

Circumference_steps = 200
Delta_theta = (2*math.pi)/Circumference_steps

if e < 1:
    theta_range = int(Circumference_steps)
    theta_range = int(Circumference_steps)
    Left_extreme_distance = (-C+miu/H**2)**-1
    Left_extreme_velocity = (2*(E + miu/Left_extreme_distance))**0.5
    
    if Right_extreme_distance <= Left_extreme_distance:
        Pericenter_distance = Right_extreme_distance
        Pericenter_velocity = Right_extreme_velocity
        Apocenter_distance = Left_extreme_distance
        Apocenter_velocity = Left_extreme_velocity
    
    else:
        Pericenter_distance = Left_extreme_distance
        Pericenter_velocity = Left_extreme_velocity
        Apocenter_distance = Right_extreme_distance
        Apocenter_velocity = Right_extreme_velocity

    a = (Pericenter_distance + Apocenter_distance)/2
    Orbital_Period = 2*math.pi*(a**3/miu)**0.5
    Delta_V = Pericenter_velocity

    print("\n--------Initial orbit -------\n")
    print("Pericenter velocity: ", round(Pericenter_velocity/1000,4), "km/s")
    print("Pericenter distance: ", round(Pericenter_distance/1000,4), "km")
    print("Apocenter velocity: ", round(Apocenter_velocity/1000,4), "km/s")
    print("Apocenter distance: ", round(Apocenter_distance/1000,1), "km")
    print("Excentricity: ", round(e,4))
    print("Semi major axis",  round(a/1000,1), "km")
    print("Orbital_Period :", round(Orbital_Period/3600,4), " h")
    print("Delta V from surface and rest to pericenter: ", round(Delta_V/1000,4), " km/s")
    print("\n----------------------------\n")
    Delta_V = 0

else:
    theta_range = int(Circumference_steps * 0.2)
    print("The orbit is not eliptical")
    
    
#--- 3D PLOT

Orbit_X_lists = []
Orbit_Y_lists = []
Orbit_Z_lists = []

theta_list,X_list,Y_list,Z_list = [],[],[],[]
    
for i in range(theta_range+1):
    theta = Delta_theta * i
    r = (C*math.cos(theta)+miu/H**2)**-1
    X = r * math.cos(theta)
    Y = r * math.sin(theta) * math.cos(Inclination)
    Z = r * math.sin(theta) * math.sin(Inclination)
    
    theta_list.append(theta)
    X_list.append(X)
    Y_list.append(Y)
    Z_list.append(Z)
    
Orbit_X_lists.append(X_list)
Orbit_Y_lists.append(Y_list)
Orbit_Z_lists.append(Z_list)


#---3D-MANEUVERS
Maneuvers_done = False

if e < 1:
    Maneuvers= str(input("Would you like to make a maneuver? \n 1- Yes \n 2- No \n Choose: "))

while e < 1 and Maneuvers == "1":
    Maneuvers_done = True
    Type= str(input("Type: \n 1- Left extreme delta V \n 2- Right extreme delta V \n 3- Inclination change \n Choose: "))
    
    if Type == "1":
        
        #------------------------------------------------ Reset lists and save values
        X_list,Y_list,Z_list = [],[],[]
        
        Velocity_before_maneuver = Left_extreme_velocity
        Velocity_to_perfect_circle = (miu/Left_extreme_distance)**0.5
        print("The required velocity to a perfeclty circular orbit is [km/s]: ",round(Velocity_to_perfect_circle/1000,3))
        
        #------------------------------------------------ Input For new velocity
        print("Your current left extreme velocity is :", round(Left_extreme_velocity/1000,3), " km/s")
        Left_extreme_velocity= 1000*float(input("Insert your new left extreme velocity [km/s]: "))
        
        #------------------------------------------------ New orbital parameters from left
        E = 1/2*Left_extreme_velocity**2 - miu/Left_extreme_distance
        H = Left_extreme_distance * Left_extreme_velocity
        C = 1/Left_extreme_distance *(1- miu/(Left_extreme_distance * Left_extreme_velocity**2))
        e = C*H**2/miu
        
        #------------------------------------------------ Setting interval for plotting
        if e < 1:
            theta_range = int(Circumference_steps)
        else:
            theta_range = int(Circumference_steps * 0.2)
        
        #------------------------------------------------ Creating the data for plotting
        for i in range(theta_range+1):
            theta = Delta_theta * i
            theta = Delta_theta * i
            r = (C*math.cos(theta)+miu/H**2)**-1
            X = - r * math.cos(theta)
            Y = r * math.sin(theta) * math.cos(Inclination)
            Z = r * math.sin(theta) * math.sin(Inclination)
            
            #------------------------------------------------ Storing the coordinates in lists same length
            theta_list.append(theta)
            X_list.append(X)
            Y_list.append(Y)
            Z_list.append(Z)
        
        #------------------------------------------------ Find new right extreme for next maneuver
        if e < 1:
            Right_extreme_distance = (-C+miu/H**2)**-1
            Right_extreme_velocity = (2*(E + miu/Right_extreme_distance))**0.5
            print("Left extreme distance: [km]", round(Left_extreme_distance/1000))
            print("Right extreme distance [km]: ", round(Right_extreme_distance/1000))
            print("Right extreme velocity [km/s]: ", round(Right_extreme_velocity/1000,3))
        
        #------------------------------------------------ Add the work required to the maneuver
        Delta_V += abs(Left_extreme_velocity - Velocity_before_maneuver)
        
        print("The excentricity of your orbit is :", round(e,4))
        
        print("\n----------------------------\n")
        
    if Type == "2":
        
        #------------------------------------------------ Reset lists and save values
        X_list,Y_list,Z_list = [],[],[]
        
        Velocity_before_maneuver = Right_extreme_velocity
        Velocity_to_perfect_circle = (miu/Right_extreme_distance)**0.5
        print("The required velocity to a perfeclty circular orbit is [km/s]: ",round(Velocity_to_perfect_circle/1000,3))
        
        #------------------------------------------------ Input For new velocity
        print("Your current right extreme velocity is :", round(Right_extreme_velocity/1000,3), " km/s")
        Right_extreme_velocity= 1000*float(input("Insert your new right extreme velocity [km/s]: "))
        
        #------------------------------------------------ New orbital parameters from right
        E = 1/2*Right_extreme_velocity**2 - miu/Right_extreme_distance
        H = Right_extreme_distance * Right_extreme_velocity
        C = 1/Right_extreme_distance *(1- miu/(Right_extreme_distance * Right_extreme_velocity**2))
        e = C*H**2/miu
        
        #------------------------------------------------ Setting interval for plotting
        if e < 1:
            theta_range = int(Circumference_steps)
        else:
            theta_range = int(Circumference_steps * 0.2)
        
        #------------------------------------------------ Creating the data for plotting
        for i in range(theta_range+1):
            theta = Delta_theta * i
            r = (C*math.cos(theta)+miu/H**2)**-1
            X = r * math.cos(theta)
            Y = r * math.sin(theta) * math.cos(Inclination)
            Z = r * math.sin(theta) * math.sin(Inclination)
            
            #------------------------------------------------ Storing the coordinates in lists same length
            theta_list.append(theta)
            X_list.append(X)
            Y_list.append(Y)
            Z_list.append(Z)
        
        #------------------------------------------------ Find new right extreme for next maneuver
        if e < 1:
            Left_extreme_distance = (-C+miu/H**2)**-1
            Left_extreme_velocity = (2*(E + miu/Left_extreme_distance))**0.5
            print("Right extreme distance [km]: ", round(Right_extreme_distance/1000))
            print("Left extreme distance [km]: ", round(Left_extreme_distance/1000))
            print("Left extreme velocity [km/s]: ", round(Left_extreme_velocity/1000,3))
        
        #------------------------------------------------ Add the work required to the maneuver
        Delta_V += abs(Right_extreme_velocity - Velocity_before_maneuver)
        print("The excentricity of your orbit is :", round(e,4))
        
        print("\n----------------------------\n")
        
    if Type == "3":
        X_list,Y_list,Z_list = [],[],[]
        
        Inclination_before_maneuver = Inclination
        Inclination = math.radians(float(input("Incert your new inclination in degrees: ")))
        
        E = 1/2*Right_extreme_velocity**2 - miu/Right_extreme_distance
        H = Right_extreme_distance * Right_extreme_velocity
        C = 1/Right_extreme_distance *(1- miu/(Right_extreme_distance * Right_extreme_velocity**2))
        e = C*H**2/miu
        
        for i in range(theta_range+1):
            theta = Delta_theta * i
            r = (C*math.cos(theta)+miu/H**2)**-1
            X = r * math.cos(theta)
            Y = r * math.sin(theta) * math.cos(Inclination)
            Z = r * math.sin(theta) * math.sin(Inclination)
            
            theta_list.append(theta)
            X_list.append(X)
            Y_list.append(Y)
            Z_list.append(Z)
        
        Location = str(input("Where would you like to make the maneuver? \n 1- Left extreme \n 2- Right extreme \n Choose: "))
        
        if Location == "1":
            
            Delta_V += (2*(Left_extreme_velocity**2)*(1-math.cos(abs(Inclination-Inclination_before_maneuver))))**0.5
            
            if Left_extreme_velocity > Right_extreme_velocity:
                print("Are you smartn't? ")

        if Location == "2":

            Delta_V += (2*(Right_extreme_velocity**2)*(1-math.cos(abs(Inclination-Inclination_before_maneuver))))**0.5
            
            if Right_extreme_velocity > Left_extreme_velocity:
                print("Why you doing this? (╯°□°)╯︵ ┻━┻ ")

        print("\n----------------------------\n")
        
    if Type == "3" or Type == "2" or Type == "1":
        Orbit_X_lists.append(X_list)
        Orbit_Y_lists.append(Y_list)
        Orbit_Z_lists.append(Z_list)
        
        if e < 1:
            Maneuvers= str(input("Would you like to make more maneuvers? \n 1- Yes \n 2- No \n Choose: "))
        if e >= 1:
            print("The orbit is not eliptical")
    
    else:
        print("\n Please choose one of the options next time >>:( \n")

#--- EARTH SPHERE

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Make data
u = np.linspace(0, 2 * math.pi, Circumference_steps)
v = np.linspace(0, math.pi, Circumference_steps)
x = r_e * np.outer(np.cos(u), np.sin(v)) 
y = r_e * np.outer(np.sin(u), np.sin(v))
z = r_e * np.outer(np.ones(np.size(u)), np.cos(v))
# outer function make a linear combination of the np arrays "u" and "v". So there will be several arrays of x values for a fixed u and all the v values


ax.plot_surface(x, y, z,color="orange",zorder = 1)

#--- FINAL ORBIT STATS
if Maneuvers_done == True and e<1:
    if Right_extreme_distance <= Left_extreme_distance:
        Pericenter_distance = Right_extreme_distance
        Pericenter_velocity = Right_extreme_velocity
        Apocenter_distance = Left_extreme_distance
        Apocenter_velocity = Left_extreme_velocity
    
    else:
        Pericenter_distance = Left_extreme_distance
        Pericenter_velocity = Left_extreme_velocity
        Apocenter_distance = Right_extreme_distance
        Apocenter_velocity = Right_extreme_velocity
    
    
    E = 1/2*Right_extreme_velocity**2 - miu/Right_extreme_distance
    H = Right_extreme_distance * Right_extreme_velocity
    C = 1/Right_extreme_distance *(1- miu/(Right_extreme_distance * Right_extreme_velocity**2))
    e = C*H**2/miu
    
    a = (Pericenter_distance + Apocenter_distance)/2
    Orbital_Period = 2*math.pi*(a**3/miu)**0.5

    print("\n--------Final orbit -------\n")
    print("Pericenter velocity: ", round(Pericenter_velocity/1000,4), "km/s")
    print("Pericenter distance: ", round(Pericenter_distance/1000,1), "km")
    print("Apocenter velocity: ", round(Apocenter_velocity/1000,4), "km/s")
    print("Apocenter distance: ", round(Apocenter_distance/1000,1), "km")
    print("Excentricity: ", round(e,4))
    print("Semi major axis",  round(a/1000,1), "km")
    print("Inclination :",  math.degrees(Inclination), "Degrees")
    print("Orbital_Period :", round(Orbital_Period/3600,4), " h")
    print("Net Delta V from manuevers: ", round(Delta_V/1000,4), " km/s")
    print("\n----------------------------\n")

#--- OBJECT POSITION VS TIME

if e< 1:
    n = (2*math.pi)/Orbital_Period
    Tracking = str(input("Would you like to know the postion vs time? \n 1- Yes \n 2- No \n Choose: "))
    
    if Tracking == "1":
        Tracking_type = str(input("Would you like to know the time or postion? \n 1- Time \n 2- Position \n Choose: "))
        
        if Tracking_type == "1":
            theta_tracking = math.radians(float(input("Insert the true anomaly in degrees from right extreme: ")))
            
            Complete_orbits  = (theta_tracking / (2*math.pi)) - (theta_tracking / (2*math.pi))%1
            theta_tracking = (theta_tracking / (2*math.pi))%1 * 2*math.pi
            
            
            if abs(theta_tracking) >= math.pi:
                
                E_anomaly = 2*math.pi + math.atan(math.tan(theta_tracking/2)/((1+e)/(1-e))**0.5)*2
                
                print(E_anomaly)
                
            if abs(theta_tracking) < math.pi:
                
                E_anomaly = math.atan(math.tan(theta_tracking/2)/((1+e)/(1-e))**0.5)*2
            
            time_tracking = (E_anomaly - e*math.sin(E_anomaly))/n + Complete_orbits * Orbital_Period
            
            print("Time interval from pericenter in hours: ", round(time_tracking/3600,4))
            
                
            
        if Tracking_type == "2":
            time_tracking = float(input("Insert the time from right extreme in hours: "))*3600
            
            Complete_orbits  = (time_tracking / Orbital_Period) - (time_tracking / Orbital_Period)%1
            time_tracking = (time_tracking / Orbital_Period)%1 * Orbital_Period
            
            E_anomaly = 0
            Iteration_count = int(0)
            
            while abs(E_anomaly-e*math.sin(E_anomaly) - time_tracking * n) > abs(time_tracking * n)/1000:
                
                E_anomaly_attempt = e*math.sin(E_anomaly) + time_tracking * n
                
                E_anomaly = E_anomaly_attempt
                
                Iteration_count += 1
                
            print("Number of iterations: ", Iteration_count)
            
            if abs(E_anomaly) > math.pi:
                
                theta_tracking = 2*math.pi - math.atan(((1+e)/(1-e))**0.5*math.tan(E_anomaly/2))*2
                
                
            if abs(E_anomaly) < math.pi:
                
                theta_tracking = math.atan(((1+e)/(1-e))**0.5*math.tan(E_anomaly/2))*2
            
            print("Change in true anomaly right extreme in degrees: ", round(math.degrees(theta_tracking),4))
            

        r_pt = (C*math.cos(theta_tracking)+miu/H**2)**-1
        X_pt = r_pt * math.cos(theta_tracking)
        Y_pt = r_pt * math.sin(theta_tracking) * math.cos(Inclination)
        Z_pt = r_pt * math.sin(theta_tracking) * math.sin(Inclination)
        
        ax.scatter(X_pt, Y_pt, Z_pt, label= "Object position", color ="black",zorder= 25)

#--- ORBIT PLOTS

for i in range(len(Orbit_X_lists)):
    name = "Orbit: " + str(i+1)
    ax.plot(Orbit_X_lists[i], Orbit_Y_lists[i], Orbit_Z_lists[i], label= name, zorder=5)
    
answer = str(input("The orbit will be plotted, Would you like to hide the axis? \n 1- Yes \n 2- No \n Answer:"))
if answer == "1":
    plt.axis('off')
    plt.grid(b=None)

ax.set_aspect('equal')
ax.legend()

plt.show()
