
def Solar_Torque(A_bus,r_b,r_a,h,theta,CR):
    import math
    import matplotlib.pyplot as plt
    import numpy as np
    from sympy import sin, cos, tan, pi

    CR = 2
    GscoverR2 = 590 # /m^2
    c = 3*10**6 # m/s

    A_bus = 2.045 # m^2
    r_b = 0.480 # m
    r_a = 1.04938 # m
    h = 1.050 # m

    def cross_prod(a, b):
        result = [a[1]*b[2] - a[2]*b[1],
                a[2]*b[0] - a[0]*b[2],
                a[0]*b[1] - a[1]*b[0]]

        return result

    #%% Rotating Force unit vector 

    n_F2 = []
    for i in theta:
        if  0 < i < pi/2:
            for i in theta:
                n_F2.append([sin(i), 0, cos(i)])
        elif pi/2 < i < pi:
            for i in theta:
                n_F2.append([sin(i), 0, -cos(i)])
        elif pi < i < 3/2*pi:
            for i in theta:
                n_F2.append([sin(i), 0, -cos(i)])
        elif 3/2*pi < i < 2*pi:
            for i in theta:
                n_F2.append([sin(i), 0, cos(i)])

    n_F = np.array(n_F2)


    # %% Bus
    phiQ2 = []
    phiQ3 = []
    phiQ4 = []
    for i in theta:
        phiQ2.append(np.pi-i)
        phiQ3.append(i-np.pi)
        phiQ4.append(2*np.pi-i)
                
    r_bus = np.array([0.028367, -0.015365, -0.232788])

    A_bus2 = []
    r_bus2 = []
    for i in theta:
        if 0 < i < pi/2 or pi < i < 3/2*pi:
            for i in theta:
                y = r_a - r_b/cos(i) - h*sin(i) + r_b*tan(i)*sin(i)
                x = y / sin(i)
                if A_bus*sin(i) - 2*r_b* y >= 0:
                    A_bus2.append(A_bus*sin(i) - 2*r_b* y)
                    x_COA = -x/2 # delta z in catia
                    deltaRbus = np.array([0, 0, x_COA])
                    r_bus2.append(np.add(r_bus, deltaRbus))
                else:
                    A_bus2.append(0)
                    r_bus2.append(r_bus)
            
        elif pi/2 < i < pi:
            for i in phiQ2:
                y = r_a - r_b/cos(i) - h*sin(i) + r_b*tan(i)*sin(i)
                x = y / sin(i)
                if A_bus*sin(i) - 2*r_b* y >= 0:
                    A_bus2.append(A_bus*sin(i) - 2*r_b* y)
                    x_COA = -x/2 # delta z in catia
                    deltaRbus = np.array([0, 0, x_COA])
                    r_bus2.append(np.add(r_bus, deltaRbus))
                else:
                    A_bus2.append(0)
                    r_bus2.append(r_bus)
            
        elif pi < i < 3/2*pi:
            for i in phiQ3:
                y = r_a - r_b/cos(i) - h*sin(i) + r_b*tan(i)*sin(i)
                x = y / sin(i)
                if A_bus*sin(i) - 2*r_b* y >= 0:
                    A_bus2.append(A_bus*sin(i) - 2*r_b* y)
                    x_COA = -x/2 # delta z in catia
                    deltaRbus = np.array([0, 0, x_COA])
                    r_bus2.append(np.add(r_bus, deltaRbus))
                else:
                    A_bus2.append(0)
                    r_bus2.append(r_bus)
            
        elif 3/2*pi < i < 2*pi:
            for i in phiQ4:
                y = r_a - r_b/cos(i) - h*sin(i) + r_b*tan(i)*sin(i)
                x = y / sin(i)
                if A_bus*sin(i) - 2*r_b* y >= 0:
                    A_bus2.append(A_bus*sin(i) - 2*r_b* y)
                    x_COA = -x/2 # delta z in catia
                    deltaRbus = np.array([0, 0, x_COA])
                    r_bus2.append(np.add(r_bus, deltaRbus))
                else:
                    A_bus2.append(0)
                    r_bus2.append(r_bus)

    F_bus = []
    for i in range(len(theta)):
        F_bus.append(n_F[i] * CR * GscoverR2 / c * A_bus2[i])


    tau_bus = []
    for i in range(len(theta)):
        tau_bus.append(cross_prod(r_bus2[i], F_bus[i]))

    #%% Antenna

    A_ant = 3.47 # m^2
    r_ant = np.array([0.028367, -0.015365, 1.827212])
    F_ant0 = np.array(CR * GscoverR2 / c * A_ant)

    F_ant = []
    for i in range(len(theta)):
        F_ant.append(n_F[i] * F_ant0)

    tau_ant = []
    for i in range(len(theta)):
        tau_ant.append(cross_prod(r_ant, F_ant[i]))

    #%% Arrays

    A_array = 2.015 # m^2
    r_array1 = np.array([0.028367, 1.765072, -0.344661])
    r_array2 = np.array([0.028367, -1.795802, -0.344661])

    F_array0 = np.array(CR * GscoverR2 / c * A_array)

    F_array = []
    for i in range(len(theta)):
        F_array.append(n_F[i] * F_array0)

    tau_array1 = []
    tau_array2 = []
    for i in range(len(theta)):
        tau_array1.append(cross_prod(r_array1, F_array[i]))
        tau_array2.append(cross_prod(r_array2, F_array[i]))

    #%% Total torque per axis

    # THIS IS PART OF A VERY BIG COMMENT
    # THIS IS PART OF A VERY BIG COMMENT
    # THIS IS PART OF A VERY BIG COMMENT
    # THIS IS PART OF A VERY BIG COMMENT
    # THIS IS PART OF A VERY BIG COMMENT

    tau_total_x = []
    tau_total_y = []
    tau_total_z = []
    ST_vector_list = []
    for i in range(len(theta)):
        ST_vector_list.append([(tau_array1[i][0] + tau_array2[i][0] + tau_bus[i][0] + tau_ant[i][0]),
                              (tau_array1[i][1] + tau_array2[i][1] + tau_bus[i][1] + tau_ant[i][1]),
                              (tau_array1[i][2] + tau_array2[i][2] + tau_bus[i][2] + tau_ant[i][2])])
        
    # THIS IS PART OF A VERY BIG COMMENT
    # THIS IS PART OF A VERY BIG COMMENT
    # THIS IS PART OF A VERY BIG COMMENT
    # THIS IS PART OF A VERY BIG COMMENT
    # THIS IS PART OF A VERY BIG COMMENT
    # THIS IS PART OF A VERY BIG COMMENT
        
    return ST_vector_list