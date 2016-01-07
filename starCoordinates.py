import csv, math, re
import matplotlib.pyplot as plot
import numpy as np
from matplotlib.patches import Circle
from numpy import linalg as la

########################################################################################
## Attributes: MPG;Cylinders;Displacement;Horsepower;Weight;Acceleration;Model;Origin ##
########################################################################################
def csv2mat(filepath, removed_keys):
    data_mat = []
    with open(filepath, 'rb') as csvfile:
        items = csv.reader(csvfile)
        for item in items:
            item = item and item[0]
            item = item.split(';')[1:]
            data_mat.append([str2float(value) for value in item])
    data_mat = data_mat[1:]
    return data_mat

################################################################
## Input  ==> (Dict) data_dict                                ##
## Output ==> (List) data_mat[size_of_dataSet * num_of_attrs] ## 
################################################################
#def dict2mat(data_dict):
#    mat = []
#    for key in data_dict:
#        vector = []
#        vector.append(key)
#        attrs = data_dict[key]
#        for attr in attrs:
#            value = str2float(attrs[attr])
#            vector.append(value)
#        mat.append(vector)
#    return mat

##############################################################
## Usage: Convert any data with Type string into Type float ##
##############################################################
def str2float(string, maxLength=5):
    try:
        res = float(string)
    except ValueError:
        res = 0
        string = string.lower()
        re.sub('[^a-z0-9]', '', string)
        if len(string) < maxLength:
            maxLength = len(string)
        for i in range(maxLength):
            res += (ord(string[i])-ord('a')+1)*pow(27.0, maxLength-1-i) 
    return res

########################################
## Usage: Regulization in each column ##
########################################
def dataRegulization(data_mat):
    num_of_row = len(data_mat)
    num_of_col = len(data_mat[0])
    data_mat_reg = [[0 for col in range(num_of_col)] for row in range(num_of_row)] 
    col_max = []
    col_min = []
    for i in range(num_of_col):
        col_max.append(max([row[i] for row in data_mat]))
        col_min.append(min([row[i] for row in data_mat]))    
        col_range = col_max[i] - col_min[i] + 0.0001
        for j in range(num_of_row):
            data_mat_reg[j][i] = (data_mat[j][i]-col_min[i])/col_range
    return data_mat_reg 

#################################################################
## Usage: Convert n-Dimensional space to Cartesian Coordinates ##
#################################################################
def nD2cc(num_of_attrs):
    theta = 2*np.pi / num_of_attrs
    cc_cos = ''
    cc_sin = ''
    for i in range(num_of_attrs):
        cc_cos += str(np.cos(theta*(1+i))) + ' '
        cc_sin += str(np.sin(theta*(1+i))) + ' '
    np_cc = np.matrix(cc_cos[:-1]+';'+cc_sin[:-1])

    return np_cc

##################################################################
## Usage: Display the scatter figure to show data in 2D from nD ##
##################################################################
def figure(np_coord, np_data, attrs):
    ## scatters
    np_osc = np.dot(np_coord, np_data.T)
    coord_x = np.array(np_coord[0,:]).tolist()
    coord_y = np.array(np_coord[1,:]).tolist()
    osc_x = np.array(np_osc[0,:]).tolist()
    osc_x = osc_x[0]
    osc_y = np.array(np_osc[1,:]).tolist()
    osc_y = osc_y[0]
    fig = plot.figure()
    ax = fig.add_subplot(111)
    ax.scatter(osc_x, osc_y, color='#5CACEE')
    r_max = 0.0
    dims = np_coord.size / 2
    for i in range(dims):
        ## lines 
        ax.plot([np_coord[0,i],0], [np_coord[1,i],0], linewidth=1.5, color='#FF6A6A')
        ## points
        ax.plot(np_coord[0,i], np_coord[1,i], 'o', color='#7FFFD4')
        ## annotations
        ax.annotate(attrs[i], xy=(np_coord[0,i], np_coord[1,i]))
        ## circles 
        r = math.sqrt(np_coord[0,i]*np_coord[0,i]+np_coord[1,i]*np_coord[1,i])
        if r > r_max:
            r_max = r
        circle = Circle((0.0, 0.0), radius=r, fill=False, color='#B5B5B5')
        ax.add_patch(circle)
    ## axis
    ax.axis([-1.2, 1.2, -1.2, 1.2])
    ax.axis('equal')
    return [ osc_x, osc_y]

#################################
## Usage: Display coords in 2D ##
#################################
def coordsFigure(np_coord, attrs=[]):
    dims = np_coord.size / 2
    fig = plot.figure()
    ax = fig.add_subplot(111)
    for i in range(dims):
        ax.plot([np_coord[0,i],0], [np_coord[1,i],0], linewidth=1.5, color='#FF6A6A')
        ax.plot(np_coord[0,i], np_coord[1,i], 'o', color='#7FFFD4')
        try:
            ax.annotate(attrs[i], xy=(np_coord[0,i], np_coord[1,i]))
        except:
            pass
        r = math.sqrt(np_coord[0,i]**2+np_coord[1,i]**2)
        circle = Circle((0.0, 0.0), radius=r, fill=False, color='#B5B5B5')
        ax.add_patch(circle)
    ax.axis([-1.2, 1.2, -1.2, 1.2])
    ax.axis('equal')
    return 0
      

##############################################
## Usage: Make star coordinates orthgraphic ##
##############################################
def osc(sc_now, Type='em', fix=[], rad=[], direct=[], iterate=9999):
    sc_ori = sc_now.copy()
    num_of_attrs = sc_now.size / 2
    sc_next = np.matrix([[0 for i in range(num_of_attrs)] for j in range(2)], dtype='double') 
    if Type == 'Leh_re':
        coord_x = sc_now[0,:]
        coord_y = sc_now[1,:]
        coord_x = coord_x / la.norm(coord_x)
        coord_y = coord_y - np.sum(np.multiply(coord_y, coord_x))*coord_x
        sc_next[0,:] = coord_x
        sc_next[1,:] = coord_y / la.norm(coord_y)
    elif Type == 're':
        coord_x = sc_now[0,:]
        coord_y = sc_now[1,:]
        print sc_now
        xreg = 0
        yreg = 0
        for i in fix:
            xreg += coord_x[0,i] ** 2;
            yreg += coord_y[0,i] ** 2;
        co_xreg = math.sqrt(1-xreg)
        co_yreg = math.sqrt(1-yreg)
        for i in range(num_of_attrs):
            if i not in fix:
               coord_x[0,i] *= xreg
               coord_y[0,i] *= yreg
        print coord_x, coord_y
        coord_y = coord_y - np.sum(np.multiply(coord_y, coord_x))*coord_x
        sc_next[0,:] = coord_x
        sc_next[1,:] = coord_y / la.norm(coord_y)
        print sc_next[:,7]
    elif Type == 'em':
        ## Backtrack Line Search
        sc_next = sc_now
        np_delta = np.matrix([[0 for i in range(num_of_attrs)] for j in range(2)], dtype='double')
        [cost, xx, yy, xy] = costFunc(sc_now)
        step = 1
        alpha = 0.25
        beta = 0.8
        error = 0.001
        while (cost > error) and (iterate > 0):
            rad_delta_sum = 0.0
            direct_delta_sum = 0.0
            for i in rad:
                rad_delta_sum += (math.sqrt(sc_next[0,i]**2+sc_next[1,i]**2)-math.sqrt(sc_ori[0,i]**2+sc_ori[1,i]**2))**2
            for i in direct:
                direct_delta_sum += ((sc_next[0,i]*sc_ori[0,i]+sc_next[1,i]*sc_ori[1,i])/(math.sqrt(sc_next[0,i]**2+sc_next[1,i]**2)*math.sqrt(sc_ori[0,i]**2+sc_ori[1,i]**2))-1)**2 
            for i in range(num_of_attrs):
                x = sc_next[0,i]
                y = sc_next[1,i]
                sr = math.sqrt(x**2+y**2)
                x0 = sc_ori[0,i]
                y0 = sc_ori[1,i]
                sr0 = math.sqrt(x0**2+y0**2)

                np_delta[0,i] = x * (xx - 1) + 0.5 * xy * x
                np_delta[1,i] = y * (yy - 1) + 0.5 * xy * y

                if i in rad:
                    np_delta[0,i] *= 0.3
                    np_delta[0,i] += 0.7 * rad_delta_sum * x * (1 - sr0 / sr)
                    np_delta[1,i] *= 0.3
                    np_delta[1,i] = 0.7 * rad_delta_sum * y * (1 - sr0 / sr)
                if i in direct:
                    np_delta[0,i] *= 0.3
                    np_delta[0,i] += 0.7 * direct_delta_sum*((x*x0+y*y0)/(sr*sr0)-1)/sr0*(x0*y**2-y0*x*y)/sr**3
                    np_delta[1,i] *= 0.3
                    np_delta[1,i] += 0.7 * direct_delta_sum*((x*x0+y*y0)/(sr*sr0)-1)/sr0*(x0*y**2-y0*x*y)/sr**3
                if i in fix:
                    np_delta[0,i] = 0
                    np_delta[1,i] = 0

            cost_gen = cost
            while cost_gen > cost - alpha*step*(np.sum(np.multiply(np_delta, np_delta))):
                step = beta * step 
                [cost_gen, xx_gen, yy_gen, xy_gen] = costFunc(sc_next-step*np_delta)
            sc_next = sc_next - step * np_delta
            [cost, xx, yy, xy] = [cost_gen, xx_gen, yy_gen, xy_gen]
            iterate = iterate - 1
    return sc_next

## First Radius Next Direciton
def FRND(sc_now, sc_prev, fix):
    num_of_attrs = sc_now.size / 2 
    fx_prev = sc_prev[0,fix]
    fy_prev = sc_prev[1,fix]
    fx = sc_now[0,fix]
    fy = sc_now[1,fix]
    
    theta_prev = math.atan2(fy_prev, fx_prev)
    scale_prev = math.sqrt((fx_prev**2+fy_prev**2))
    theta_now = math.atan2(fy, fx)
    scale_now = math.sqrt((fx**2+fy**2))
    sc_sca = sc_now.copy()
    sc_sca[0,fix] = math.cos(theta_prev) * scale_now  
    sc_sca[1,fix] = math.sin(theta_prev) * scale_now
    osc_sca = osc(sc_sca, 'em', [fix])

    sc_rad = osc_sca.copy()
    sc_rad[0,fix] = math.cos(theta_now) * scale_now
    sc_rad[1,fix] = math.sin(theta_now) * scale_now
    osc_rad = osc(sc_rad, 'em', [fix])
    return [osc_sca, osc_rad]


def costFunc(np_sc):
    num_of_attrs = np_sc.size / 2 
    sum_of_xx = np.sum(np.multiply(np_sc[0,:], np_sc[0,:]))
    sum_of_yy = np.sum(np.multiply(np_sc[1,:], np_sc[1,:]))
    sum_of_xy = np.sum(np.multiply(np_sc[0,:], np_sc[1,:]))
    cost = (sum_of_xx - 1) * (sum_of_xx - 1) + (sum_of_yy - 1) * (sum_of_yy - 1) + sum_of_xy * sum_of_xy
    return [cost, sum_of_xx, sum_of_yy, sum_of_xy]

def coordMove(np_ori, num, x, y):
    np_alpha = np_ori.copy()
    for n in range(len(num)):
        np_alpha[0,num[n]] = x[n]
        np_alpha[1,num[n]] = y[n]
    return np_alpha

def REandEM_demo():
    coords = np.matrix([[0.66, 0.91, 0.38, -0.77, -1.17], [0.83, 0.05, -0.52, -0.09, 1.08]])
    re = osc(coords, 'Leh_re')
    em = osc(coords, 'em')
    coordsFigure(coords,[])
    coordsFigure(re,[])
    coordsFigure(em,[])
    diff(re, em)
    return 0

def diff(co,co_):
    dims = co.size / 2
    fig = plot.figure()
    ax = fig.add_subplot(111)
    for i in range(dims):
        ax.plot([co[0,i],0], [co[1,i],0], linewidth=1.5, color='#FF6A6A')
        ax.plot([co_[0,i],0], [co_[1,i],0], linewidth=1.5, color='#FF6A6A')

        ax.plot([co[0,i],co_[0,i]], [co[1,i],co_[1,i]], linewidth=1.5, color='#FF6A6A', linestyle="--")

        ax.plot(co[0,i], co[1,i], 'o', color='#7FFFD4')
        ax.plot(co_[0,i], co_[1,i], 'o', color='#F98614')
        try:
            ax.annotate(attrs[i], xy=(co[0,i], co[1,i]))
            ax.annotate(attrs[i], xy=(co_[0,i], co_[1,i]))
        except:
            pass
        r = math.sqrt(co[0,i]**2+co[1,i]**2)
        r_ = math.sqrt(co_[0,i]**2+co_[1,i]**2)
        circle = Circle((0.0, 0.0), radius=r, fill=False, color='#B5B5B5')
        circle_ = Circle((0.0, 0.0), radius=r_, fill=False, color='#C1D2D0')
        ax.add_patch(circle)
        ax.add_patch(circle_)
    ax.axis([-1.2, 1.2, -1.2, 1.2])
    ax.axis('equal')
    return 0 
def conditional_demo():
    coords = nD2cc(7)
    coordsFigure(coords)
    coords = coordMove(coords, [1,3,5], [-0.3,-0.7,0.8], [0.7, 0.4, -0.8])
    coordsFigure(coords)
    coords_fix = osc(coords, 'em', fix=[1,4])
    coordsFigure(coords_fix)
    coords_rad = osc(coords, 'em', rad=[1,4])
    coordsFigure(coords_rad)
    coords_direct = osc(coords, 'em', direct=[1,4])
    coordsFigure(coords_direct)
    return 0
def fig_point(ax, coords):
    return 0

def interaction_demo():
    num_of_attrs = 5
    delta_x = 0.4
    delta_y = 0.3
    iterate = 5
    coords = nD2cc(num_of_attrs)
    fig = plot.figure()
    ax = fig.add_subplot(111)
    for i in range(num_of_attrs):
        ax.plot(coords[0,i], coords[1,i], 'o', color='#7FFFD4')
        ax.plot([coords[0,i], 0], [coords[1,i], 0], linewidth=1.5, color='#FF6A6A')
        r = math.sqrt(coords[0,i]**2+coords[1,i]**2)
        circle = Circle((0.0, 0.0), radius=r, fill=False, color='#B5B5B5')
        ax.add_patch(circle)
    coords[0,0] += delta_x
    coords[1,0] += delta_y
    ax.plot([coords[0,0], coords[0,0]-delta_x], [coords[1,0], coords[1,0]-delta_y], linewidth=1.5, color='#FF6A6A', linestyle='--')
    ax.plot(coords[0,0], coords[1,0], 'o', color='#F98614')
    ax.plot([coords[0,0], 0], [coords[1,0], 0], linewidth=1.5, color='#FF6A6A')
    r = math.sqrt(coords[0,0]**2+coords[1,0]**2)
    circle = Circle((0.0, 0.0), radius=r, fill=False, color='#C1D2D0')
    ax.add_patch(circle)
    ax.axis([-1.2, 1.2, -1.2, 1.2])
    ax.axis('equal')
    ## reconditioning 
    coords_re = osc(coords, 'Leh_re')
    diff(coords, coords_re)
    ## energy minimisation
    coords_ori = coords.copy()
    fig_em = plot.figure()
    ax_em = fig_em.add_subplot(111)
    for i in range(num_of_attrs):
        ax_em.plot(coords_ori[0,i], coords_ori[1,i], 'o', color='#7FFFD4')
        ax_em.plot([coords_ori[0,i], 0], [coords_ori[1,i], 0], linewidth=1.5, color='#FF6A6A')
        r = math.sqrt(coords_ori[0,i]**2+coords_ori[1,i]**2)
        circle = Circle((0.0, 0.0), radius=r, fill=False, color='#B5B5B5')
        ax_em.add_patch(circle)
    for it in range(iterate):
        coords_em = osc(coords_ori, 'em', [0], iterate=iterate)
        for i in range(num_of_attrs):
            ax_em.plot(coords_em[0,i], coords_em[1,i], 'o', color='#F98614')
            ax_em.plot([coords_em[0,i], 0], [coords_em[1,i], 0], linewidth=1.5, color='#FF6A6A')
            ax_em.plot([coords_ori[0,i], coords_em[0,i]], [coords_ori[1,i], coords_em[1,i]], linewidth=1.5, color='#FF6A6A', linestyle="--")
        coords_ori = coords_em
    for i in range(num_of_attrs):
        r = math.sqrt(coords_em[0,i]**2+coords_em[1,i]**2)
        circle = Circle((0.0, 0.0), radius=r, fill=False, color='#C1D2D0')
        ax_em.add_patch(circle)
    ax_em.axis([-1.2, 1.2, -1.2, 1.2])
    ax_em.axis('equal')
    return 0

if __name__ == '__main__':
    
    ## Demo of RE and EM Algorithms 
    #REandEM_demo = REandEM_demo() 
    
    ## Demo of Interaction on RE and EM Algorithms
    #interaction_demo()

    ## Demo of Conditional Interaction on EM Algorithms
    conditional_demo()
    data_mat = csv2mat('./cars.csv', 'Car')
    attrs = ['MGP', 'Cylinders', 'Displacement', 'HorsePower', 'Weight', 'Acceleration', 'Model', 'Origin']
    #data_mat = csv2mat('./test.csv', 'Test')
    data_mat_reg = dataRegulization(data_mat)
    np_mat = np.matrix(data_mat_reg)
    num_of_attrs = len(data_mat[0])
    ## nD -> 2D
    np_cc = nD2cc(num_of_attrs)
    ## origin coords
    ##coordsFigure(np_ori, attrs)
    #print costFunc(np_gen)
    ##figure(np_gen, np_mat, attrs)
    #np_gen = osc(np_cc, 'Leh_re')
    #print costFunc(np_gen)
    #figure(np_gen, np_mat, attrs)
    #np_gen = osc(np_cc, 're', [7])
    #figure(np_gen, np_mat, attrs)
    ##[osc_rad, osc_sca] = FRND(np_cc,np_ori, 7)
    ##figure(osc_rad, np_mat, attrs)
    plot.show()
