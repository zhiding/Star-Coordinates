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
    osc_y = np.array(np_osc[1,:]).tolist()
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
    ax.axis('equal')
    ax.axis('scaled')
    return 0


##############################################
## Usage: Make star coordinates orthgraphic ##
##############################################
def osc(sc_now, Type='em', fix=[]):
    sc_ori = sc_now
    num_of_attrs = sc_now.size / 2
    sc_next = np.matrix([[0 for i in range(num_of_attrs)] for j in range(2)], dtype='double') 
    if Type == 'Leh_re':
        coord_x = sc_now[0,:]
        coord_y = sc_now[1,:]
        coord_x = coord_x / la.norm(coord_x)
        coord_y = coord_y - np.sum(np.multiply(coord_y, coord_x))*coord_x
        sc_next[0,:] = coord_x
        sc_next[1,:] = coord_y / la.norm(coord_y)
        print sc_next[:,7]
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
        error = 0.000001
        iterate = 9999
        while (cost > error) and (iterate > 0):
            for i in range(num_of_attrs):
                if i in fix:
                    np_delta[0,i] = 0
                    np_delta[1,i] = 0
                else:
                    np_delta[0,i] = 4 * sc_next[0,i] * (xx - 1) + 2 * xy * sc_next[1,i]
                    np_delta[1,i] = 4 * sc_next[1,i] * (yy - 1) + 2 * xy * sc_next[0,i]
            cost_gen = cost
            while cost_gen > cost - alpha*step*(np.sum(np.multiply(np_delta, np_delta))):
                step = beta * step 
                [cost_gen, xx_gen, yy_gen, xy_gen] = costFunc(sc_next-step*np_delta)
            sc_next = sc_next - step * np_delta
            [cost, xx, yy, xy] = [cost_gen, xx_gen, yy_gen,  xy_gen]
            iterate = iterate - 1
    return sc_next

## First Radius Next Direciton
def FRND(sc_now, sc_prev, fix):
    num_of_attrs = sc_now.size / 2 
    fx = sc_now[0,fix]
    fy = sc_now[1,fix]
    theta = math.atan2(fy, fx)
    scale = math.sqrt((fx**2+fy**2) / (sc_prev[0,fix]**2+sc_prev[1,fix]**2))
    #sc_rad = sc_now.copy()
    #sc_rad[0,fix] = math.cos(theta)
    #sc_rad[1,fix] = math.sin(theta)
    #sc_sca = sc_rad.copy()
    #sc_sca[0,fix] = sc_rad[0,fix] * scale
    #sc_sca[1,fix] = sc_rad[1,fix] * scale
    sc_sca = sc_now.copy()
    sc_sca[0,fix] = sc_now[0,fix] * scale
    sc_sca[1,fix] = sc_now[1,fix] * scale
    sc_rad = sc_sca.copy()
    sc_rad[0,fix] = math.cos(theta)
    sc_rad[1,fix] = math.sin(theta)
    osc_rad = osc(sc_rad, 'em', [fix])
    osc_sca = osc(sc_sca, 'em', [fix])
    return [osc_rad, osc_sca]


def costFunc(np_sc):
    num_of_attrs = np_sc.size / 2 
    sum_of_xx = np.sum(np.multiply(np_sc[0,:], np_sc[0,:]))
    sum_of_yy = np.sum(np.multiply(np_sc[1,:], np_sc[1,:]))
    sum_of_xy = np.sum(np.multiply(np_sc[0,:], np_sc[1,:]))
    cost = (sum_of_xx - 1) * (sum_of_xx - 1) + (sum_of_yy - 1) * (sum_of_yy - 1) + sum_of_xy * sum_of_xy
    return [cost, sum_of_xx, sum_of_yy, sum_of_xy]

if __name__ == '__main__':
    data_mat = csv2mat('./cars.csv', 'Car')
    attrs = ['MGP', 'Cylinders', 'Displacement', 'HorsePower', 'Weight', 'Acceleration', 'Model', 'Origin']
    #data_mat = csv2mat('./test.csv', 'Test')
    data_mat_reg = dataRegulization(data_mat)
    np_mat = np.matrix(data_mat_reg)
    num_of_attrs = len(data_mat[0])
    np_cc = nD2cc(num_of_attrs)
    #np_ccc = np.matrix([[0 for i in range(num_of_attrs)] for j in range(2)])
    #np_ccc = np_cc
    np_ccc = np_cc.copy()
    #figure(np_cc, np_mat, attrs)
    np_cc[0,-1] = 0.3
    np_cc[1,-1] = 0.4
    #print 'np_ccc',np_ccc
    #np_cc[0,0] = 0.5*np_cc[0,0]
    #np_cc[1,0] = 0.5*np_cc[1,0]
    #figure(np_cc, np_mat, attrs)
    np_gen = osc(np_cc, 'em', [7])
    #print costFunc(np_gen)
    figure(np_gen, np_mat, attrs)
    #np_gen = osc(np_cc, 'Leh_re')
    #print costFunc(np_gen)
    #figure(np_gen, np_mat, attrs)
    #np_gen = osc(np_cc, 're', [7])
    #figure(np_gen, np_mat, attrs)
    [osc_rad, osc_sca] = FRND(np_cc,np_ccc, 7)
    figure(osc_rad, np_mat, attrs)
    figure(osc_sca, np_mat, attrs)
    plot.show()
