import csv, math, re
import numpy as np
import matplotlib.pyplot as plot
from matplotlib.patches import Circle

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
def figure(np_coord, np_data):
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
        ax.plot(np_coord[0,i], np_coord[1,i], 'o', color='#7FFFD4')
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
def osc(np_ori, method='rc'):
    num_of_attrs = np_ori.size / 2
    if method == 'em':
        np_delta = np.matrix([[0 for i in range(num_of_attrs)] for j in range(2)])
        [cost_ori, xx_ori, yy_ori, xy_ori] = cost(np_ori)
        step = 0.01
        error = 0.0001
        for i in range(num_of_attrs):
            np_delta[0,i] = 4 * np_ori[0,i] * (xx_ori - 1) + 2 * xy_ori * np_ori[1,i]
            np_delta[1,i] = 4 * np_ori[1,i] * (yy_ori - 1) + 2 * xy_ori * np_ori[0,i]
        np_gen = np_ori - step * np_delta
        [cost_gen, xx_gen, yy_gen, xy_gen] = cost(np_gen)
        print [cost_ori, xx_ori, yy_ori, xy_ori]
        print [cost_gen, xx_gen, yy_gen, xy_gen]
        print np_ori
        print np_gen
    else:
        pass
    return np_gen

def cost(np_sc):
    num_of_attrs = np_sc.size / 2 
    sum_of_xx = 0
    sum_of_yy = 0
    sum_of_xy = 0
    for i in range(num_of_attrs):
        sum_of_xx += np_sc[0,i] * np_sc[0,i]
        sum_of_yy += np_sc[1,i] * np_sc[1,i]
        sum_of_xy += np_sc[0,i] * np_sc[1,i]
    cost = (sum_of_xx - 1) * (sum_of_xx - 1) + (sum_of_yy - 1) * (sum_of_yy - 1) + sum_of_xy * sum_of_xy
    return [cost, sum_of_xx, sum_of_yy, sum_of_xy]


if __name__ == '__main__':
    data_mat = csv2mat('./cars.csv', 'Car')
    #data_mat = csv2mat('./test.csv', 'Test')
    data_mat_reg = dataRegulization(data_mat)
    np_mat = np.matrix(data_mat_reg)
    num_of_attrs = len(data_mat[0])
    np_cc = nD2cc(num_of_attrs)
    figure(np_cc, np_mat)
    np_cc[:,num_of_attrs-1] *= 5
    figure(np_cc, np_mat)
    np_gen = osc(np_cc, 'em')
    figure(np_gen, np_mat)
    plot.show()
