import csv, re
import numpy as np
import matplotlib.pyplot as plot

########################################################################################
## Attributes: MPG;Cylinders;Displacement;Horsepower;Weight;Acceleration;Model;Origin ##
########################################################################################
def csv2dict(filepath, primary_key):
    data_dict = {}
    with open(filepath, 'rb') as csvfile:
        items = csv.DictReader(csvfile, delimiter=';')
        item_id = 0
        for item in items:
            data_dict[str(item_id)] = item
            item_id += 1
    return data_dict

################################################################
## Input  ==> (Dict) data_dict                                ##
## Output ==> (List) data_mat[size_of_dataSet * num_of_attrs] ## 
################################################################
def dict2mat(data_dict):
    mat = []
    for key in data_dict:
        vector = []
        vector.append(key)
        attrs = data_dict[key]
        for attr in attrs:
            value = str2float(attrs[attr])
            vector.append(value)
        mat.append(vector)
    return mat

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

##############################################
## Descript: 1. Remove colmun with data_id  ##
##           2. Regulization in each column ##
##############################################
def dataRegulization(data_mat):
    data_mat_wfc = [row[1:] for row in data_mat]
    num_of_row = len(data_mat_wfc)
    num_of_col = len(data_mat_wfc[0])
    data_mat_reg = [[0 for col in range(num_of_col)] for row in range(num_of_row)] 
    col_max = []
    col_min = []
    for i in range(num_of_col):
        col_max.append(max([row[i] for row in data_mat_wfc]))
        col_min.append(min([row[i] for row in data_mat_wfc]))    
        col_range = col_max[i] - col_min[i]
        for j in range(num_of_row):
            data_mat_reg[j][i] = (data_mat_wfc[j][i]-col_min[i])/(col_max[i]-col_min[i])
    return data_mat_reg 

#################################################################
## Usage: Convert n-Dimensional space to Cartesian Coordinates ##
#################################################################
def nD2cc(num_of_attrs):
    theta = 2*np.pi / num_of_attrs
    cc_cos = ''
    cc_sin = ''
    for i in range(num_of_attrs):
        cc_cos += str(np.cos(theta*i)) + ' '
        cc_sin += str(np.sin(theta*i)) + ' '
    np_cc = np.matrix(cc_cos[:-1]+';'+cc_sin[:-1])
    return np_cc

if __name__ == '__main__':
    data_dict = csv2dict('./cars.csv', 'Car')
    data_mat = dict2mat(data_dict)
    data_mat_reg = dataRegulization(data_mat)
    np_mat = np.matrix(data_mat_reg)

    num_of_attrs = len(data_mat[0])-1
    np_cc = nD2cc(num_of_attrs)
    sc_mat = np.dot(np_cc, np_mat.T)
    sc_mat_len = sc_mat.size/sc_mat.ndim
    sc_array = np.array(sc_mat).reshape(-1,).tolist()
    sc_dx = sc_array[:sc_mat_len]
    sc_dy = sc_array[sc_mat_len:]
    cc_array = np.array(np_cc).reshape(-1,).tolist()
    cc_x = cc_array[:num_of_attrs]
    cc_y = cc_array[num_of_attrs:]
    figure = plot.figure(1)
    plot.scatter(sc_dy ,sc_dx)
    for i in range(num_of_attrs):
        plot.plot([cc_y[i],0],[cc_x[i],0])
    plot.show()

