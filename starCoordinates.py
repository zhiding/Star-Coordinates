import csv, re
import numpy as np

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

def coordinateInit(num_of_attrs):
    theta = 2*np.pi / num_of_attrs
    mat_str = ''
    for i in range(num_of_attrs):
        mat_str += str(np.cos(theta*i)) + ' ' + str(np.sin(theta*i)) + ';'
    mat = np.matrix(mat_str[:-1])
    return mat

#########
## Input: (List) dataSet[size_of_dataSet * (num_of_attrs-1)]
#########
def matStats(data_mat, opt):
    mat = np.matrix(data_mat[1:][:])
    return mat
    
def nD2sc(unit_vecs, data):
    proj = np.zeros((1,2))
    d_max = np.amax(data)
    d_min = np.amin(data)
    return 0

if __name__ == '__main__':
    data_dict = csv2dict('./cars.csv', 'Car')
    data_mat = dict2mat(data_dict)
    mat = matStats(data_mat, '')
    print mat[0:2]
