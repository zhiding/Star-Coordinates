import csv, math, re
import matplotlib.pyplot as plot
import numpy as np
from matplotlib.patches import Circle
from numpy import linalg as la
import random

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
        cc_cos += str(np.cos(theta*(i))) + ' '
        cc_sin += str(np.sin(theta*(i))) + ' '
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
        #ax.annotate(attrs[i], xy=(np_coord[0,i], np_coord[1,i]))
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
def osc(sc_now, Type='em', fix=[], rad=[], direct=[], iterate=50):
    sc_ori = sc_now.copy()
    num_of_attrs = sc_now.size / 2
    sc_next = np.matrix([[0 for i in range(num_of_attrs)] for j in range(2)], dtype='double') 
    if Type == 're':
        coord_x = sc_now[0,:]
        coord_y = sc_now[1,:]
        coord_x = coord_x / la.norm(coord_x)
        coord_y = coord_y - np.sum(np.multiply(coord_y, coord_x))*coord_x
        sc_next[0,:] = coord_x
        sc_next[1,:] = coord_y / la.norm(coord_y)
    elif Type == 'em':
        ## Backtrack Line Search
        sc_next = sc_now.copy()
        np_delta = np.matrix([[0 for i in range(num_of_attrs)] for j in range(2)], dtype='double')
        [cost, xx, yy, xy] = costFunc(sc_next)
        step = 1
        alpha = 0.25
        beta = 0.8
        error = 0.000001
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

                np_delta[0,i] = 4 * x * (xx - 1) + 2 * xy * y
                np_delta[1,i] = 4 * y * (yy - 1) + 2 * xy * x

                if i in rad:
                    np_delta[0,i] *= 0.3
                    np_delta[1,i] *= 0.3
                    np_delta[0,i] += 0.7 * 4 * rad_delta_sum * x * (1 - sr0 / sr)
                    np_delta[1,i] += 0.7 * 4 * rad_delta_sum * y * (1 - sr0 / sr)
                if i in direct:
                    np_delta[0,i] *= 0.3
                    np_delta[1,i] *= 0.3
                    np_delta[0,i] += 0.7 * 4 * direct_delta_sum*((x*x0+y*y0)/(sr*sr0)-1)/sr0*(x0*y**2-y0*x*y)/sr**3
                    np_delta[1,i] += 0.7 * 4 * direct_delta_sum*((x*x0+y*y0)/(sr*sr0)-1)/sr0*(y0*x**2-x0*x*y)/sr**3
                if i in fix:
                    np_delta[0,i] = 0
                    np_delta[1,i] = 0
            cost_gen = costFunc(sc_next-step*np_delta)[0]
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
    cost = (sum_of_xx - 1) ** 2 + (sum_of_yy - 1) ** 2 + sum_of_xy ** 2
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
    coords = coordMove(coords, [1,3,5], [-0.2, -0.4, 0.8], [0.4, -0.5, -0.8])
    coordsFigure(coords)

    coords_fix = osc(coords, 'em', fix=[1,6])
    fig_fix = plot.figure()
    ax_fix = fig_fix.add_subplot(111)
    add_point(ax_fix, coords_fix, color='#F98614')
    add_point(ax_fix, coords)
    add_p2p(ax_fix, coords, coords_fix)
    add_o2p(ax_fix, coords)
    add_o2p(ax_fix, coords_fix)
    add_circle(ax_fix, coords)
    add_circle(ax_fix, coords_fix)
    ax_fix.axis([-1.2, 1.2, -1.2, 1.2])
    ax_fix.axis('equal')

    coords_rad = osc(coords, 'em', rad=[1,6])
    fig_rad = plot.figure()
    ax_rad = fig_rad.add_subplot(111)
    add_point(ax_rad, coords_rad, color='#F98614')
    add_point(ax_rad, coords)
    add_p2p(ax_rad, coords, coords_rad)
    add_o2p(ax_rad, coords)
    add_o2p(ax_rad, coords_rad)
    add_circle(ax_rad, coords)
    add_circle(ax_rad, coords_rad)
    ax_rad.axis([-1.2, 1.2, -1.2, 1.2])
    ax_rad.axis('equal')

    coords_direct = osc(coords, 'em', direct=[1,6])
    fig_direct = plot.figure()
    ax_direct = fig_direct.add_subplot(111)
    add_point(ax_direct, coords_direct, color='#F98614')
    add_point(ax_direct, coords)
    add_p2p(ax_direct, coords, coords_direct)
    add_o2p(ax_direct, coords)
    add_o2p(ax_direct, coords_direct)
    add_circle(ax_direct, coords)
    add_circle(ax_direct, coords_fix)
    ax_direct.axis([-1.2, 1.2, -1.2, 1.2])
    ax_direct.axis('equal')
    return 0

def add_point(ax, coords, color='#7FFFD4', marker='o', markersize=5):
    dims = coords.size / 2
    for i in range(dims):
        ax.plot(coords[0,i], coords[1,i], marker=marker, markersize=markersize, color=color)
    return 0

def add_circle(ax, coords):
    dims = coords.size / 2
    for i in range(dims):
        r = math.sqrt(coords[0,i]**2+coords[1,i]**2)
        circle = Circle((0.0, 0.0), radius=r, fill=False, color='#B5B5B5')
        ax.add_patch(circle)
    return 0

def add_o2p(ax, co):
    dims = co.size / 2 
    for i in range(dims):
        ax.plot([co[0,i], 0], [co[1,i], 0], linewidth=1.5, color='#FF6A6A')
    return 0

def add_p2p(ax, co, co_):
    dims = co.size / 2 
    for i in range(dims):
        ax.plot([co[0,i], co_[0,i]], [co[1,i], co_[1,i]], linewidth=1.5, color='#FF6A6A', linestyle='--')
    return 0

def interaction_demo():
    num_of_attrs = 5
    delta_x = 0.2
    delta_y = 0.1
    iterate = 5
    coords = nD2cc(num_of_attrs)
    fig = plot.figure()
    ax = fig.add_subplot(111)
    add_point(ax, coords, '#7FFFD4')
    add_o2p(ax, coords)
    add_circle(ax, coords)
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
    add_point(ax_em, coords_ori, '#7FFFD4')
    add_o2p(ax_em, coords_ori)
    add_circle(ax_em, coords_ori)
    for it in range(iterate):
        coords_em = osc(coords_ori, 'em', [0], iterate=iterate)
        add_point(ax_em, coords_em, '#F98614')
        add_o2p(ax_em, coords_em)
        add_p2p(ax_em, coords_em, coords_ori)
        coords_ori = coords_em
    add_circle(ax_em, coords_em)
    ax_em.axis([-1.2, 1.2, -1.2, 1.2])
    ax_em.axis('equal')
    return 0

def morphing(start, end, k=100, method='simple', alg='em'):
    fig = plot.figure()
    ax = fig.add_subplot(111)
    med = start.copy()
    for i in range(k+1):
        if method == 'simple':
            med = start * (1-1.0*i/k) + end * 1.0*i/k
        elif method == 'blending':
            med = start * (1-1.0*i/k) + end * 1.0*i/k
            if alg == 're':
                med = osc(med, 're')
            else:
                med = osc(med, 'em')
        elif method == 'stepwise':
            med = med * (1-1.0*i/k) + end * 1.0*i/k
            if alg == 're':
                med = osc(med, 're')
            else:
                med = osc(med, 'em')
        add_point(ax, med)
        add_o2p(ax, med)
    ax.axis('equal')
    ax.axis('scaled')
    return 0

def morphing_demo(k=100, method='simple', alg='em'):
    start = nD2cc(5)
    end = coordMove(start, [1, 2, 4], [-0.3, -0.4, -0.5], [0.6, -0.3, -0.2])
    start = osc(start, 'em')
    end = osc(end, 'em', fix=[1, 4])

    #fig_ex = plot.figure()
    #ax_ex = fig_ex.add_subplot(111)
    #add_point(ax_ex, start)
    #add_point(ax_ex, end, color='#F98614')
    #add_o2p(ax_ex, start)
    #add_o2p(ax_ex, end)
    #add_circle(ax_ex, start)
    #ax_ex.axis('equal')
    #ax_ex.axis('scaled')

    fig = plot.figure()
    ax = fig.add_subplot(111)
    med = start.copy()
    for i in range(k+1):
        if method == 'simple':
            med = start * (1-1.0*i/k) + end * 1.0*i/k
        elif method == 'blending':
            med = start * (1-1.0*i/k) + end * 1.0*i/k
            if alg == 're':
                med = osc(med, 're')
            else:
                med = osc(med, 'em')
        elif method == 'stepwise':
            med = med * (1-1.0*i/k) + end * 1.0*i/k
            if alg == 're':
                med = osc(med, 're')
            else:
                med = osc(med, 'em')
        add_point(ax, med)
        add_o2p(ax, med)
    ax.axis('equal')
    ax.axis('scaled')
    return 0

def genCirclenum(ox, oy, r_max, n=2000):
    circle = [[],[]]
    for i in range(n):
        sign_x = 2 * int(random.random() > 0.5) - 1
        sign_y = 2 * int(random.random() < 0.5) - 1
        r = random.random()*r_max
        x = random.random()*r*sign_x-ox
        y = math.sqrt(r**2-(x+ox)**2)*sign_y-oy
        circle[0].append(x)
        circle[1].append(y)
    circle = np.matrix(circle)
    return circle 

def application_demo():
    n = 500
    data1 = genCirclenum(-0.6, 0.6, 0.1, n)
    data2 = genCirclenum(0.1, 0.3, 0.1, n)
    data3 = genCirclenum(-0.1, -0.1, 0.1, n)

    start = osc(nD2cc(5), 're')
    end = start.copy()
    end[:,3] = 0.5*end[:,3]+1.2*end[:,4]

    fig_s = plot.figure()
    ax_s = fig_s.add_subplot(111)
    add_point(ax_s, start)
    add_point(ax_s, end[:,3], color='#F98614')
    add_o2p(ax_s, start)
    add_o2p(ax_s, end[:,3])
    add_p2p(ax_s, start[:,3], end[:,3])
    add_point(ax_s, data1, color='#FFCC5C')
    add_point(ax_s, data2, color='#66CCFF')
    add_point(ax_s, data3, color='#68C4AF')
    add_circle(ax_s, start)
    ax_s.axis('equal')
    ax_s.axis('scaled')

    inv = start.T*(start*start.T)**(-1)
    end1 = osc(end, 'em', fix=[3])
    end2 = osc(end, 'em', rad=[3])
    end3 = osc(end, 'em', direct=[3])
    ends = [end1, end2, end3]
    for e in ends:
        A = e*inv
        cos = 0.5 * (A[0,0] + A[1,1])
        A = np.matrix([[cos, -math.sqrt(1-cos**2)],[math.sqrt(1-cos**2), cos]])
        d1 = A*data1 
        d2 = A*data2 
        d3 = A*data3 
        fig = plot.figure()
        ax = fig.add_subplot(111)
        add_point(ax, e)
        add_o2p(ax, e)
        add_point(ax, d1, color='#FFCC5C')
        add_point(ax, d2, color='#66CCFF')
        add_point(ax, d3, color='#68C4AF')
        add_circle(ax, e)
        ax.axis('equal')
        ax.axis('scaled')
    k = 100
    e = end1
    fig_m = plot.figure()
    ax_m = fig_m.add_subplot(111)
    med = start.copy()
    for i in range(k+1):
        print i
        med = med * (1-1.0*i/k) + e * 1.0*i/k
        med = osc(med, 'em')
        mA = med*inv
        cos = 0.5*(A[0,0]+A[1,1])
        A = np.matrix([[cos, -math.sqrt(1-cos**2)],[math.sqrt(1-cos**2), cos]])
        d1 = mA*data1 
        d2 = mA*data2 
        d3 = mA*data3 
        add_point(ax_m , med)
        add_o2p(ax_m, med)
        add_point(ax_m, d1, color='#FFCC5C', marker='.')
        add_point(ax_m, d2, color='#66CCFF', marker='.')
        add_point(ax_m, d3, color='#68C4AF', marker='.')
    ax_m.axis('equal')
    ax_m.axis('scaled')
    return 0

if __name__ == '__main__':
    
    ## Demo of RE and EM Algorithms 
    #REandEM_demo = REandEM_demo() 
    
    ## Demo of Interaction on RE and EM Algorithms
    #interaction_demo()

    ## Demo of Conditional Interaction on EM Algorithms
    #conditional_demo()
    
    ## Demo of Morphing
    #morphing_demo(method='simple')
    #morphing_demo(method='blending', alg='re') 
    #morphing_demo(method='stepwise', alg='re')

    ## Demo of Applications 
    #application_demo()
    
    ## Demo of DataSet
    data_mat = csv2mat('./cars.csv', 'Car')
    #attrs = ['MGP', 'Cylinders', 'Displacement', 'HorsePower', 'Weight', 'Acceleration', 'Model', 'Origin']
    attrs = []
    #data_mat = csv2mat('./test.csv', 'Test')
    data_mat_reg = dataRegulization(data_mat)
    np_mat = np.matrix(data_mat_reg)
    num_of_attrs = len(data_mat[0])
    ## nD -> 2D
    np_ori = nD2cc(num_of_attrs)
    figure(np_ori, np_mat, attrs)
    np_gen = np_ori.copy()
    np_gen = osc(np_ori, 'em', fix=[7])
    figure(np_gen, np_mat, attrs)
    np_gen[:,7] = np_gen[:,7]+[[0.4],[0.2]]
    np_gen = osc(np_gen, 'em', fix=[7])
    figure(np_gen, np_mat, attrs)
    np_gen = osc(np_ori, 'em', fix=[2,4])
    figure(np_gen, np_mat, attrs)
    np_gen[:,7] = np_gen[:,7]+[[0.4],[0.2]]
    np_gen = osc(np_gen, 'em', fix=[7,2])
    figure(np_gen, np_mat, attrs)
    #np_gen = osc(np_cc, 're', [7])
    #figure(np_gen, np_mat, attrs)
    ##[osc_rad, osc_sca] = FRND(np_cc,np_ori, 7)
    ##figure(osc_rad, np_mat, attrs)
    plot.show()
