import numpy as np
import re
import shapefile
import matplotlib.pyplot as plt
import glob


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]

def read_shp(fl):
    sf = shapefile.Reader(fl)
    print 'number of shapes imported:',len(sf.shapes())
    print ' '
    print 'geometry attributes in each shape:'
    for name in dir(sf.shape()):
        if not name.startswith('__'):
           print name
    return sf


def plot_single_shape(sf):
    plt.figure()
    ax = plt.axes()
    ax.set_aspect('equal')
    shape_ex = sf.shape(1)
    x_lon = np.zeros((len(shape_ex.points),1))
    y_lat = np.zeros((len(shape_ex.points),1))
    for ip in range(len(shape_ex.points)):
        x_lon[ip] = shape_ex.points[ip][0]
        y_lat[ip] = shape_ex.points[ip][1]

    plt.plot(x_lon,y_lat,'k') 
    plt.show()


def plot_all_shapes(sf):
    plt.figure()
    ax = plt.axes()
    ax.set_aspect('equal')
    for shape in list(sf.iterShapes()):
        x_lon = np.zeros((len(shape.points),1))
        y_lat = np.zeros((len(shape.points),1))
        for ip in range(len(shape.points)):
            x_lon[ip] = shape.points[ip][0]
            y_lat[ip] = shape.points[ip][1]
        plt.plot(x_lon,y_lat) 
    plt.show()



def plot_all_shapes_and_parts(sf, fl_num):
    colors = ['black', 'red', 'yellow',
              'orange', 'blue', 'cyan',
              'green', 'limegreen', 'lightgrey',
              'grey', 'm', 'olive',
              'maroon', 'lightcoral', 'aqua']

    ax = plt.axes() # add the axes
    ax.set_aspect('equal')

    for shape in list(sf.iterShapes()):
        npoints=len(shape.points) # total points
        nparts = len(shape.parts) # total parts

        if nparts == 1:
            x_lon = np.zeros((len(shape.points),1))
            y_lat = np.zeros((len(shape.points),1))
            for ip in range(len(shape.points)):
                x_lon[ip] = shape.points[ip][0]
                y_lat[ip] = shape.points[ip][1]
            plt.plot(x_lon,y_lat, c=colors[int(fl_num)-1]) 

        else: # loop over parts of each shape, plot separately
            for ip in range(nparts): # loop over parts, plot separately
                i0=shape.parts[ip]
                if ip < nparts-1:
                   i1 = shape.parts[ip+1]-1
                else:
                   i1 = npoints

                seg=shape.points[i0:i1+1]
                x_lon = np.zeros((len(seg),1))
                y_lat = np.zeros((len(seg),1))
                for ip in range(len(seg)):
                    x_lon[ip] = seg[ip][0]
                    y_lat[ip] = seg[ip][1]

                plt.plot(x_lon,y_lat, c=colors[int(fl_num)-1]) 
    #plt.savefig('../vector_plots/event_%s.pdf' % fl_num)

if __name__ == '__main__':
    '''
    fls = glob.glob('../data/FSR_DATA/spread_vectors/SVs_Event_*.shp')
    fls = sorted(fls, key=natural_sort_key)
    plt.figure()
    for fl in fls:
        fl_num = fl.split('_')[-1].strip('.shp')
        print fl_num
        sf = shapefile.Reader(fl)

        plot_all_shapes_and_parts(sf, fl_num)
    plt.savefig('../vector_plots/all_events.pdf')
    '''

    fl = '/Users/mcurrie/GitRepos/CAwind/data/FSR_DATA/fire_fronts_timed/' \
        + 'Event 12/Active_Flame_Front_(14_09).shp'

    sf = shapefile.Reader(fl)
    shape = list(sf.iterShapes())[0]
    x_lon = np.zeros((len(shape.points),1))
    y_lat = np.zeros((len(shape.points),1))
    for ip in range(len(shape.points)):
        x_lon[ip] = shape.points[ip][0]
        y_lat[ip] = shape.points[ip][1]

    print x_lon, y_lat
