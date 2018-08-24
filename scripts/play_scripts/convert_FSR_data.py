import numpy as np
import matplotlib.pyplot as plt
import shapefile
import glob

if __name__ == '__main__':

    data_dir = '/Users/mcurrie/GitRepos/CAwind/data/FSR_DATA/'
    '''
    front_events = glob.glob(data_dir + 'fire_fronts_timed/Event*')
    for event_dir in front_events:
        event = event_dir.split(' ')[-1]
        print 'Working on event', event
        active_flame_fronts = glob.glob(event_dir + '/Active_Flame_Front*shp')
        for front in active_flame_fronts:
            date = front.split(')')[0].split('(')[-1]
            print 'date', date
            sf = shapefile.Reader(front)
            shape = list(sf.iterShapes())[0]

            with open('../converted_data/event%s_front%s.txt'
                      % (event, date), 'wb') as fl:
                fl.write('#x_lon    y_lat\n')
                for ip in range(len(shape.points)):
                    x_lon = shape.points[ip][0]
                    y_lat = shape.points[ip][1]
                    print x_lon, y_lat
                    fl.write('%.04f    %.04f\n' % (x_lon, y_lat))
    '''
    spread_vectors = glob.glob(data_dir + 'spread_vectors/*Event*shp')
    plt.figure()
    for vector_fl in spread_vectors:
        event = vector_fl.split('.')[0].split('_')[-1]

        print 'Working on event', event
        sf = shapefile.Reader(vector_fl)
        shapes = list(sf.iterShapes())
        for shape in shapes:
            for ip in range(len(shape.points)):
                x_lon = shape.points[ip][0]
                y_lat = shape.points[ip][1]
                print x_lon, y_lat

            x_lon_start = shape.points[0][0]
            x_lon_stop = shape.points[1][0]
            y_lat_start = shape.points[0][1]
            y_lat_stop = shape.points[1][1]
            dx = x_lon_stop - x_lon_start
            dy = y_lat_stop - y_lat_start


            plt.arrow(x_lon_start, y_lat_start, dx, dy,
                      head_width=20, head_length=20)
            plt.scatter([x_lon_start, x_lon_stop], [y_lat_start, y_lat_stop],
                        s=[100,20])
            with open('../converted_data/SV_event%s.txt' % event, 'wb') as fl:
                fl.write('#x_lon    #y_lat\n')
                fl.write('%.05f    %.05f\n' % (x_lon_start, y_lat_start))
                fl.write('%.05f    %.05f\n' % (x_lon_stop, y_lat_stop))

    plt.show()
