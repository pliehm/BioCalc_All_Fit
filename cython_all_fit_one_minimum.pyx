#cython: boundscheck=False
#cython: profile=False
#cython: cdivision = True

import numpy as np
cimport numpy as np
cimport cython
import matplotlib.pyplot as plt
import time as time

# open file to write
#result = open("result.txt",'w')

# define types for cython, not sure if this is really necessary, but I found it in one of the examples I used to improve the speed
DTYPE1 = np.uint16
ctypedef np.uint16_t DTYPE_t
DTYPE2 = np.float
ctypedef np.float DTYPE_t2
DTYPE3 = np.uint32
ctypedef np.uint32_t DTYPE_t3

# fit simulated and experimental minima to calculate thickness

# define a faster function to calculate the absolute value of a number
cdef float _abs(float a): return a if a>=0 else -a

#plt.ion()

###################################################
# Function to fit thickness with thickness limits #
###################################################

# almost the same function as above, the only difference is that only thicknesses within the thickness_limit are considered for the fitting, this improves speed and avoids jumps between cavities
@cython.boundscheck(False)
cdef Fit_limit(np.ndarray[DTYPE_t, ndim=1] thickness,np.ndarray[DTYPE_t3, ndim=1] array_thickness_len_pos,np.ndarray[DTYPE_t, ndim=1] array_length_block, np.ndarray[double,ndim=1] exp_waves, float tolerance,np.ndarray[double,ndim=1] sim_minima_blocks, current_index, thickness_list, thickness_limit):

# cython type definition of the used objects
 
    cdef list sim_min_waves = [[],[]] # dummy values to have two lists
    cdef unsigned short i, k, min_thickness_i, max_thickness_i
    cdef unsigned short L_exp_waves = len(exp_waves) # too not use the len() function too often
    cdef float summe, summe_temp=0
    cdef unsigned int counter = 0
    cdef unsigned int position, len_block 
    cdef unsigned int breaker = 0
    cdef list diff_list = []
    cdef float exp_wave

    if L_exp_waves > 0: # first test if there is a minimum
       # for waves in s_waves_arrays:
        min_thickness_i = current_index-thickness_limit #
        max_thickness_i = current_index+thickness_limit #
        for i in range(min_thickness_i,max_thickness_i): # do the following calculations for every simulated thickness
            summe = 0
            summe_temp = 0
            position = array_thickness_len_pos[i]
            len_block = array_length_block[i]

            exp_wave = exp_waves[0]
            if len_block >0:
                for k in xrange(len_block): # can I just use the last minimum here? Do I actually have to????, hm, less points fitted....
                    #diff_list.append(_abs(sim_minima_blocks[position+k]-exp_wave)) # i thin this is slow
                    summe_temp = _abs(sim_minima_blocks[position+k]-exp_wave)
                    if k == 0:
                        summe = summe_temp
                    if summe_temp<summe and summe != 0:
                        summe = summe_temp

                
                #summe = _abs(sim_minima_blocks[position+len_block-1]-exp_waves[0])
            else:
                summe = 1000
           #print summe
            if summe<= tolerance:
                sim_min_waves[0].append(thickness[i])
                sim_min_waves[1].append(summe)



# return the thickness with minimum value
            #result.write(str(sim_min_waves)+"\n")
        if  len(sim_min_waves[0])>1 and (min(sim_min_waves[1]) < tolerance):
            #print sim_min_waves[0][sim_min_waves[1].index(min(sim_min_waves[1]))], thickness_list.index(sim_min_waves[0][sim_min_waves[1].index(min(sim_min_waves[1]))])
            # print i
            # print sim_min_waves[0]

            return sim_min_waves[0][sim_min_waves[1].index(min(sim_min_waves[1]))]#, thickness_list.index(sim_min_waves[0][sim_min_waves[1].index(min(sim_min_waves[1]))]), sim_min_waves
        else: 
            return 0

    else:
        return 0



#########################################################
# function to get minima of the intensity profile array #
#########################################################

cdef list peakdetect(np.ndarray[DTYPE_t, ndim=1] y_axis, list x_axis = None, unsigned short lookahead_min=5, unsigned short lookahead_max=3, unsigned short delta = 0,unsigned short enhanced_resolution=10,unsigned short average_window=5):
    
    # define output container
    #cdef list max_peaks=[]
    cdef list min_peaks = []
    cdef  list min_peaks_new = []
    cdef list dump = [] # used to pop the first hit which almost always is false
    cdef list y_axis_list, x_axis_list = []  # convert array to list, min() is faster for list
    # check input data --> this makes the algorithm 5 times slower
    #x_axis, y_axis = _datacheck_peakdetect(x_axis, y_axis) 
    
    # store data length for later use
    cdef unsigned int length 

    #perform some checks
    #if lookahead < 1:
    #    raise ValueError, "Lookahead must be '1' or above in value"
    #if not (np.isscalar(delta) and delta >= 0):
    #    raise ValueError, "delta must be a positive number"

    #maxima and minima candidates are temporarily stored in 
    #mx and mn respectively
    cdef float mn = 100000
    cdef float mx = -100000
    cdef unsigned short index, window_len 
    cdef float x, y 
    cdef np.ndarray[np.double_t, ndim=1] y_interp, x_interp,s, y_temp,w

#


    if enhanced_resolution != 1:

        ################################################
        ### Interpolate data from 1nm to 0.1nm steps ###
        ################################################

        x_interp = np.linspace(x_axis[0],x_axis[-1],num=int((len(x_axis)-1)*enhanced_resolution+1))

        y_interp = np.interp(x_interp,x_axis,y_axis)

        #print x_interp
        #print y_interp

        #y_axis = y_interp
        #x_axis = x_interp



        #########################################
        ### Smooth (by averaging) the profile ###
        #########################################

        #plt.plot(x_axis,y_axis)

        window_len = average_window*enhanced_resolution # this needs to be checked, maybe it can be smaller or has to be larger

        s=np.r_[abs(2*y_interp[0]-y_interp[window_len:1:-1]), y_interp, abs(2*y_interp[-1]-y_interp[-1:-window_len:-1])]

        #print s
        w = np.ones(window_len,'d')/window_len
        y_temp = np.convolve(w, s, mode='same')

        y_axis_list = y_temp[window_len-1:-window_len+1].tolist()

        x_axis= x_interp.tolist()


    else:
        y_axis_list = y_axis.tolist()
        #x_axis = x_axis_t

    #y_axis_list = y_axis#.tolist()
    length = len(y_axis_list)

    #plt.figure(1)
    #plt.clf()
    #plt.plot(x_axis,y_axis_list)
    #plt.ylim((1000,9000))
    #print len(x_axis), len(y_axis)
    #print x_axis
    #print y_axis
    for index, (x,y) in enumerate(zip(x_axis[:-lookahead_min], y_axis_list[:-lookahead_min])):
        #print index
        
        if y > mx:
            mx = y
            mxpos = x

        if y < mn:
            mn = y
            mnpos = x

        #### look for max ####
        
        if y < mx-delta and mx != 100000:
            #Maxima peak candidate found
            # lool ahead in signal to ensure that this is a peak and not jitter
            if max(y_axis_list[index:index+lookahead_max]) < mx:
                #max_peaks.append([mxpos, mx])
                dump.append(True)
                #set algorithm to only find minima now
                
                mx = 100000
                mn = 100000
                if index+lookahead_min >= length:
                    #end is within lookahead no more peaks can be found
                    break
                continue

        #### look for min ####    
        
        if y > mn+delta and mn != -100000:
            #Minima peak candidate found
            # look ahead in signal to ensure that this is a peak and not jitter
            if min(y_axis_list[index:index+lookahead_min]) > mn:
                min_peaks.append(mnpos)
                dump.append(False)
                #set algorithm to only find maximum now
                mn = -100000
                mx = -100000
                if index+lookahead_min >= length:
                    #end is within lookahead no more peaks can be found
                    break


    #Remove the false hit on the first value of the y_axis
    if len(dump)>0:
        if not dump[0]:
            min_peaks.pop(0)
        else:
            pass    
    
        #no peaks were found, should the function return empty lists?

    # take only the last minimum
    #plt.vlines(min_peaks,min(y_axis)-100,max(y_axis))
    if len(min_peaks)>0:
        min_peaks = [min_peaks[-1]]

    #result.write(str(min_peaks)+"\n")


    #plt.vlines(min_peaks,min(y_axis)-100,max(y_axis), color='r')
    #plt.pause(0.00001)

    #print min_peaks
    return min_peaks # maybe this can actually be a number and not a list?

    

####################################################################################
# Main function which calls the minima finding algorithm and the thickness fitting #
####################################################################################


# the following parameters are passed to the function:
# start wavelength, end wavelength, all images, list of thickness/blocklength/position,list of waves, tolerance, lookahead_min, lookahead_max, delta, delta variations, minima blocks, thickness limits in use?, thickness limit, 
def c_Fit_Pixel(unsigned int start,unsigned int ende, np.ndarray[DTYPE_t, ndim=3] data, list thickness_len_pos, list waves, float tolerance, unsigned short lookahead_min,unsigned short lookahead_max, unsigned short delta, unsigned short delta_vary, list list_minima_blocks, use_thickness_limits, unsigned int thickness_limit, unsigned short area_avrg, unsigned short init_guess, unsigned short enhanced_resolution, unsigned short average_window):

    ########################################
    # definition of all the variable types #
    ########################################

    cdef unsigned int Image_width = len(data[0][0])
    cdef np.ndarray[DTYPE_t, ndim=2] thickness_ready = np.zeros((ende-start,Image_width),np.uint16 )
    cdef unsigned short column, row, column_c, row_c
    cdef np.ndarray[DTYPE_t, ndim=1] intensity
    cdef np.ndarray[double,ndim=1] minima_exp
    cdef unsigned int counter=start, 
    cdef np.ndarray[double,ndim=1] sim_minima_blocks
    cdef np.ndarray[DTYPE_t,ndim=1] array_length_block, thickness
    cdef np.ndarray[DTYPE_t3,ndim=1] array_thickness_len_pos
    cdef unsigned int current_thickness = 0, last_index = 0, limit_counter = 0 # choose b = 10 to use interp + smoothing of data
    cdef float last_thickness = 0
    cdef list a = [] # dummy list
    cdef list thickness_list = []
    cdef list init_guess_pos = []
    # TEST
    cdef np.ndarray[double, ndim=2] minima_ready = np.zeros((ende-start,Image_width),np.float)

    #cdef unisgned short init_guess_thickness = init_guess[2]
    #cdef np.ndarray[double, ndim=1] x_interp = np.linspace(waves[0],waves[-1],num=int((len(waves)-1)/0.1+1))
    #cdef np.ndarray[double, ndim=1] y_interp = np.interp(x_interp,waves,data[:,0,0])

    # make a list which contains all the thicknesses (maybe one could also use the array below?)
    for i in range(len(thickness_len_pos)):
        thickness_list.append(int(thickness_len_pos[i][0]))

    # make an array which contains all minima, but not seperated as before
    for block in list_minima_blocks:
        for i in range(len(block)):
            a.append(block[i])
    sim_minima_blocks = np.array(a,dtype=np.float)

    # make arrays of positions, lengths, thickness
    array_thickness_len_pos = np.array(zip(*thickness_len_pos)[2],dtype=np.uint32)
    array_length_block = np.array(zip(*thickness_len_pos)[1],dtype=np.uint16)
    thickness = np.array(zip(*thickness_len_pos)[0],dtype=np.uint16)

    print 'x ', len(data) 
    print 'y ', len(data[0])
    print 'z ', len(data[0][0])

    ######################################### 
    # do calculations with thickness limits #
    #########################################
    #plt.figure(2)
    if use_thickness_limits:
        #print "using thickness limit: ", thickness_limit

        for row in range(len(data[0])):
            print counter
            # plt.figure(2)
            # plt.clf()
            # plt.imshow(thickness_ready, aspect='auto')
            # plt.clim(8300,8850)
            # plt.colorbar()
            # plt.pause(0.0001)
            counter+=1
            for column in xrange(Image_width):
                #print column
                last_thickness = 0
                if column == 0 and row ==0:
                     last_thickness =init_guess

                last_index = 0
                current_thickness = 0
                limit_counter = 0

                # write loop to consider the area around the current pixel

                
                if column !=0:
                    if thickness_ready[row][column-1] != 0:
                        last_thickness = thickness_ready[row][column-1]


                if last_thickness == 0:
                    if column >1 and row > 0:
                        if thickness_ready[row-1][column-1] != 0:
                            last_thickness = thickness_ready[row-1][column-1]
                if last_thickness == 0:

                    # iterate over rows, distance from current row is 0 to area_avrg
                    for row_c in range(area_avrg+1):
                        # if the considered row is not in the range, stop iteration over rows
                        if row  < row_c:
                            #print "stopped 1"
                            break

                        # block to consider values with subtracted column values
                        
                        # iterate over column between 1 to area_avrg to not sonsider the same column as current 
                        for column_c in range(1,area_avrg+1):
                            # if the considered column is below the data range, stop loop
                            if column < column_c:
                                #print "stopped 2"
                                break

                            # check if the considered value is not zero, if so add the value to a sum and increase the count for later averaging
                            if thickness_ready[row-row_c][column-column_c] != 0:
                                last_thickness+=  thickness_ready[row-row_c][column-column_c]
                                limit_counter += 1

                        # block to consider values with added column values and same column
                        for column_c in range(area_avrg+1):
                            # if the considered row is the current row, stop, because there are no values
                            if (row_c == 0):
                                #print "stopped 3"
                                break
                            # if the considered column is above the data range, stop loop
                            if (column + column_c) >= Image_width:
                                #print "stopped 4"
                                break
                            # check if the considered value is not zero, if so add the value to a sum and increase the count for later averaging
                            if thickness_ready[row-row_c][column+column_c] != 0:
                                last_thickness+=  thickness_ready[row-row_c][column+column_c]
                                limit_counter += 1

                            # # check if other values in the area were found
                if limit_counter != 0:
                    #print "limit counter is:", limit_counter
                    # calculate average of the area
                    last_thickness = last_thickness/float(limit_counter)

                if last_thickness == 0:
                    last_thickness = init_guess





                # # no previous thickness found
                # else:

                #     last_thickness = init_guess
                #print column, last_thickness
                # if the thickness in the area is in the thickness list, search for the index of that thickness and store it
                if last_thickness > (thickness_list[0] + 2*thickness_limit):
                    last_index = last_thickness - thickness_list[0]# thickness_list.index(int(last_thickness))

                # print current_thickness, last_thickness, last_index
                # last_index = thickness_list.index(int(last_thickness))

                # get array with the intensity profile for the current pixel
                intensity = data[:,row, column]
                
                # find the minima in the profile

                minima_exp = np.array(peakdetect(intensity, waves, lookahead_min*enhanced_resolution,lookahead_max*enhanced_resolution, delta,enhanced_resolution,average_window),dtype=np.float)

                if len(minima_exp)>0:
                    minima_ready[row][column] = minima_exp[0]
                else:
                    minima_ready[row][column] = np.nan
                # start calculations with limits
                #if (last_thickness != 0) and (last_index > 0):

                # only do the calculation if a minimum was found
                if len(minima_exp) != 0: 

                    current_thickness = (Fit_limit(thickness,array_thickness_len_pos, array_length_block, minima_exp,tolerance,sim_minima_blocks,last_index,thickness_list, thickness_limit))
                        #if current_thickness != 0:
                    thickness_ready[row][column] = current_thickness
                else:
                    thickness_ready[row][column] = 0
                    
                #plt.imshow(thickness_ready)

                # if current_thickness != 0:
                #    last_thickness, last_index = current_thickness, current_index
                # print current_thickness
                # if current_thickness != 0:
                #     init_guess = current_thickness

                #plt.figure(2)
                #plt.clf()
                #plt.imshow(thickness_ready, aspect='auto')
                #plt.clim(7500,8000)
                #plt.colorbar()
                #plt.pause(0.0001)
    # feed the fitted thicknesses back to the main program
    #result.close()
    #plt.ioff()
    np.savetxt("minima_map.txt",minima_ready)
    return thickness_ready#, sim_min_waves