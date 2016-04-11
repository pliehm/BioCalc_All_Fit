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
cdef Fit_limit(np.ndarray[DTYPE_t, ndim=1] thickness,np.ndarray[DTYPE_t3, ndim=1] array_thickness_len_pos,np.ndarray[DTYPE_t, ndim=1] array_length_block, float exp_wave, float tolerance,np.ndarray[double,ndim=1] sim_minima_blocks, unsigned int current_index, thickness_list, unsigned int thickness_limit):

# cython type definition of the used objects
 
    cdef list sim_min_waves = [[],[]] # dummy values to have two lists
    cdef unsigned short i, k, l, min_thickness_i, max_thickness_i
    cdef float summe, summe_temp=0, summe_min = 1000
    cdef unsigned int counter = 0
    cdef unsigned int position, len_block 
    cdef unsigned int thickness_temp, breaker = 0
    cdef double value = 0



   
   # for waves in s_waves_arrays:
    min_thickness_i = current_index-thickness_limit #
    max_thickness_i = current_index+thickness_limit #
    for i in range(min_thickness_i,max_thickness_i): # do the following calculations for every simulated thickness
        summe = 0
        summe_temp = 0
        position = array_thickness_len_pos[i]
        len_block = array_length_block[i]


        if len_block >0:

            for k in xrange(len_block): # can I just use the last minimum here? Do I actually have to????, hm, less points fitted....
                value = sim_minima_blocks[position+k]
                summe_temp = _abs(value-exp_wave) #_abs(sim_minima_blocks[position+k]-exp_wave)

                if k == 0:
                    summe = summe_temp
                if summe_temp<summe and summe != 0:
                    summe = summe_temp


            #summe = _abs(sim_minima_blocks[position+len_block-1]-exp_waves[0])
        else:
            summe = 1000
       #print summe
        if summe<= tolerance:
            if summe < summe_min:
                summe_min = summe
                thickness_temp = thickness[i]
            # sim_min_waves[0].append(thickness[i])
            # sim_min_waves[1].append(summe)



# return the thickness with minimum value
        #result.write(str(sim_min_waves)+"\n")


    # if  len(sim_min_waves[0])>1:

    #     summe_temp = min(sim_min_waves[1])
    #     if summe_temp < tolerance:
    #     #print sim_min_waves[0][sim_min_waves[1].index(min(sim_min_waves[1]))], thickness_list.index(sim_min_waves[0][sim_min_waves[1].index(min(sim_min_waves[1]))])
    #     # print i
    #     # print sim_min_waves[0]


    #         return sim_min_waves[0][sim_min_waves[1].index(summe_temp)]
    #     else:
    #         return 0

    if summe_min != 1000:
        return thickness_temp
    else: 
        return 0





#########################################################
# function to get minima of the intensity profile array #
#########################################################

cdef tuple peakdetect(np.ndarray[DTYPE_t, ndim=1] y_axis, list x_axis = None, unsigned short lookahead_min=5, unsigned short lookahead_max=3, unsigned short delta = 0,unsigned short average_window=5, double minima_guess=0, unsigned int index_guess = 0):
    
    # define output container
    #cdef list max_peaks=[]
    cdef list min_peaks = []
    cdef  list min_peaks_new = []
    cdef list dump = [] # used to pop the first hit which almost always is false
    cdef list y_axis_list, x_axis_list # convert array to list, min() is faster for list
    # check input data --> this makes the algorithm 5 times slower
    #x_axis, y_axis = _datacheck_peakdetect(x_axis, y_axis) 
    
    # store data length for later use
    cdef unsigned int length, x_min, x_i, min_index

    cdef np.ndarray[DTYPE_t, ndim=1] y_temp 
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
    cdef float x, y, y_min, min_peak,wave_max, mnpos,mxpos


    fail = 0

    length = np.size(y_axis)

    if minima_guess != 0:

        if (lookahead_min-1)< index_guess < (length-lookahead_min): # length could be a variable as it should always be the same
            y_temp = y_axis[(index_guess-lookahead_min):(index_guess+lookahead_min+1)] 
            #y_axis_list = y_temp.tolist()
            y_min = 1000000
            x_min = 0 
            for x_i in xrange(len(y_temp)):
                if y_temp[x_i] < y_min:
                    y_min = y_temp[x_i]
                    x_min = x_i 
            #print x_min
            # check that the minimum found is not the last or first value in the array
            if ((lookahead_min-1) < (x_min-lookahead_min+index_guess) < (length-lookahead_min)) and (x_min != 0):
                # save found minimum in list
                min_peaks.append([x_min-lookahead_min+index_guess,x_axis[x_min-lookahead_min+index_guess]])
                
                # check if minimum is in first half of the range
                if (x_min-lookahead_min+index_guess) < (length/2.0):
                    
                    index_guess = length-lookahead_min
                    y_temp = y_axis[(index_guess-lookahead_min):(index_guess+lookahead_min+1)] # array faster?
                    #y_axis_list = y_temp.tolist()
                    y_min = 1000000
                    x_min = 0 
                    for x_i in xrange(len(y_temp)):
                        if y_temp[x_i] < y_min:
                            y_min = y_temp[x_i]
                            x_min = x_i
                    if ((lookahead_min-1) < (x_min-lookahead_min+index_guess) < (length-lookahead_min)) and (x_min != 0):
                        min_peaks.append([x_min-lookahead_min+index_guess,x_axis[x_min-lookahead_min+index_guess]])  
                    min_peaks.sort()        
            else:
                fail = 1
        else: 
            #print "fail"
            fail = 1
    else:
        fail = 1


    if fail == 1:
        y_axis_list = y_axis.tolist()
        length = len(y_axis_list)
        #print "failed"

        # plt.figure(1)
        # plt.clf()
        # plt.plot(x_axis,y_axis_list)
        # plt.ylim((1000,9000))
        #print len(x_axis), len(y_axis)
        #print x_axis
        #print y_axis
        #for index, (x,y) in enumerate(zip(x_axis[:-lookahead_min], y_axis_list[:-lookahead_min])):
        #for index in xrange(length-lookahead_min):

        wave_max = x_axis[-1]-5
        for index in xrange(length):
            y = y_axis_list[index]
            x = x_axis[index]
            # plt.vlines(x_axis[index],5000,9000,color='g') 
            #print index
            
            if y > mx:
                mx = y
                mxpos = x

            if y < mn:
                mn = y
                mnpos = x

            if y < mx-delta and mx != 100000:
                #Maxima peak candidate found
                #max_peaks.append([mxpos, mx])
                dump.append(True)
                #set algorithm to only find minima now
                
                mx = 100000
                mn = 100000
                continue    


            #### look for min ####    
            
            if y > mn+delta and mn != -100000:
                #Minima peak candidate found
                # look ahead in signal to ensure that this is a peak and not jitter
                if mnpos < wave_max:
                    min_peaks.append([index,mnpos])
                    dump.append(False)
                    #set algorithm to only find maximum now
                    mn = -100000
                    mx = -100000

        #Remove the false hit on the first value of the y_axis
        if len(dump)>0:
            if not dump[0]:
                min_peaks.pop(0)
            else:
                pass    
        
            #no peaks were found, should the function return empty lists?

    # take only the last minimum
    # b = np.asarray(min_peaks)
    # plt.vlines(b[:,1],5000,9000)
    # plt.vlines(685,5000,9000)
    if len(min_peaks)>0:
        min_index = min_peaks[-1][0]
        min_peak = min_peaks[-1][1]
    else:
        min_peak = 0
        min_index = 0
        
    #result.write(str(min_peaks)+"\n")


    # plt.vlines(min_peak,5000,9000, color='r')
    # plt.pause(0.00001)

    #print min_peaks
    return min_peak, min_index 

    

####################################################################################
# Main function which calls the minima finding algorithm and the thickness fitting #
####################################################################################


# the following parameters are passed to the function:
# start wavelength, end wavelength, all images, list of thickness/blocklength/position,list of waves, tolerance, lookahead_min, lookahead_max, delta, delta variations, minima blocks, thickness limits in use?, thickness limit, 
def c_Fit_Pixel(unsigned int start,unsigned int ende, np.ndarray[DTYPE_t, ndim=3] data, list thickness_len_pos, list waves, float tolerance, unsigned short lookahead_min,unsigned short lookahead_max, unsigned short delta, unsigned short delta_vary, list list_minima_blocks, use_thickness_limits, unsigned int thickness_limit, unsigned short area_avrg, unsigned short init_guess, unsigned short average_window):

    ########################################
    # definition of all the variable types #
    ########################################

    cdef unsigned int Image_width = len(data[0][0])
    cdef np.ndarray[DTYPE_t, ndim=2] thickness_ready = np.zeros((ende-start,Image_width),np.uint16 )
    cdef unsigned short column, row, column_c, row_c
    cdef float minima_exp
    cdef unsigned int counter=start, 
    cdef np.ndarray[double,ndim=1] sim_minima_blocks
    cdef np.ndarray[DTYPE_t,ndim=1] array_length_block, thickness
    cdef np.ndarray[DTYPE_t3,ndim=1] array_thickness_len_pos
    cdef unsigned int current_thickness = 0, current_index = 0, last_index = 0, limit_counter = 0# choose b = 10 to use interp + smoothing of data
    cdef float last_thickness = 0, temp_value = 0
    cdef list a = [] # dummy list
    cdef list thickness_list = []
    cdef list init_guess_pos = []
    # TEST
    cdef np.ndarray[double, ndim=2] minima_ready = np.zeros((ende-start,Image_width),np.float)
    cdef unsigned int min_index, min_thickness, min_thickness_p_2_limit
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

    # calculate values which will be needed later in each loop
    min_thickness = thickness_list[0]
    min_thickness_p_2_limit = min_thickness + 2*thickness_limit

    print 'x ', len(data[0][0])
    print 'y ', len(data[0])
    print 'z ', len(data) 

    print "lookahead_min; ", lookahead_min
    ######################################### 
    # do calculations with thickness limits #
    #########################################
    #plt.figure(2)

    if use_thickness_limits:
        #print "using thickness limit: ", thickness_limit

        for row in xrange(len(data[0])):
            #print counter
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
                limit_counter = 0

                # write loop to consider the area around the current pixel

                if column !=0:
                    if thickness_ready[row,column-1] != 0:
                        last_thickness = thickness_ready[row,column-1]


                if last_thickness == 0:
                    if column >1 and row > 0:
                        if thickness_ready[row-1,column-1] != 0:
                            last_thickness = thickness_ready[row-1,column-1]

                # this will actually only run if the two cases before did not yield a thickness, maybe it should run always?
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
                            if thickness_ready[row-row_c,column-column_c] != 0:
                                last_thickness+=  thickness_ready[row-row_c,column-column_c]
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
                            if thickness_ready[row-row_c,column+column_c] != 0:
                                last_thickness+=  thickness_ready[row-row_c,column+column_c]
                                limit_counter += 1

                            # # check if other values in the area were found

                
                if limit_counter != 0:
                    #print "limit counter is:", limit_counter
                    # calculate average of the area
                    last_thickness = last_thickness/float(limit_counter)

                if last_thickness == 0:
                    last_thickness = init_guess

                # if the thickness in the area is in the thickness list, search for the index of that thickness and store it

                if last_thickness > (min_thickness_p_2_limit):
                    last_index = int(last_thickness - min_thickness)

                # find the minima in the profile
                if column>0:
                    minima_exp, min_index = peakdetect(data[:,row,column], waves, lookahead_min,lookahead_max, delta,average_window, minima_ready[row,column-1], min_index)
                    #minima_exp, min_index = peakdetect(intensity, waves, lookahead_min,lookahead_max, delta,average_window)
                else:
                    minima_exp, min_index = peakdetect(data[:,row,column], waves, lookahead_min,lookahead_max, delta,average_window)
                    
                # only do the calculation if a minimum was found
                if minima_exp != 0: 
                    # add point to the current minima map
                    minima_ready[row,column] = minima_exp
                    # calculate thickness
                    current_thickness = (Fit_limit(thickness,array_thickness_len_pos, array_length_block, minima_exp,tolerance,sim_minima_blocks,last_index,thickness_list, thickness_limit))
                            #if current_thickness != 0:

                    # write thickness to array
                    thickness_ready[row,column] = current_thickness

                # in case no minimum was found: assign the array a 0 or nan
                else:
                    thickness_ready[row,column] = 0
                    minima_ready[row,column] = np.nan
                    
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
    #np.savetxt("minima_map.txt",minima_ready)
    return thickness_ready#, sim_min_waves