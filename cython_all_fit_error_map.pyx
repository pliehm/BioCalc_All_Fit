#cython: boundscheck=False
#cython: profile=False
#cython: cdivision = True

# import all the modules needed, make sure that numpy gets imported with cimport as well
import numpy as np
cimport numpy as np
cimport cython
import matplotlib.pyplot as plt
import time as time




# define types for cython, not sure if this is really necessary, but I found it in one of the examples I used to improve the speed
DTYPE1 = np.uint16
ctypedef np.uint16_t DTYPE_t
DTYPE2 = np.float
ctypedef np.float_t DTYPE_t2
DTYPE3 = np.uint32
ctypedef np.uint32_t DTYPE_t3

# --> datatype np.float is the same as double


# it is important to define as many variables as possible with cdef (numbers, arrays, lists), this is when cython really speeds up the code!


# define a faster function to calculate the absolute value of a number
cdef float _abs(float a): return a if a>=0 else -a

# switch on for plotting
#plt.ion()

######################################################
# Function to fit thickness without thickness limits #
######################################################

# function to find the best simulated thickness which fits to the measured set of minima, without limits
@cython.boundscheck(False)
cdef Fit(np.ndarray[DTYPE_t, ndim=1] thickness,np.ndarray[DTYPE_t3, ndim=1] array_thickness_len_pos,np.ndarray[DTYPE_t, ndim=1] array_length_block, np.ndarray[DTYPE_t2,ndim=1] exp_waves, float tolerance,np.ndarray[DTYPE_t2,ndim=1] all_minima,current_index,thickness_list):

# cython type definition of the used objects
 
    cdef list sim_min_waves = [[1000],[1000]] # dummy values to have two lists
    cdef unsigned short i, k,min_thickness_i, max_thickness_i
    cdef unsigned short L_exp_waves = len(exp_waves) # too not use the len() function too often
    cdef float summe=0
    cdef unsigned int counter = 0
    cdef unsigned int position, len_block 
    cdef unsigned int breaker = 0

    # first test if there are more than two minima
    if L_exp_waves > 2: 
        # do the following calculations for every simulated thickness
        for i in range(len(array_thickness_len_pos)):
            position = array_thickness_len_pos[i] 
            len_block = array_length_block[i]


            # the fitting will be run for different szenarios: same number of simulated and fitted minima, more sim minima, less sim minima. This yields more fitted values as the wavelength boundaries are better taken into account.

            # case for equal minimas (exp, sim)

            if array_length_block[i] == L_exp_waves: 
                breaker = 0
                summe=0
                # perform something like least-square with every exp-wavelength 
                for k in xrange(L_exp_waves): 
                    summe+=_abs(all_minima[position+k]-exp_waves[k])
                    if summe/L_exp_waves > tolerance:
                        breaker = 1
                        break
                if breaker == 1:
                    continue
                # append the thickness and error to a list
                sim_min_waves[0].append(thickness[i])
                sim_min_waves[1].append(summe/float(L_exp_waves))
                continue

            # do the same if number of exp and sim is not  equal
            if array_length_block[i] == (L_exp_waves + 1):
                breaker = 0
                summe=0
                # check if the first elements (exp and sim) or the last tow are not matching
                #if _abs(all_minima[position] - exp_waves[0]) > _abs(all_minima[position+thickness_len_pos[i][1]-1]-exp_waves[-1]):
                if _abs(all_minima[position] - exp_waves[0]) > _abs(all_minima[position+len_block-1]-exp_waves[-1]):
                    for k in xrange(L_exp_waves):
                        summe+= _abs(all_minima[position+k+1]-exp_waves[k])
                        if summe/(L_exp_waves + 1) > tolerance:
                            breaker = 1
                            break
                    if breaker == 1:
                        continue        
                    sim_min_waves[0].append(thickness[i])
                    sim_min_waves[1].append(summe/float(L_exp_waves+1))
                    continue
                else:
                    for k in xrange(L_exp_waves):
                        summe+= _abs(all_minima[position+k]-exp_waves[k])
                        if summe/L_exp_waves > tolerance:
                            breaker =1
                            break
                    if breaker == 1:
                        continue
                    sim_min_waves[0].append(thickness[i])
                    sim_min_waves[1].append(summe/float(L_exp_waves))
                    continue

            if array_length_block[i] == (L_exp_waves - 1):
                breaker = 0
                summe=0
                #sim_waves_part_part = sim_waves_part[2]
                if _abs(all_minima[position] - exp_waves[0]) > _abs(all_minima[position+len_block-1]-exp_waves[-1]):
                    for k in xrange(array_length_block[i]):
                        summe+= _abs(all_minima[position+k]-exp_waves[k+1])
                        if summe/(L_exp_waves-1) > tolerance:
                            breaker =1
                            break
                    if breaker ==1:
                        continue
                    sim_min_waves[0].append(thickness[i])
                    sim_min_waves[1].append(summe/float(L_exp_waves))
                    continue
                else:
                    for k in xrange(array_length_block[i]):
                        summe+= _abs(all_minima[position+k]-exp_waves[k])
                        if summe/(L_exp_waves-1) > tolerance:
                            breaker =1
                            break
                    if breaker ==1:
                        continue
                    sim_min_waves[0].append(thickness[i])
                    sim_min_waves[1].append(summe/float(L_exp_waves))
                    continue

            # if array_length_block[i] == (L_exp_waves + 2):
            #     breaker = 0
            #     summe=0
            #     #sim_waves_part_part = sim_waves_part[2]

            #     for k in xrange(array_length_block[i]-2):
            #         summe+= abs(all_minima[position+k+1]-exp_waves[k+1])
            #         if summe/(L_exp_waves-1) > tolerance:
            #             breaker = 0
            #             break
            #     if breaker == 1:
            #         continue
            #     sim_min_waves[0].append(thickness[i])
            #     sim_min_waves[1].append(summe/float(L_exp_waves))
            #     continue

# return the thickness with minimum value
        if  len(sim_min_waves[0])>1 and (min(sim_min_waves[1]) < tolerance):
            return sim_min_waves[0][sim_min_waves[1].index(min(sim_min_waves[1]))], thickness_list.index(sim_min_waves[0][sim_min_waves[1].index(min(sim_min_waves[1]))]), sim_min_waves
        else: 
            return 0, 0,0

    else:
        return 0, 0,0

###################################################
# Function to fit thickness with thickness limits #
###################################################

# almost the same function as above, the only difference is that only thicknesses within the thickness_limit are considered for the fitting, this improves speed and avoids jumps between cavities
@cython.boundscheck(False)
cdef Fit_limit(np.ndarray[DTYPE_t, ndim=1] thickness,np.ndarray[DTYPE_t3, ndim=1] array_thickness_len_pos,np.ndarray[DTYPE_t, ndim=1] array_length_block, np.ndarray[DTYPE_t2,ndim=1] exp_waves, float tolerance,np.ndarray[DTYPE_t2,ndim=1] all_minima, current_index, thickness_list, thickness_limit):

# cython type definition of the used objects
 
    cdef list sim_min_waves = [[1000],[1000]] # dummy values to have two lists
    cdef unsigned short i, k, min_thickness_i, max_thickness_i
    cdef unsigned short L_exp_waves = len(exp_waves) # too not use the len() function too often
    cdef float summe=0
    cdef unsigned int counter = 0
    cdef unsigned int position, len_block 
    cdef unsigned int breaker = 0

    if L_exp_waves > 2: # first test if there are more than two minima
       # for waves in s_waves_arrays:
        min_thickness_i = current_index-thickness_limit #
        max_thickness_i = current_index+thickness_limit #
        for i in range(min_thickness_i,max_thickness_i): # do the following calculations for every simulated thickness
            position = array_thickness_len_pos[i]
            len_block = array_length_block[i]
            if array_length_block[i] == L_exp_waves: # case for equal minimas (exp, sim)
                breaker = 0
                summe=0
                # perform something like least-square with every exp-wavelength 
                for k in xrange(L_exp_waves): 
                    summe+=_abs(all_minima[position+k]-exp_waves[k])
                    if summe/L_exp_waves > tolerance:
                        breaker = 1
                        break
                # append the thickness and error to a list
                if breaker == 1:
                    continue
                sim_min_waves[0].append(thickness[i])
                sim_min_waves[1].append(summe/float(L_exp_waves))
                continue

            # do the same if number of exp and sim is not  equal
            if array_length_block[i] == (L_exp_waves + 1):
                breaker=0
                summe=0
                # check if the first elements (exp and sim) or the last tow are not matching
                #if _abs(all_minima[position] - exp_waves[0]) > _abs(all_minima[position+thickness_len_pos[i][1]-1]-exp_waves[-1]):
                if _abs(all_minima[position] - exp_waves[0]) > _abs(all_minima[position+len_block-1]-exp_waves[-1]):
                    for k in xrange(L_exp_waves):
                        summe+= _abs(all_minima[position+k+1]-exp_waves[k])
                        if summe/L_exp_waves > tolerance:
                            breaker = 1
                            break
                    if breaker == 1:
                        continue
                    sim_min_waves[0].append(thickness[i])
                    sim_min_waves[1].append(summe/float(L_exp_waves))
                    continue
                else:
                    for k in xrange(L_exp_waves):
                        summe+= _abs(all_minima[position+k]-exp_waves[k])
                        if summe/L_exp_waves > tolerance:
                            breaker = 1
                            break
                    if breaker == 1:
                        continue
                    sim_min_waves[0].append(thickness[i])
                    sim_min_waves[1].append(summe/float(L_exp_waves))
                    continue

            if array_length_block[i] == (L_exp_waves - 1):
                breaker = 0
                summe=0
                #sim_waves_part_part = sim_waves_part[2]
                if _abs(all_minima[position] - exp_waves[0]) > _abs(all_minima[position+len_block-1]-exp_waves[-1]):
                    for k in xrange(array_length_block[i]):
                        summe+= _abs(all_minima[position+k]-exp_waves[k+1])
                        if summe/(L_exp_waves-1) > tolerance:
                            breaker = 0
                            break
                    if breaker == 1:
                        continue
                    sim_min_waves[0].append(thickness[i])
                    sim_min_waves[1].append(summe/float(L_exp_waves))
                    continue
                else:
                    for k in xrange(array_length_block[i]):
                        summe+= _abs(all_minima[position+k]-exp_waves[k])
                        if summe/(L_exp_waves-1) > tolerance:
                            breaker = 1
                            break
                    if breaker == 1:
                        continue
                    sim_min_waves[0].append(thickness[i])
                    sim_min_waves[1].append(summe/float(L_exp_waves))
                    continue

            # if array_length_block[i] == (L_exp_waves + 2):
            #     breaker = 0
            #     summe=0
            #     #sim_waves_part_part = sim_waves_part[2]

            #     for k in xrange(array_length_block[i]-2):
            #         summe+= abs(all_minima[position+k+1]-exp_waves[k+1])
            #         if summe/(L_exp_waves-1) > tolerance:
            #             breaker = 0
            #             break
            #     if breaker == 1:
            #         continue
            #     sim_min_waves[0].append(thickness[i])
            #     sim_min_waves[1].append(summe/float(L_exp_waves))
            #     continue
                
# return the thickness with minimum value
        if  len(sim_min_waves[0])>1 and (min(sim_min_waves[1]) < tolerance):
            return sim_min_waves[0][sim_min_waves[1].index(min(sim_min_waves[1]))], thickness_list.index(sim_min_waves[0][sim_min_waves[1].index(min(sim_min_waves[1]))]), sim_min_waves
        else: 
            return 0, 0,0

    else:
        return 0, 0,0

#########################################################
# function to get minima of the intensity profile array #
#########################################################

cdef list peakdetect(np.ndarray[DTYPE_t, ndim=1] y_axis, list x_axis = None, unsigned short lookahead_min=5, unsigned short lookahead_max=3, unsigned short delta = 0,unsigned short enhance_resolution=10, unsigned short average_window=5):    
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
    #    raise ValueError, "Lookahead must enhance_resolutione '1' or above in value"
    #if not (np.isscalar(delta) and delta >= 0):
    #    raise ValueError, "delta must be a positive number"

    #maxima and minima candidates are temporarily stored in 
    #mx and mn respectively
    cdef float mn = 100000
    cdef float mx = -100000
    cdef unsigned short index, window_len 
    cdef float x, y 
    # why am I using np.double_t here? could it not be np.float?
    cdef np.ndarray[np.double_t, ndim=1] y_interp, x_interp,s, y_temp,w 

#


    if enhance_resolution != 1:

        ################################################
        ### Interpolate data from 1nm to 0.1nm steps ###
        ################################################

        # interpolate x
        x_interp = np.linspace(x_axis[0],x_axis[-1],num=int((len(x_axis)-1)*enhance_resolution+1))

        # interpolate intensity
        y_interp = np.interp(x_interp,x_axis,y_axis)

        #########################################
        ### Smooth (by averaging) the profile ###
        #########################################

        # define window for averaging / smoothing
        window_len = average_window*enhance_resolution # this needs to be checked, maybe it can be smaller or has to be larger

        # create extrapolated array to have better averaging at the boundary
        s=np.r_[abs(2*y_interp[0]-y_interp[window_len:1:-1]), y_interp, abs(2*y_interp[-1]-y_interp[-1:-window_len:-1])]

        # normalised array for convolution (window, e.g. [0.1,0.1,0.1,...,0.1])
        w = np.ones(window_len,'d')/window_len
        
        # convolve signal with window
        y_temp = np.convolve(w, s, mode='same')

        # only take the original range and convert to list (minima finding is faster with list)
        y_axis_list = y_temp[window_len-1:-window_len+1].tolist()

        # convert to list because it's faster
        x_axis= x_interp.tolist()

        #plt.figure(1)
        #plt.clf()
        #plt.plot(x_axis,y_axis)
        #plt.ylim((1000,9000))

    else:
        y_axis_list = y_axis.tolist()
        #x_axis = x_axis_t

    #y_axis_list = y_axis#.tolist()
    length = len(y_axis_list)
    #print len(x_axis), len(y_axis)
    #print x_axis
    #print y_axis
    # plt.figure(1)
    # plt.clf()
    # plt.plot(x_axis,y_axis_list)
    # plt.ylim((1000,9000))
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

    # plt.vlines(min_peaks,min(y_axis)-100,max(y_axis))
    # plt.pause(0.00001)
    # if len(min_peaks)>0:
    #     min_peaks = [min_peaks[-1]]

    #result.write(str(min_peaks)+"\n")


    #plt.vlines(min_peaks,min(y_axis)-100,max(y_axis), color='r')

    #print min_peaks
    return min_peaks

    

####################################################################################
# Main function which calls the minima finding algorithm and the thickness fitting #
####################################################################################


# the following parameters are passed to the function:
def c_Fit_Pixel(unsigned int start,unsigned int ende, np.ndarray[DTYPE_t, ndim=3] all_images, list thickness_len_pos, list waves, float tolerance, unsigned short lookahead_min,unsigned short lookahead_max, unsigned short delta, unsigned short delta_vary, list list_all_minima_blocks, use_thickness_limits, unsigned int thickness_limit, unsigned short area_avrg, unsigned short enhance_resolution, unsigned short average_window, error_map_path):

    ########################################
    # definition of all the variable types #
    ########################################


    cdef unsigned int Image_width = len(all_images[0][0])
    cdef np.ndarray[DTYPE_t, ndim=2] result = np.zeros((ende-start,Image_width),np.uint16 )
    cdef np.ndarray[double, ndim = 2] error_map = np.zeros((ende-start,Image_width), dtype=np.float)
    cdef unsigned short column, row, column_c, row_c
    cdef np.ndarray[DTYPE_t, ndim=1] intensity
    cdef np.ndarray[DTYPE_t2,ndim=1] minima_exp
    cdef unsigned int counter=start 
    # create different arrays to speed up the search for the right thickness
    cdef np.ndarray[DTYPE_t2,ndim=1] all_minima
    cdef np.ndarray[DTYPE_t,ndim=1] array_length_block, thickness
    cdef np.ndarray[DTYPE_t3,ndim=1] array_thickness_len_pos

    cdef unsigned int current_thickness = 0, current_index = 0, last_index = 0, limit_counter = 0
    cdef float last_thickness = 0
    cdef list a = [] # dummy list
    cdef list thickness_list = []
    cdef float t_dfit, t_dfit_s, t_

    # make an array which contains all minima, but not seperated as before
    for block in list_all_minima_blocks:
        for i in range(len(block)):
            a.append(block[i])
    all_minima = np.array(a,dtype=np.float)

    # make arrays of positions, lengths, thickness
    array_thickness_len_pos = np.array(zip(*thickness_len_pos)[2],dtype=np.uint32)
    array_length_block = np.array(zip(*thickness_len_pos)[1],dtype=np.uint16)
    thickness = np.array(zip(*thickness_len_pos)[0],dtype=np.uint16)

    # make a list which contains all the thicknesses --> kind of reduntant because I alrady have the thickness array, but I want to use the index function in the list, but the speed of the array. Most likely there is a better solution than using two variables for this.
    for i in range(len(thickness_len_pos)):
        thickness_list.append(int(thickness_len_pos[i][0]))


    # print to the user the simulation dimensions
    print 'x ', len(all_images) 
    print 'y ', len(all_images[0])
    print 'z ', len(all_images[0][0])


    ######################################### 
    # do calculations with thickness limits #
    #########################################

    # should be used as a standard

    if use_thickness_limits:
        #print "using thickness limit: ", thickness_limit

        # perform the follwing loop for all rows
        for row in range(len(all_images[0])):
            print counter
            counter+=1
            # perform the following loop for all coloumns
            for column in xrange(Image_width):
                # set parameters for averaging to 0
                last_thickness = 0
                last_index = 0
                current_thickness = 0
                current_index = 0
                limit_counter = 0

                # write loop to consider the area around the current pixel
                # this basically looks at values around the current pixel --> builds sum for all non-zero pixels and devides by the amount of considered pixels to estimate the average thickness as a guess 

                # iterate over rows, distance from current row is 0 to area_avrg
                for row_c in range(area_avrg+1):
                    # if the considered row is not in the range, stop iteration over rows
                    if row  < row_c:
                        break

                    # block to consider values with subtracted column values
                    
                    # iterate over column between 1 to area_avrg to not sonsider the same column as current 
                    for column_c in range(1,area_avrg+1):
                        # if the considered column is below the data range, stop loop
                        if column < column_c:
                            break

                        # check if the considered value is not zero, if so add the value to a sum and increase the count for later averaging
                        if result[row-row_c][column-column_c] != 0:
                            last_thickness+=  result[row-row_c][column-column_c]
                            limit_counter += 1

                    # block to consider values with added column values and same column
                    for column_c in range(area_avrg+1):
                        # if the considered row is the current row, stop, because there are no values
                        if (row_c == 0):
                            break
                        # if the considered column is above the data range, stop loop
                        if (column + column_c) >= Image_width:
                            break
                        # check if the considered value is not zero, if so add the value to a sum and increase the count for later averaging
                        if result[row-row_c][column+column_c] != 0:
                            last_thickness+=  result[row-row_c][column+column_c]
                            limit_counter += 1

                # check if other values in the area were found
                if limit_counter != 0:
                    # calculate average of the area
                    last_thickness = last_thickness/float(limit_counter)

                # if the thickness in the area is in the thickness list, search for the index of that thickness and store it --> this will be used to guess the thickness
                if last_thickness > (thickness_list[0] + 2*thickness_limit):
                    last_index =  last_thickness - thickness_list[0]# thickness_list.index(int(last_thickness))

                # get array with the intensity profile for the current pixel
                intensity = all_images[:,row, column]
                # find the minima in the profile
                minima_exp = np.array(peakdetect(intensity, waves, lookahead_min*enhance_resolution,lookahead_max*enhance_resolution, delta,enhance_resolution,average_window),dtype=np.float)

                # start calculations with limits

                # if a guess for the thickness has been found:
                if (last_thickness != 0) and (last_index > 0):
                    # call the thickness fitting function
                    current_thickness, current_index, sim_min_waves = (Fit_limit(thickness,array_thickness_len_pos, array_length_block, minima_exp,tolerance,all_minima,last_index,thickness_list, thickness_limit))
                    # if a thickness was fitted, save it
                    if current_thickness != 0:
                        result[row][column]=current_thickness
                        try:
                            len(sim_min_waves)
                            error_map[row][column] = min(sim_min_waves[1])
                        except TypeError:
                            pass  
                    
                    # if no thickness could be fitted, vary delta (fitting parameter in the minima algorithm) and try again
                    if current_thickness == 0: 
                        for new_delta in range(1,5):
                            # test with "+" new_delta*delta_vary
                            minima_exp = np.array(peakdetect(intensity, waves, lookahead_min*enhance_resolution,lookahead_max*enhance_resolution, delta + new_delta*delta_vary,enhance_resolution,average_window),dtype=np.float)
                            current_thickness, current_index, sim_min_waves = (Fit_limit(thickness,array_thickness_len_pos, array_length_block, minima_exp,tolerance,all_minima, last_index, thickness_list,thickness_limit))
                            
                            if current_thickness != 0:
                                result[row][column] = current_thickness
                                try:
                                    len(sim_min_waves)
                                    error_map[row][column] = min(sim_min_waves[1])
                                except TypeError:
                                    pass
                                break
                            # test with "-" new_delta*delta_vary
                            minima_exp = np.array(peakdetect(intensity, waves, lookahead_min*enhance_resolution,lookahead_max*enhance_resolution, delta - new_delta*delta_vary,enhance_resolution,average_window),dtype=np.float)
                            current_thickness, current_index, sim_min_waves = (Fit_limit(thickness,array_thickness_len_pos, array_length_block, minima_exp,tolerance,all_minima, last_index,thickness_list,thickness_limit))
                
                            if current_thickness != 0:
                                result[row][column] = current_thickness
                                try:
                                    len(sim_min_waves)
                                    error_map[row][column] = min(sim_min_waves[1])
                                except TypeError:
                                    pass
                                break



                # if there is still no thickness fitted, try to fit it without thickness limits
                if current_thickness == 0:
                    current_thickness, current_index, sim_min_waves = (Fit(thickness,array_thickness_len_pos, array_length_block, minima_exp,tolerance,all_minima, current_index,thickness_list))
                    result[row][column] = current_thickness
                    try:
                        len(sim_min_waves)
                        error_map[row][column] = min(sim_min_waves[1])
                    except TypeError:
                        pass

                # if no thickness could be fitted, vary delta and try again (still without limits)
                if current_thickness == 0: 
                    for new_delta in range(1,5):
                        minima_exp = np.array(peakdetect(intensity, waves, lookahead_min*enhance_resolution,lookahead_max*enhance_resolution, delta + new_delta*delta_vary,enhance_resolution,average_window),dtype=np.float)
                        current_thickness, current_index, sim_min_waves = (Fit(thickness,array_thickness_len_pos, array_length_block, minima_exp,tolerance,all_minima, current_index, thickness_list))
                        
                        if current_thickness != 0:
                            result[row][column] = current_thickness
                            try:
                                len(sim_min_waves)
                                error_map[row][column] = min(sim_min_waves[1])
                            except TypeError:
                                pass
                            break
                        minima_exp = np.array(peakdetect(intensity, waves, lookahead_min*enhance_resolution,lookahead_max*enhance_resolution, delta - new_delta*delta_vary,enhance_resolution,average_window),dtype=np.float)
                        current_thickness, current_index, sim_min_waves = (Fit(thickness,array_thickness_len_pos, array_length_block, minima_exp,tolerance,all_minima, current_index,thickness_list))
                        
                        if current_thickness != 0:
                            result[row][column] = current_thickness
                            try:
                                len(sim_min_waves)
                                error_map[row][column] = min(sim_min_waves[1])
                            except TypeError:
                                pass
                            break

    ################################################
    # Do the calculations without thickness limits #
    ################################################


    # for comments, see the section before, it's basically the same, just another function is called
    else:
        #print "using no thickness limits"
        for row in range(len(all_images[0])):
            print counter
            counter+=1
            for column in xrange(Image_width):
                intensity = all_images[:,row, column]
                minima_exp = np.array(peakdetect(intensity, waves, lookahead_min*enhance_resolution,lookahead_max*enhance_resolution, delta,enhance_resolution,average_window),dtype=np.float)
                current_thickness, current_index, sim_min_waves = (Fit(thickness,array_thickness_len_pos, array_length_block, minima_exp,tolerance,all_minima,current_index,thickness_list))
                result[row][column] = current_thickness
                try:
                    len(sim_min_waves)
                    error_map[row][column] = min(sim_min_waves[1])
                except TypeError:
                    pass   

                # if no thickness was fitted, try again with other delta
                if current_thickness == 0: 
                    for new_delta in range(1,5):
                        minima_exp = np.array(peakdetect(intensity, waves, lookahead_min*enhance_resolution,lookahead_max*enhance_resolution, delta + new_delta*delta_vary,enhance_resolution,average_window),dtype=np.float)
                        current_thickness, current_index, sim_min_waves  = (Fit(thickness,array_thickness_len_pos, array_length_block, minima_exp,tolerance,all_minima, current_index, thickness_list))
                        
                        if current_thickness != 0:
                            result[row][column] = current_thickness
                            try:
                                len(sim_min_waves)
                                error_map[row][column] = min(sim_min_waves[1])
                            except TypeError:
                                pass
                            break
                        minima_exp = np.array(peakdetect(intensity, waves, lookahead_min*enhance_resolution,lookahead_max*enhance_resolution, delta - new_delta*delta_vary,enhance_resolution,average_window),dtype=np.float)
                        current_thickness, current_index, sim_min_waves = (Fit(thickness,array_thickness_len_pos, array_length_block, minima_exp,tolerance,all_minima,current_index,thickness_list))
                        
                        if current_thickness != 0:
                            result[row][column] = current_thickness
                            try:
                                len(sim_min_waves)
                                error_map[row][column] = min(sim_min_waves[1])
                            except TypeError:
                                pass
                            break

    # write error map to file
    file_name = error_map_path + "_error_map" + time.strftime("_%Y%m%d_%H%M%S")+'.txt'
    error_map[error_map==0]= np.nan
    np.savetxt(file_name,error_map,fmt='%0.5f')
    # feed the fitted thicknesses back to the main program, with sim_min_waves one can access the error, this can be done better --> in progress
    return result#, sim_min_waves