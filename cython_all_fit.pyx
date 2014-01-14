#cython: boundscheck=False
#cython: profile=False
#cython: cdivision = True
import numpy as np
cimport numpy as np
cimport cython


DTYPE1 = np.uint16
ctypedef np.uint16_t DTYPE_t
DTYPE2 = np.float
ctypedef np.float DTYPE_t2


# fit simulated and experimental minima to calculate thickness
# problem: speed!! --> one main problem is the access of the large list sim_waves --> can this be written as an array?

cdef float _abs(float a): return a if a>=0 else -a

@cython.boundscheck(False)
cdef Fit(np.ndarray[DTYPE_t, ndim=1] thickness,np.ndarray[DTYPE_t, ndim=1] array_thickness_pos,np.ndarray[DTYPE_t, ndim=1] array_length_block, np.ndarray[double,ndim=1] exp_waves, unsigned short tolerance,np.ndarray[double,ndim=1] sim_wave_blocks_array,current_index,thickness_list):

# cython type definition of the used objects
 
    cdef list sim_min_waves = [[1000],[1000]] # dummy values to have two lists
    cdef unsigned short i, k,min_thickness_i, max_thickness_i
    cdef unsigned short L_exp_waves = len(exp_waves) # too not use the len() function too often
    cdef float summe=0
    cdef unsigned int counter = 0
    cdef unsigned int position, len_block 
    cdef unsigned int breaker = 0

    if L_exp_waves > 2: # first test if there are more than two minima
        for i in range(len(array_thickness_pos)): # do the following calculations for every simulated thickness
            position = array_thickness_pos[i] 
            len_block = array_length_block[i]

            if array_length_block[i] == L_exp_waves: # case for equal minimas (exp, sim)
                breaker = 0
                summe=0
                # perform something like least-square with every exp-wavelength 
                for k in xrange(L_exp_waves): 
                    summe+=_abs(sim_wave_blocks_array[position+k]-exp_waves[k])
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
                #if _abs(sim_wave_blocks_array[position] - exp_waves[0]) > _abs(sim_wave_blocks_array[position+thickness_pos[i][1]-1]-exp_waves[-1]):
                if _abs(sim_wave_blocks_array[position] - exp_waves[0]) > _abs(sim_wave_blocks_array[position+len_block-1]-exp_waves[-1]):
                    for k in xrange(L_exp_waves):
                        summe+= _abs(sim_wave_blocks_array[position+k+1]-exp_waves[k])
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
                        summe+= _abs(sim_wave_blocks_array[position+k]-exp_waves[k])
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
                if _abs(sim_wave_blocks_array[position] - exp_waves[0]) > _abs(sim_wave_blocks_array[position+len_block-1]-exp_waves[-1]):
                    for k in xrange(array_length_block[i]):
                        summe+= _abs(sim_wave_blocks_array[position+k]-exp_waves[k+1])
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
                        summe+= _abs(sim_wave_blocks_array[position+k]-exp_waves[k])
                        if summe/(L_exp_waves-1) > tolerance:
                            breaker =1
                            break
                    if breaker ==1:
                        continue
                    sim_min_waves[0].append(thickness[i])
                    sim_min_waves[1].append(summe/float(L_exp_waves))
                    continue

# return the thickness with minimum value
        if  len(sim_min_waves[0])>1 and (min(sim_min_waves[1]) < tolerance):
            return sim_min_waves[0][sim_min_waves[1].index(min(sim_min_waves[1]))], thickness_list.index(sim_min_waves[0][sim_min_waves[1].index(min(sim_min_waves[1]))])
        else: 
            return 0, 0

    else:
        return 0, 0

@cython.boundscheck(False)
cdef Fit_2(np.ndarray[DTYPE_t, ndim=1] thickness,np.ndarray[DTYPE_t, ndim=1] array_thickness_pos,np.ndarray[DTYPE_t, ndim=1] array_length_block, np.ndarray[double,ndim=1] exp_waves, unsigned short tolerance,np.ndarray[double,ndim=1] sim_wave_blocks_array, current_index, thickness_list, thickness_limit):

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
            position = array_thickness_pos[i]
            len_block = array_length_block[i]
            if array_length_block[i] == L_exp_waves: # case for equal minimas (exp, sim)
                breaker = 0
                summe=0
                # perform something like least-square with every exp-wavelength 
                for k in xrange(L_exp_waves): 
                    summe+=_abs(sim_wave_blocks_array[position+k]-exp_waves[k])
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
                #if _abs(sim_wave_blocks_array[position] - exp_waves[0]) > _abs(sim_wave_blocks_array[position+thickness_pos[i][1]-1]-exp_waves[-1]):
                if _abs(sim_wave_blocks_array[position] - exp_waves[0]) > _abs(sim_wave_blocks_array[position+len_block-1]-exp_waves[-1]):
                    for k in xrange(L_exp_waves):
                        summe+= _abs(sim_wave_blocks_array[position+k+1]-exp_waves[k])
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
                        summe+= _abs(sim_wave_blocks_array[position+k]-exp_waves[k])
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
                if _abs(sim_wave_blocks_array[position] - exp_waves[0]) > _abs(sim_wave_blocks_array[position+len_block-1]-exp_waves[-1]):
                    for k in xrange(array_length_block[i]):
                        summe+= _abs(sim_wave_blocks_array[position+k]-exp_waves[k+1])
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
                        summe+= _abs(sim_wave_blocks_array[position+k]-exp_waves[k])
                        if summe/(L_exp_waves-1) > tolerance:
                            breaker = 1
                            break
                    if breaker == 1:
                        continue
                    sim_min_waves[0].append(thickness[i])
                    sim_min_waves[1].append(summe/float(L_exp_waves))
                    continue

# return the thickness with minimum value
        if  len(sim_min_waves[0])>1 and (min(sim_min_waves[1]) < tolerance):
            return sim_min_waves[0][sim_min_waves[1].index(min(sim_min_waves[1]))], thickness_list.index(sim_min_waves[0][sim_min_waves[1].index(min(sim_min_waves[1]))])
        else: 
            return 0, 0

    else:
        return 0, 0

# function to get minima of one array


cdef list peakdetect(y_axis, x_axis = None, unsigned short lookahead_min=5, unsigned short lookahead_max=3, unsigned short delta = 0):
    
    # define output container
    #cdef list max_peaks=[]
    cdef list min_peaks = []
    cdef list dump = [] # used to pop the first hit which almost always is false
    cdef list y_axis_list = y_axis.tolist() # convert array to list, min() is faster for list

    # check input data --> this makes the algorithm 5 times slower
    #x_axis, y_axis = _datacheck_peakdetect(x_axis, y_axis) 
    
    # store data length for later use
    cdef unsigned int length = len(y_axis)

    #perform some checks
    #if lookahead < 1:
    #    raise ValueError, "Lookahead must be '1' or above in value"
    #if not (np.isscalar(delta) and delta >= 0):
    #    raise ValueError, "delta must be a positive number"

    #maxima and minima candidates are temporarily stored in 
    #mx and mn respectively
    cdef int mn = 1000
    cdef int mx = -1000
    cdef unsigned short x,y, index

    for index, (x,y) in enumerate(zip(x_axis[:-lookahead_min], y_axis[:-lookahead_min])):
        
        if y > mx:
            mx = y
            mxpos = x

        if y < mn:
            mn = y
            mnpos = x

        #### look for max ####
        
        if y < mx-delta and mx != 1000:
            #Maxima peak candidate found
            # lool ahead in signal to ensure that this is a peak and not jitter
            if max(y_axis_list[index:index+lookahead_max]) < mx:
                #max_peaks.append([mxpos, mx])
                dump.append(True)
                #set algorithm to only find minima now
                
                mx = 1000
                mn = 1000
                if index+lookahead_min >= length:
                    #end is within lookahead no more peaks can be found
                    break
                continue

        #### look for min ####    
        
        if y > mn+delta and mn != -1000:
            #Minima peak candidate found
            # look ahead in signal to ensure that this is a peak and not jitter
            if min(y_axis_list[index:index+lookahead_min]) > mn:
                min_peaks.append(mnpos)
                dump.append(False)
                #set algorithm to only find maximum now
                mn = -1000
                mx = -1000
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
    
    return min_peaks

    
# find reflection minima for every pixel

def c_Fit_Pixel(unsigned int start,unsigned int ende, np.ndarray[DTYPE_t, ndim=3] data, list thickness_pos, list waves, unsigned short tolerance, unsigned short lookahead_min,unsigned short lookahead_max, unsigned short delta, list sim_wave_blocks_list, use_thickness_limits, unsigned int thickness_limit):
    cdef np.ndarray[DTYPE_t, ndim=2] thickness_ready = np.zeros((ende-start,1280),np.uint16 )
    cdef unsigned short spalte, zeile
    cdef np.ndarray[DTYPE_t, ndim=1] intensity
    cdef np.ndarray[double,ndim=1] minima_exp
    cdef unsigned int counter=start, 
    cdef np.ndarray[double,ndim=1] sim_wave_blocks_array
    cdef np.ndarray[DTYPE_t,ndim=1] array_thickness_pos, array_length_block, thickness
    cdef unsigned int current_thickness = 0, current_index = 0, last_thickness = 0, last_index = 0
    # build another sim_waves_m list with different list class
    cdef list a = [] # dummy list
    cdef list thickness_list = []

    for i in range(len(thickness_pos)):
        thickness_list.append(int(thickness_pos[i][0]))

    for block in sim_wave_blocks_list:
        for i in range(len(block)):
            a.append(block[i])
    sim_wave_blocks_array = np.array(a,dtype=np.float)

    array_thickness_pos = np.array(zip(*thickness_pos)[2],dtype=np.uint16)
    array_length_block = np.array(zip(*thickness_pos)[1],dtype=np.uint16)
    thickness = np.array(zip(*thickness_pos)[0],dtype=np.uint16)

    if use_thickness_limits:
        #print "using thickness limit: ", thickness_limit

        for zeile in range(start, ende):
            print counter
            counter+=1
            for spalte in xrange(1280):
                intensity = data[:,zeile, spalte]
                minima_exp = np.array(peakdetect(intensity, waves, lookahead_min,lookahead_max, delta),dtype=np.float)
                if last_thickness != 0:
                    current_thickness, current_index = (Fit_2(thickness,array_thickness_pos, array_length_block, minima_exp,tolerance,sim_wave_blocks_array,last_index,thickness_list, thickness_limit))
                    if current_thickness != 0:
                        thickness_ready[zeile-start][spalte]=current_thickness
                        last_thickness, last_index = current_thickness, current_index
                    
                    if current_thickness == 0: 
                        for new_delta in range(1,5):
                            minima_exp = np.array(peakdetect(intensity, waves, lookahead_min,lookahead_max, delta + new_delta),dtype=np.float)
                            current_thickness, current_index = (Fit_2(thickness,array_thickness_pos, array_length_block, minima_exp,tolerance,sim_wave_blocks_array, last_index, thickness_list,thickness_limit))
                            
                            if current_thickness != 0:
                                thickness_ready[zeile-start][spalte] = current_thickness
                                break
                            minima_exp = np.array(peakdetect(intensity, waves, lookahead_min,lookahead_max, delta - new_delta), dtype=np.float)
                            current_thickness, current_index = (Fit_2(thickness,array_thickness_pos, array_length_block, minima_exp,tolerance,sim_wave_blocks_array, last_index,thickness_list,thickness_limit))
                
                            if current_thickness != 0:
                                thickness_ready[zeile-start][spalte] = current_thickness
                                break




                if current_thickness == 0:
                    current_thickness, current_index = (Fit(thickness,array_thickness_pos, array_length_block, minima_exp,tolerance,sim_wave_blocks_array, current_index,thickness_list))
                    thickness_ready[zeile-start][spalte] = current_thickness

                if current_thickness == 0: 
                    for new_delta in range(1,5):
                        minima_exp = np.array(peakdetect(intensity, waves, lookahead_min,lookahead_max, delta + new_delta),dtype=np.float)
                        current_thickness, current_index = (Fit(thickness,array_thickness_pos, array_length_block, minima_exp,tolerance,sim_wave_blocks_array, current_index, thickness_list))
                        
                        if current_thickness != 0:
                            thickness_ready[zeile-start][spalte] = current_thickness
                            break
                        minima_exp = np.array(peakdetect(intensity, waves, lookahead_min,lookahead_max, delta - new_delta), dtype=np.float)
                        current_thickness, current_index = (Fit(thickness,array_thickness_pos, array_length_block, minima_exp,tolerance,sim_wave_blocks_array, current_index,thickness_list))
                        
                        if current_thickness != 0:
                            thickness_ready[zeile-start][spalte] = current_thickness
                            break
                if current_thickness != 0:
                    last_thickness, last_index = current_thickness, current_index

    else:
        #print "using no thickness limits"
        for zeile in range(start, ende):
            print counter
            counter+=1
            for spalte in xrange(1280):
                intensity = data[:,zeile, spalte]
                minima_exp = np.array(peakdetect(intensity, waves, lookahead_min,lookahead_max, delta),dtype=np.float)
                current_thickness, current_index = (Fit(thickness,array_thickness_pos, array_length_block, minima_exp,tolerance,sim_wave_blocks_array,current_index,thickness_list))
                thickness_ready[zeile-start][spalte] = current_thickness   

                if current_thickness == 0: # if no thickness was fitted, try again 
                    for new_delta in range(1,5):
                        minima_exp = np.array(peakdetect(intensity, waves, lookahead_min,lookahead_max, delta + new_delta),dtype=np.float)
                        current_thickness, current_index = (Fit(thickness,array_thickness_pos, array_length_block, minima_exp,tolerance,sim_wave_blocks_array, current_index, thickness_list))
                        
                        if current_thickness != 0:
                            thickness_ready[zeile-start][spalte] = current_thickness
                            break
                        minima_exp = np.array(peakdetect(intensity, waves, lookahead_min,lookahead_max, delta - new_delta), dtype=np.float)
                        current_thickness, current_index = (Fit(thickness,array_thickness_pos, array_length_block, minima_exp,tolerance,sim_wave_blocks_array,current_index,thickness_list))
                        
                        if current_thickness != 0:
                            thickness_ready[zeile-start][spalte] = current_thickness
                            break


    return thickness_ready