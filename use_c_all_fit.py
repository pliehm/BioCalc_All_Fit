##############################################################################
### Script to calculate the cavity thicknes of a Biosensor for every pixel ###
##############################################################################

#############
### INPUT ###
#############


# enter folder with data, no subfolders allowed
data_folder = 'H9'


# chose wavelength range and step-width

wave_start = 550    # [nm]
wave_end = 750      # [nm]


# enter average deviation of experiment to simulation in nanometer, "1" is a good value to start

tolerance = 1

# define parameters for minima detection  

lookahead_min = 5 # something like peak width for the minima
delta = 7    # something like peak height

# chose elastomer thickness range , the smaller the range the faster the program. If you are not sure, just take d_min = 1000, d_max = 19000

d_min= 5000   # [nm]
d_max= 11000 # [nm]

use_thickness_limits = True # Enter "True" if you want to do calculation with thickness limits and "False" if not. I recommend starting with "False"

thickness_limit = 50 # [nm] enter the thickness limit (if thickness was found, next on will be: last_thickness +- thickness_limit)


# parameters for printing
# color map is calculated like (mean_thickness - color_min, mean_thickness + color_max) 

color_min = 500
color_max = 500

# enter True if you want to enable this smoothing
x_y_smooth = False
# enter sigma for the gaussian smoothing
x_y_sigma = 0.1

lambda_smooth = False
lambda_sigma = 1

############################
### END of INPUT SECTION ###
############################



#############################
#### start of the program ###
#############################

import cython_all_fit as Fit # import self-written cython code
import numpy as np
import time
import os 
import Image as im
from scipy import ndimage
import multiprocessing as mp
import matplotlib.pyplot as plt

t_a_start = time.time() # start timer for runtime measurement

if __name__ == '__main__':
    print os.listdir(data_folder)
    folder_list = os.listdir(data_folder)
    for folder in folder_list:
        # enter number of cpu cores, this has to be an integer number!
        # number of physical cores is a good start, but you can try with a larger as well

        multi_p = True   # True for multiprocessing, False for single core (Windows)
        cores = 4

        
        # enter name of simulation_file

        sim_file = 'Sim_0.5Cr_10Au_50SiO2_Elastomer_RT601_15Au_500_760nm.txt'

        lookahead_max = lookahead_min-1 # for the maxima --> should not be larger than lookahead_min

        # make wavelength list

        wave_step = 1       # [nm]
        
        waves=[]

        waves=[wave_start + i*wave_step for i in xrange((wave_end-wave_start)/wave_step + 1)]

        ## read image data 
        dateien=os.listdir(data_folder+'/'+folder)
        dateien.sort()
        
        #generates an empty array --> image grey values 
        alle=np.zeros(((wave_end-wave_start)/wave_step + 1,1024,1280),np.uint16)
        alle_x_y_smooth = np.zeros(((wave_end-wave_start)/wave_step + 1,1024,1280),np.uint16)
        # define function to convert the image-string to an array
        def image2array(Img):
            newArr= np.fromstring(Img.tostring(),np.uint8)
            newArr= np.reshape(newArr, (1024,1280))
            return newArr

        # read every image in folder and check if it is in the wavelength range --> 
        # write grey values into array
        
        counter=0
        print 'reading images from folder: ', folder
        if x_y_smooth == True:
            print 'smoothing of every wavelength image with sigma = ', x_y_sigma, ' before thickness calculation'
        else:
            print 'no smoothing of the wavelength images'    
        for i in xrange(len(dateien)):
            if dateien[i][-5:]=='.tiff':
                if int(dateien[i][:3]) >= wave_start and int(dateien[i][:3]) <= wave_end:
                    #print dateien[i]
                    #print counter
                    Img=im.open(data_folder + '/'+folder + '/' + dateien[i]).convert('L')
                    alle[counter]=image2array(Img)
                    # smoothing x-y direction
                    if x_y_smooth == True:
                        Img_s = ndimage.gaussian_filter(Img, sigma=x_y_sigma)
                        alle_x_y_smooth[counter] = image2array(Img_s)
                    counter+= 1

##################################
##### Section for smoothing ######
##################################

        if x_y_smooth == True:
            alle = alle_x_y_smooth.copy()

        if lambda_smooth == True:
            print 'smoothing over wavelength with sigma = ', lambda_sigma
            for zeile in range(1024):
                #print zeile
                for spalte in range(1280):
                    alle[:,zeile,spalte] = ndimage.gaussian_filter1d(alle[:,zeile,spalte],lambda_sigma)

#####################################
##### END section for smoothing #####
#####################################


        # read simulation file 
        print 'read simulated data'

        p= open(sim_file,'r')

        string=p.read()

        p.close()

        positioncounter = 0

        sim_waves=[]
        wave_block=[]
        s_waves_arrays = []
        position = 0

        for thisline in string.split('\n'):
            if ('\t' in thisline) == False and len(thisline) != 0:
                thickness = thisline
            if ('\t' in thisline) == True and int(thickness) >= d_min and int(thickness) <= d_max:
                positioncounter == 0
                for word in thisline.split('\t'):
                    if len(word)<6 and float(word) >= wave_start + lookahead_min and float(word)<= wave_end - lookahead_min: # use only minimas which are in the wave-range + lookahead_min
                        wave_block.append(float(word))
            if len(thisline) == 0 and int(thickness) >= d_min and int(thickness) <= d_max:
                sim_waves.append([np.uint16(thickness),np.uint16(len(wave_block)),np.uint16(position)]) # calculate length of the waveblock since it will be needed later
                s_waves_arrays.append(np.array(wave_block,dtype=np.float))
                position += len(wave_block)
                wave_block=[]

        print 'perform the calculations'
        t1 = time.time()
        # define queue for the multiprocessing

        if multi_p == True:

            def put_into_queue(start,ende,que,alle, sim_waves, waves, tolerance, lookahead_min, lookahead_max, delta, use_thickness_limits, thickness_limit):

                que.put(Fit.c_Fit_Pixel(start,ende,alle, sim_waves, waves, tolerance, lookahead_min, lookahead_max, delta,s_waves_arrays, use_thickness_limits, thickness_limit)) # calls the C-Fit-function
                #print 'Schlange ist satt'

            
            
            # devide the rows by the core-number --> to split it equally, assing the rest to the last process
            Zeile_Teil = 1024/cores
            Zeile_Rest = 1024%cores

            # start multiprocessing with queues

            Prozesse = []
            Queues = []

            for i in range(cores):
                Queues.append(mp.Queue())

            for i in range(cores):
                if i < cores-1:
                    Prozesse.append(mp.Process(target=put_into_queue,args=(i*Zeile_Teil,(i+1)*Zeile_Teil,Queues[i],alle[:,(i*Zeile_Teil):((i+1)*Zeile_Teil),:], sim_waves, waves, tolerance, lookahead_min, lookahead_max, delta, use_thickness_limits, thickness_limit)))
                if i == cores-1:
                    Prozesse.append(mp.Process(target=put_into_queue,args=(i*Zeile_Teil,(i+1)*Zeile_Teil+Zeile_Rest,Queues[i],alle[:,(i*Zeile_Teil):((i+1)*Zeile_Teil),:], sim_waves, waves, tolerance, lookahead_min, lookahead_max, delta, use_thickness_limits, thickness_limit)))
            for i in range(cores):
                Prozesse[i].start()
                

            # initialise array for thicknesses
            dicke = np.ndarray((0,1280),dtype=np.uint16)

            for i in range(cores):
                #print 'queuet', i
                dicke = np.append(dicke,Queues[i].get(),axis=0)

            for i in range(cores):
                #print 'joint', i
                Prozesse[i].join()

        if multi_p == False:
            start = 0
            ende = 1024

            dicke = Fit.c_Fit_Pixel(start,ende,alle, sim_waves, waves, tolerance, lookahead_min, lookahead_max, delta,s_waves_arrays, use_thickness_limits, thickness_limit)
        t2 = time.time()


        print t2-t1, 'seconds just for the calculation'

        # count not fitted values

        not_fitted = 1024*1280 - np.count_nonzero(dicke)
        not_fitted_percent = 100.0/(1024*1280)*not_fitted
        print 'not fitted values',not_fitted
        print 'in percent:', not_fitted_percent

        print 'write data to file'
        # use numpy function to save array to file, '0' and not '-' used for missing values
        HEADER = time.strftime("%d.%m.%Y at %H:%M:%S")+'\n' + 'folder with data = ' + folder + '\n' + 'simulation file = ' + sim_file + '\n' + 'wave_start = '+str(wave_start) + '\n' + 'wave_end = ' + str(wave_end) + '\n' + 'lookahead_min = ' + str(lookahead_min) + '\n'  + 'lookahead_max = ' + str(lookahead_max) + '\n' + 'delta = ' + str(delta) + ' delta was varied +-5'+ '\n' + 'tolerance = ' + str(tolerance) + '\n' + 'thickness limits used: ' + str(use_thickness_limits) + '\n' + 'thickness limits: ' + str(thickness_limit) + '\n' +  'not fitted values: ' + str(not_fitted) + ', percentage of whole image: ' + str(not_fitted_percent)  + '\n'
        if x_y_smooth == True:
            HEADER+= 'x_y_smoothing done with sigma = ' + str(x_y_sigma) + '\n'
        if lambda_smooth == True:
            HEADER+= 'lambda smoothing done with sigma = ' + str(lambda_sigma) + '\n'

        HEADER+= '\n'

        if x_y_smooth == True and lambda_smooth == False:
            file_name = folder + time.strftime("_%Y%m%d_%H%M%S")+'_x_y_sigma_' + str(x_y_sigma) + '_smoothed.txt'
            
        if lambda_smooth == True and x_y_smooth == False:
            file_name = folder + time.strftime("_%Y%m%d_%H%M%S")+'_lambda_sigma_' + str(lambda_sigma) + '_smoothed.txt'

        if lambda_smooth == True and x_y_smooth == True:
            file_name = folder + time.strftime("_%Y%m%d_%H%M%S")+'_x_y_sigma_' + str(x_y_sigma) + '_lambda_sigma_' +str(lambda_sigma) + '_smoothed.txt'

        if lambda_smooth == False and x_y_smooth == False: 
            file_name = folder + time.strftime("_%Y%m%d_%H%M%S")+'.txt'

        np.savetxt(file_name,dicke,fmt='%d',header=HEADER )

        ### script to replace a certain string or string combination in a file

        #t_replace_1 = time.time()

        # do it twice because 0 0 would not be - - but - 0 since the empty character before the second 0 has already been read
        for i in range(2):
            p = open(file_name,'r')

            string = p.read()

            p.close()

            string = string.replace(' 0 ', ' - ')
            string = string.replace('\n'+'0 ', '\n'+'- ')
            string = string.replace(' 0'+'\n', ' -'+'\n')

            p = open(file_name,'w')

            p.write(string)

            p.close()
        #print (time.time()-t_replace_1), ' seconds for replaceing 0 with -'
        print (time.time()-t_a_start), ' seconds for the whole program'

        # plot datta
        
        plt.figure(folder)
        plt.imshow(dicke)
        plt.clim(dicke.mean()-color_min,dicke.mean()+color_max)
        plt.colorbar()
        #plt.show()


