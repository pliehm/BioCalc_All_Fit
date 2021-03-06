##############################################################################
### Script to calculate the cavity thicknes of a Biosensor for every pixel ###
##############################################################################

#############
### INPUT ###
#############


# enter all the folders which contain subfolders with the image files, 
# the folder with the image files MUST be in another folder (its name has to be entered here)
# e.g. you have 5 folders with images: 1,2,3,4,5 --> all these 5 folders have to be e.g. in a folder
# named "data". "data" is what you would enter in the list below. You can enter more then one folder 
# in this list (e.g. for 3 different long-time measurements). But you can also run different instances of the program to make use of multiple cores.

data = ['data']

# enter name of simulation_file, copy and paste the file name of the
# simulation file corresponding to your layer structure

sim_file = '2015_10_26_standard.txt'

# chose wavelength range for calculation
wave_start = 650    # [nm]
wave_end = 690     # [nm]

# choose if the images are all in one tiff stack or in separate files (separate files was how we used to have it, tiff stack is new but will make data transfer faster)

tiff_stack = True

# enter wavelength range of stack (only needed if tiff_stack = True)

stack_wave_start = 550 # [nm]
stack_wave_end = 750 # [nm]

# enter a value to apply binning to run the calculation faster
binning = 1


# chose elastomer thickness range , the smaller the range the faster the program. If you are not sure, just take d_min = 1000, d_max = 19000

d_min= 6000 # [nm]
d_max= 10000 # [nm]

####################
# Advanced options #
####################

# one minimum fit option --> only the last fitted minimum in the wavelength range will be considered and fitted

one_minimum_fit = True

# guess the thickness at a certain position. THIS HAS TO BE A LIST!!! e.g. [7500] or [7500,7800,7400]
# the values have to correspond to the folders
init_guess = [7667] 


# enter average allowed deviation of experiment to simulation in nanometer, "1" is a good value to start

tolerance = 1

# define parameters for minima detection  

# something like peak width for the minima, 5 is a good value
# should be larger for thinner cavities and smaller for bigger cavities (maybe FSR/4)

lookahead_min = 5 

# set parameter for how many waves shall be interpolated b = 1, --> none, b=10 --> 0.1nm steps
# 5 is a good value, 10 does not improve the result that much
enhance_resolution = 5

# average window, parameter which defines the window length for the smoothing: window_len = enhanced_resolution*average_window --> 1 equals 1 nm window

average_window = 5

 # Enter "True" if you want to do calculation with thickness limits and "False" if not. I recommend starting with "True" if you see artifacts try to use "False"
use_thickness_limits = True

# [nm] enter the thickness limit (if thickness was found, next on will be: last_thickness +- thickness_limit) --> this has to be smaller than 230nm/2, better 230nm/4
thickness_limit = 50 

# this number defines how many pixel are considerd for an average to guess the new thickness, e.g.: 2 means that all pixels in a range of 2 rows above, two columns to the left and the right are averaged, that makes 12px, 1 --> 4px, 2 --> 12px, 3--> 24px
area_avrg = 2 



# check what your experimental data for each stack is, multiply by the enhance_resolution factor, multiply with number of 
# simulations you are planning to run in parallel, multiply with 2. For example: data = 550MB, enhanced_resolution = 5 
# --> 550*5 = 2750*2 = 5500MB. This is how much memory one calculation occupies, if you want to run 4 in parallel: 
# 5500*4 = 22000 MB = 22 GB You would need 22 GB of RAM for this!!
# If this is more (or close to) the amount of RAM you have, chose "save_ram = True". This will be a bit slower but more memory efficient.

save_ram = False

############################
### END of INPUT SECTION ###
############################

#########################################
# Developer parameters, DONT CHANGE!!!! #
#########################################

plot_error = False # standard is False

#############################
#### start of the program ###
#############################

# load all the python modules needed, they should all be part of the anaconda python distribution
# import self-written cython code
import cython_all_fit as Fit
import numpy as np
import scipy as sp
import time
import os 
from PIL import Image as im
import matplotlib.pyplot as plt
from skimage.io import imread
from skimage import transform

 # only import other modules if needed
if one_minimum_fit == True:
    import cython_all_fit_one_minimum as Fit_one

if plot_error == True:
    import cython_all_fit_error_map  as Fit_error

# change version of the release here which will be included in the results files
version = 'BioCalc 2.3'

t_a_start = time.time() # start timer for runtime measurement




# define function for interpolation

def self_interpolate(raw,enhance_resolution):

    Image_width = raw.shape[2]
    Image_height = raw.shape[1]
    Image_depth = raw.shape[0]

    arr_interp = np.zeros(((Image_depth-1)*enhance_resolution+1,Image_height,Image_width),np.uint16)


    interp_counter = 0
    stack_counter = 0
    raw = raw.astype(np.float32)
    t1= time.time()
    for i in range(len(arr_interp)):
        if interp_counter == 0 and stack_counter==0:
            #print "c"
            arr_interp[i] = raw[stack_counter]
            temp_arr = (raw[stack_counter+1]-raw[stack_counter])/enhance_resolution
            interp_counter+=1
            continue
        if interp_counter != enhance_resolution:
            #print "b"
            arr_interp[i] = raw[stack_counter]+interp_counter*temp_arr
            interp_counter += 1
            continue

        if interp_counter == enhance_resolution:
            stack_counter += 1
            if stack_counter != Image_depth-1:
                arr_interp[i] = raw[stack_counter]
                temp_arr = (raw[stack_counter+1]-raw[stack_counter])/enhance_resolution
            
                interp_counter = 1
                continue
            else:
                arr_interp[i] = raw[stack_counter]
                continue

    return arr_interp

#######################################################################
# Read the simulation file and get the minima for the thickness range #
#######################################################################


# read simulation file 
print 'read simulated data'

# open file
p= open('Simulated_minima/' + sim_file,'r')

# read the whole content of the file into a string
string=p.read()

# close the file
p.close()

# list which will contain, thickness, length of a set of wavelengths, starting position of the block --> this is all useful to have to improve speed
thickness_len_pos=[]

# current set of minima for one thickness (temporary list)
minima_block=[]

# list which contains all minima_blocks arrays one after another --> in combination with the list thickness_len_pos one has all the data needed
list_all_minima_blocks = []

# variable which defines at which position a set of minima starts
position = 0

# Start reading the file

# read every line of the string
for thisline in string.split('\n'):
    # check if the current line contains a thickness
    if ('\t' in thisline) == False and len(thisline) != 0:
        thickness = thisline
    # check if this is a line with no thickness, but one is still in the thickness range
    if ('\t' in thisline) == True and int(thickness) >= d_min and int(thickness) <= d_max:
        # split line into words
        for word in thisline.split('\t'):
            # check word if it is a wavelength, only append it if its in the wavelength range +- lookahead, len(word)<7 assures that word is a wavelength and not the reflectivity value
            if len(word)<7 and float(word) >= wave_start + lookahead_min and float(word)<= wave_end - lookahead_min: 
                # append the wavelength to the current block of wavelengths
                minima_block.append(float(word))

    # check if the line is empty and inside the thickness range
    if len(thisline) == 0 and int(thickness) >= d_min and int(thickness) <= d_max:

        # append thickness, length of waveblock, position of block to a list, convert to 16 or 32 bit, that is important to assure that the values are not corrupted because they are outside the right data range
        thickness_len_pos.append([np.uint16(thickness),np.uint16(len(minima_block)),np.uint32(position)]) # calculate length of the waveblock since it will be needed later

        # append waveblock as an array (faster data handling later on)
        list_all_minima_blocks.append(np.array(minima_block,dtype=np.float))

        # update the current starting position of the next waveblock
        position += len(minima_block)
        # clear the waveblock to write new wavelengths into it
        minima_block=[]

#print thickness_len_pos[200]
#print list_all_minima_blocks[10]
#quit()

# counter which will be index for the data list to iterate folders and initial guesses at the same time
folder_counter = 0
# for every folder with images in the data folder do:
for data_folder in data:


    # make a list of the wavelength images in the folder
    folder_list = os.listdir(data_folder)

    # sort the list of folders
    folder_list.sort()

    # for all folders in the list do
    for folder in folder_list:

        # for the maxima --> should not be larger than lookahead_min
        lookahead_max = lookahead_min-1 

        # make wavelength list

        # This should not be changed to another step size
        wave_step = 1.0/enhance_resolution      
        
        # create an empty list which will later contain all wavelengths, e.g. 550,551,552,...,750
        waves=[]

        # write wavelengths into the list
        waves=[wave_start + i*wave_step for i in xrange(int((wave_end-wave_start)/wave_step + 1))]

        ###################
        # read image data #
        ###################

        # make and sort list of file in the folder
        files = os.listdir(data_folder+'/'+folder)
        files.sort()
        

        # check if the images should be read from a tiff stack or from single images
        if tiff_stack == False:
            # get size, bit-depth of the images
            for i in range(len(files)):
                # only consider files with the ending tiff, tif
                if files[i][-5:]=='.tiff' or files[i][-4:]=='.tif': 
                    Img=im.open(data_folder + '/'+folder + '/' + files[i])
                    break # stop the loop if an image was found


            # calculate dimensions of the images, this will be needed later on
            Image_width = Img.size[0]/binning
            Image_height = Img.size[1]/binning

            # get colour and bit depth of image, this is important to know the range of the values (e.g. 8-bit is 0-255, 16-bit is 0-65535)
            Image_mode = Img.mode 
            if Image_mode == 'RGB' or Image_mode == 'P' or Image_mode == 'L':
                Image_bit = '8'
            elif Image_mode == 'I;16' or Image_mode == 'I;12' or Image_mode=='I;16B':
                Image_bit = '16'
            else:
                print 'Image mode is:', Img.mode
                Image_bit = '16B'

            # create an empty array which will contain all image data 
            all_images = np.zeros(((wave_end-wave_start)/wave_step + 1,Image_height,Image_width),np.uint16)
           

            # read every image in folder and check if it is in the wavelength range --> write grey values into array
            
            # set a counter to check how many images have been processed
            counter=0

            print 'reading images from folder: ', folder

            # start iterating over the files
            if enhance_resolution == 1: 
                for i in xrange(len(files)):
                    # only consider files with the ending tiff, tif
                    if files[i][-5:]=='.tiff' or files[i][-4:]=='.tif':
                        # check if the current file is in the wavelength range the user specified at the beginning
                        if float(files[i][:-4]) >= wave_start and float(files[i][:-4]) <= wave_end:
                            print files[i]
                            #print counter

                            # check if its 8-bit, convert, load
                            if Image_bit == '8':
                                Img = im.open(data_folder + '/'+folder + '/' + files[i])
                                Img = Img.convert('L')
                                all_images[counter]=transform.rescale(np.asarray(Img),1.0/binning,preserve_range=True).round().astype(np.uint16)

                            # read all other formats with imread from skimage
                            else:
                                Img = data_folder + '/'+folder + '/' + files[i]
                                all_images[counter]=transform.rescale(imread(Img, as_grey=True),1.0/binning,preserve_range=True).round().astype(np.uint16)

                            counter+= 1
            else:
                sys.exit("Sorry, but you cannot use the enhanced resolution function without using a tiff stack. Please use Imagej to create a tiff stack from the single files.")

        if tiff_stack == True:
            print 'reading images from folder: ', folder

            # find filename for stack --> only on file per folder allowed
            for i in xrange(len(files)):
                if files[i][-5:]=='.tiff' or files[i][-4:]=='.tif':
                    stack_name = files[i]
                    break

            # use imread to load multipage tiff
            stack = imread(data_folder + '/'+folder + '/' + stack_name)

            # calculate dimensions of the images, this will be needed later on
            Image_width = int(np.round(len(stack[0][0])/float(binning)))
            Image_height = int(np.round(len(stack[0])/float(binning)))

            # create an empty array which will contain all image data 
            #all_images = np.zeros(((wave_end-wave_start)/wave_step + 1,Image_height,Image_width),np.uint16)

            # write image stack to numpy array
            if binning != 1:
                binned_stack = np.zeros(((wave_end-wave_start)/(wave_step*enhance_resolution) + 1,Image_height,Image_width),np.uint16)
                for i in range(len(stack[wave_start-stack_wave_start:wave_end-stack_wave_start+1])):
                    binned_stack[i]=transform.rescale(stack[i+(wave_start-stack_wave_start)],1.0/binning,preserve_range=True).round().astype(np.uint16)

                # change name of stack
                stack = binned_stack
                #all_images = sp.ndimage.uniform_filter(sp.ndimage.interpolation.zoom(all_images, [((len(all_images)-1.0)*enhance_resolution+1)/len(all_images),1,1], order = 1),(enhance_resolution**2,0,0))
            

            if enhance_resolution != 1:
                # time interpolation
                test_time_1 = time.time()
                
                # if there is binning the iteration over the raw image stack has already been done
                if binning != 1:
                    all_images = self_interpolate(stack,enhance_resolution)
                else:
                    all_images = self_interpolate(stack[wave_start-stack_wave_start:wave_end-stack_wave_start+1],enhance_resolution)

                test_time_2 = time.time()
                print "interpolation takes: ", test_time_2-test_time_1

                # more memory efficient but also much more time consuming, but necessary to avoid running out of memory
                if save_ram == True:
                    for y in xrange(len(all_images[0])):
                        #print y
                        for x in xrange(len(all_images[0][0])):
                            all_images[:,y,x] = sp.ndimage.uniform_filter1d(all_images[:,y,x],size=average_window*enhance_resolution)
                else:
                    all_images = sp.ndimage.uniform_filter1d(all_images,axis=0,size=average_window*enhance_resolution)
                print "smoothing takes: ", time.time()-test_time_2
                
            # case for no resolution enhancement
            else:
                if binning != 1:
                    all_images = stack
                else: 
                    all_images = stack[wave_start-stack_wave_start:wave_end-stack_wave_start+1]
                    
            #sys.exit("smoothing")
            #quit()
            # bit depth
            Image_bit = '16'



###############################################################################
## Find new delta (minima finding algorithm) to deal with different dynamics ##
###############################################################################
           
        # use different delta scaling for 8,12,16 bit, delta something like the peak height of the minima, it differs of course significantly for 8-bit images and 16-bit images, just because of the different range


        # the proper delta shall be a 10th of the mean value of all images
        new_delta = int(all_images.mean()/10)
        # in case delta is larger than 7, take it
        if new_delta > 7:
            delta = new_delta

        # for 16 or 12 bit images do the following
        if Image_bit == '16':
            # the proper delta shall be a 10th of the mean value of all images
            new_delta = int(all_images.mean()/10)
            # in case delta is larger than 7, take it
            if new_delta > 7:
                delta = new_delta

        # for 8-bit images just take a delta of 7 as this is sufficient
        if Image_bit == '8':
            delta = 7    

        # TEST TEST TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # since the data is smoothed there is less noise and the delta can be chosen smaller
        # does this work for all values for enhanced_resolution equally??
        if enhance_resolution != 1 and delta>(4*7):
            delta = delta/4
        
        #delta = 100

        # calculate how much the delta should be varied in case no value is found
        # this has been found empirically to lead to a good fit of the minima
        delta_vary = int(delta/5-delta/20)
        print 'delta = ', delta

#######################
## Start calculations #
#######################

        
        # start row
        start = 0
        # last row
        ende = Image_height

        print 'perform the calculations'
        # get start time for simulation to check later how long it took
        t1 = time.time()

        if one_minimum_fit == True:
            print "Using one minimum algorithm" 
            if use_thickness_limits == True:
            
                # call the external one minimum cython/c++ function with all the parameters
                result = Fit_one.c_Fit_Pixel(start,ende,all_images, thickness_len_pos, waves, tolerance, lookahead_min*enhance_resolution, lookahead_max*enhance_resolution, delta,delta_vary,list_all_minima_blocks, use_thickness_limits, thickness_limit,area_avrg,init_guess[folder_counter])
            else:
                print "\nERROR: You have to set: use_thickness_limits == True\n"
                quit()
        elif plot_error == True:
            print "An error map will be plotted"
            error_map_path = data_folder + '/' + folder 
            result = Fit_error.c_Fit_Pixel(start,ende,all_images, thickness_len_pos, waves, tolerance, lookahead_min*enhance_resolution, lookahead_max*enhance_resolution, delta,delta_vary,list_all_minima_blocks, use_thickness_limits, thickness_limit,area_avrg,error_map_path)
        else:
            print "Standard calculation"
            # call the external cython/c++ function with all the parameters
            result = Fit.c_Fit_Pixel(start,ende,all_images, thickness_len_pos, waves, tolerance, lookahead_min*enhance_resolution, lookahead_max*enhance_resolution, delta,delta_vary,list_all_minima_blocks, use_thickness_limits, thickness_limit,area_avrg)#[0]
       
        t2 = time.time()



        print t2-t1, 'seconds just for the calculation'



        ###########################
        # count not fitted values #
        ###########################

        not_fitted = Image_height*Image_width - np.count_nonzero(result)
        not_fitted_percent = 100.0/(Image_height*Image_width)*not_fitted
        print 'not fitted values',not_fitted
        print 'in percent:', not_fitted_percent

        ######################
        # Write data to file #
        ######################

        print 'write data to file'
        
        # generate a header with all parameters
        HEADER = time.strftime('Version = ' + version + '\n' + "%d.%m.%Y at %H:%M:%S")+'\n' + 'folder with data = ' + folder + '\n' + 'simulation file = ' + sim_file + '\n' + 'wave_start = '+str(wave_start) + '\n' + 'wave_end = ' + str(wave_end) + '\n' + 'binning = ' + str(binning) + '\n' +  'lookahead_min = ' + str(lookahead_min) + '\n'  + 'lookahead_max = ' + str(lookahead_max) + '\n' + 'delta = ' + str(delta) + ' delta was varied +-'+str(delta_vary*5)+ '\n' + 'tolerance = ' + str(tolerance) + '\n' + 'enhance_resolution = ' + str(enhance_resolution) + '\n' +'average_window = '+ str(average_window) + '\n' + 'one_minimum_fit = ' + str(one_minimum_fit) + '\n' +'init_guess = ' + str(init_guess[folder_counter]) + '\n' +'thickness limits used: ' + str(use_thickness_limits) + '\n' + 'thickness limits: ' + str(thickness_limit) + '\n' + 'area_avrg = ' +str(area_avrg) + '\n' + 'not fitted values: ' + str(not_fitted) + ', percentage of whole image: ' + str(not_fitted_percent)  + '\n'

        HEADER+= '\n'

        # generate a useful filename
        file_name = data_folder + '_' + folder + time.strftime("_%Y%m%d_%H%M%S")+'.txt'

        # better way, no need to replace stuff afterwards:
        #convert result to float to insert nan
        result = result.astype(np.float)
        # replace "0" with "nan"
        result[result==0] = np.nan
        # use numpy function to save array to file
        # The Header is currently not written into the file, just uncomment the last part and remove one of the brackets
        np.savetxt(data_folder + '/' + file_name,result,fmt='%0.0f')#,header=HEADER )


        print (time.time()-t_a_start), ' seconds for the whole program'


        ##############
        # plot datta #
        ##############

        # parameters for printing
        # color map is calculated like (mean_thickness - color_min, mean_thickness + color_max) 

        # color_min = 500
        # color_max = 500
        # # make a new figure
        # plt.figure(folder)
        # # create plot of the results
        # plt.imshow(result)
        # # set the color scale to the limits provided
        # plt.clim(np.nanmean(result)-color_min,np.nanmean(result)+color_max)
        # # plot a color bar
        # plt.colorbar()

        # remove "#" to show the plot after the calculation
        #plt.show()


    # write parameter file instead of header

    parameter_file = data_folder +'/' + "Calculation_Parameters_" + data_folder + ".txt"
    p = open(parameter_file,'w')

    p.write(HEADER)
    p.close()

    # increase folder_counter
    folder_counter += 1
