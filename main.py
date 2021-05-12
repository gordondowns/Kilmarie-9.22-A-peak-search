import numpy as np
import AmorphousProfileCharacterization
from glob import glob
import pandas as pd
from time import time

image_directory = "images_and_CSVs/"

# Prepare IMA lists to get mineral chemistry
IMA_list = pd.read_csv('IMA_mineral_list_with_chemistries.csv')
IMA_mineral_names = IMA_list['Mineral Name (plain)']
IMA_mineral_names = list(IMA_mineral_names)
IMA_chemistries = IMA_list['IMA Chemistry (concise)']
IMA_chemistries = list(IMA_chemistries)

# Get list of difdata. Difdata are essentially files that contain preprocessed CIF data.
dif_directory_name = 'possible_mars_difs/'
dif_file_paths = glob(dif_directory_name+'*')
dif_file_names = [file_name.split('\\')[-1][:-4] for file_name in dif_file_paths]

# Initialize loop values
t0 = time()
results = []
live_list = open('temp/_live_list_'+str(t0).split('.')[0]+'.csv','w')

for d in range(len(dif_file_paths)):
    if d % 100 == 0:
        print(d,'/',len(dif_file_names))
    try:
        # Instantiate the xrd analysis object
        # This pattern is offset-corrected Kilmarie (sub 113 microns), ignoring data below 10 degrees two theta
        apc = AmorphousProfileCharacterization.TopLevel("phases/2389_KM32_EDA2389_01-60_nSsub113rvm.img150501.dif",twotheta_ranges=[(10.0,52.0)],print_warnings=False)
        apc.AddDif(dif_file_paths[d],refine_cell_parameters=False)
        
        # Check for a d-spacing between 8-11 A with intensity at least 20% of the highest peak's intensity.
        # If there isn't a d-spacing in the target region with sufficiently high intensity, discard this phase.
        skip_file = True
        for i in range(len(apc.difs[0].d_spacings)):
            if apc.difs[0].d_spacings[i] > 8.0 and apc.difs[0].d_spacings[i] < 11.0 and apc.difs[0].peak_intensities[i] > 20.0:
                skip_file = False
                break
        if skip_file:
            continue
        
        apc.AddDif('phases/difs/Andesine0001052.txt',refine_scaling=False,refine_cell_parameters=False,refine_FWHM=False,scaling=0.011057592067424713,FWHM=0.3729429121706677,FWHM_bounds=(0.02, 20.0),cell_parameters=[8.144720797464297, 12.886996967587267, 7.119394030276795, 93.3540045435841, 116.27707405953119, 90.11922703272907],cell_parameter_bounds=[(7.3611, 8.996900000000002), (11.592, 14.168000000000003), (6.4008, 7.823200000000001), (84.096, 102.784), (104.589, 127.831), (81.20700000000001, 99.25300000000001)])
        apc.AddDif('phases/difs/Anhydrite0005117.txt',refine_scaling=False,refine_cell_parameters=False,refine_FWHM=False,scaling=0.032020587527150424,FWHM=0.4732777939596451,FWHM_bounds=(0.02, 20.0),cell_parameters=[7.0011358057016055, 6.989218056554633, 6.2392146609485195, 90.0, 90.0, 90.0],cell_parameter_bounds=[(6.2937, 7.692300000000001), (6.2955000000000005, 7.694500000000001), (5.6205, 6.8695), (81.0, 99.00000000000001), (81.0, 99.00000000000001), (81.0, 99.00000000000001)])
        apc.AddDif('phases/difs/Bassanite0013868.txt',refine_scaling=False,refine_cell_parameters=False,refine_FWHM=False,scaling=0.00651524020708094,FWHM=0.5119682944694989,FWHM_bounds=(0.02, 20.0),cell_parameters=[12.03645329449536, 6.967468308599716, 12.756734885011403, 90.0, 90.27482675300189, 90.0],cell_parameter_bounds=[(10.82853, 13.234870000000003), (6.23421, 7.6195900000000005), (11.40408, 13.938320000000003), (81.0, 99.00000000000001), (81.243, 99.297), (81.0, 99.00000000000001)])
        apc.AddDif('phases/difs/Enstatite0010680.txt',refine_scaling=False,refine_cell_parameters=False,refine_FWHM=False,scaling=0.004531047068740678,FWHM=0.45706684789515356,FWHM_bounds=(0.02, 20.0),cell_parameters=[18.429727388551076, 8.790084619628933, 5.222134964943968, 90.0, 90.0, 90.0],cell_parameter_bounds=[(16.389000000000003, 20.031000000000002), (7.9308, 9.693200000000001), (4.6602, 5.6958), (81.0, 99.00000000000001), (81.0, 99.00000000000001), (81.0, 99.00000000000001)])
        apc.AddDif('phases/difs/Siderite0017561.txt',refine_scaling=False,refine_cell_parameters=False,refine_FWHM=False,scaling=0.022444076715808906,FWHM=0.44256560173412995,FWHM_bounds=(0.02, 20.0),cell_parameters=[4.693052058491399, 4.693052058491399, 15.342228424754369, 90.0, 90.0, 120.0],cell_parameter_bounds=[(4.2246, 5.1634), (4.2246, 5.1634), (13.887, 16.973000000000003), (81.0, 99.00000000000001), (81.0, 99.00000000000001), (108.0, 132.0)])
        apc.AddDif('phases/difs/Sanidine0010740.txt',refine_scaling=False,refine_cell_parameters=False,refine_FWHM=False,scaling=0.010774060955917486,FWHM=0.7336085257053988,FWHM_bounds=(0.02, 20.0),cell_parameters=[8.670700390352078, 12.983396190567678, 7.112458785565398, 90.0, 116.00347703783208, 90.0],cell_parameter_bounds=[(7.6941, 9.4039), (11.725200000000001, 14.330800000000002), (6.4692, 7.9068000000000005), (81.0, 99.00000000000001), (104.41799999999999, 127.622), (81.0, 99.00000000000001)])
        apc.AddDif('phases/difs/goethite.dif',refine_scaling=False,refine_cell_parameters=False,refine_FWHM=False,scaling=0.0026325800623211272,FWHM=0.3866085922977733,FWHM_bounds=(0.02, 20.0),cell_parameters=[4.561298557804432, 9.953572396655263, 3.0254571430105157, 90.0, 90.0, 90.0],cell_parameter_bounds=[(4.13811, 5.057690000000001), (8.955900000000002, 10.946100000000001), (2.71602, 3.31958), (81.0, 99.00000000000001), (81.0, 99.00000000000001), (81.0, 99.00000000000001)])
        apc.AddDif('phases/difs/Hematite0000143.txt',refine_scaling=False,refine_cell_parameters=False,refine_FWHM=False,scaling=0.0032790016997333964,FWHM=0.45801522501103553,FWHM_bounds=(0.02, 20.0),cell_parameters=[5.011502715962327, 5.011502715962327, 13.813392230140012, 90.0, 90.0, 120.0],cell_parameter_bounds=[(4.5342, 5.541800000000001), (4.5342, 5.541800000000001), (12.3948, 15.149200000000002), (81.0, 99.00000000000001), (81.0, 99.00000000000001), (108.0, 132.0)])
        apc.AddDif('phases/difs/Quartz0006362.txt',refine_scaling=False,refine_cell_parameters=False,refine_FWHM=False,scaling=0.005023536732570089,FWHM=0.5413841877233957,FWHM_bounds=(0.02, 20.0),cell_parameters=[5.072361082702881, 5.072361082702881, 5.473712866793275, 90.0, 90.0, 120.0],cell_parameter_bounds=[(4.4223300000000005, 5.405070000000001), (4.4223300000000005, 5.405070000000001), (4.86423, 5.945170000000001), (81.0, 99.00000000000001), (81.0, 99.00000000000001), (108.0, 132.0)])
        apc.AddDif('phases/difs/Cristobalite0001629.txt',refine_scaling=False,refine_cell_parameters=False,refine_FWHM=False,scaling=0.004819069312631354,FWHM=0.4149135698249561,FWHM_bounds=(0.02, 20.0),cell_parameters=[5.0010706354413434, 5.0010706354413434, 6.92388771645422, 90.0, 90.0, 90.0],cell_parameter_bounds=[(4.474530000000001, 5.468870000000001), (4.474530000000001, 5.468870000000001), (6.23007, 7.61453), (81.0, 99.00000000000001), (81.0, 99.00000000000001), (81.0, 99.00000000000001)])
        apc.AddDif('phases/difs/Magnetite0000945.txt',refine_scaling=False,refine_cell_parameters=False,refine_FWHM=False,scaling=0.0,FWHM=0.3758482507169937,FWHM_bounds=(0.02, 20.0))
        apc.AddProfile('phases/profiles/10A_Bent_Otay_rvm1.mdi',refine_scaling=False,scaling=0.03716610456530104)
        apc.AddProfile('phases/profiles/10A_Sapn_GriffithPark1_rvm1.mdi',refine_scaling=False,scaling=0.023107937756938715)
        apc.AddProfile('phases/profiles/10A_Nont_Garfield_1_rvm1.mdi',refine_scaling=False,scaling=0.018548501285925414)
        apc.AddProfile('phases/profiles/Opal-A_AUSOPL2_LT0150_RIR1_500Frame_AIR-ka20171112.mdi',refine_scaling=False,scaling=0.010541327388249021)
        apc.AddProfile('phases/profiles/Ferrihydrite-TNA31-230-ka20150324-edit_sub-mylar-bkg.txt',refine_scaling=False,scaling=0.02030341373275215)
        apc.AddProfile('phases/profiles/13A_591_EmptyCell.mdi',refine_scaling=False,scaling=0.0041664041744977455)
        apc.SetBackgroundGuess([4.067933831846612])
        apc.refine_background = False
        # Norm for above fit is 13.704499270541868

        # First run: fit scaling and FWHM of new phase and freeze everything else 
        apc.DoOptimization(print_start_and_finish=False)

        # If error was not reduced by more than 0.0045, discard this phase
        if apc.GetNorm() > 13.7:
            continue

        # Second run: fit scaling, FWHM, and cell parameters of new phase and freeze everything else
        apc.difs[0].refine_cell_parameters = True
        apc.DoOptimization(print_start_and_finish=False)

        # Third run: full refinement
        apc.refine_background = True
        for dif in apc.difs:
            dif.refine_scaling = True
            dif.refine_FWHM = True
            dif.refine_cell_parameters = True
        for profile in apc.profiles:
            profile.refine_scaling = True
        apc.DoOptimization(print_start_and_finish=False)

        # Plot results
        apc.Plot(image_file_path=image_directory+dif_file_names[d]+'.png')
        # Save CSV of computed patterns
        apc.OutputCSV(image_directory+dif_file_names[d]+'.csv')

        # Prepare results string
        cur_mineral_name = dif_file_names[d][:-7]
        cur_chemistry = IMA_chemistries[IMA_mineral_names.index(cur_mineral_name)].replace(',',';')
        cur_cell_parameters = str(apc.difs[0].cell_parameters_optimized).replace(',','')
        cur_cell_parameters_delta = str(np.array(apc.difs[0].cell_parameters_optimized) - np.array(apc.difs[0].cell_parameters_guess)).replace(',','')
        r = str(apc.GetNorm()) +','+ dif_file_names[d] +','+ cur_chemistry +','+ cur_cell_parameters_delta +','+ cur_cell_parameters \
            +','+ str(apc.difs[0].scaling_optimized) +','+str(apc.difs[0].FWHM_optimized)
        
        # Print result to console and write to temporary results file
        print(r)
        live_list.write(r+'\n')
        live_list.flush()
        results.append(r)

    except:
        print('error with file:',dif_file_names[d])
        continue

live_list.close()

# Print elapsed time
print('-----')
minutes = (time() - t0)/60.0
hours,minutes = divmod(minutes,60)
print(hours,'hours,',minutes,'minutes')

# Sort results by norm and write to final results file
results.sort()
f = open('search_results.csv','w')
f.write('Sum of squares error,Mineral name and AMCSD code,Chemistry,Change in cell parameters,Original cell parameters,Scaling,FWHM\n')
f.write("\n".join(results))
f.close()

# After running this file, interpret the results by hand:
# First, look at search_results.csv.
# Discard phases if their cell parameters had to vary too much to produce a good fit
#  (a good place to find the acceptable range of cell parameters for a given mineral is https://rruff.info/ima/).
# Discard phases if they have too much of an element not commonly found on Mars.
# For the remaining phases, look at their plots or CSVs.
# Discard phases that did not fit the 9.22 A peak well.
# Discard phases that have significant peaks in places where the Kilmarie pattern does not.