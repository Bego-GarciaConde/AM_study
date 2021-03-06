
import numpy as np
        ####################################
        #             Settings             #
        ####################################
 

r_max = 5
delta_vir =         333   # WHICH OVERDENSITY ?

#angular_momentum_ref = np.array([-5.718599521275628e+29,-1.0252942573706002e+29,2.570348003882171e+29])
#angular_momentum_ref =np.array([-1.850951505301476e+29,-3.150082206215513e+28,8.734402694190449e+28])
#angular_momentum_ref =np.array([-2.740155533109571e+29,-4.545889508503972e+28,1.1010570705289006e+29])
angular_momentum_ref =np.array([-5.7202837196079984e+29,-1.0243667905866969e+29,2.569850545344377e+29])

#angular_momentum_ref_dm =np.array([-3.5550889734836453e+28,3.147127171767552e+27,5.100549273333895e+28])
angular_momentum_ref_dm =np.array([-6.8751390926297565e+28,1.2874646084830964e+28,5.7977225441833435e+28])
#angular_momentum_ref_dm =np.array([-2.8516996436835407e+28,-4.6817751385933117e+26,4.340760159051236e+28])
#angular_momentum_ref = np.array([-1.63535369e+26, -1.90795637e+27, -9.11225245e+27])
#angular_momentum_ref = np.array([-1.6337010768551683e+28,-1.137774486844715e+26,2.782572704217265e+28])
        ####################################
        #      Data and saving paths      #
        ####################################


#path_guardado_snapshots = "/media/temp/bego/snapshots_resim/" #path of csv saving
path_csv = "/media/temp/bego/snapshots_resim/"
path_save_data = "results/"    #Where to save/load csv?
path_datos = "/home/bego/GARROTXA/datos_GARROTXA_resim/"    #Where to save/load centers and Rvir?
path_save_angular_momentum_figs = "/media/temp2/bego/darkmatter_2Rvir/Angular_momentum_figs/"



        ####################################
        #      Snapshots to analyse        #
        ####################################

#label of the snapshots a=0.xxx
snapshots_analysis = [520,523,525, 527,530,532,535, 537,539,541,
543, 545,547, 550, 553, 555,557, 560, 563, 565, 567,570,573, 577, 580,
583, 585,587,590, 592,594,596,598,600,
602,604, 608, 610, 612, 614, 616, 618, 620, 622, 624, 626, 
629, 630, 632, 634, 636, 639, 640, 642, 644, 646, 648, 650, 652, 654, 656, 658, 660, 662, 
664, 666, 668,670, 672, 674, 676, 679, 681, 682, 684, 687, 689,
690, 692, 694, 698, 704, 706, 708,711, 712,714, 716, 718, 720, 
722, 724, 726, 728, 731, 732, 734, 736, 739, 740, 742, 744, 746, 748, 751,752,
 755, 756, 758, 761,763, 764, 766, 768, 770, 772, 774, 776, 778,780, 
782, 784, 786, 788, 790, 792, 794, 797, 798, 802, 805, 806, 808, 810, 812, 814, 816,
 818, 820, 822, 824, 826, 828, 
830, 832, 834, 836, 839, 840, 842, 844, 846, 848, 850,
853, 855, 856, 858, 860, 862, 864, 867, 870, 872, 875, 877, 879, 881, 883, 884, 886, 888,
890, 892, 894, 898, 900, 902, 904, 907, 908, 910, 912, 915, 
916, 918, 921, 922, 924, 927, 929, 
930, 932, 934, 937, 939, 941,942, 944, 946, 948, 950, 952, 954,956, 
958, 961, 963, 965, 966, 968, 970, 972, 974, 976, 979,
 980, 982, 984, 989, 990, 993, 994, 996, 999] 
#snapshots_analysis = [610]