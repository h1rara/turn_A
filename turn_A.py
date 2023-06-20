import os
from astropy.table import Table, vstack
from astropy.io import fits
import astropy.units as u
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

from astropy import wcs
import warnings
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy



for i in range(18):

    hdu_1 = fits.open('deg_snr_fb'+str(i)+'_1_clusters_asgn.fits')[0]
    a = hdu_1.data
    hd = hdu_1.header


    data1 = fits.getdata('deg_snr_fb'+str(i)+'_1_catalog.fits', 1)
    cat= Table(data1)

    if i == 17:

    #identify high and low latitude clouds

        lo_edge = a[:,0:2,:]
        lo_edge_list = list(np.unique(lo_edge))

        up_edge = a[:,a.shape[1]-2:a.shape[1],:]
        up_edge_list = list(np.unique(up_edge))


    #identify left edge

        l_edge = a[:,:,0:2]
        l_edge_list = list(np.unique(l_edge))


    #identify right edge

        r_edge = a[:,:,a.shape[2]-2:a.shape[2]]
        r_edge_list = list(np.unique(r_edge))


    #full list to remove

        remove_list = l_edge_list + r_edge_list + lo_edge_list + up_edge_list


    #fix catalog

        for g in range(len(cat)-1, -1, -1):
            if cat['_idx'][g] in remove_list:
                cat.remove_row(g)

        cat.write('clean_cat_region'+str(i)+'.fits', format='fits', overwrite=True)


    #fix cube

        hdu_1 = fits.open('deg_snr_fb17_1_clusters_asgn.fits')[0]
        a = hdu_1.data
        hd = hdu_1.header


        for k in remove_list:

            a[a==k] = -1.0

        fits.writeto('clean_cohrs_17.fits',a, hd, overwrite=True)



    #for all others

    else:

    #identify high and low latitude clouds

        lo_edge = a[:,0:2,:]
        lo_edge_list = list(np.unique(lo_edge))

        up_edge = a[:,a.shape[1]-2:a.shape[1],:]
        up_edge_list = list(np.unique(up_edge))


    #everything that toucehs the middle line

        p_A = int(a.shape[2]/2)
        A = a[:,:, p_A]
        list_A = list(np.unique(A))

    #every thing between the left edge and the middle line

        AB = a[:,:, 0:p_A]
        list_AB = list(np.unique(AB))

    #everything in AB that does not touches the middle line

        left_block_list = list(set(list_AB) - set(list_A))


    #identify right edge

        r_edge = a[:,:,a.shape[2]-2:a.shape[2]]
        r_edge_list = list(np.unique(r_edge))

    #list to remove

        remove_list = left_block_list + r_edge_list + lo_edge_list + up_edge_list


        #fix catalog

        for g in range(len(cat)-1, -1, -1):
            if cat['_idx'][g] in remove_list:
                cat.remove_row(g)

        cat.write('clean_cat_region'+str(i)+'.fits', format='fits', overwrite=True)


    #fix cube

        hdu_1 = fits.open('deg_snr_fb'+str(i)+'_1_clusters_asgn.fits')[0]
        a = hdu_1.data
        hd = hdu_1.header


        for k in remove_list:

            a[a==k] = -1.0

        fits.writeto('clean_cohrs_'+str(i)+'.fits',a, hd, overwrite=True)


















###join catalogues


data1 = fits.getdata('clean_cat_region0.fits', 1)
cat_0= Table(data1)
cat_0['REGION'] = [0]*len(cat_0)

data1 = fits.getdata('clean_cat_region1.fits', 1)
cat_1= Table(data1)
cat_1['REGION'] = [1]*len(cat_1)


data1 = fits.getdata('clean_cat_region2.fits', 1)
cat_2= Table(data1)
cat_2['REGION'] = [2]*len(cat_2)

data1 = fits.getdata('clean_cat_region3.fits', 1)
cat_3= Table(data1)
cat_3['REGION'] = [3]*len(cat_3)


data1 = fits.getdata('clean_cat_region4.fits', 1)
cat_4= Table(data1)
cat_4['REGION'] = [4]*len(cat_4)


data1 = fits.getdata('clean_cat_region5.fits', 1)
cat_5= Table(data1)
cat_5['REGION'] = [5]*len(cat_5)


data1 = fits.getdata('clean_cat_region6.fits', 1)
cat_6= Table(data1)
cat_6['REGION'] = [6]*len(cat_6)


data1 = fits.getdata('clean_cat_region7.fits', 1)
cat_7= Table(data1)
cat_7['REGION'] = [7]*len(cat_7)


data1 = fits.getdata('clean_cat_region8.fits', 1)
cat_8= Table(data1)
cat_8['REGION'] = [8]*len(cat_8)


data1 = fits.getdata('clean_cat_region9.fits', 1)
cat_9= Table(data1)
cat_9['REGION'] = [9]*len(cat_9)


data1 = fits.getdata('clean_cat_region10.fits', 1)
cat_10= Table(data1)
cat_10['REGION'] = [10]*len(cat_10)

data1 = fits.getdata('clean_cat_region11.fits', 1)
cat_11= Table(data1)
cat_11['REGION'] = [11]*len(cat_11)

data1 = fits.getdata('clean_cat_region12.fits', 1)
cat_12= Table(data1)
cat_12['REGION'] = [12]*len(cat_12)

data1 = fits.getdata('clean_cat_region13.fits', 1)
cat_13= Table(data1)
cat_13['REGION'] = [13]*len(cat_13)

data1 = fits.getdata('clean_cat_region14.fits', 1)
cat_14= Table(data1)
cat_14['REGION'] = [14]*len(cat_14)

data1 = fits.getdata('clean_cat_region15.fits', 1)
cat_15= Table(data1)
cat_15['REGION'] = [15]*len(cat_15)

data1 = fits.getdata('clean_cat_region16.fits', 1)
cat_16= Table(data1)
cat_16['REGION'] = [16]*len(cat_16)


data1 = fits.getdata('clean_cat_region17.fits', 1)
cat_17= Table(data1)
cat_17['REGION'] = [17]*len(cat_17)




t_list = (cat_1, cat_2, cat_3, cat_4, cat_5, cat_6, cat_7, cat_8, cat_9, cat_10, cat_11, cat_12,  cat_13, cat_14, cat_15, cat_16, cat_17)



T = cat_0
for i in (range(17)):
    #print(t_list[i])
    T = vstack([T, t_list[i]])


galid= []
for i in range(len(T)):
    galid.append(str(T['REGION'][i]) +'_' + str(T['_idx'][i]))

T['gal_id'] = galid


T.write('full_cohrs_cat.fits', format='fits', overwrite = True)












