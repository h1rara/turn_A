import os
from astropy.table import Table, vstack, unique, hstack
from astropy.io import fits, ascii
import astropy.units as u
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import statistics
from astropy import wcs
import warnings
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy
import copy
import math
import cmath
import random



os.chdir('/home/siwa/data/COHRS/')

#prepare catalogues

for i in range(18):


    hdu_1 = fits.open('deg_snr_fb'+str(i)+'_3_clusters_asgn.fits')[0]
    a = hdu_1.data
    hd = hdu_1.header

  #  hdu_1 = fits.open('deg_snr_fb'+str(i)+'_A_leaf_asgn.fits')[0]
   # l = hdu_1.data
    #hd = hdu_1.header

    data1 = fits.getdata('deg_snr_fb'+str(i)+'_3_catalog.fits', 1)
    cat= Table(data1)

    b = np.unique(a)

    for j in range(len(cat)-1, -1, -1):
        print(i, j)

        if cat['_idx'][j] not in b:
            cat.remove_row(j)

    cat.write('ready_cat_region_'+str(i)+'.fits', format='fits', overwrite=True)



#remove clouds touching the edges of the field of observation

for i in range(18):

    hdu_1 = fits.open('deg_snr_fb'+str(i)+'_3_clusters_asgn.fits')[0]
    a = hdu_1.data
    hd = hdu_1.header

  #  hdu_1 = fits.open('deg_snr_fb'+str(i)+'_A_leaf_asgn.fits')[0]
   # l = hdu_1.data
    #hd = hdu_1.header



    data1 = fits.getdata('ready_cat_region_'+str(i)+'.fits', 1)
    cat= Table(data1)

    if i == 16:

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
            if float(cat['_idx'][g])in remove_list:
                cat.remove_row(g)

        cat.write('clean_cat_region_'+str(i)+'.fits', format='fits', overwrite=True)


    #fix cube


        for k in remove_list:

            a[a==k] = -1.0

        fits.writeto('clean_cohrs_'+str(i)+'.fits',a, hd, overwrite=True)


    #for the last region remove every thing on the right instead of the left, since it is alreasdy contained in the previous region

    elif i == 17:

    #identify high and low latitude clouds

        lo_edge = a[:,0:2,:]
        lo_edge_list = list(np.unique(lo_edge))

        up_edge = a[:,a.shape[1]-2:a.shape[1],:]
        up_edge_list = list(np.unique(up_edge))


    #everything that toucehs the middle line

        p_A = int(a.shape[2]/2)
        A = a[:,:, p_A]
        list_A = list(np.unique(A))

    #every thing between the right edge and the middle line

        AB = a[:,:, p_A:a.shape[2]]
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

        cat.write('clean_cat_region_'+str(i)+'.fits', format='fits', overwrite=True)


    #fix cube

        for k in remove_list:

            a[a==k] = -1.0

        fits.writeto('clean_cohrs_'+str(i)+'.fits',a, hd, overwrite=True)


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

        cat.write('clean_cat_region_'+str(i)+'.fits', format='fits', overwrite=True)


        #fix cube

        for k in remove_list:

            a[a==k] = -1.0

        fits.writeto('clean_cohrs_'+str(i)+'.fits',a, hd, overwrite=True)



#####OK!




#remove snr < 10 and add dimensions

for i in range(18):


data = fits.getdata('asu.fits', 1)
ag_table = Table(data)
mask = (ag_table['_Glon']>27.0)
ag_table = ag_table[mask]
mask = (ag_table['_Glon']<47.0)
ag_table = ag_table[mask]

    print(i)

    #open catalog
    data1 = fits.getdata('clean_cat_region_'+str(i)+'.fits', 1)
    cat= Table(data1)


    #open cluster assignments
    hdu_1 = fits.open('clean_cohrs_'+str(i)+'.fits')[0]
    a = hdu_1.data
    hd = hdu_1.header

    #open fb map
    hdu_1 = fits.open('deg_snr_fb'+str(i)+'.fits')[0]
    fb = hdu_1.data
    hd_fb = hdu_1.header



data = fits.getdata('asu.fits', 1)
ag_table = Table(data)
mask = (ag_table['_Glon']>27.0)
ag_table = ag_table[mask]
mask = (ag_table['_Glon']<47.0)
ag_table = ag_table[mask]

    #fix nan
    #a[np.isnan(a)]= 0.0
    fb[np.isnan(fb)]= 0.0

    ll = []
    lb = []
    lv = []


    #for each entry in the region catalog
    for j in range(len(cat)-1, -1, -1):

        print(i, j, cat['_idx'][j])

        #isolate each cloud
        ix = float(cat['_idx'][j])

        b = copy.deepcopy(a)
        b[b!=ix]=-1.0
        b[b==ix]=1.0
        b[b<0.0]=0.0


       # c = np.max(em, axis = 0)
       # plt.imshow(c)
       # plt.show()

        #emission signature
        em = np.multiply(b, fb)
        print(np.max(em))


        if np.max(em) < 10.0:


            #remove entry from cat
            cat.remove_row(j)
            #remove cloud from cube
            a[a==ix]=-1.0

        else:

            #find dimensions
            voxels = np.where(b==1)

            if len(voxels[2]) == 0:
                lmin = 0
                lmax = 0

            else:
                lmin = min(voxels[2])
                lmax = max(voxels[2])


            if len(voxels[1]) == 0:
                bmin = 0
                bmax = 0

            else:
                bmin = min(voxels[1])
                bmax = max(voxels[1])

            if len(voxels[0]) == 0:
                vmin = 0
                vmax = 0

            else:
                vmin = min(voxels[0])
                vmax = max(voxels[0])


            l = float(lmax - lmin)
            b = float(bmax - bmin)
            v = float(vmax - vmin)


            ll.append(l)
            lb.append(b)
            lv.append(v)

    #add columns
    cat['l_size'] = ll
    cat['b_size'] = lb
    cat['v_size'] = lv


    #write new cube and catalog
    cat.write('final_cat_region'+str(i)+'.fits', format='fits', overwrite=True)
    fits.writeto('final_cohrs_'+str(i)+'.fits',a, hd, overwrite=True)






#join catalogs

data1 = fits.getdata('final_cat_region0.fits',

data = fits.getdata('asu.fits', 1)
ag_table = Table(data)
mask = (ag_table['_Glon']>27.0)
ag_table = ag_table[mask]
mask = (ag_table['_Glon']<47.0)
ag_table = ag_table[mask]
1)
cat_0= Table(data1)
cat_0['REGION'] = [0]*len(cat_0)

data1 = fits.getdata('final_cat_region1.fits', 1)
cat_1= Table(data1)
cat_1['REGION'] = [1]*len(cat_1)


data1 = fits.getdata('final_cat_region2.fits', 1)
cat_2= Table(data1)
cat_2['REGION'] = [2]*len(cat_2)

data1 = fits.getdata('final_cat_region3.fits', 1)
cat_3= Table(data1)
cat_3['REGION'] = [3]*len(cat_3)


data1 = fits.getdata('final_cat_region4.fits', 1)
cat_4= Table(data1)
cat_4['REGION'] = [4]*len(cat_4)


data1 = fits.getdata('final_cat_region5.fits', 1)
cat_5= Table(data1)
cat_5['REGION'] = [5]*len(cat_5)


data1 = fits.getdata('final_cat_region6.fits', 1)
cat_6= Table(data1)
cat_6['REGION'] = [6]*len(cat_6)


data1 = fits.getdata('final_cat_region7.fits', 1)
cat_7= Table(data1)
cat_7['REGION'] = [7]*len(cat_7)


data1 = fits.getdata('final_cat_region8.fits', 1)
cat_8= Table(data1)
cat_8['REGION'] = [8]*len(cat_8)


data1 = fits.getdata('final_cat_region9.fits', 1)
cat_9= Table(data1)
cat_9['REGION'] = [9]*len(cat_9)
https://www.youtube.com/watch?v=DieG5aetBL0

data1 = fits.getdata('final_cat_region10.fits', 1)
cat_10= Table(data1)
cat_10['REGION'] = [10]*len(cat_10)

data1 = fits.getdata('final_cat_region11.fits', 1)
cat_11= Table(data1)
cat_11['REGION'] = [11]*len(cat_11)

data1 = fits.getdata('final_cat_region12.fits', 1)
cat_12= Table(data1)
cat_12['REGION'] = [12]*len(cat_12)

data1 = fits.getdata('final_cat_region13.fits', 1)
cat_13= Table(data1)
cat_13['REGION'] = [13]*len(cat_13)

data1 = fits.getdata('final_cat_region14.fits', 1)
cat_14= Table(data1)
cat_14['REGION'] = [14]*len(cat_14)

data1 = fits.getdata('final_cat_region15.fits', 1)
cat_15= Table(data1)
cat_15['REGION'] = [15]*len(cat_15)

data1 = fits.getdata('final_cat_region16.fits', 1)
cat_16= Table(data1)
cat_16['REGION'] = [16]*len(cat_16)

T.write('full_cohrs_cat.fits', format='fits', overwrite = True)

data1 = fits.getdata('final_cat_region17.fits', 1)
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

gid = []

for i in range(len(T)):
    gid.append(str(T['gal_id'][i]))

T['gal_id'] = gid


T.write('full_cohrs_cat.fits', format='fits', overwrite = True)

https://www.youtube.com/watch?v=DieG5aetBL0


####remove clouds that are smaller than 4x4x2

#open catalog
data1 = fits.getdata('full_cohrs_cat.fits', 1)
cat= Table(data1)

for i in range(len(cat)-1, -1, -1):
    if cat['l_size'][i] < 4 :
        cat.remove_row(i)



for i in range(len(cat)-1, -1, -1):
    if cat['b_size'][i] < 4 :
        cat.remove_row(i)


for i in range(len(cat)-1, -1, -1):
    if cat['v_size'][i] < 2 :
        cat.remove_row(i)


cat.write('full_cohrs_cat.fits', format='fits', overwrite = True)






#make clean cubes

data1 = fits.getdata('full_cohrs_cat.fits', 1)
cat= Table(data1)

for i in range(18):
    #creat a list of clouds in a given region

    mask = (cat['REGION']==i)
    cat2 = cat[mask]

    clouds_to_keep = list(cat2['_idx'])

    #open the corresponding cube

    hdu = fits.open('final_cohrs_'+str(i)+'.fits')[0]
    a = hdu.data
    hd = hdu.header

    #clouds in cube
    clouds_in_cube = list(np.unique(a))

    #remove the clouds that are not in clouds_to_keep


    for j in clouds_in_cube:

            if j not in clouds_to_keep:
                a[a==j]=-1.0

    #save cube
    fits.writeto('reduced_cohrs_'+str(i)+'.fits',a, hd, overwrite=True)






##functions


def get_label (image, hd, l, b, v):
    'Given an image and a set of galactic coordinates, it returns the value at the voxel containing those coordinates'



#we need to convert coordinated to pixels

# Read in the header keywords
#v
    crpix3=hd["CRPIX3"]
    crval3=hd["CRVAL3"]
    cdelt3=hd["CDELT3"]
    naxis3=hd["NAXIS3"]

#l
    crpix1=hd["CRPIX1"]
    crval1=hd["CRVAL1"]
    cdelt1=hd["CDELT1"]
    naxis1=hd["NAXIS1"]

#b
    crpix2=hd["CRPIX2"]
    crval2=hd["CRVAL2"]
    cdelt2=hd["CDELT2"]
    naxis2=hd["NAXIS2"]


#add velocity in grid coordinates

   # vpix=int(((v-crval3)/cdelt3)+(crpix3-1))

    x=[None]*naxis3

    for i in range(0, naxis3):
                x[i]=((i-(crpix3-1))*cdelt3+crval3)

    X=Table()

    X['x']=x

    xl=[None]*naxis1

    for i in range(0, naxis1):
                xl[i]=((i-(crpix1-1))*cdelt1+crval1)

    XL=Table()
    XL['l']=xl

    xb=[None]*naxis2

    for i in range(0, naxis2):
                xb[i]=((i-(crpix2-1))*cdelt2+crval2)

    XB=Table()
    XB['b']=xb


#output label


    lpix=(np.abs(XL['l']-l)).argmin()
    bpix=(np.abs(XB['b']-b)).argmin()
    vpix=(np.abs(X['x']-v)).argmin()

    label=image[vpix,bpix,lpix]

    return label




# #prepare catalog, add peaks
#
#
#
# data1 = fits.getdata('solenoidal.fits', 1)
# cat= Table(data1)
#
#
# data1 = fits.getdata('mc552021.fits', 1)
# scimes_cat= Table(data1)
#
# peakx = []
# peaky = []
# peakz = []
# galid = list(cat['gal_id'])
# scimes_galid  = list(scimes_cat['gal_id'])
#
# for d in galid:
#     print(d)
#
#     ind = scimes_galid.index(d)
#     peakx.append(scimes_cat['l_peak'][ind])
#     peaky.append(scimes_cat['b_peak'][ind])
#     peakz.append(scimes_cat['v_peak'][ind])
#
#
# cat['peak_x'] = peakx
# cat['peak_y'] = peaky
# cat['peak_z'] = peakz
#
#
# cat.write('solenoidal.fits', format= 'fits', overwrite = True)
#
#



#scimes open CHIMPS maps

hdu_1 = fits.open('reduced_cohrs_0.fits')[0]
chimps0 = hdu_1.data
hd_0 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_1.fits')[0]
chimps1 = hdu_1.data
hd_1 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_2.fits')[0]
chimps2 = hdu_1.data
hd_2 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_3.fits')[0]
chimps3 = hdu_1.data
hd_3 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_4.fits')[0]
chimps4 = hdu_1.data
hd_4 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_5.fits')[0]
chimps5 = hdu_1.data
hd_5 = hdu_1.header


hdu_1 = fits.open('reduced_cohrs_6.fits')[0]
chimps6 = hdu_1.data
hd_6 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_7.fits')[0]
chimps7 = hdu_1.data
hd_7 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_8.fits')[0]
chimps8 = hdu_1.data
hd_8 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_9.fits')[0]
chimps9 = hdu_1.data
hd_9 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_10.fits')[0]
chimps10 = hdu_1.data
hd_10 = hdu_1.header


hdu_1 = fits.open('reduced_cohrs_11.fits')[0]
chimps11 = hdu_1.data
hd_11 = hdu_1.header


hdu_1 = fits.open('reduced_cohrs_12.fits')[0]
chimps12 = hdu_1.data
hd_12 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_13.fits')[0]
chimps13 = hdu_1.data
hd_13 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_14.fits')[0]
chimps14 = hdu_1.data
hd_14 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_15.fits')[0]
chimps15 = hdu_1.data
hd_15 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_16.fits')[0]
chimps16 = hdu_1.data
hd_16 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_17.fits')[0]
chimps17 = hdu_1.data
hd_17 = hdu_1.header


chimps_data = (chimps0, chimps1, chimps2, chimps3, chimps4, chimps5, chimps6, chimps7, chimps8, chimps9, chimps10, chimps11, chimps12, chimps13, chimps14, chimps15, chimps16, chimps17)
chimps_hd = (hd_0, hd_1, hd_2, hd_3, hd_4, hd_5, hd_6, hd_7, hd_8, hd_9, hd_10, hd_11, hd_12, hd_13, hd_14, hd_15, hd_16, hd_17)


data1 = fits.getdata('solenoidal.fits', 1)
cat= Table(data1)

cohrs_galid = []


#####OK


for i in range(len(cat)):

    print(i)

    x = cat['peak_x'][i]
    y = cat['peak_y'][i]
    z = cat['peak_z'][i] #/1000 #km/s

    idx = []

    for g in range(18):

        idx = get_label(chimps_data[g], chimps_hd[g], x, y, z)

        if idx != -1.0:
            cohrs_galid.append(str(g)+'_'+str(idx))

            break


        if g == 17: #the loop wasnt broken before, so idx = -1.0 in all regions
            cohrs_galid.append('NA')




cat['cohrs_gal_id'] = cohrs_galid

cat.write('solenoidal_with_cohrs_galid.fits', format='fits', overwrite=True)



uni = unique(cat, keys= ['cohrs_gal_id']) #full table

uni2 = unique(cat, keys= ['cohrs_gal_id'])['cohrs_gal_id'] #column

#remove 'NA' entry (last of the list)

uni = uni[0:len(uni)-1]
uni2 = uni2[0:len(uni2)-1]


dist = []
#galcen = []

for i in uni2:
    print(i)

    mask = (cat['cohrs_gal_id'] == i)
    cat2 = cat[mask]

    dist.append(statistics.mode(cat2['gistance']))
#    galcen.append(statistics.mode(cat2['galcen_distance']))


#uni['cohrs_galcen_distance'] = galcen
uni['cohrs_distance'] = dist

uni.write('cohrs_distances_from_chimps.fits', format='fits', overwrite=True)






####Distances with Atlasgal

#start with fc (same idea as with CHIMPS)

#scimes open CHIMPS maps

hdu_1 = fits.open('reduced_cohrs_0.fits')[0]
chimps0 = hdu_1.data
hd_0 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_1.fits')[0]
chimps1 = hdu_1.data
hd_1 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_2.fits')[0]
chimps2 = hdu_1.data
hd_2 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_3.fits')[0]
chimps3 = hdu_1.data
hd_3 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_4.fits')[0]
chimps4 = hdu_1.data
hd_4 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_5.fits')[0]
chimps5 = hdu_1.data
hd_5 = hdu_1.header


hdu_1 = fits.open('reduced_cohrs_6.fits')[0]
chimps6 = hdu_1.data
hd_6 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_7.fits')[0]
chimps7 = hdu_1.data
hd_7 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_8.fits')[0]
chimps8 = hdu_1.data
hd_8 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_9.fits')[0]
chimps9 = hdu_1.data
hd_9 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_10.fits')[0]
chimps10 = hdu_1.data
hd_10 = hdu_1.header


hdu_1 = fits.open('reduced_cohrs_11.fits')[0]
chimps11 = hdu_1.data
hd_11 = hdu_1.header


hdu_1 = fits.open('reduced_cohrs_12.fits')[0]
chimps12 = hdu_1.data
hd_12 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_13.fits')[0]
chimps13 = hdu_1.data
hd_13 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_14.fits')[0]
chimps14 = hdu_1.data
hd_14 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_15.fits')[0]
chimps15 = hdu_1.data
hd_15 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_16.fits')[0]
chimps16 = hdu_1.data
hd_16 = hdu_1.header

hdu_1 = fits.open('reduced_cohrs_17.fits')[0]
chimps17 = hdu_1.data
hd_17 = hdu_1.header


chimps_data = (chimps0, chimps1, chimps2, chimps3, chimps4, chimps5, chimps6, chimps7, chimps8, chimps9, chimps10, chimps11, chimps12, chimps13, chimps14, chimps15, chimps16, chimps17)
chimps_hd = (hd_0, hd_1, hd_2, hd_3, hd_4, hd_5, hd_6, hd_7, hd_8, hd_9, hd_10, hd_11, hd_12, hd_13, hd_14, hd_15, hd_16, hd_17)


#Fix Atlasgal catalog

#hi-Gal table



data = fits.getdata('asu.fits', 1)
ag_table = Table(data)
mask = (ag_table['_Glon']>27.0)
ag_table = ag_table[mask]
mask = (ag_table['_Glon']<47.0)
ag_table = ag_table[mask]


cat = ag_table

#check radial velocities
for i in range(len(cat)-1, -1, -1):
    print(i)
    if np.isnan(cat['VLSR'][i]):
        cat.remove_row(i)


cat.write('ag_chimps_cat.fits', format = 'fits', overwrite = True)

data1 = fits.getdata('ag_chimps_cat.fits', 1)
cat= Table(data1)

cohrs_galid = []

for i in range(len(cat)):

    print(i)

    x = cat['_Glon'][i]
    y = cat['_Glat'][i]
    z = cat['VLSR'][i] #/1000 #km/s

    idx = []

    for g in range(18):

        idx = get_label(chimps_data[g], chimps_hd[g], x, y, z)

        if idx != -1.0:
            cohrs_galid.append(str(g)+'_'+str(idx))

            break


        if g == 17: #the loop wasnt broken before, so idx = -1.0 in all regions
            cohrs_galid.append('NA')




cat['cohrs_gal_id'] = cohrs_galid



uni = unique(cat, keys= ['cohrs_gal_id']) #full table

uni2 = unique(cat, keys= ['cohrs_gal_id'])['cohrs_gal_id'] #column

#remove 'NA' entry (last of the list)

uni = uni[0:len(uni)-1]
uni2 = uni2[0:len(uni2)-1]


dist = []
#galcen = []

for i in uni2:
    print(i)

    mask = (cat['cohrs_gal_id'] == i)
    cat2 = cat[mask]

    dist.append(statistics.mode(cat2['Dist']))
#    galcen.append(statistics.mode(cat2['galcen_distance']))


#uni['cohrs_galcen_distance'] = galcen
uni['cohrs_distance'] = dist

uni.write('cohrs_distances_from_ag.fits', format='fits', overwrite=True)





######Make distance assignments



#open COHRS catalogue and remove all entries that are in uni



data1 = fits.getdata('full_cohrs_cat.fits', 1)
fc= Table(data1)


f = copy.deepcopy(fc)



#cohrs distance from chimps - priority

data1 = fits.getdata('cohrs_distances_from_chimps.fits', 1)
cc= Table(data1)

galid_cc = list(cc['cohrs_gal_id'])


#cohrs sources without chimps distances
for i in range(len(fc)-1, -1, -1):
    if fc['gal_id'][i] in galid_cc:
        fc.remove_row(i)

 #cohrs sources with CHIMPS distances
for i in range(len(f)-1, -1, -1):
     if f['gal_id'][i] not in galid:
        f.remove_row(i)




data1 = fits.getdata('cohrs_distances_from_ag.fits', 1)
ac= Table(data1)

galid_ag = list(ac['cohrs_gal_id'])


#cohrs sources without AG distances
for i in range(len(fc)-1, -1, -1):
    if fc['gal_id'][i] in galid_ag:
        fc.remove_row(i)




distance = [0]*len(f)

for i in range(len(f)):
    if f['gal_id'][i] in galid_ag:
        ind = galid_ag.index(f['gal_id'][i])
        distance[i] = ac['cohrs_distance'][ind]

    if f['gal_id'][i] in galid_cc:
        ind = galid_cc.index(f['gal_id'][i])
        distance[i] = cc['cohrs_distance'][ind]

f['distance'] = distance


f.write('cohrs_cat.fits', format='fits', overwrite = True)



###Galcen distance and distances with rotation curve


#open cohrs_cat.fits
data1 = fits.getdata('cohrs_cat.fits', 1)
t= Table(data1)


def omega (l, b, v):

    l = math.radians(l)
    b = math.radians(b)

    K = 1 + v/(220*math.sin(l)*math.cos(b))

    return K


# def f(a1, c, a3, K):
#
#     return a1*Y**c + a3/Y - K

a1, a2, a3 = 1.0077, 0.0394, 0.0071

R0, omega0 = 8.5, 220



omegas = []
Ks = []


for i in range(len(t)):

    K= omega(t['x_cen'][i], t['y_cen'][i], t['v_cen'][i]/1000)
    print(K)
    Ks.append(K)



guess1 = []


for i in range(len(t)):

    if t['v_cen'][i] < 17.4:
        guess1.append(2.0)

    elif t['v_cen'][i] > 113.4:
        guess1.append(6.0)

    else:
        guess1.append(4.0)


gcd2 = []


from scipy.optimize import fsolve

for i in range(len(t)):

    K = Ks[i]
    f=lambda x: a1*x**a2 -K*x + a3
    gcd2.append(fsolve(f, guess1[i]))




GCD_2 = []




for i in range(len(gcd2)):

    GCD_2.append(gcd2[i][0]*8.5)


t['galcen_dist'] = GCD_2

t.write('cohrs_cat.fits', format='fits', overwrite = True)


##now heliocentric distances

#open cohrs_cat.fits
data1 = fits.getdata('cohrs_cat.fits', 1)
t= Table(data1)

dist = []

for i in range(len(t)):

    if t['distance'][i] == 0.0:
        a = (math.cos(math.radians(t['y_cen'][i])))**2
        b = -2*R0*math.cos(math.radians(t['y_cen'][i]))*math.cos(math.radians(t['x_cen'][i]))
        c = R0**2 - t['galcen_dist'][i]**2

        d = (b**2) - (4*a*c)

        if d < 0.0:
            dist.append(R0*math.cos(math.radians(t['x_cen'][i]))/math.cos(math.radians(t['y_cen'][i])))
            print('complex')

        elif d ==0:
            dist.append((-b-cmath.sqrt(d))/(2*a))


        elif d > 0:
            print('else')
            sol1 = np.real((-b-cmath.sqrt(d))/(2*a))
            sol2 = np.real((-b+cmath.sqrt(d))/(2*a))

            if sol1  > sol2:
                near = sol2
                far = sol1
                h =  R0*math.sin(math.radians(t['y_cen'][i]))

                #check condition
                if h > 0.120:
                    dist.append(near)

                else:
                    dist.append(far)

            else:
                near = sol1
                far = sol2
                h =  R0*math.sin(math.radians(t['y_cen'][i]))

                #check condition
                if h > 0.120:
                    dist.append(near)

                else:
                    dist.append(far)


    else:
        dist.append(t['distance'][i])



t['distance2'] = dist

t.write('cohrs_cat.fits', format='fits', overwrite = True)





