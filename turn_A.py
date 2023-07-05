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
import copy




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




data1 = fits.getdata('final_cat_region0.fits', 1)
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







#add distances



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

hdu_1 = fits.open('final_cohrs_0.fits')[0]
chimps0 = hdu_1.data
hd_0 = hdu_1.header

hdu_1 = fits.open('final_cohrs_1.fits')[0]
chimps1 = hdu_1.data
hd_1 = hdu_1.header

hdu_1 = fits.open('final_cohrs_2.fits')[0]
chimps2 = hdu_1.data
hd_2 = hdu_1.header

hdu_1 = fits.open('final_cohrs_3.fits')[0]
chimps3 = hdu_1.data
hd_3 = hdu_1.header

hdu_1 = fits.open('final_cohrs_4.fits')[0]
chimps4 = hdu_1.data
hd_4 = hdu_1.header

hdu_1 = fits.open('final_cohrs_5.fits')[0]
chimps5 = hdu_1.data
hd_5 = hdu_1.header


hdu_1 = fits.open('final_cohrs_6.fits')[0]
chimps6 = hdu_1.data
hd_6 = hdu_1.header

hdu_1 = fits.open('final_cohrs_7.fits')[0]
chimps7 = hdu_1.data
hd_7 = hdu_1.header

hdu_1 = fits.open('final_cohrs_8.fits')[0]
chimps8 = hdu_1.data
hd_8 = hdu_1.header

hdu_1 = fits.open('final_cohrs_9.fits')[0]
chimps9 = hdu_1.data
hd_9 = hdu_1.header

hdu_1 = fits.open('final_cohrs_10.fits')[0]
chimps10 = hdu_1.data
hd_10 = hdu_1.header


hdu_1 = fits.open('final_cohrs_11.fits')[0]
chimps11 = hdu_1.data
hd_11 = hdu_1.header


hdu_1 = fits.open('final_cohrs_12.fits')[0]
chimps12 = hdu_1.data
hd_12 = hdu_1.header

hdu_1 = fits.open('final_cohrs_13.fits')[0]
chimps13 = hdu_1.data
hd_13 = hdu_1.header

hdu_1 = fits.open('final_cohrs_14.fits')[0]
chimps14 = hdu_1.data
hd_14 = hdu_1.header

hdu_1 = fits.open('final_cohrs_15.fits')[0]
chimps15 = hdu_1.data
hd_15 = hdu_1.header

hdu_1 = fits.open('final_cohrs_16.fits')[0]
chimps16 = hdu_1.data
hd_16 = hdu_1.header

hdu_1 = fits.open('final_cohrs_17.fits')[0]
chimps17 = hdu_1.data
hd_17 = hdu_1.header


chimps_data = (chimps0, chimps1, chimps2, chimps3, chimps4, chimps5, chimps6, chimps7, chimps8, chimps9, chimps10, chimps11, chimps12, chimps13, chimps14, chimps15, chimps16, chimps17)
chimps_hd = (hd_0, hd_1, hd_2, hd_3, hd_4, hd_5, hd_6, hd_7, hd_8, hd_9, hd_10, hd_11, hd_12, hd_13, hd_14, hd_15, hd_16, hd_17)


data1 = fits.getdata('solenoidal.fits', 1)
cat= Table(data1)

cohrs_galid = []


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
galcen = []

for i in uni2:

    mask = (cat['cohrs_gal_id'] == i)
    cat2 = cat[mask]

    dist.append(statistics.mode(cat2['gistance']))
    galcen.append(statistics.mode(cat2['galcen_distance']))


uni['cohrs_galcen_distance'] = galcen
uni['cohrs_distance'] = dist

uni.write('cohrs_distances_from_chimps.fits', format='fits', overwrite=True)


#open COHRS catalogue and remove all entries that are in uni



data1 = fits.getdata('full_cohrs_cat.fits', 1)
fc= Table(data1)


f = fc


data1 = fits.getdata('cohrs_distances_from_chimps.fits', 1)
cc= Table(data1)

galid = list(cc['cohrs_gal_id'])


for i in range(len(fc)-1, -1, -1):
    if fc['gal_id'][i] in galid:
        fc.remove_row(i)


galid = galid[0:(len(galid)-2)]

for i in range(len(f)-1, -1, -1):
    if f['gal_id'][i] not in galid:
        f.remove_row(i)




# estimate Reid distances for the sources that are left

p=[0.5]*len(fc)
   # z3['p']=p

z3 = fc


ascii.write([z3['gal_id'], z3['x_cen'], z3['y_cen'], z3['v_cen']/1000, p, z3['gal_id']], 'sources_info.inp', names=['id','l', 'b', 'v', 'p', 'extra'], overwrite=True)



#read in Reid distances

ta1=Table.read('output.dat', format='ascii')

dist_col=ta1['col4']

dist_col = [0.0 if math.isnan(y) else y for y in dist_col]

z3['distance']=dist_col


dist2 = []

for i in range(len(z3)):

    upper_limit = z3['distance'][i] - 1.0
    lower_limit = z3['distance'][i] - 1.0
    dist2.append(random() * (upper_limit - lower_limit) + lower_limit)

z3['distance2'] = dist2



f['distance'] = cc['cohrs_distance']


f['distance2'] = 0*len(f)

T = vstack([f, ze])

T.write('full_cohrs_with_distances', format='fits', overwrite = True)













