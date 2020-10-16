#!/usr/bin/env python
import numpy as np
import math
from netCDF4 import Dataset
from pylab import *
#import pdb
import netCDF4 as nc

#based on a simple Lehmer "random" number generator
#where the final result is scaled from 0 to 1
def lrng(count,m,a,initval):

        r = m % a
        q = m // a

        z=np.zeros(count+1)
        z[0]=initval

        for i in range(0,count):
                t = a * (z[i] % q) - r * (z[i] // q)
                if (t > 0):
                        z[i+1]=t
                else:
                        z[i+1]=t+m
        #z = z/m

        #normalize within [0,1]
        z=z-np.min(z)
        z=z/(np.max(z))

        return z

def Create_iceberg_restart_file(Number_of_bergs, lon,lat,thickness,width,mass,mass_scaling,iceberg_num,Ice_geometry_source,static_berg):

	print 'Writing iceberg restart files, with ' , Number_of_bergs  , 'icebergs..'
	# To copy the global attributes of the netCDF file

	#Input and output files
	#Create Empty restart file. This is later read so that the attributes can be used.
	Empty_restart_filename='output_files/Empty_icebergs.res.nc'
	create_empty_iceberg_restart_file(Empty_restart_filename)
	#Empty_restart_filename='input_files/icebergs.res.nc'

	#Read empty restart file
	f=Dataset(Empty_restart_filename,'r') # r is for read only
	#Write a new restart file
	g=Dataset('output_files/' + Ice_geometry_source + '_icebergs.res.nc','w', format='NETCDF3_CLASSIC') # w if for creating a file

	for attname in f.ncattrs():
		    setattr(g,attname,getattr(f,attname))


	# To copy the dimension of the netCDF file
	for dimname,dim in f.dimensions.iteritems():
		# if you want to make changes in the dimensions of the new file
		# you should add your own conditions here before the creation of the dimension.
		#g.createDimension(dimname,len(dim))
		g.createDimension(dimname,Number_of_bergs)

	# To copy the variables of the netCDF file

	for varname,ncvar in f.variables.iteritems():
		# if you want to make changes in the variables of the new file
		# you should add your own conditions here before the creation of the variable.
		var = g.createVariable(varname,ncvar.dtype,ncvar.dimensions)
		#Proceed to copy the variable attributes
		for attname in ncvar.ncattrs():
			setattr(var,attname,getattr(ncvar,attname))
		#Finally copy the variable data to the new created variable
		#var[:] = ncvar[0]  #I commented out this line because it was causing errors. I'm not sure if it is needed.

		if varname=='i':
			var[:]=Number_of_bergs

		if varname=='iceberg_num':
			for j in range(Number_of_bergs):
				#var[j]=j+1
				var[j]=iceberg_num[j]

		if varname=='uvel' or varname=='vvel' or varname=='uvel_old' or varname=='vvel_old' or varname=='axn' or varname=='ayn'\
		or varname=='bxn' or varname=='byn' or  varname=='halo_berg' or varname=='heat_density' or varname=='lon_old' or varname=='lat_old' \
		or varname=='mass_of_bits' or varname=='start_mass' or  varname=='start_day' or varname=='start_year' or varname=='start_lon' \
		or varname=='start_lat' or varname=='start_mass' or  varname=='start_day' or varname=='start_year' or varname=='start_lon' or varname=='lat_old':\
			var[:]=0

		if varname=='mass_scaling':
			var[:]=mass_scaling

		if varname=='thickness':
			for j in range(Number_of_bergs):
				var[j]=thickness[j]

		if varname=='mass':
			for j in range(Number_of_bergs):
				var[j]=mass[j]

		if varname=='width'  or varname=='length':
			for j in range(Number_of_bergs):
				var[j]=width[j]

		if varname=='lon':
			for j in range(Number_of_bergs):
				var[j]=lon[j]

		if varname=='lat':
			for j in range(Number_of_bergs):
				var[j]=lat[j]

		if varname=='static_berg':
			for j in range(Number_of_bergs):
				var[j]=static_berg[j]


	f.close()
	g.close()


def create_empty_iceberg_restart_file(Empty_restart_filename):

	f = Dataset(Empty_restart_filename,'w', format='NETCDF3_CLASSIC')

	i=f.createDimension('i', None)
	lon=f.createVariable('i','i')

	lon=f.createVariable('lon','d',('i'))
	lon.long_name = "longitude" ;
	lon.units = "degrees_E" ;
	lon.checksum = "               0" ;

	lat=f.createVariable('lat','d',('i'))
	lat.long_name = "latitude" ;
	lat.units = "degrees_N" ;
	lat.checksum = "               0" ;

	uvel=f.createVariable('uvel','d',('i'))
	uvel.long_name = "zonal velocity" ;
	uvel.units = "m/s" ;
	uvel.checksum = "               0" ;

	vvel=f.createVariable('vvel','d',('i'))
	vvel.long_name = "meridional velocity" ;
	vvel.units = "m/s" ;
	vvel.checksum = "               0" ;

	mass=f.createVariable('mass','d',('i'))
	mass.long_name = "mass" ;
	mass.units = "kg" ;
	mass.checksum = "               0" ;

	axn=f.createVariable('axn','d',('i'))
	axn.long_name = "explicit zonal acceleration" ;
	axn.units = "m/s^2" ;
	axn.checksum = "               0" ;

	ayn=f.createVariable('ayn','d',('i'))
	ayn.long_name = "explicit meridional acceleration" ;
	ayn.units = "m/s^2" ;
	ayn.checksum = "               0" ;

	bxn=f.createVariable('bxn','d',('i'))
	bxn.long_name = "inplicit zonal acceleration" ;
	bxn.units = "m/s^2" ;
	bxn.checksum = "               0" ;

	byn=f.createVariable('byn','d',('i'))
	byn.long_name = "implicit meridional acceleration" ;
	byn.units = "m/s^2" ;
	byn.checksum = "               0" ;

	ine=f.createVariable('ine','i',('i'))
	ine.long_name = "i index" ;
	ine.units = "none" ;
	ine.packing = 0 ;
	ine.checksum = "               0" ;

	jne=f.createVariable('jne','i',('i'))
	jne.long_name = "j index" ;
	jne.units = "none" ;
	jne.packing = 0 ;
	jne.checksum = "               0" ;

	thickness=f.createVariable('thickness','d',('i'))
	thickness.long_name = "thickness" ;
	thickness.units = "m" ;
	thickness.checksum = "               0" ;

	width=f.createVariable('width','d',('i'))
	width.long_name = "width" ;
	width.units = "m" ;
	width.checksum = "               0" ;

	length=f.createVariable('length','d',('i'))
	length.long_name = "length" ;
	length.units = "m" ;
	length.checksum = "               0" ;

	start_lon=f.createVariable('start_lon','d',('i'))
	start_lon.long_name = "longitude of calving location" ;
	start_lon.units = "degrees_E" ;
	start_lon.checksum = "               0" ;

	start_lat=f.createVariable('start_lat','d',('i'))
	start_lat.long_name = "latitude of calving location" ;
	start_lat.units = "degrees_N" ;
	start_lat.checksum = "               0" ;

	start_year=f.createVariable('start_year','i',('i'))
	start_year.long_name = "calendar year of calving event" ;
	start_year.units = "years" ;
	start_year.packing = 0 ;
	start_year.checksum = "               0" ;

	iceberg_num=f.createVariable('iceberg_num','i',('i'))
	iceberg_num.long_name = "identification of the iceberg" ;
	iceberg_num.units = "dimensionless" ;
	iceberg_num.packing = 0 ;
	iceberg_num.checksum = "               0" ;

	start_day=f.createVariable('start_day','d',('i'))
	start_day.long_name = "year day of calving event" ;
	start_day.units = "days" ;
	start_day.checksum = "               0" ;

	start_mass=f.createVariable('start_mass','d',('i'))
	start_mass.long_name = "initial mass of calving berg" ;
	start_mass.units = "kg" ;
	start_mass.checksum = "               0" ;

	mass_scaling=f.createVariable('mass_scaling','d',('i'))
	mass_scaling.long_name = "scaling factor for mass of calving berg" ;
	mass_scaling.units = "none" ;
	mass_scaling.checksum = "               0" ;

	mass_of_bits=f.createVariable('mass_of_bits','d',('i'))
	mass_of_bits.long_name = "mass of bergy bits" ;
	mass_of_bits.units = "kg" ;
	mass_of_bits.checksum = "               0" ;

	heat_density=f.createVariable('heat_density','d',('i'))
	heat_density.long_name = "heat density" ;
	heat_density.units = "J/kg" ;
	heat_density.checksum = "               0" ;

	halo_berg=f.createVariable('halo_berg','d',('i'))
	halo_berg.long_name = "halo_berg" ;
	halo_berg.units = "dimensionless" ;
	halo_berg.checksum = "               0" ;

	static_berg=f.createVariable('static_berg','d',('i'))
	static_berg.long_name = "static_berg" ;
	static_berg.units = "dimensionless" ;
	static_berg.checksum = "               0" ;

	f.sync()
	f.close()



#-----------------------------------#
#               Main                #
#-----------------------------------#

#ice params
rho_ice=850

#grid parameters (m)
grdres=5000
grdxmin=0; grdxmax=1000e3
grdymin=0; grdymax=1000e3

#--- berg locs: ---
nbergs_x=10 #num of conglomerate bergs in x-dir
nbergs_y=10 #num of conglomerate bergs in y-dir
nbergs=nbergs_x*nbergs_y
start=50e3; end=950e3; #m
#center coords of each (rectangular) conglomerate berg (CB=conglomerate berg)
CBxc=np.linspace(start,end,num=nbergs_x); CByc=np.linspace(start,end,num=nbergs_y)
#CBxc=np.tile(CBxc,nbergs_y); CByc=np.tile(CByc,nbergs_x)
CBxc,CByc=np.meshgrid(CBxc,CByc)
CBxc=CBxc.flatten()
CByc=CByc.flatten()

#reduced version
# mask = abs(500e3-CByc)
# CBxc=CBxc[mask<100e3]
# CByc=CByc[mask<100e3]
# mask = abs(500e3-CBxc)
# CBxc=CBxc[mask<300e3]
# CByc=CByc[mask<300e3]
# nbergs=len(CBxc)

#reduced version 2
#CBxc=CBxc[CByc<100e3]
#CByc=CByc[CByc<100e3]
#nbergs=len(CBxc)

#--- Lehmer pseudo-random number generator params ---
m=2147483647; a=16807l
iv=a #first value

#--- side lengths for the rectangular bergs ---
m=2147483647; a=16807 #Lehmer pseudo-random number generator params
rn=lrng(nbergs+1,m,a,iv) #the pseudo-random numbers
minn=5e3; maxn=60e3 #bounds of side lengths, in m (diag should be < length of a processor)
start=0; end=nbergs #slice of random number array with which to work
CBxl=(maxn-minn)*rn[start:end]+minn #array of berg lengths in the x-dir, adjusted for bounds
CByl=rn[start+1:end+1] #array of berg lengths in the y-dir
CByl[CByl<0.33]=0.33 #keep the aspect ratio somewhat realistic
CBxl=np.around(CBxl)
CByl=CByl*CBxl
CByl=np.around(CByl)
CByl[CByl<minn]=minn
CByl=CByl.astype(int); CBxl=CBxl.astype(int)
#flip CBxl and CByl for half of them, so that CBxl is not always > CByl
yl2=CByl.copy(); xl2=CBxl.copy()
CByl[start:end:2]=xl2[start:end:2]; CBxl[start:end:2]=yl2[start:end:2]

#--- xmax, xmin, ymax, ymin for each CB --
CBxmin=CBxc-(0.5*CBxl); CBxmax=CBxc+(0.5*CBxl)
CBymin=CByc-(0.5*CByl); CBymax=CByc+(0.5*CByl)
#this shouldn't happen:
CBxmin[CBxmin<grdxmin]=grdxmin; CBxmax[CBxmax>grdxmax]=grdxmax
CBymin[CBxmin<grdymin]=grdymin; CBymax[CBymax>grdymax]=grdymax
CBxc=0.5*(CBxmin+CBxmax); CByc=0.5*(CBymin+CBymax)

#--- min and max thicknesses for each CB ---
minn=50; maxn=350
#minn=250; maxn=250
iv=rn[-1]; rn=lrng(nbergs+1,m,a,iv)
h1=(maxn-minn)*rn[start:end]+minn
h2=(maxn-minn)*rn[start+1:end+1]+minn
CBhmax=np.maximum(h1,h2); CBhmin=np.minimum(h1,h2)

#--- berg element radii ---
R_frac=0.45
maxn=(np.sqrt(3)/2.)*(R_frac*grdres) #S is < 0.5 grid res
minn=maxn/2 # smaller radii are allowed
iv=rn[-1]; rn=lrng(nbergs+1,m,a,iv)
CBrad=(maxn-minn)*rn[start:end]+minn

print('maxn,minn',maxn,minn)

#-- create berg elements for each conglomerate berg --
# CBxc,CByc=center coords for each conglom berg
# CByl,CBxl=side lengths for each conglom berg
# CBhmax,CBhmin=max and min thicknesses on each berg
# CBrad=radius of each berg element
# CBnx,CBny=number of bergs in each dir
# CBxmin,CBxmax,CBymin,CBymax

berg_x=[]; berg_y=[]
berg_id=[]; berg_static=[]
berg_width=[]; berg_bonds=[]
berg_h=[]; berg_mass_scaling=[]
berg_mass=[]

#berg_CBid=[];
berg_count=0

for i in range(nbergs):
        x_start=CBxmin[i]+(CBrad[i]*2./np.sqrt(3))

        #pdb.set_trace()
        if x_start>CBxmax[i]:
                x_start=CBxmax[i]
        y_start0=CBymin[i]+CBrad[i]
        if y_start0>CBymax[i]:
                y_start0=CBymax[i]
        element_area=(3.*np.sqrt(3.)/2.)*((4./3.)*(CBrad[i])**2)
        #H (thickness) is a linear function of a berg element's position from the
        #center of the CB. H=Hmax at the center of the CB and H=Hmin at a corner of the
        #CB, where the dist of a corner from the center is:
        cdistb=np.sqrt((CBxmin[i]-CBxc[i])**2+(CBymin[i]-CByc[i])**2)
        j=0
        x_val=x_start
        #berg_count_start=berg_count
        while x_val<=CBxmax[i]:
                y_start=y_start0+((j%2)*CBrad[i])
                k=0
                y_val=y_start
                while y_val<CBymax[i]:
                        berg_count=berg_count+1
                        berg_id.append(berg_count)
                        berg_x.append(x_val)
                        berg_y.append(y_val)
                        berg_width.append(sqrt(element_area))
                        #dist of berg elem from center of CB
                        bdistc=np.sqrt((x_val-CBxc[i])**2+(y_val-CByc[i])**2)
                        bh=CBhmin[i]*bdistc/cdistb + CBhmax[i]*(1-bdistc/cdistb)
                        berg_h.append(bh) #thickness
                        berg_mass_scaling.append(1)
                        berg_mass.append(bh*rho_ice*element_area)
                        berg_static.append(0)
                        #berg_CBid.append(i)
                        k=k+1
                        y_val=y_start+(2*k*CBrad[i])
                j=j+1
                x_val=x_start+(np.sqrt(3)*CBrad[i]*j)

print('Number of bergs',berg_count)

#pdb.set_trace()

#Create iceberg restart file
Ice_geometry_source='Generic'
Create_iceberg_restart_file(berg_count,berg_x,berg_y,berg_h,berg_width,berg_mass,berg_mass_scaling,berg_id,Ice_geometry_source,berg_static)
