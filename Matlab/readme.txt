The matlab code is split into four functions. 
So command lines have to be used to run those functions, otherwise, the results will not appear in the workspace.

The python code is for reading the GW sky localisation information fits file,
and writing the information into a txt file.

Following are the steps to read a fits file of a gravitational signal triggers 
and to generate a follow-up observing strategy for that gravitational signal:
1. Run the python script ReadingFits to read the GW information from the interested fits file.
Depending on the way of ordering (ring or nested), you might need to change the setting of the function hp.pix2ang
in the code.
One possible source of fits files to test the code is the following link:
http://www.ligo.org/scientists/first2years/
if you prefer to download a txt file that is ready for testing the code, you can go to this link:
https://mega.nz/#!1dY1wLhS!bxpETr6FyAHwOVCR7dlJOljr4D8VKc_g7-imBkwnKiE

28700.txt is a simulated GW event taken from (Singer, et al, 2014). The ID of the event is 28700.
This size of 90% credible region is ~300 deg^2. The resolution of this file is 9, which corresponds to more than
3 millions healpix points on the skymap.


2. Run Pem: 
This function may take more than one hour to run depending on the number of samples of tau (observation times).
For one telescope, this function will only need to be run once if the number of samples is the same. 
So saving the data after a run is recommended.

3. Run Greedy:
This function may take ~40 minutes to run depending on the number of healpix points on the skymap. 
For the same event and the same telescope, the result would be the same if the number of fields is the same.
So saving after a run is also recommended.

4. Run main:
This function will call function Solver, which will return the detection probability, and the time allocation
given that the total observation time is distributed equally across all the observed fields 

5. rotx; roty; rotz:
These are three functions that rotate the coordinate system. Usually there is no need to run them individiually. 

Using the following command lines will generate an observing strategy for a telescope with the following properties 
for simulated GW event 28700 assumming the event at a distance equal to 100 Mpc:

telescopes:

1.sensitivity: 21 mag in R band in 180 seconds;
2.radius of the aperture: 0.5m;
3.fov: 23X23 arc minc;

the command lines:
1. [tau,prob]= Pem(0.5,21,180,100e6,30e6,61,[],[],[]);
2. [PGW,fields_location]=Greedy(23/60*23/60,100,'28700.txt',9,0.9);
3. [tprob,time_allocation]= Main(6,0,1,PGW,tau,prob);





	HEALPix Pixel Information
Res	NSide	NPixels		Mean Spacing　  		  Area
				　　　　(deg)    		(sterad)
0	1	12		   58.6323	     1.0471976 X 10+00
1	2	48		   29.3162	     2.6179939 X 10-01
2	4	192		   14.6581	     6.5449847 X 10-02
3	8	768		    7.3290	     1.6362462 X 10-02
4	16	3072		    3.6645	     4.0906154 X 10-03
5	32	12288		    1.8323	     1.0226539 X 10-03
6	64	49152		    0.9161	     2.5566346 X 10-04
7	128	196608		    0.4581	     6.3915866 X 10-05
8	256	786432		    0.2290	     1.5978967 X 10-05
9	512	3145728		    0.1145	     3.9947416 X 10-06
10	1024	12582912	    0.0573	     9.9868541 X 10-07
