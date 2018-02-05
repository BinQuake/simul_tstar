c23456789a123456789b123456789c123456789d123456789e123456789f123456789g12
c  Program SIMUL2014 (started January 2014)
c   Has transverse mercator in for coordinate conversion.  
c    (AK already had TM, now NZ has NZTM2000, and TM general option.)
c   This has S-P residuals related to both dvp and dvpvs, changes
c   to ttmder
c   This includes teleseismic events
c  Program SIMUL2010 (started May 2009; mostly done Sept 2010 now compiles on both sun an linux)
c
c   Updated to allow for "flexible" gridding - master nodes and linked
c   nodes that have their values linked and partial derivatives combined,
c   by C. Thurber.
c
c   Modified for Q inversion 15-March-1998 following Andreas Rietbrock's version.
c   The input velocity file has Vp then Q.
c   The input data is tstar with the hypocentre from the 3-D Vp
c   solution (which must be already done for the model region).
c   The input hypocentres are treated as shots and so are not
c   perturbed in the Q inversion.
c
c   Then (when iuseq=1), the second part of the model array
c   is used to store Vp*Qp (not Vp/Vs).
c   The solution is ONLY for Q, not for Vp though.  
c   THUS the solution arrays that would otherwise be related
c   to the first part of the model array (Vp) are used.
c   The fixed points, damping and dvmax otherwise related to 
c   Vp are related to Q or V*Q for iuseq=1.  
c   Fixed node indices for Q are as for Vp otherwise (ie: z indices 2 to nz-1).
c   As this is slightly confusing to the user, qdamp and dqmax are also 
c   included as input values.
c
c   The Q option works for Qp if the t* are P and the input velocity is Vp,
c   or for Qs if t* for S and the input velocity model was Vs and Qs.
c
c   Early simul had Vp and Vs, but Vs had poor resolution and meant 
c   poor Vs and Vp/Vs in areas of good Vp.  This version solves for
c   Vp and Vp/Vs.  The input data are P travel-times and S-P times
c   (previous versions used P travel-times and S travel-times.  See
c   Thurber (ref 1 below) for discussion of solution for Vp/Vs
c   using S-P data.
c   This version also allows the user to vary the weighting, in the 
c   hypocenter solution, of the S-P data relative to the P data. (see wtsp)
c   Vp and Vp/Vs are input and are used to compute Vs at each iteration.
c   The S velocities are on the same grid as the p velocities and are
c   stored in additional nz layer positions.  
c   Uses up to 'maxpar' velocity nodes (and station parameters), 
c   but only invert for up to 'mxpari'.
c   'maxpar' and 'mxpari' are set in ' simul2014_common.inc'.
c   Has the option of outputing the raypaths.
c   Has Thurber's pseudo-bending ray-tracing from his SIMUL3M version
c   Allows less curvature for initial arcuate rays below "moho"
c   (see i3d input parameter).
c
c   Corresponding code authors:
c        Prof. Clifford H. Thurber
c        University of Wisconsin-Madison
c        email:  thurber@geology.wisc.edu
c
c        Donna Eberhart-Phillips, Ph.D.
c        GNS Science
c        email:  d.eberhart-phillips@gns.cri.nz
c
c  If you are using this program please read and reference some of the following
c     chapters by Thurber and Eberhart-Phillips:
c
c  1)  Thurber, C. H., and D. Eberhart-Phillips, Local earthquake tomography with
c            flexible gridding, Computers and Geoscience, v. 25, p. 809-818,
c            1999.
c  2)  Thurber, C. H., Local earthquake tomography: velocities and Vp/Vs - theory,
c            in Seismic Tomography: Theory and Practice, edited by
c            H. M. Iyer and K. Hirahara, 1993. 
c  3)  Eberhart-Phillips, D., Local earthquake tomography: earthquake source regions,
c            in Seismic Tomography: Theory and Practice, edited by
c            H. M. Iyer and K. Hirahara, 1993. 
c  4)  Thurber, C. H., Earthquake locations and three-dimensional crustal structure
c             in the Coyote Lake area, central California, J. Geophys. Res., 
c             v. 88, p. 8226-8236, 1983.
c  5)  For Q inversion:
c      Rietbrock, A., P-wave attenuation structure in the fault area of the 1995
c            Kobe earthquake, J. Geophys. Res., 106, 4141-4154, 2001.
c  6)  For rdt,edt:
c      Eberhart-Phillips, D., Reyners, M., 2012. Imaging the Hikurangi plate 
c            interface region with improved local-earthquake tomography. Geophys. 
c             J. Int., 190, 1221-1242, DOI: 10.1111/j.1365-246X.2012.05553.x.
c
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c  USGS disclaimer statement:
c  ALTHOUGH THIS PROGRAM HAS BEEN USED BY THE USGS, NO WARRANTY,
c  EXPRESSED OR IMPLIED, IS MADE BY THE USGS OR THE UNITED
c  STATES GOVERNMENT AS TO THE ACCURACY AND FUNCTIONING OF THE
c  PROGRAM AND RELATED PROGRAM MATERIAL NOR SHALL THE FACT OF
c  DISTRIBUTION CONSTITUTE ANY SUCH WARRANTY, AND NO 
c  RESPONSIBILITY IS ASSUMED BY THE USGS IN CONNECTION THEREWITH.
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c
c There are many output files possible with various formats that
c were given names by Haslinger
cfh suggestion: incorporate  file-naming for input and output:
c      input:     file 01: CNTL
c                 file 02: STNS
c                 file 03: MOD
c                 file 04: EQKS
c                 file 07: SHOT
c                 file 08: BLAS
c                 file 10: TELE
c                 file 09: CLUS
c      output:    file 12: summary
c                 file 13: finalsmpout
c                 file 15: ev'xxxx'.rp
c                 file 16: output
c                 file 17: resol.out
c                 file 18: nodes.out
c                 file 19: ttdiff.out
c                 file 20: residuals
c                 file 21: rdtres.out
c                 file 22: newstns
c                 file 23: velomod.out
c                 file 24: tteq.out
c                 file 25: vl_layer.out
c                 file 26: pseudo
c                 file 27: ttsht.out
c                 file 28: ttbls.out
c                 file 29: ttcev.out
c                 file 30: ceres.out
c                 file 31: edtres.out
c                 file 32: hypo.gmt
c                 file 43: DWSALL
c                 file 44: DWSALL.GRD  DO NOT NEED 
c                 file 33: hyp_prt.out
c                 file 34: hypo71list
c                 file 36: itersum
c                 file 45: covarid.out
c                 file 53: vlxyxltln.out
c                 file 59: resDRE.out
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c  Written by Cliff Thurber  as part of his PhD thesis.
c  Modified by W. Prothero to include station delays.
c  Parts of this were programmed by Steve Taylor of LLL.
c  Obtained from Prothero&Thurber in 1983 and subsequently
c  modified by Donna Eberhart-Phillips, U.S. Geological Survey.
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c  Input data file list:
c     File 01 - Control Parameters  (read in subroutine INPUT1)
c      line 1- (free format)
c        neqs - number of earthquakes
c        nsht  - number of shots
c        nbls  - number of blasts with known location, but unknown
c                origin time
c        wtsht - weighting of shots relative to quake weighting
c        kout - output control parameter
c               value    files created
c                 0      16,36
c                 1      16,13,36
c                 2      16,13,22,23,24,36
c                 3      16,13,22,23,24,34,36
c                 4      16,13,22,23,24,25,34,36
c                 5      16,12,13,22,23,24,25,34,36
c        kout2 - Printout control parameter
c                0 = full printout including station residuals and
c                    location steps
c                1 = printout station residuals
c                2 = printout location steps
c                3 = don't printout location steps or station residuals
c                4 = as 3 above, also don't printout stations in input2
c                5 = as 4. Also printout 1/(diag. covariance) to for016 
c                    and for045, if ires.gt.0.
c        kout3 - Yet another output control parameter
c                0 = Don't output raypath points or tt differences
c                1 = Output raypath points to file 15, for all raypaths
c                    Print out travel-time differences between ART and
c                    pseudo-bending to File 19.
c                    This is useful for making plots of raypaths.
c                    (For instance, if you want to test a range of pseudo-
c                    bending parameters (xfac,tlim,nitpb) )
c                    **NOTE that this option should not be used regularly,
c                    ** but only to check a few selected events, since
c                    ** it creates a lot of output.
c
c      line 2 - (free format)
c        iuserdt = 1 to use rdt (receiver differential times)
c              iuserdt=1 for eq only, =2 for rdt on shots and blasts also. (25-jan-13)
c        edisrdt = max distance between receivers for rdt
c        facrayr = factor times receiver distance to allow partials along raypaths
c        mxiqrdt = maximum iqual to use ttime for rdt
c        wtrdtsht = relative weight of rdt data
c        res4 = resedt must be less than res4
c      line 3 - (free format)
c        iuseclu = 1 to use edt data for groups of earthquakes
c        nclu = number of clusters to read in from file 9
c        edisedt = max distance between eq for edt.
c        facraye =  factor times eq distance to allow partials along raypaths
c        mxiqedt = maximum iqual to use ttime for edt
c        nstaepair = maximum number of edt for specific eq pair
c        wtepdt = wt factor for eq-pair dt, like wtsp or wtrdtsht
c        wrmscemx = max wrms for ce event, deleted if greater
c      line 4 - (free format)
c        nitloc - max of iterations for hypocenter location.
c        wtsp   - for hypocenter solution, weight of S-P residual 
c                 relative to P residual (ie:wtsp=1.0 gives equal wt,
c                 wtsp<1 downweights S-P)
c        eigtol - SVD cutoff in hypocentral adjustments. If smallest
c                 eigenvalue in geiger's matrix is < eigtol, the depth
c                 is not adjusted, and a message is printed.
c        rmscut - value for rms residual below which hypocentral adjustments
c                 are terminated.
c        zmin - minimum hypocenter depth
c        dxmax - maximum horizontal hypocentral adjustment allowed in
c                each hypocenter iteration.
c
c        rderr - estimate of reading error, used to estimate hypocenter error
c                used for adding error to synthetic times, not for t*
c        ercof - for hypoinverse-like error calculations. Set > 0 and
c                < 1 if you want to include rms.res in hypocenter error
c                estimate. (sigsq=rderr**2 + ercof * rms**2)
c      line 5 - (free format)
c        nhitct - of observations for a parameter to be included
c                 in the inversion. This uses the variable khit.
c        dvpmx - maximum P-velocity adjustment allowed per iteration. 
c        dvsmx - maximum S-velocity adjustment allowed per iteration.
c        corrvpmax = max corrected vp allowed
c        corrvpmin = min corrected vp allowed
c      line 6 - (free format)
c        idmp -  set to 1 to recalculate the damping value for succeeding
c                iterations.  set to 0 to have constant damping
c        vdamp - damping parameter used in velocity inversion.
c                 vdamp(1)=damping for p-velocity
c                 vdamp(2)=damping for Vp/Vs
c                 vdamp(3)=damping for station delays
c        stepl - (km) used for calculation of partial derivatives along
c                the raypath.
c
c      line 7 - (free format)
c        ires - set to 1 to compute the resolution and print diagonal elements,  
c             also prints out 1/(diagonal covariance)
c           2 to print full resolution to file 17 (recomputed on each
c              iteration.)
c           3 to calculate full resolution on 1st iteration only.
c           0 no resolution calculations.
c        i3d - flag for using pseudo-bending:
c              -1 = only use straight path 
c              0=no pseudo-bending
c              1=use in forward problem to compute velocity partial derivatives
c              2=also use in earthquake location subroutine
c              3=pseudo-bending; also use less curvature for initial arcuate
c              3=pseudo-bending; also test more squashed paths (3rd root arc)
c                that may be better initial paths for long ray paths (moho-like).
c                rays below "moho".  Assumes last z grid (k=nz-1) is "moho".
c        nitmax - max of iterations of the velocity inversion-hypocenter
c                 relocation loop.
c                 For locations only, set nitmax=0.
c                 For creating synthetic data, set nitmax= -1, and
c                   use file007 for hypocenters and stations.
c        snrmct - cutoff value for solution norm. Program will stop
c                 iterating if the solution norm is less than snrmct.
c        ihomo - flag for using 1-d starting model (1=yes,0=no)
c               if flagged, on first ihomo iterations do 2-d art
c        rmstop - overall rms residual for termination of program.
c        ifixl - number of velocity inversion steps to keep hypocenters
c                fixed at start of progressive inversion
c
c      line 8 - (free format)
c        delt1, delt2 - distance weighting factors. The weight is 1 for
c               x<delt1, but tapers linearly to 0 between delt1 and delt2.
c                Note that this is epicentral distance in simulps13
c        res1,res2 - same pattern as above, but for residual weighting.
c        res3 = allowing some high residual data during location iterations 
c             downweighting(linear) 0 to 98% res1 to res2, 98 to 100% res2 to res3
c
c      line 9 - (free format)
c        ndip - of rotation angles of the plane of the ray, which will
c               be computed in the exhaustive search for the fastest time.
c        iskip - of rotation angles which will be skipped.
c                ndip=9, iskip=3 will give a vertical plane, and 2 swung
c                at angles of 22.5 degrees on each side.
c                ndip=9, iskip=4 will give only the vertical plane, and should
c                be used for 1 dim models, to save computer time.
c           *** Most of the computer time is spent in the raytracing. Careful
c               selection of the starting model and use of iskip can save
c               a lot.
c         scale1  - set scale1 to the step length for the travel
c                           time computation. Set no larger than the grid
c                           spacing.
c         scale2 - scale for the number of paths tried in the raytracing.
c                  Cliff uses a value of 1, and this seems ok, but would
c                  need to be tested in detail. If scale2 is smaller, the
c                  number of paths increases, and the computation time also
c                  goes up.
c
c      line 10 - (free format)
c        xfac - Convergence enhancement factor for pseudo-bending
c               Cliff suggests 1.2 to 1.5
c        tlim - Cutoff value for travel time difference to terminate
c               iteration. (0.0005 to 0.002 s)
c               Make lower tlim for t* (iuseq=1) than for t.
c        nitpb - maximum permitted iterations for pseudo-bending (5 to 10)
c                nitpb(1) for shorter raypaths
c                nitpb(2) for raypaths > delt1
c              
c      line 11 - (free format)
c        kttfor - flag for travel-time data format
c               1=cnsp, 2=phs, 3=6-letter-sta phs(Japanese), 4=stn5 and net (ncedc)
c        iusep - flag to tell whether or not to use P arrivals (and invert for
c                P velocities).  0 = no, 1 = yes.
c        iuses - flag to tell whether or not to use S arrivals (and invert for
c                S velocities).  0 = no, 1 = yes.
c        invdel - flag to control inclusion of station delays in the inversion.
c                 =0 to not include stn delays.
c                 =1 to include stn delays.
c        iuse2t = 1 to allow 2 origin times (PANDA array with different master clock)
c
c      line 12 - (free format)
c        iuseq - flag to tell whether inversion is for velocity or Q
c                0 = velocity,  1 = Q 
c        dqmax - maximum Q adjustment allowed per iteration
c        qdamp - damping parameter for Q inversion
c                actually damps solution for Q*V
c        qmin -  minimum Q value for solution (prior versions just had 1.0)
c        qrmax = maximum Q estimate allowed to discard data
c               using simple slant path with t* data
c
c     File 02 - Station Data  (read in subroutine INPUT2)
c        line 1 - (free format)
c         ltdo, oltm, lndo, olnm, rota, nzco, cmerid - sets origin and rotation of the
c              coordinate system. Choose the origin at the lower right corner
c              of the region.
c              Y points in North direction.
c              X points West, then orig lat should be negative for east in simul
c              rota - angle of rotation counterclockwise (degrees). This is
c                     used to rotate the entire coordinate system.
c              nzco - select coordinate transformation
c                nzco= 0 use origin latitude for central meridian
c                nzco= 1 New Zealand (lat S and long E), use NZTM2000
c                nzco= 2 Alaska, use state plane zone AK4 for central Alaska
c                nzco= 3 use input central meridian (this could be selected from UTM)
c                nzco= 4 use NZ map grid as in previous simul versions. distances are
c                        similar to NZTM2000, but some rotation especially offshore.
c                nzco= 5 use short distance conversion as in Thurber original simul
c              cmerid - central meridian for coordinate conversion, used for nzco=3
c                       Note UTM system has 6 deg wide zones, so TM okay for about 700km wide
c
c         line 2 - (free format)
c         nsts - number of stations in station list to follow.
c
c         station list: see documentation.
c
c     File 03 - velocity model (read in subroutine INPUT3)
c         This contains the number of nodes in the x, y , and z directions.
c         Then the velocity model. See the documentation for a complete
c         description.  In simulps12 this has Vp, followed by Vp/Vs.
c         For iuseq=1, has Vp, followed by Qp.
c         Specify which nodes (if any) you want to fix velocity for:
c           2  3  4   means x2(2nd column), y3(3rd row), z4(4th layer)
c         then linked nodes
c           master i j k. linking type 1=linear, 2=constant.  slaves i j k 
c
c     File 04 - travel-time data for earthquakes. (read in subroutine INPUT4)
c         Note that S data should be made into S-P data.
c         Use program 'convert6' to convert hypo71 summary and phase data
c         to this format.
c         (should still work with old 'convert3' files also, as well as
c         Michelini's format)
c         kttfor in input1 allows phs or Japanese format
c               1=cnsp, 2=phs, 3=6-letter-sta phs(Japanese)
c
c     File 07 - travel-time data for shots
c         For iuseq=1, read in t* data as shots.
c         (For iuseq=1, neqs and nbls should be zero.)
c
c     File 08 - travel-time data for blasts
c
c     File 09 - travel-time data for clusters
c         Eq groups start "BEGIN CLUSTER", end "END CLUSTER"
c
c  OUTPUT FILES:
c     File 16 - printed output, send to line printer
c     File 12 - hypocenters from each iteration, hypoinverse format with error ellipse
c     File 13 - Final hypocenters in hypo71 summary card format
c     File 15 - Written when kout3=1, A separate file for each event
c         containing all the raypath points, useful for plotting
c     File 17 - Full resolution matrix.  Created if ires.ge.2.
c     File 18 - Pointer from full nodes to inversion nodes.
c          Output when fixed nodes are used.
c     File 19 - Written when kout3=1, Contains the difference in travel 
c          time between ART and pseudo-bending
c     File 20 - Station residual output (created for kout2=0 or 1)
c     File 21 - Residuals for receiver-pair differential times 
c     File 22 - Station data with new P and S delays.  This
c          can be used as input to future runs. NOT created if
c          invdel=0 (no station delays inverted for).
c     File 23 - Final velocity model.  This can be used as input
c          to future runs.
c     File 24 - Earthquake travel-time data for new hypocenters.
c          This can be used as input in future runs.
c     File 25 - New velocities in station format so can be plotted 
c          with qplot.
c     File 26 - List of observations that used maximum allowed 
c          pseudo-bending iterations (nitpb)
c     File 28 - Blast travel-time data with new origin times.
c     File 29 - Revised travel-time data for shots
c     File 30 - Cluster Eq: Station residual output (for kout2=0 or 1)
c     File 31 - Residuals for earthquake-pair differential times
c     File 33 - Final hypocenters as in simul printout with x,y,z
c     File 34 - Output file similar to HYPO71 listing file, which can
c          be used as input to FPFIT fault plane solution program.
c          Only for earthquakes and shots since called from Loceqk.
c          (Note that DIST is the hypocentral distance, whereas HYPO71
c          outputs the epicentral distance.)  Written on last (nitmax) iteration.
c     File 36 - Summary output file that contains key solution
c          statistics.  This is useful when doing numerous
c          damping runs and wanting to compare variances.
c     File 43 - Derivative weighted sum from ttmder for all nodes, including
c          fixed and linked.  So this is NOT the DWS actually used in the
c          inversion, but rather an approximation of the data distribution
c          over all velocity nodes.  This is useful for deciding how to set
c          fixed and linked nodes.
c          Only created for nitmax=1 (since use this on preliminary test runs). 
c     File 44 - DWS in print format
c     File 45 - 1/diagonal elements of covariance matrix.  In same
c          format as velocity model input.  Created if ires.ge.2.
c     File 53 - velocity points in Table format with x y z Latitude longitude
c     File 59 - diagonal resolution element in table format
c   HINTS:
c       1. Invert for a one-dimensional model first. Use nx=3, ny=3, nz=
c          whatever you want in the final model. Set ndip=9, iskip=4.
c          Or, better yet, use VELEST, locate your events with the VELEST
c          model, then fix the earthquake locations on the 1st iteration.
c
c       2. Run the inversion on calculated data. This will tell you what
c          kind of averaging is inherent in the inversion process, since
c          the result will probably be different from your original input
c          model.
c
c       3. When the program crashes due to divide by 0, etc, the cause can
c          most often be traced to errors in the setup of the velocity
c          model. No part of a ray must reach over halfway to an outer node.
c          Try increasing the distance of the outer node.
c
c  ARRAY DIMENSIONING:
c      11-mar-94 These notes are now in the simul2014_common.inc file.
c      PLEASE READ that file to understand the array dimensions
c      and how many parameters you are allowed to invert for.

c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c  some additional changes by f.haslinger, eth zurich 
c  15-june-98
c  - incorporate changes (from e.kissling) to deal with N/S , E/W
c    problem. requires correct coordinate qualifier (N/S or E/W) in
c    station and event coordinates and uses cns, cew to get correct
c    values (N, W are positive; S, E are negative - U.S centered world...)
c 
c  Record of changes I have made in program (DMEP).
c  22-feb-17 Use 2*qrmax to check for bad teleseismic tstar
c 21-sep-16 outhit had extra 0.0 row for dws, after 436 continue
c 08-feb-16 c corrected input9 for reading tstar for nz data (not kttfor=4)
c  20-Jul-15 Add kttfor=4 for phase data that has 5-letter station and net (ncedc)
c  09-jun-15 changed OUTLIS to have format okay for 6-letter station (kttfor=3)
c  31-mar-15 have vlxyxltln.out file 53 be like table with vp and vpvs on same line
c      Have file 59 be diag resol elem in velocity file format, resDRE.out
c  05-Feb-14 Have synthetic theoretical without error written to hypolist file
c     file 34: hypo71list
c  09-Dec-14 Added error to synthetic times (not t*), use rderr and iquality
c     changes to OUTEND, added RNORMAL
c  22-apr-14 have S-P related to both dvp and dvpvs, changes to ttmder, ttmderce
c  25-jan-13 iuserdt=1 for eq only, =2 for rdt on shots and blasts also.
c  11-jun-12 put in facrayr and facraye instead of 1.5
c  13-feb-12 Change to getcedt to get S-P obs before P for a given event pair
c  21-oct-09 have rdlta array to have slant distance
c  24-jan-06 Added search of more squashed arc paths in rayweb to help
c     with longer ray paths.  Added subroutine cmdsvm.
c  17-dec-03 Put in fix to BEND from Cliff Thurber
c  05-jul-00  Print out to file 38, a list of event,station,tt for use in Woodward program
c             if kout2=0 and nitmax=0
c  21-Sep-99  Now allow 2 origin times.  For the NZ Panda experiment, we have
c     data from both the permanent network and a telemetered temporary array which
c     has a bad master clock.  Thus we solve for 2 origin times (real time and
c     temp master clock).  Many changes.  Input added iuse2t and iclock (sta)
c     so that all travel-time to certain stations are used to solve for seco2.
c  09-Feb-99 Output file 44 which has layers, y-grids and x-grids
c     of dws. This is useful format for scanning to decide links.
c  18-Jan-99 For nzco=1, use New Zealand Map Grid coordinates for
c     all lat-lon to cartesian conversions.
c     Note that rota is about 2 deg different for S. Island
c     northing compared to short distance conversions
c     Also note that with the ETH now using S and N latitude,
c     the rota must be plus 180 deg for NZ from simulps13.
c  05-Jan-99 Have hit(DWS) and hitall (DWS all nodes) now include
c      weight factor.  Have combined weight (as used for inversion)
c      now saved to array wtcomb(maxobs,maxev).
c  9-dec-98  Allow S-P data input when char2,3 of phase remark are 'Sp'
c 11-mar-94 Changes throughout program to have common blocks
c    in a separate file that is included: simul2014_common.inc
c  03-jun-93 Change to FORWRD.  Corrected way synthetic data was computed for
c      S-P observations.
c  5-jan-93 Changes to VELADJ to output sol norm and damping
c      for both vp and vp/vs.
c      Added "wtsp" to allow user to change relative weighting of
c      S-P observations in hypocentral solution.  Changes to WTHYP.
c  30-dec-92 Made changes so that vpvs array is used directly. (Cliff
c       had solved for vp/vs but then used the vp/vs perturbation in
c       veladj to perturb the associated vs element.)  Now vp/vs is 
c       input instead of vs, and vs is calculated from vp and vp/vs
c       in input3 and at end of veladj.
c  5-nov-92 Now have option to create synthetic travel-time data.
c           Use nitmax= -1, use file007 for input hypocenters and
c           stations.  Calculated travel-times will be output to
c           file028.  Made changes to Main and Outend.
c  24-feb-92 Added "bld" to input3.  Bld is factor for setting up velocity
c      interpolation arrays (1.0 or 0.1). Now velocity nodes can be defined to 1/10th km
c      if desired.  Changes to BLDMAP, INPUT3, INTMAP, OUTEND.
c  7-sep-90 Write out 1/covariance(diag) if kout2.eq.5 and ires.gt.0
c      to for016 and for045; changes for OUTEND, RESCOV.
c  12-jul-90 Put in a second linear residual weighting, 98% of downweighting
c      is done res1 to res2, last 2% of downweighting is done 
c      res2 to res3.  This is useful in case you happen to get a poor
c      hypocenter with mostly very high residuals.  If res3=res2,
c      100% of downweighting is done res1 to res2, as before.
c      Remove "inew" from MAIN since it's never used.
c      In MAIN, if nitmax.eq.0, skip parsep for blasts & shots also.
c  11-apr-90 Allow ires=3 to only compute resolution on 1st iteration.
c      Saves cpu and is appropriate when idmp=1(damping increases,
c      resolution decreases on later iterations.)
c      Changes to MAIN.
c  4-oct-89 Can now invert for station corrections at selected stations
c      by fixing delays (to 0) at other stations.  Changes to common
c      OBSERVE, routines INPUT2, INPUT3, TTMDER, VELADJ.
c  22-sep-87 Made a variety of changes to make the program work properly
c     when station delays are included in the inversion (also for the
c     case when all velocity gridpoints are fixed).   
c
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c
c  Main Program
c
c  declaration statements:
      character*26 dash
c
c  common block variables:
      include 'simul2014_common.inc'
c
      dash='--------------------------'
c
c
c   open output files
        open(unit=16,file='output')
        open(unit=56,file='TEST_PRINT')
        rewind (16)
        open(unit=36,file='itersum')
        rewind (36)
cfhdmep
        open(unit=18,file='nodes.out')
c    file to check pseudo-bending
        open(unit=26,file='pseudo')
        rewind (26)
c
c  input routines
c
c  input control parameters
      open(unit=01,status='old',file='CNTL',form='formatted')
      rewind (01)
      call input1
      close(01)
c  open resolution file if ires=2
      if(ires.ge.2) open(unit=17,file='resol.out')
c  open residual output file if kout2=0,1
      open(unit=20,file='residuals')
      if((kout2.le.1).and.(iuserdt.gt.0)) open(21,file='rdtres.out')
      if(iuseclu.gt.0) open(30,file='ceres.out')
      if(iuseclu.gt.0) open(31,file='edtres.out')
c  open file to write out event,station,tt list if kout2=0 and
c  nitmax=0
      if((kout2.eq.0).and.(nitmax.eq.0)) then
        open(unit=38,form='formatted')
        write(38,3801)
 3801   format(' Event Latitude Longitude   Depth  O-T  Sta ',
     2  'Latitude Longitude Depth Delay  Travel-Time')
      endif
c  open files to write raypath points, tt differences to if kout3=1
      if(kout3.eq.0) goto 70
      open(unit=19,file='ttdiff.out')
      write(19,1901)
 1901 format('  ne  stn  delta  fstime  ttime   tdif')
c  initializations
   70 call strt(0)
      nit=0
c  input station list, set up center of coordinates, calculate
c  cartesian coordinates
      open(unit=02,status='old',file='STNS',form='formatted')
      rewind (02)
      call input2
      close(02)
c  input medium model
      open(unit=03,status='old',file='MOD',form='formatted')
      rewind (03)
      call input3
      close(03)
      call outadj(nit,0,0)
      if(neqs.eq.0) goto 71
      open(unit=04,status='old',file='EQKS',form='formatted')
   71 if(nsht.eq.0) goto 72
      open(unit=07,status='old',file='SHOT',form='formatted')
   72 if(nbls.eq.0) goto 62
      open(unit=08,status='old',file='BLAS',form='formatted')
   62 if(ntel.eq.0) goto 73
      open(unit=10,status='old',file='TELE',form='formatted')
   73 if(iuseclu.eq.0) goto 74
      open(unit=09,status='old',file='CLUS',form='formatted')
c
   74 if(kout.lt.3) goto 75
      open(unit=34,file='hypo71list')
      rewind(34)
      if(kout.lt.4) goto 75
c       open file for hypocenter summary card output on each iteration
      open(unit=12,file='summary')
        rewind (12)
c
   75 istop=0
      istop1=0
c
      neb=neqs+nbls
      nevt=neb+nsht
c
c  iterative inversion loop
c
    1 continue
c
      netemp=neqs
      nbtemp=nbls
      nc=0
c  if flagged, on first ifixl iterations, treat all quakes as blasts
      if (ifixl.le.nit) go to 110
      nbls=nbls+neqs
      neqs=0
  110 continue
c
      istemp=iskip
      ndtemp=ndip
c  if flagged, on first ihomo iterations do 2-d art
      if (ihomo.le.nit) go to 111
      iskip = 0
      ndip = 1
  111 continue
c
c
      write(16,1005) nit
 1005 format(///,' iteration step',i3,'; hypocenter adjustments')
      if (neqs.eq.0) go to 11
c
      nrdtsh=0
c  loop over all earthquakes
      ne=1
    9 continue
c  input observations of event
      if((nit.eq.0).or.(kout2.eq.0)) write(16,130)
     2 dash,dash,dash,dash,dash
  130 format(1x,5a26)
      if(ne.gt.neqs) goto 10
      if (nit.eq.0) call input4(ne,4)
c  locate individual earthquake to reduce residuals
      call loceqk(ne,nit,nwr)
      if((nit.eq.0).and.(nwr.lt.4)) goto 9
      if((nit.gt.0).and.(nwr.lt.4)) goto 9998
c  skip if only locating
      if((nitmax.le.0).and.(kout3.eq.0).and.(kout2.gt.1)) goto 208
c  do forward problem and calculate partial derivatives
      call forwrd(ne,nit)
      if((nitmax.le.0).and.(kout2.gt.1)) goto 208
c  create receiver-pair shot
      if(iuserdt.gt.0) then
        if(nit.eq.0) call getrdtsht(ne)
        if(kobsrdt(ne).eq.0) goto 207
        call forwrdrdt(ne,nit)
        if(kobsrdt(ne).gt.0) nrdtsh=nrdtsh+1
      endif
c  perform parameter separation
  207 call parsep(nc,ne,1,nwr)
c  if last iteration, write out station residuals
  208 if((nit.eq.nitmax).and.(kout2.lt.2)) call outres(ne)
      if((nit.eq.nitmax).and.(kout.ge.3)) call outlis(ne)

c  skip this event if not enough good readings,
c        stop if not first iteration
      if((nit.gt.0).and.(nwr.lt.4)) goto 9998
c  for receiver-pair shot, get derivatives and set up matrix
      if((iuserdt.eq.0).or.(nwr.lt.4)) goto 8
      if(kobsrdt(ne).eq.0) goto 8
      call medder(nc,ne,2,nwrrdt)
      call parsep(nc,ne,2,nwrrdt)
      if((nit.eq.nitmax).and.(kout2.lt.2)) call outresrdt(ne)
    8 if(nwr.ge.4) ne=ne+1
      if(ne.le.neqs) goto 9
   10 continue
c
c  loop over all blasts
   11 if(nbls.eq.0) goto 15
      write(16,1613) nbls
 1613 format(/,2x,'Following ',i5,' Events are Blasts with Unknown',
     2 ' Origin Time')
      nb=1
   12 ne=nb+neqs
   13 infile=8
      if((nit.eq.0).or.(kout2.eq.0)) write(16,130)
     2 dash,dash,dash,dash,dash
      if((ifixl.gt.nit).and.(ne.le.netemp)) infile=4
      if(nit.eq.0) then
        call input4(ne,infile)
        if(ne.gt.neb) goto 14
        if((ne.gt.netemp).and.(infile.eq.4)) goto 13
      endif
      call loceqk(ne,nit,nwr)
c should have 2 observations for blast
      if((nit.eq.0).and.(nwr.lt.2)) goto 13
      if((nit.gt.0).and.(nwr.lt.2)) goto 9998
      if((nitmax.le.0).and.(kout3.eq.0).and.(kout2.gt.1)) goto 308
      call forwrd(ne,nit)
      if((nitmax.le.0).and.(kout2.gt.1)) goto 308
c  create receiver-pair shot
      if(iuserdt.eq.2) then
        if(nit.eq.0) call getrdtsht(ne)
        if(kobsrdt(ne).eq.0) goto 307
        call forwrdrdt(ne,nit)
        if(kobsrdt(ne).gt.0) nrdtsh=nrdtsh+1
      endif
  307 call parsep(nc,ne,1,nwr)
  308 if((nit.eq.nitmax).and.(kout2.lt.2)) call outres(ne)
      if((nit.eq.nitmax).and.(kout.ge.3)) call outlis(ne)
c  for receiver-pair shot, get derivatives and set up matrix
      if(iuserdt.ne.2) goto 50
      if(kobsrdt(ne).eq.0) goto 50
      call medder(nc,ne,2,nwrrdt)
      call parsep(nc,ne,2,nwrrdt)
   50 if(nwr.ge.2) nb=nb+1
      if(nb.le.nbls) goto 12
   14 continue
c
c  loop over all shots
   15 continue
      if (nsht.eq.0) go to 25
      write(16,1614) nsht
 1614 format(/,2x,'Following ',i5,' Events are Shots with Known',
     2 ' Origin Time')
      ns=1
c  input observations of shot
   16 ne=ns+neqs+nbls
      if((nit.eq.0).or.(kout2.eq.0)) write(16,130)
     2 dash,dash,dash,dash,dash
      if (nit.gt.0) goto 17
      call input4(ne,7)
      if(ns.gt.nsht) goto 20
c  do forward problem and calculate partial derivatives
      if(nit.eq.0) then
        jfl=0
      else
        jfl=2
        if(ihomo.eq.nit) jfl=1
      end if
   17 call forwrd(ne,nit)
c  create receiver-pair shot
      if(iuserdt.eq.2) then
        if(nit.eq.0) call getrdtsht(ne)
        if(kobsrdt(ne).eq.0) goto 407
        call forwrdrdt(ne,nit)
        if(kobsrdt(ne).gt.0) nrdtsh=nrdtsh+1
      endif
c  add medium derivatives from shot to medium matrix
  407 call medder(nc,ne,1,nwr)
      if((nit.eq.0).and.(nwr.lt.1)) goto 16
      if((nit.gt.0).and.(nwr.lt.1)) goto 9998
c** Change for synthetic data, nitmax= -1
      if((nitmax.le.0).and.(kout2.gt.1)) goto 408
c  put derivatives into g matrix
      call parsep(nc,ne,1,nwr)
  408 ns=ns+1
      if((nit.eq.nitmax).and.(kout2.lt.2)) call outres(ne)
      if((nit.eq.nitmax).and.(kout.ge.3)) call outlis(ne)
      if(nitmax.eq.-1) call outlis(ne)
c  for receiver-pair shot, get derivatives and set up matrix
      if((nitmax.eq.-1).or.(iuserdt.ne.2)) goto 18
      if(kobsrdt(ne).eq.0) goto 18
      call medder(nc,ne,2,nwrrdt)
      call parsep(nc,ne,2,nwrrdt)
   18 if(ns.le.nsht) goto 16
   20 continue
c
   25 continue
      nebs=neqs+nbls+nsht
c  loop over all teleseismic events
      if(ntel.eq.0) go to 35
      write(16,1615) ntel
 1615 format(/,2x,'Following ',i5,' Events are Teleseismic, use path',
     2 ' from piercing point, solve for Tele path delay')
      nt=1
c  input observations of tele
   26 ne=nt+nebs
      if((nit.eq.0).or.(kout2.eq.0)) write(16,130)
     2 dash,dash,dash,dash,dash
      if (nit.gt.0) goto 27
      call input10(ne,10)
      if(nt.gt.ntel) goto 20
c  do forward problem and calculate partial derivatives
      if(nit.eq.0) then
        jfl=0
      else
        jfl=2
        if(ihomo.eq.nit) jfl=1
      end if
   27 call forwrd(ne,nit)
c  add medium derivatives from shot to medium matrix
      call medder(nc,ne,5,nwr)
      if((nit.eq.0).and.(nwr.lt.ntobmin)) goto 26
      if((nit.gt.0).and.(nwr.lt.ntobmin)) goto 9998
c** Change for synthetic data, nitmax= -1
      if((nitmax.le.0).and.(kout2.gt.1)) goto 409
c  put derivatives into g matrix
      call parsep(nc,ne,1,nwr)
  409 nt=nt+1
      if((nit.eq.nitmax).and.(kout2.lt.2)) call outres(ne)
      if(nt.le.ntel) goto 26
   30 continue
c
   35 continue
      if(iuseclu.eq.0) goto 80
c  loop over all clusters
      nc=1
   55 continue
c  input observations of cluster
      if((nit.eq.0).or.(kout2.eq.0)) write(16,130)
     2 dash,dash,dash,dash,dash
      if(nc.gt.nclu) goto 76
      if (nit.eq.0) call input9(nc,9)
      if((nit.gt.0).and.(iuseq.eq.1)) goto 61
c  locate earthquake cluster, including cedt; getcedt for nit=0
      call locclu(nc,nit,nwrce)
      nec=ncev(nc)
      necpar=nec*4
      if(nwrce.lt.necpar) goto 55
c loop through earthquakes in cluster 
   61 do 65 ne=1,nec
c  skip if only locating
        if((nitmax.le.0).and.(kout3.eq.0).and.(kout2.gt.1)) goto 218
c  do forward problem and calculate partial derivatives
        call forwrdce(nc,ne,nit)
        if((nitmax.le.0).and.(kout2.gt.1)) goto 218
c  put derivatives into g matrix
        call medder(nc,ne,3,nwre)
        call parsep(nc,ne,3,nwre)
c  if last iteration, write out station residuals
  218   if((nit.eq.nitmax).and.(kout2.lt.2)) call outresce(nc,ne)
c don't bother with outlis for cluster eq
c        if((nit.eq.nitmax).and.(kout.ge.3)) call outlis(ne)
   65 continue
c  for earthquake-pair diff-times, get derivatives and set up matrix
      if(kobsedt(nc).eq.0) goto 68
CCC NOT SURE OF OUTRES  WITH NITMAX=0 -CHECK LATER- dmep 15feb10
      if(nitmax.le.0) goto 220
      call forwrdedt(nc,nit,kobsok)
      if(kobsok.eq.0) then
        write(16,1668) nc,kobsedt(nc)
 1668   format('*** after ttmderce Cluster:',i4,', kobsedt=',i6,
     2   '; ALL nnodej=0 : NO EDT For Medium-Derivatives ***')
        kobsedt(nc)=0
        goto 68
      endif
      call medder(nc,ne,4,nwredt)
      call parsep(nc,ne,4,nwredt)
  220 if((nit.eq.nitmax).and.(kout2.lt.2)) call outresedt(nc)
   68 if(nwrce.ge.necpar) nc=nc+1
      if(nc.le.nclu) goto 55
   76 continue
c
c** Change for synthetic data, nitmax= -1
      if(nitmax.le.0) goto 9999
c Print out dws for all nodes if nitmax=1
   80 if(nitmax.eq.1) call outhit(nit)
c  continue iterating or terminate?
      call decide(istop,nit)
c     if (nit.gt.0) call decide(istop)
        if(istop1.gt.2) goto 9999  !If backed-up twice, don't adjust again,end
        if(istop.eq.2) goto 28     !If variance ratio bad, backup
        if(istop1.eq.2) goto 9999  !If backed-up once and now okay, end
      if(nit.ge.nitmax) goto 9999  !As in CT version, do hyp last
c
 900  continue
      nit=nit+1
      write(16,1015) nit
 1015 format(///,' iteration step',i3,'; simultaneous inversion')
c  invert for velocity model adjustments
      call veladj(nit)
      if(nit.eq.99) goto 9999
c  output results of iteration
      call outadj(nit,istop,istop1)
c  compute resolution?
      if(ires.eq.0) goto 901
      if((ires.ne.3).or.(nit.eq.1)) call rescov
  901 if(istop.eq.1) go to 9999
      if(istop.eq.0) goto 34
c
   28   write(16,902)
        write(36,902)
  902   format(' ******* Variance Ratio is less than Critical Ratio',
     2 /,12x,' ******* Backup Parameters Halfway *******',/)
        call velbku(istop1)
        istop1=istop1+istop
        call outadj(nit,istop,istop1)
c
   34 continue
c
c  restore neqs,nbls,iskip,ndip
      neqs=netemp
      nbls=nbtemp
      iskip=istemp
      ndip=ndtemp
c
      call strt(nit)
      rewind(26)
c
      go to 1
c
 9998 continue
      write(16,1620)
 1620 format(//,'  ********** STOP **********')
 9999 continue
c      WRITE(6,1621) nit,istop
      write(16,1621) nit,istop
 1621 format('* * * *',/,'  Starting outend, nit=',i3,' istop=',i2)
c
      call outend
c
      close(04)
      close(07)
      close(08)
      close(09)
      close(12)
      close(16)
      close(26)
      close(36)
      close(38)
      if(ires.eq.2) close(17)
      if(kout3.eq.0) stop
      close(19)
      stop
c***** end of main program *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine avsd(cnull,x,nx,sd,av,devtot)
c  program to find the average and standard deviation of a 
c  list of numbers.  Those with value cnull are not included.
c
      real x(15000)
      sum=0
      i=0
      do 50 ix=1,nx
         if(x(ix).eq.cnull) goto 50
         sum=sum+x(ix)
         i=i+1
   50 continue
      nx1=nx
      nx=i
      if(nx.eq.0) goto 800
      av=sum/float(nx)
      devtot=0
      do 260 i=1,nx1
         if(x(i).eq.cnull) goto 260 
         dev=x(i)-av
         devtot=devtot+dev*dev
  260 continue
      sd=sqrt(devtot/float(nx))
c     write(6,620) nx,av,sd
  620 format(' nx=',i5,', average=',e11.4,', sd=',e11.4)
      return
  800 continue
      sd=0.00
      devtot=0.00
      av=0.00
      nx=0.00
      return
c ***** end of subroutine avsd *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine aztoa(x1,x2,y1,y2,z1,z2,xr,yr,azim,tkofan)
c  this subroutine computes the azimuth and take-off-angle
c  for an individual observation
c  (called from Loceqk)
c  pr is station, p1 and p2 are hypocenter and adjoining point
c  on the raypath
c
      common/shortd/ xltkm,xlnkm,rota,nzco,xlt,xln,snr,csr,
     2 onorth,oeast,rotadg,cmerid,rlatdo,rlondo
      parameter (drad=1.7453292d-02)
      parameter (pi=3.1415926536)
      parameter (twopi=6.2831853072)
c
c  Azimuth
c
c  start cht 1998
c
      xd=x2-x1
      yd=y2-y1
c
c  end cht 1998
c
      xda=abs(xd)
      yda=abs(yd)
      phi=atan(xda/yda)
c  compute correct azimuth depending on quadrant
      if(xd.ge.0.0) then
        if(yd.ge.0.0) then
          theta=twopi-phi
        else
         theta=pi+phi
        endif
      else
        if(yd.ge.0.0) then
          theta=phi
        else
          theta=pi-phi
        endif
      endif
c  rotate back to real north, convert to degrees
      azim=(theta-rota)/drad
      if(azim.lt.0.0) azim=azim+360.0
c
c  Take-off-angle
      xd=x2-x1
      yd=y2-y1
      r=sqrt(xd*xd+yd*yd)
      zd=z2-z1
      zda=abs(zd)
      phi=atan(r/zda)
      if(zd.lt.0.0) phi=pi-phi
      tkofan=phi/drad
c
      return
c ***** end of subroutine aztoa *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine bend(isp,xfac)
c*****this routine perturbs the initial path in the direction
c      of the normal to the ray path tangent at each point
c      by the optimal distance r
c
      common/pathm/x(260),y(260),z(260),v(260),vq(260),tra,qtra,n,nn
      common/temp/xtemp(260),ytemp(260),ztemp(260),rtemp(260),ttemp(260)
c
c ***
      xtemp(1)=x(1)
      ytemp(1)=y(1)
      ztemp(1)=z(1)
c ***
      do 200 k=2,nn
c
         kk=k-1
         kkk=k+1
c
c*****compute the normal direction of maximum gradient of velocity
c
         dx=x(kkk)-xtemp(kk)
         dy=y(kkk)-ytemp(kk)
         dz=z(kkk)-ztemp(kk)
         dn=dx*dx+dy*dy+dz*dz
         ddn=sqrt(dn)
         rdx=dx/ddn
         rdy=dy/ddn
         rdz=dz/ddn
c
         xk=0.5*dx+xtemp(kk)
         yk=0.5*dy+ytemp(kk)
         zk=0.5*dz+ztemp(kk)
c ***
         call vel3eft(isp,xk,yk,zk,vk)
         call veld(isp,xk,yk,zk,vx,vy,vz)
c
c ***
         vrd=vx*rdx+vy*rdy+vz*rdz
         rvx=vx-vrd*rdx
         rvy=vy-vrd*rdy
         rvz=vz-vrd*rdz
c
         rvs=sqrt(rvx*rvx+rvy*rvy+rvz*rvz)
         if(rvs.eq.0.0) goto 200
         rvx=rvx/rvs
         rvy=rvy/rvs
         rvz=rvz/rvs
c
c*****compute the optimal distance r
          rcur=vk/rvs
c  BUG FIX OCTOBER 2003
          if (rcur*rcur.gt.0.25*dn) then
          rtemp(k)=rcur-sqrt(rcur*rcur-0.25*dn)
          else
          rtemp(k)=rcur
          endif
c
c*****compute the new points and distance of perturbations
c
         xxk=xk+rvx*rtemp(k)
         yyk=yk+rvy*rtemp(k)
         zzk=zk+rvz*rtemp(k)
c
c  convergence enhancement
         xxk=xfac*(xxk-x(k))+x(k)
         yyk=xfac*(yyk-y(k))+y(k)
         zzk=xfac*(zzk-z(k))+z(k)
c
         ttemp(k)=sqrt((x(k)-xxk)**2+(y(k)-yyk)**2+(z(k)-zzk)**2)
         xtemp(k)=xxk
         ytemp(k)=yyk
         ztemp(k)=zzk
         call vel3eft(isp,xxk,yyk,zzk,vk)
         v(k)=vk
         call vel3eft(1,xxk,yyk,zzk,vkk)
         vq(k)=vkk
c ***
200   continue
c
      return
c ***** end of subroutine bend *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine bldmap
c  common block variables:
      include 'simul2014_common.inc'
c
c     array size limits
c
c     write(6,400)
c 400 format(' subroutine bldmap')
      xl=bld-xn(1)
      ixmax=(xn(nx)+xl)/bld
      yl=bld-yn(1)
      iymax=(yn(ny)+yl)/bld
      zl=bld-zn(1)
      izmax=(zn(nz)+zl)/bld
c     write(6,402)ixmax,iymax,izmax
c 402 format(' array sizes: ',3i5)
c
c  Check for array size overflow
      if(ixmax.gt.ixkms.or.iymax.gt.iykms.or.izmax.gt.izkms)goto 330
      ix=1
      do 10 i=1,ixmax
c
         ix1=ix+1
c
         xnow=float(i)*bld-xl
         if (xnow.ge.xn(ix1)) ix=ix1
c
         ixloc(i)=ix
   10 continue
c  Fill remainder of array with zeroes.
      do 12 i=ixmax,ixkms
         ixloc(i)=0
   12 continue
c
c
      iy=1
      do 15 i=1,iymax
c
         iy1=iy+1
c
         ynow=float(i)*bld-yl
         if (ynow.ge.yn(iy1)) iy=iy1
c
         iyloc(i)=iy
   15 continue
c
c  Fill rest of array with zeroes.
      do 17 i=iymax,iykms
         iyloc(i)=0
 17   continue
c
      iz=1
      do 20 i=1,izmax
c
         iz1=iz+1
c
         znow=float(i)*bld-zl
         if (znow.ge.zn(iz1)) iz=iz1
c
         izloc(i)=iz
   20 continue
c
c  Fill remainder of array with zeroes.
      do 22 i=izmax,izkms
         izloc(i)=0
  22  continue
      return
 330   continue
      write(16,331)ixkms,iykms,izkms
 331  format(' ***** error in array size in common/locate/',/,
     *' maximum map dimensions (km)=',/,' x=',i5,' y=',i5,' z=',i5)
      write(16,332)ixmax,iymax,izmax
  332 format(' Actual map size (km): ',/,' x=',i5,' y=',i5,' z=',i5)
      stop
c***** end of subroutine bldmap *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine cmpdpv(xe,ye,ze,xr,yr,zr,scale2,ndip,dipvec)
c
c  parameters
      real xe,ye,ze,xr,yr,zr,scale2,dipvec(3,9)
c
      integer ndip
c  local variables
      real dx,dy,dz,xh1,yh1,zh1,xh2,yh2,zh2,size,xv,yv,zv,
     *     rescal,x451,y451,z451,x452,y452,z452
c
      integer nv
c
      dx=xr-xe
      dy=yr-ye
      dz=zr-ze
c
c  near-vertical vector
      xv=-dx*dz
      yv=-dy*dz
      zv=dx*dx+dy*dy
c  rescale vector to length scale2
      size=sqrt(xv*xv+yv*yv+zv*zv)
      rescal=scale2/size
c
      xv=xv*rescal
      yv=yv*rescal
      zv=zv*rescal
c
c  store this vector
      nv=(ndip+1)/2
      dipvec(1,nv)=xv
      dipvec(2,nv)=yv
      dipvec(3,nv)=zv
c
      if (ndip.eq.1) return
c
c  horizontal vectors
      xh1=dy
      yh1=-dx
      zh1=0.0
      xh2=-dy
      yh2=dx
      zh2=0.0
c  rescale the vectors to length scale2
      size=sqrt(xh1*xh1+yh1*yh1)
      rescal=scale2/size
c
      xh1=xh1*rescal
      yh1=yh1*rescal
      xh2=xh2*rescal
      yh2=yh2*rescal
c
c  store these two vectors
      dipvec(1,1)=xh1
      dipvec(2,1)=yh1
      dipvec(3,1)=zh1
c
      dipvec(1,ndip)=xh2
      dipvec(2,ndip)=yh2
      dipvec(3,ndip)=zh2
c
      if (ndip.eq.3) return
c
c  determine two 45 degree dip vectors
      rescal=0.7071068
c
      n1=(1+nv)/2
      n2=(nv+ndip)/2
c
      x451=(xh1+xv)*rescal
      y451=(yh1+yv)*rescal
      z451=(zh1+zv)*rescal
c
      x452=(xh2+xv)*rescal
      y452=(yh2+yv)*rescal
      z452=(zh2+zv)*rescal
c
      dipvec(1,n1)=x451
      dipvec(2,n1)=y451
      dipvec(3,n1)=z451
c
      dipvec(1,n2)=x452
      dipvec(2,n2)=y452
      dipvec(3,n2)=z452
c
      if (ndip.eq.5) return
c
c  determine four 22.5 degree dip vectors
      rescal=0.5411961
c
      dipvec(1,2)=(xh1+x451)*rescal
      dipvec(2,2)=(yh1+y451)*rescal
      dipvec(3,2)=(zh1+z451)*rescal
c
      dipvec(1,4)=(x451+xv)*rescal
      dipvec(2,4)=(y451+yv)*rescal
      dipvec(3,4)=(z451+zv)*rescal
c
      dipvec(1,6)=(xv+x452)*rescal
      dipvec(2,6)=(yv+y452)*rescal
      dipvec(3,6)=(zv+z452)*rescal
c
      dipvec(1,8)=(x452+xh2)*rescal
      dipvec(2,8)=(y452+yh2)*rescal
      dipvec(3,8)=(z452+zh2)*rescal
c
c***** end of subroutine cmpdpv *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine cmpdsv(ndip,iskip,ns,dipvec,disvec)
c
c  parameters
      real dipvec(3,9),disvec(780,9)
c
      integer ndip,iskip,ns
c  local variables
      real darc(257)
c
      integer inc,narc,ndp,np,n,nd
c  coefficients of standard arc
c
       data darc/
     *0.0000000,0.0171173,0.0342346,0.0509603,0.0676860,0.0840283,
     *0.1003707,0.1163377,0.1323047,0.1479038,0.1635029,0.1787412,
     *0.1939796,0.2088641,0.2237485,0.2382856,0.2528226,0.2670183,
     *0.2812141,0.2950745,0.3089350,0.3224657,0.3359963,0.3492027,
     *0.3624091,0.3752962,0.3881833,0.4007561,0.4133289,0.4255921,
     *0.4378553,0.4498133,0.4617713,0.4734285,0.4850857,0.4964461,
     *0.5078065,0.5188740,0.5299417,0.5407202,0.5514988,0.5619919,
     *0.5724850,0.5826960,0.5929070,0.6028394,0.6127717,0.6224285,
     *0.6320853,0.6414696,0.6508538,0.6599685,0.6690831,0.6779310,
     *0.6867787,0.6953623,0.7039459,0.7122679,0.7205898,0.7286526,
     *0.7367154,0.7445213,0.7523272,0.7598785,0.7674298,0.7747287,
     *0.7820274,0.7890757,0.7961241,0.8029239,0.8097238,0.8162769,
     *0.8228301,0.8291384,0.8354468,0.8415120,0.8475771,0.8534007,
     *0.8592244,0.8648080,0.8703916,0.8757367,0.8810817,0.8861896,
     *0.8912975,0.8961695,0.9010416,0.9056790,0.9103164,0.9147204,
     *0.9191245,0.9232962,0.9274679,0.9314083,0.9353487,0.9390589,
     *0.9427691,0.9462499,0.9497307,0.9529830,0.9562353,0.9592599,
     *0.9622845,0.9650821,0.9678797,0.9704511,0.9730224,0.9753681,
     *0.9777138,0.9798344,0.9819550,0.9838510,0.9857470,0.9874189,
     *0.9890908,0.9905390,0.9919872,0.9932120,0.9944367,0.9954384,
     *0.9964401,0.9972190,0.9979979,0.9985541,0.9991102,0.9994439,
     *0.9997776,0.9998888,1.0000000,0.9998888,0.9997776,0.9994439,
     *0.9991102,0.9985541,0.9979979,0.9972190,0.9964401,0.9954384,
     *0.9944367,0.9932120,0.9919872,0.9905390,0.9890908,0.9874189,
     *0.9857470,0.9838510,0.9819550,0.9798344,0.9777138,0.9753681,
     *0.9730224,0.9704511,0.9678797,0.9650821,0.9622845,0.9592599,
     *0.9562353,0.9529830,0.9497307,0.9462499,0.9427691,0.9390589,
     *0.9353487,0.9314083,0.9274679,0.9232962,0.9191245,0.9147204,
     *0.9103164,0.9056790,0.9010416,0.8961695,0.8912975,0.8861896,
     *0.8810817,0.8757367,0.8703916,0.8648080,0.8592244,0.8534007,
     *0.8475771,0.8415120,0.8354468,0.8291384,0.8228301,0.8162769,
     *0.8097238,0.8029239,0.7961241,0.7890757,0.7820274,0.7747287,
     *0.7674298,0.7598785,0.7523272,0.7445213,0.7367154,0.7286526,
     *0.7205898,0.7122679,0.7039459,0.6953623,0.6867787,0.6779310,
     *0.6690831,0.6599685,0.6508538,0.6414696,0.6320853,0.6224285,
     *0.6127717,0.6028394,0.5929070,0.5826960,0.5724850,0.5619919,
     *0.5514988,0.5407202,0.5299417,0.5188740,0.5078065,0.4964461,
     *0.4850857,0.4734285,0.4617713,0.4498133,0.4378553,0.4255921,
     *0.4133289,0.4007561,0.3881833,0.3752962,0.3624091,0.3492027,
     *0.3359963,0.3224657,0.3089350,0.2950745,0.2812141,0.2670183,
     *0.2528226,0.2382856,0.2237485,0.2088641,0.1939796,0.1787412,
     *0.1635029,0.1479038,0.1323047,0.1163377,0.1003707,0.0840283,
     *0.0676860,0.0509603,0.0342346,0.0171173,0.0000000/
      inc=256/ns
        ndip1=1+iskip
      ndip2=ndip-iskip
c
c  loop over dips
      do 30 ndp=ndip1,ndip2
         narc=1
         nd=3
c  loop over points on the path (skip first and last)
         do 20 np=2,ns
            narc=narc+inc
c
c  loop over x,y,z
            do 10 n=1,3
               nd=nd+1
               disvec(nd,ndp)=darc(narc)*dipvec(n,ndp)
c
   10       continue
   20    continue
   30 continue
c
c***** end of subroutine cmpdsv *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine cmdsvm(ndip,iskip,ns,dipvec,disvcm)
c
c  This subroutine uses the 3rd root of the arc from cmpdsv
c  This allows for searching (in rayweb) through flatter paths
c  that may be more appropriate for long ray paths.
c
c  parameters
      real dipvec(3,9),disvcm(780,9)
c
      integer ndip,iskip,ns
c  local variables
      real darc3rt(257)
c
      integer inc,narc,ndp,np,n,nd
c  coefficients of standard arc
c
       data darc3rt/
     *0.0000000,0.2577182,0.3247046,0.3707467,0.4075363,0.4380012,
     *0.4647317,0.4881727,0.5095558,0.5288426,0.5468168,0.5633024,
     *0.5788757,0.5933185,0.6070904,0.6199632,0.6323225,0.6439424,
     *0.6551574,0.6657491,0.6760140,0.6857426,0.6952028,0.7041943,
     *0.7129620,0.7213146,0.7294781,0.7372702,0.7449011,0.7521963,
     *0.7593527,0.7662035,0.7729338,0.7793840,0.7857291,0.7918155,
     *0.7978099,0.8035643,0.8092375,0.8146871,0.8200648,0.8252332,
     *0.8303376,0.8352453,0.8400959,0.8447610,0.8493752,0.8538138,
     *0.8582067,0.8624330,0.8666182,0.8706449,0.8746347,0.8784732,
     *0.8822783,0.8859388,0.8895692,0.8930610,0.8965256,0.8998570,
     *0.9031639,0.9063426,0.9094990,0.9125319,0.9155447,0.9184381,
     *0.9213133,0.9240729,0.9268162,0.9294474,0.9320638,0.9345714,
     *0.9370657,0.9394543,0.9418309,0.9441046,0.9463673,0.9485298,
     *0.9506826,0.9527375,0.9547835,0.9567339,0.9586765,0.9605255,
     *0.9623674,0.9641177,0.9658617,0.9675159,0.9691644,0.9707248,
     *0.9722802,0.9737490,0.9752133,0.9765925,0.9779677,0.9792591,
     *0.9805471,0.9817523,0.9829547,0.9840754,0.9851936,0.9862313,
     *0.9872667,0.9882225,0.9891765,0.9900517,0.9909254,0.9917210,
     *0.9925154,0.9932324,0.9939485,0.9945877,0.9952263,0.9957886,
     *0.9963503,0.9968364,0.9973219,0.9977322,0.9981421,0.9984772,
     *0.9988120,0.9990721,0.9993322,0.9995178,0.9997033,0.9998146,
     *0.9999259,0.9999629,1.0000000,0.9999629,0.9999259,0.9998146,
     *0.9997033,0.9995178,0.9993322,0.9990721,0.9988120,0.9984772,
     *0.9981421,0.9977322,0.9973219,0.9968364,0.9963503,0.9957886,
     *0.9952263,0.9945877,0.9939485,0.9932324,0.9925154,0.9917210,
     *0.9909254,0.9900517,0.9891765,0.9882225,0.9872667,0.9862313,
     *0.9851936,0.9840754,0.9829547,0.9817523,0.9805471,0.9792591,
     *0.9779677,0.9765925,0.9752133,0.9737490,0.9722802,0.9707248,
     *0.9691644,0.9675159,0.9658617,0.9641177,0.9623674,0.9605255,
     *0.9586765,0.9567339,0.9547835,0.9527375,0.9506826,0.9485298,
     *0.9463673,0.9441046,0.9418309,0.9394543,0.9370657,0.9345714,
     *0.9320638,0.9294474,0.9268162,0.9240729,0.9213133,0.9184381,
     *0.9155447,0.9125319,0.9094990,0.9063426,0.9031639,0.8998570,
     *0.8965256,0.8930610,0.8895692,0.8859388,0.8822783,0.8784732,
     *0.8746347,0.8706449,0.8666182,0.8624330,0.8582067,0.8538138,
     *0.8493752,0.8447610,0.8400959,0.8352453,0.8303376,0.8252332,
     *0.8200648,0.8146871,0.8092375,0.8035643,0.7978099,0.7918155,
     *0.7857291,0.7793840,0.7729338,0.7662035,0.7593527,0.7521963,
     *0.7449011,0.7372702,0.7294781,0.7213146,0.7129620,0.7041943,
     *0.6952028,0.6857426,0.6760140,0.6657491,0.6551574,0.6439424,
     *0.6323225,0.6199632,0.6070904,0.5933185,0.5788757,0.5633024,
     *0.5468168,0.5288426,0.5095558,0.4881727,0.4647317,0.4380012,
     *0.4075363,0.3707467,0.3247046,0.2577182,0.0000000/
      inc=256/ns
        ndip1=1+iskip
      ndip2=ndip-iskip
c
c  loop over dips
      do 30 ndp=ndip1,ndip2
         narc=1
         nd=3
c  loop over points on the path (skip first and last)
         do 20 np=2,ns
            narc=narc+inc
c
c  loop over x,y,z
            do 10 n=1,3
               nd=nd+1
               disvcm(nd,ndp)=darc3rt(narc)*dipvec(n,ndp)
c
   10       continue
   20    continue
   30 continue
c
c***** end of subroutine cmdsvm *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine cmpsep(path,pthsep,ns)
c
c  parameters
      real path(780),pthsep(260)
c
      integer ns
c  local variables
      integer nx,ny,nz,nx1,ny1,nz1
c
      nx=-2
c  loop over pairs of points in one set of stored vectors
      do 10 n=1,ns
         nx=nx+3
         ny=nx+1
         nz=nx+2
         nx1=nx+3
         ny1=nx+4
         nz1=nx+5
c
         pthsep(n)=sqrt((path(nx1)-path(nx))**2
     *              +(path(ny1)-path(ny))**2
     *              +(path(nz1)-path(nz))**2)
c
   10 continue
c
c***** end of subroutine cmpsep *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine curvdr(isp,nc,ndp,dtt,tt1,tt2,disvec,pthsep,
     *  strpth,trpth1)
c  This subroutine computes the difference in time between the nc
c  and the nc+1 curve with dip = ndp. dtt is the tt difference,
c  tt1 is the tt of the nc curve, and tt2 is the tt of the nc+1 curve.
      dimension disvec(780,9),pthsep(260),strpth(780),trpth1(780)
      common/raytr/trpath(780,9),npt,ns
      ncv=nc
      call curvtm(isp,ncv,ndp,tt,disvec,pthsep,strpth,trpth1)
      tt1=tt
      nc1=nc+1
      call curvtm(isp,nc1,ndp,tt,disvec,pthsep,strpth,trpth1)
      tt2=tt
      dtt=tt2-tt1
      return
c***** end of subroutine curvdr *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine curvtm(isp,nc,ndp,ttm,disvec,pthsep,strpth,trpth1)
c  This computes the travel time for curve nc at dip ndp. It is
c  used in the faster search programmed by Prothero.
      dimension disvec(780,9),pthsep(260),strpth(780),trpth1(780)
      common/raytr/trpath(780,9),npt,ns
c  loop to determine points along one path
      npt2=npt-2
      do 45 np=1,npt2
         n1=3*np+1
         n3=n1+2
         do 44 nn=n1,n3
            trpath(nn,ndp)=nc*disvec(nn,ndp)+strpth(nn)
            trpth1(nn)=trpath(nn,ndp)
  44     continue
  45  continue
c  set up pthsep array for travel time calculations
      call cmpsep(trpth1,pthsep,ns)
c  compute travel time along the path
      call ttime(isp,ns,npt,trpth1,pthsep,tt)
      ttm=tt
      return
c***** end of subroutine curvtm *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine decide(istop,nit)
c  common block variables:
      include 'simul2014_common.inc'
c
      istop=0
      nobtt=nobt+nordt+noedt
      write(16,1618) nobtt,nobt,nordt,noedt,ssqr
 1618 format('CHECK:nobtt,nobt,nordt,noedt',4i8,', ssqr=',f12.0)
      rms=sqrt(ssqr/float(nobtt))
      dvar=ssqrw/wnobt
      rmsw=sqrt(dvar)
c  CKECK PRINT
      write(16,1610) rms,rmsw,dvar
      write(36,1610) rms,rmsw,dvar
 1610 format(//,' unweighted rms=',f8.5,'; weighted rms=',f8.5,
     2 ' data var.(ssqrw/wnobt ie:rms**2)=',f14.5)
      dvarp=ssqrwp/wnobtp
      if(iuses.eq.1) then
        write(16,1611) dvarp
        write(36,1611) dvarp
 1611   format(50x,'P data var.=',f14.6)
      else
        dvars=ssqrws/wnobts
        write(16,1612) dvarp,dvars
        write(36,1612) dvarp,dvars
 1612   format(40x,'P data var.=',f14.6,'  S-P data var.=',f14.6)
      endif
      if (rmsw.lt.rmstop) then
        istop=1
        write(16,1621) rmsw,rmstop,istop
 1621   format('* * * *',/,'  decide rmsw',g10.3,' < rmstop ',g10.3,
     2  ' istop=',i2,/,'* * * *')
      endif
c
c  f-test
      mbl0=0
c  mbl1 was wrong, did not include stations, 1-april-1983, dmep
      do 10 n=1,npar
         if(khit(n).eq.0) goto 10
         if ((hit(n).lt.hitct).or.(nfix(n).eq.1).or.(imerge(n).eq.1))
     &  go to 10
         mbl0=mbl0+1
   10 continue
      ndof=nobtt-nparhy*neqs-mbl0-nbls
      wndof=wnobt-float(nparhy*neqs+mbl0+nbls)
        write(16,1601) mbl0,mbl1
 1601 format(' subroutine decide, mbl0 = ',i6,' mbl1 = ',i6)
      var=ssqr/float(ndof)
      varw=ssqrw/wndof
      write(16,1606) ssqr,var,ssqrw,varw
 1606 format(/,'  Unweighted:  ssqr =',f12.3,'    var. (ssqr/ndof) =',
     2  f14.6,/,
     3  '  Weighted:   ssqrw =',f12.3,'  varw.(ssqrw/wndof) =',f14.6,/)
c
      if(nit.eq.0) then 
        write(36,3606) nit,ssqrw,varw
 3606   format(' Iteration:',i3,', ssqrw =',f12.3,
     *    '  varw.(ssqrw/wndof) =',f9.6)
        return
      endif
c
      call ftest(ndof,ndof1,ratio)
c
      rat=var1/var
      ratw=varw1/varw
c      write(16,6601) nit, var,ndof,var1,ndof1,rat,ratio
 6601 format(//,' f-test iteration',i3,/,
     *  ' new variance and ndof =',f11.6,i10,
     *  /,' old variance and ndof =',f11.6,i10,
     *  /,' variance ratio and critical ratio =',2f10.3)
      write(16,1605) 
 1605 format (//,'*****WEIGHTED****** use for f-test since weighted ',
     2  'throughout inversion',/)
      iwndof=nint(wndof)
      iwndf1=nint(wndof1)
      call ftest(iwndof,iwndf1,ratio)
      write(16,6602) nit,varw,wndof,varw1,wndof1,ratw,ratio
 6602 format(//,' f-test iteration',i3,/,
     *  ' new variance and ndof =',f11.6,f9.0,
     *  /,' old variance and ndof =',f11.6,f9.0,
     *  /,' variance ratio and critical ratio =',2f10.3)
      write(36,6601) nit,varw,ndof,varw1,ndof1,ratw,ratio
c
      if (ratw.lt.ratio) istop=2
c
c***** end of subroutine decide *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine delcev(ne,nc,nwav)
c  This subroutine removes an event from a cluster.  It works from locclu, after the
c  initial single event location and prior to the getcedt.
c  It moves down all the reamining events in the cluster.
c
      include 'simul2014_common.inc'
c
      integer nwav(2)
c
      write(16,1658) nc,ne,iyrmoce(ne,nc),idayce(ne,nc),ihrce(ne,nc),
     * minoce(ne,nc),secoce(ne,nc)
 1658  format(' ** WARNING:  Cluster:',i5,' DELETING EVENT',i5,
     2  ', with bad wrms: ',a4,a2,1x,a2,i2,1x,f5.2)
      ncevold=ncev(nc)
      nobsold=kobsce(ne,nc)
      do 200 me=ne+1,ncevold
        mobs=kobsce(me,nc)
        me1=me-1
        iyrmoce(me1,nc)=iyrmoce(me,nc)
        idayce(me1,nc)=idayce(me,nc)
        ihrce(me1,nc)=ihrce(me,nc)
        minoce(me1,nc)=minoce(me,nc)
        secoce(me1,nc)=secoce(me,nc)
        ltdce(me1,nc)=ltdce(me,nc)
        celtm(me1,nc)=celtm(me,nc)
        lndce(me1,nc)=lndce(me,nc)
        celnm(me1,nc)=celnm(me,nc)
        rmagce(me1,nc)=rmagce(me,nc)
          cevc(i,me1,nc)=cevc(i,me,nc)
   80   continue
        do 100 m=1,mobs
          dltac(m,me1,nc)=dltac(m,me,nc)
          rdltac(m,me1,nc)=rdltac(m,me,nc)
          iwc(m,me1,nc)=iwc(m,me,nc)
          istoc(m,me1,nc)=istoc(m,me,nc)
          secpc(m,me1,nc)=secpc(m,me,nc)
          rmkc(m,me1,nc)=rmkc(m,me,nc)
          intspc(m,me1,nc)=intspc(m,me,nc)
          wtc(m,me1,nc)=wtc(m,me,nc)
  100   continue
        kobsce(me1,nc)=kobsce(me,nc)
  200 continue
      ncev(nc)=ncev(nc)-1
      nobsc(nc)=nobsc(nc)-nobsold
      nobt=nobt-nobsold
      nobtp=nobtp-nwav(1)
      nobts=nobts-nwav(2)
      nobtce=nobtce-nobsold
      nobteq=nobteq-nobsold
c
c compute centroid of cluster
  500 continue
      ne=ncev(nc)
      do 515 i=1,3
      csum=0.0
        do 510 j=1,ne
        csum=csum+cevc(i,j,nc)
  510   continue
      ccenc(i,nc)=csum/float(ne)
  515 continue
c
c***** end of subroutine delcev *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine disto(latd,rlat,lond,rlon,x,y)
c  this routine calculates distance of station or event
c  from given coordinate origin in terms of (possibly
c  rotated) cartesian coords x and y
c  uses short distance conversion factors from setorg
c  or New Zealand Map Grid 
c  or state plane coordinates for Alaska
c
c  declaration statements:
      double precision drad,drlt
      parameter (drad=1.7453292d-02)
      parameter (drlt=9.9330647d-01)
c
c  common block variables:
      common/shortd/ xltkm,xlnkm,rota,nzco,xlt,xln,snr,csr,
     2 onorth,oeast,rotadg,cmerid,rlatdo,rlondo
      include 'simul2014_common.inc'
c
cek correction for S and E: minutes must also be taken as negative!
      plt=abs(60.*latd+rlat)
      pln=abs(60.*lond+rlon)
      if(latd.lt.0) plt=-1.*plt
      if(lond.lt.0) pln=-1.*pln
c
c      print *,'DISTO ltds,sltm:',latd,rlat,', lnds,slnm:',
c     2 lond,rlon,', plt,pln:',plt,pln
c      print *,'  xlt,xln:',xlt,xln
      if(nzco.lt.5) then
        latdp=latd
        platm=rlat
        londp=lond
        plonm=rlon
        if(nzco.eq.4) call llnzmg(latdp,platm,londp,plonm,ynorth,xeast)
        if(nzco.eq.2) call ll2spc(latdp,platm,londp,plonm,ynorth,xeast)
        if((nzco.le.1).or.(nzco.eq.3)) call tmllxy(latdp,platm,londp,
     2    plonm,ynorth,xeast,1)
c      WRITE(6,*) 'disto: latdp,platm ',latdp,platm,' londp,plonm ',
c     2 londp,plonm,' xeast,ynorth ',xeast,ynorth
      endif
c
      if(nzco.eq.5) then
c  now convert lat and lon differences to km
c  using short distance conversion
        x=pln-xln
        y=plt-xlt
c
        xlt1=atan(drlt*tan(drad*(plt+xlt)/120.))
        x=x*xlnkm*cos(xlt1)
        y=y*xltkm
        x1=x
        y1=y
c  now do rotation
        if(rota.eq.0.0) return
        ty=csr*y+snr*x
        x=csr*x-snr*y
        y=ty
c        write(16,1602) x1,y1,x,y
c 1602 format(4x,' xsd,ysd=',2f8.2,'   rot x,y=',2f8.2)
      else
c now convert using Transverse Mercator, or NZMG,
c or State Plane Coordinates for Alaska
c since SIMUL has W=positive, must subtract point from origin for easting.
c      PRINT *,'oeast,xeast ',oeast,xeast,'onorth,ynorth',onorth,ynorth
        xnz=oeast-xeast
        ynz=ynorth-onorth
        if(rota.ne.0.0) then
          y=csr*ynz+snr*xnz
          x=csr*xnz-snr*ynz
c          write(16,1601) xnz,ynz,x,y
c 1601     format(4x,' xnz,ynz=',2f8.2,'   rot x,y=',2f8.2)
        else
          x=xnz
          y=ynz
        endif
      endif
c
      return
c***** end of subroutine disto *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine fksvd (a, s, v, mmax, nmax, m, n, p, withu, withv)
c
c  common block variables:
      include 'simul2014_common.inc'
c
      common/machin/ eta,tol
c
      integer    mmax, nmax, m, n, p
      real       r, w, cs, sn, tol, f, x, eps, gf, t, y
      real       eta, h, q, z
      integer    i, j, k, l, l1, n1, np
      logical    withu, withv
      double precision a(mmax,mmax)
      dimension        s(nmax), v(nmax,nmax)
      dimension as1(mmax),as2(mmax)
      dimension  t(mmax)
c
c     ------------------------------------------------------------------
c
c     this is a translation of a cdc 6600 fortran program to ibm 360
c     fortran iv.  this subroutine uses short precision arithmetic.
c     a long precision version is available under the name 'dsvd'.
c
c     this subroutine replaces earlier subroutines with the same name,
c    689   6       &   &0&  s of a complex ar&thmetic program, published
c     as algorithm 358.  this current program is faster, more accurate
c     and less obscure in describing its capabilities.
c
c     original programmer=  r. c. singleton
c     360 version by=       j. g. lewis
c     last revision of this subroutine=  4 december 1973
c
c     ------------------------------------------------------------------
c
c     additional subroutine needed=  rotate
c
c     ------------------------------------------------------------------
c
c
c     this subroutine computes the singular value decomposition
c     of a real m*n matrix a, i.e. it computes matrices u, s, and v
c     such that
c
c                  a = u * s * vt ,
c     where
c              u is an m*n matrix and ut*u = i, (ut=transpose
c                                                    of u),
c              v is an n*n matrix and vt*v = i, (vt=transpose
c                                                    of v),
c        and   s is an n*n diagonal matrix.
c
c     description of parameters=
c
c     a = real array. a contains the matrix to be decomposed.
c         the original data are lost.  if withv=.true., then
c         the matrix u is computed and stored in the array a.
c
c     mmax = integer variable.  the number of rows in the
c            array a.
c
c     nmax = integer variable.  the number of rows in the
c            array v.
c
c     m,n = integer variables.  the number of rows and columns
c           in the matrix stored in a.  (ng=mg=100.  if it is
c           necessary to solve a larger problem, then the
c           amount of storage allocated to the array t must
c           be increased accordingly.)  if mlt n , then either
c           transpose the matrix a or add rows of zeros to
c           increase m to n.
c
c     p = integer variable.  if p'0, then columns n+1, . . . ,
c         n+p of a are assumed to contain the columns of an m*p
c         matrix b.  this matrix is multiplied by ut, and upon
c         exit, a contains in these same columns the n*p matrix
c         ut*b. (p'=0)
c
c     withu, withv = logical variables.  if withu=.true., then
c         the matrix u is computed and stored in the array a.
c         if withv=.true., then the matrix v is computed and
c         stored in the array v.
c
c     s = real array.  s(1), . . . , s(n) contain the diagonal
c         elements of the matrix s ordered so than s(i)>=s(i+1),
c         i=1, . . . , n-1.
c
c     v = real array.  v contains the matrix v.  if withu
c         and withv are not both =.true., then the actual
c         parameter corresponding to a and v may be the same.
c
c     this subroutine is a real version of a fortran subroutine
c     by businger and golub, algorithm 358=  singular value
c     decomposition of a complex matrix, comm. acm, v. 12,
c     no. 10, pp. 564-565 (oct. 1969).
c     with revisions by rc singleton, may 1972.
c     ------------------------------------------------------------------
c
c
c     VAX version:  machine constants calculated internally in strt
c     eta is the machine epsilon (relative accuracy)
c     tol is the smallest representable real divided by eta.
c
c      PRINT *,'fksvd mmax=',mmax,' nmax=',nmax,' m=',m,' n=',n,' p=',p
      np = n + p
      n1 = n + 1
c
c     householder reduction to bidiagonal form
      gf = 0.0
      eps = 0.0
      l = 1
   10 t(l) = gf
      k = l
      l = l + 1
c
c     elimination of a(i,k), i=k+1, . . . , m
      s(k) = 0.0
      z = 0.0
      do 20 i = k,m
         z = z + a(i,k)**2
   20 continue
      if (z.lt.tol) goto 50
      gf = sqrt(z)
      f = a(k,k)
      if (f.ge.0.0) gf = - gf
      s(k) = gf
      h = gf * (f - gf)
      a(k,k) = f - gf
      if (k.eq.np) goto 50
      do 45 j = l,np
         f = 0
         do 30 i = k,m
            f = f + a(i,k)*a(i,j)
   30    continue
         f = f/h
         do 40 i = k,m
            a(i,j) = a(i,j) + f*a(i,k)
   40    continue
   45 continue
c
c     elimination of a(k,j), j=k+2, . . . , n
   50 eps = amax1(eps,abs(s(k)) + abs(t(k)))
      if (k.eq.n) goto 100
      gf = 0.0
      z = 0.0
      do 60 j = l,n
         z = z + a(k,j)**2
   60 continue
      if (z.lt.tol) goto 10
      gf = sqrt(z)
      f = a(k,l)
      if (f.ge.0.0) gf = - gf
      h = gf * (f - gf)
      a(k,l) = f - gf
      do 70 j = l,n
         t(j) = a(k,j)/h
   70 continue
      do 95 i = l,m
         f = 0
         do 80 j = l,n
            f = f + a(k,j)*a(i,j)
   80    continue
      if (abs(f).lt.1.0e-25) f=0.
         do 90 j = l,n
            a(i,j) = a(i,j) + f*t(j)
   90    continue
   95 continue
c
      goto 10
c
c     tolerance for negligible elements
  100 eps = eps*eta
c
c     accumulation of transformations
      if (.not.withv) goto 160
      k = n
      goto 140
  110    if (t(l).eq.0.0) goto 140
         h = a(k,l)*t(l)
         do 135 j = l,n
            q = 0
            do 120 i = l,n
               q = q + a(k,i)*v(i,j)
  120       continue
            q = q/h
            do 130 i = l,n
               v(i,j) = v(i,j) + q*a(k,i)
  130       continue
  135    continue
  140    do 151 j = 1,n
            v(k,j) = 0
  151    continue
         v(k,k) = 1.0
         l = k
         k = k - 1
         if (k.ne.0) goto 110
c
  160 k = n
      if (.not.withu) goto 230
      gf = s(n)
      if (gf.ne.0.0) gf = 1.0/gf
      go to 210
  170    do 180 j = l,n
            a(k,j) = 0
  180    continue
         gf = s(k)
         if (gf.eq.0.0) goto 210
         h = a(k,k)*gf
         do 205 j = l,n
            q = 0
            do 190 i = l,m
               q = q + a(i,k)*a(i,j)
  190       continue
            q = q/h
            do 200 i = k,m
               a(i,j) = a(i,j) + q*a(i,k)
  200       continue
  205    continue
         gf = 1.0/gf
  210    do 220 j = k,m
            a(j,k) = a(j,k)*gf
  220    continue
         a(k,k) = a(k,k) + 1.0
         l = k
         k = k - 1
         if (k.ne.0) goto 170
c
c     qr diagonalization
      k = n
c
c     test for split
  230    l = k
  240       if (abs(t(l)).le.eps) goto 290
            l = l - 1
            if (abs(s(l)).gt.eps) goto 240
c
c     cancellation
         cs = 0.0
         sn = 1.0
         l1 = l
         l = l + 1
         do 280 i = l,k
            f = sn*t(i)
            t(i) = cs*t(i)
            if (abs(f).le.eps) goto 290
            h = s(i)
            w = sqrt(f*f + h*h)
            s(i) = w
            cs = h/w
            sn = - f/w
            do 221 ia=1,m
               as1(ia) = sngl(a(ia,l1))
               as2(ia) = sngl(a(ia,i))
  221       continue
            if (withu) call hyrot(mmax,as1, as2, cs, sn, m)
            do 222 ia=1,m
               a(ia,l1)=dble(as1(ia))
               a(ia,i)=dble(as2(ia))
  222       continue
            if (np.eq.n) goto 280
            do 270 j = n1,np
               q = a(l1,j)
               r = a(i,j)
               a(l1,j) = q*cs + r*sn
               a(i,j) = r*cs - q*sn
  270       continue
  280    continue
c
c     test for convergence
  290    w = s(k)
         if (l.eq.k) goto 360
c
c     origin shift
         x = s(l)
         y = s(k-1)
         gf = t(k-1)
         h = t(k)
         f = ((y - w)*(y + w) + (gf - h)*(gf + h))/(2.0*h*y)
         gf = sqrt(f*f + 1.0)
         if (f.lt.0.0) gf = - gf
         f = ((x - w)*(x + w) + (y/(f + gf) - h)*h)/x
c
c     qr step
         cs = 1.0
         sn = 1.0
         l1 = l + 1
         do 350 i = l1,k
            gf = t(i)
            y = s(i)
            h = sn*gf
            gf = cs*gf
            w = sqrt(h*h + f*f)
            t(i-1) = w
            cs = f/w
            sn = h/w
            f = x*cs + gf*sn
            gf = gf*cs - x*sn
            h = y*sn
            y = y*cs
            if (withv) call hyrot(nmax,v(1,i-1),v(1,i),cs,sn,n)
            w = sqrt(h*h + f*f)
            s(i-1) = w
            cs = f/w
            sn = h/w
            f = cs*gf + sn*y
            x = cs*y - sn*gf
            do 338 ia=1,m
               as1(ia) = sngl(a(ia,i-1))
               as2(ia) = sngl(a(ia,i))
  338       continue
            if (withu) call hyrot(mmax,as1, as2, cs, sn, m)
            do 339 ia=1,m
               a(ia,i-1)=dble(as1(ia))
               a(ia,i)=dble(as2(ia))
  339       continue
            if (n.eq.np) goto 350
            do 340 j = n1,np
               q = a(i-1,j)
               r = a(i,j)
               a(i-1,j) = q*cs + r*sn
               a(i,j) = r*cs - q*sn
  340       continue
  350    continue
c
         t(l) = 0.0
         t(k) = f
         s(k) = x
         goto 230
c
c     convergence
  360    if (w.ge.0.0) goto 380
         s(k) = - w
         if (.not.withv) goto 380
         do 370 j = 1,n
            v(j,k) = - v(j,k)
  370    continue
  380    k = k - 1
         if (k.ne.0) go to 230
c
c     sort singular values
      do 450 k = 1,n
         gf = -1.0
         do 390 i = k,n
            if (s(i).lt.gf) goto 390
            gf = s(i)
            j = i
  390    continue
         if (j .eq. k) go to 450
         s(j) = s(k)
         s(k) = gf
         if (.not.withv) goto 410
         do 400 i = 1,n
            q = v(i,j)
            v(i,j) = v(i,k)
            v(i,k) = q
  400    continue
  410    if (.not.withu) goto 430
         do 420 i = 1,m
            q = a(i,j)
            a(i,j) = a(i,k)
            a(i,k) = q
  420    continue
  430    if (n.eq.np) goto 450
         do 440 i = n1,np
            q = a(j,i)
            a(j,i) = a(k,i)
            a(k,i) = q
  440    continue
  450 continue
c
      return
c***** end of subroutine fksvd *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine forwrd(ne,nit)
c  this routine calculates the forward problem
c  given the stations, and initial hypocenters and
c  velocity model, the routine calculates theoretical
c  travel times in a way appropriate to the assumed
c  form of the velocity model.  three-d ray tracing
c  is performed. travel time derivatives with respect to
c  hypocentral and medium parameters are calculated
c
c  declaration statements:
      character*8 fil15
      character*1 phs(2)
c
c  common block variables:
      common/shortd/ xltkm,xlnkm,rota,nzco,xlt,xln,snr,csr,
     2 onorth,oeast,rotadg,cmerid,rlatdo,rlondo
      include 'simul2014_common.inc'
c
      data phs/'P','S'/
      zoff=99.0
      nt=ne-nebs
c  event coordinates
      if(nt.lt.1) then
        xe=evc(1,ne)
        ye=evc(2,ne)
        ze=evc(3,ne)
      endif
c  loop over all observations of this event
      nobs=kobs(ne)
      write(26,2601) ne,iyrmo(ne),iday(ne),ihr(ne),mino(ne),
     2 seco(ne),ltde(ne),eltm(ne),lnde(ne),elnm(ne)
 2601 format(3h **,1x,i3,1x,a4,a2,1x,a2,i2,f6.2,i4,i7.2,i5,i7.2,
     2 2f7.2,3f6.2,3x,3i3)
      if(kout3.eq.0) goto 20
      write(fil15(1:8),1500) ne
 1500 format('ev',i3.3,'.rp')
      open(unit=15,file=fil15)
      write(15,1800) zoff,zoff,zoff
 1800 format('RAYPATH POINTS:',t68,f10.5,/,
     2 ' lat lon z x y -z dr v',t68,f10.5,/,
     3 'format:17x,7f10.5',t68,f10.5)
   20 do 77 no=1,nobs
c for teleseismic, use piercing point for each begin obs raypath
        if(nt.ge.1) then
          xe=tepp(1,no,nt)
          ye=tepp(2,no,nt)
          ze=tepp(3,no,nt)
        endif
         call path(ne,no,xe,ye,ze,ttime)
c      PRINT *,'forwrd ne,nt,xe,ye,ze,ttime',ne,nt,xe,ye,ze,ttime
c** Change for synthetic data, nitmax= -1
         if(nitmax.eq.-1) then
            secp(no,ne)=ttime+seco(ne)
            if(nt.ge.1) secte(no,nt)=ttime
c * * Different for S-P synthetic data * *
            isp=intsp(no,ne)
            if(isp.eq.1) then
               intsp(no,ne)=0
               call path(ne,no,xe,ye,ze,ptime)
               intsp(no,ne)=1
               smptime=ttime-ptime
               secp(no,ne)=smptime
               if(nt.ge.1) secte(no,nt)=smptime
            endif
         else
            call ttmder(ne,no,1,ttime,nit,nnodej)
         endif
         if(kout3.eq.0) goto 77
c  write out raypath points
         write(15,1801) ne,no,zoff
 1801    format(' ev=',i4,' obs=',i4,t68,f10.5)
         write(15,1803) iyrmo(ne),iday(ne),ihr(ne),mino(ne),
     2   seco(ne),zoff
 1803    format(a4,a2,1x,a2,i2,1x,f5.2,t68,f10.5)
         if(kttfor.ne.3) then
           write(15,1804) stn(isto(no,ne)),phs(intsp(no,ne)+1),
     2     dlta(no,ne),zoff
         else
           write(15,1805) stn6(isto(no,ne)),phs(intsp(no,ne)+1),
     2     dlta(no,ne),zoff
         endif
 1804    format(a4,a1,f6.2,t68,f10.5)
 1805    format(a6,a1,f6.2,t68,f10.5)
         do 50 i=1,nrp(no)
            call latlon(rp(1,i,no),rp(2,i,no),lat,xlat,lon,xlon)
            xlat=float(lat)+xlat/60.0
            xlon=float(lon)+xlon/60.0
            if(nzco.eq.1) xlon=abs(xlon)
c  also print out negative z for easy plotting
            zd= -1.0 * rp(3,i,no)
            if(i.gt.1) then
              dx=rp(1,i,no)-rp(1,1,no)
              dy=rp(2,i,no)-rp(2,1,no)
              dr=sqrt((dx*dx)+(dy*dy))
            else
              dr=0.0
            endif
c print out velocity also
            call vel3eft(isp,rp(1,i,no),rp(2,i,no),rp(3,i,no),v)
            write(15,1802) xlat,xlon,rp(3,i,no),rp(1,i,no),
     2         rp(2,i,no),zd,dr,v
   50    continue
 1802    format(17x,2f11.5,f9.3,2f10.3,f9.3,f10.3,f8.3)
   77 continue
      if(kout3.eq.1) close(15)
      return
c***** end of subroutine forwrd *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine forwrdce(nc,ne,nit)
c  This version is for cluster events
c  this routine calculates the forward problem
c  given the stations, and initial hypocenters and
c  velocity model, the routine calculates theoretical
c  travel times in a way appropriate to the assumed
c  form of the velocity model.  three-d ray tracing
c  is performed. travel time derivatives with respect to
c  hypocentral and medium parameters are calculated
c
c  declaration statements:
      character*8 fil15
      character*1 phs(2)
c
c  common block variables:
      common/shortd/ xltkm,xlnkm,rota,nzco,xlt,xln,snr,csr,
     2 onorth,oeast,rotadg,cmerid,rlatdo,rlondo
      include 'simul2014_common.inc'
c
      data phs/'P','S'/
      zoff=99.0
c  event coordinates
      xe=cevc(1,ne,nc)
      ye=cevc(2,ne,nc)
      ze=cevc(3,ne,nc)
c  loop over all observations of this event
      nobs=kobsce(ne,nc)
      write(26,2601) ne,iyrmoce(ne,nc),idayce(ne,nc),ihrce(ne,nc),
     2 minoce(ne,nc),secoce(ne,nc),
     3 ltdce(ne,nc),celtm(ne,nc),lndce(ne,nc),celnm(ne,nc)
 2601 format(3h **,1x,i3,1x,a4,a2,1x,a2,i2,f6.2,i4,i7.2,i5,i7.2,
     2 2f7.2,3f6.2,3x,3i3)
      if(kout3.eq.0) goto 20
      write(fil15(1:8),1500) ne
 1500 format('ce',i3.3,'.rp')
      open(unit=15,file=fil15)
      write(15,1800) zoff,zoff,zoff
 1800 format('RAYPATH POINTS:',t68,f10.5,/,
     2 ' lat lon z x y -z dr v',t68,f10.5,/,
     3 'format:17x,7f10.5',t68,f10.5)
   20 do 77 no=1,nobs
         call pathce(nc,ne,no,xe,ye,ze,ttime)
c** Change for synthetic data, nitmax= -1
         if(nitmax.eq.-1) then
            secpc(no,ne,nc)=ttime+secoce(ne,nc)
c * * Different for S-P synthetic data * *
            isp=intspc(no,ne,nc)
            if(isp.eq.1) then
               intspc(no,ne,nc)=0
               call pathce(nc,ne,no,xe,ye,ze,ptime)
               intspc(no,ne,nc)=1
               smptime=ttime-ptime
               secpc(no,ne,nc)=smptime
            endif
         else
            call ttmderce(nc,ne,no,1,ttime,nit,nnodej)
         endif
         if(kout3.eq.0) goto 77
c  write out raypath points
         write(15,1801) ne,no,zoff
 1801    format(' ev=',i4,' obs=',i4,t68,f10.5)
         write(15,1803) iyrmo(ne),iday(ne),ihr(ne),mino(ne),
     2   secoce(ne,nc),zoff
 1803    format(a4,a2,1x,a2,i2,1x,f5.2,t68,f10.5)
         if(kttfor.ne.3) then
           write(15,1804) stn(istoc(no,ne,nc)),phs(intspc(no,ne,nc)+1),
     2     dltac(no,ne,nc),zoff
         else
           write(15,1805)stn6(istoc(no,ne,nc)),phs(intspc(no,ne,nc)+1),
     2     dltac(no,ne,nc),zoff
         endif
 1804    format(a4,a1,f6.2,t68,f10.5)
 1805    format(a6,a1,f6.2,t68,f10.5)
         do 50 i=1,nrpce(no,ne)
            call latlon(rpce(1,i,no,ne),rpce(2,i,no,ne),lat,xlat,
     2       lon,xlon)
            xlat=float(lat)+xlat/60.0
            xlon=float(lon)+xlon/60.0
            if(nzco.eq.1) xlon=abs(xlon)
c  also print out negative z for easy plotting
            zd= -1.0 * rpce(3,i,no,ne)
            if(i.gt.1) then
              dx=rpce(1,i,no,ne)-rpce(1,1,no,ne)
              dy=rpce(2,i,no,ne)-rpce(2,1,no,ne)
              dr=sqrt((dx*dx)+(dy*dy))
            else
              dr=0.0
            endif
c print out velocity also
            call vel3eft(isp,rpce(1,i,no,ne),rpce(2,i,no,ne),
     2         rpce(3,i,no,ne),v)
            write(15,1802) xlat,xlon,rpce(3,i,no,ne),rpce(1,i,no,ne),
     2         rpce(2,i,no,ne),zd,dr,v
   50    continue
 1802    format(17x,2f11.5,f9.3,2f10.3,f9.3,f10.3,f8.3)
   77 continue
      if(kout3.eq.1) close(15)
      return
c***** end of subroutine forwrdce *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine forwrdedt(nc,nit,kobsok)
c  this routine calculates the forward problem
c  Given the pair of earthquakes and velocity model 
c  the routine calculates theoretical
c  travel times in a way appropriate to the assumed
c  form of the velocity model.
c  Travel time derivatives with respect to
c  medium parameters within distance frayedt are calculated.
c  For receiver-pair differential travel-time residual,
c  the paths are already calculated.  So forwrdrdt simply
c  loops through ttmder for the rdt observations.
c
      include 'simul2014_common.inc'
      integer nnode(maxedt)
c
c  for edt obs, they are grouped by cluster rather than eq,
c   so ne is a dummy here
      ne=1
      nobs=kobsedt(nc)
      kobsok=nobs
      do 77 no=1,nobs
      call ttmderce(nc,ne,no,2,ttime,nit,nnodej)
      nnode(no)=nnodej
   77 continue
c
      if(nit.gt.0) goto 90
c Remove observation if the edt paths did not sample inversion
c nodes (eg - only fixed nodes and then nnodej=0)
c Remove observation if the resedt is too large, as we only
c want high quality for edt. resedt must be < res4 (res4 is input)
c For edt, remove by downweighting
      do 87 no=nobs,1,-1
      if((nnode(no).gt.0).and.(resedt(no).lt.res4)) goto 87
c      if(nnode(no).eq.0) goto 79
      if(nnode(no).eq.0) then
c        write(16,*) 'forwrdedt nnodej=0'
        goto 79
      endif
      if(resedt(no).ge.res4) then
c        write(16,*) resedt(no), res4,'resedt(no), res4'
      endif
   79 isp=intsped(no,nc)
      if(isp.eq.0) noedtp=noedtp-1
      if(isp.eq.1) noedts=noedts-1
      noedt=noedt-1
      kobsok=kobsok-1
      iwedt(no,nc)=8
      wtedtobs(no,nc)=0.0
   87 continue
   90 continue
      return
c***** end of subroutine forwrdedt *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine forwrdrdt(ne,nit)
c  this routine calculates the forward problem
c  Given the pair of stations and velocity model 
c  the routine calculates theoretical
c  travel times in a way appropriate to the assumed
c  form of the velocity model.
c  Travel time derivatives with respect to
c  medium parameters within distance frayrdt are calculated.
c  For receiver-pair differential travel-time residual,
c  the paths are already calculated.  So forwrdrdt simply
c  loops through ttmder for the rdt observations.
c
      include 'simul2014_common.inc'
      integer nnode(maxobs)
c
      nobs=kobsrdt(ne)
      do 77 no=1,nobs
      call ttmder(ne,no,2,ttime,nit,nnodej)
      nnode(no)=nnodej
   77 continue
c
      if(nit.gt.0) goto 90
c Remove observation if the rdt paths did not sample inversion
c nodes (eg - only fixed nodes and then nnodej=0)
c Remove observation if the resrdt is too large, as we only
c want high quality for rdt. resrdt must be < res4 (res4 is input)
      kobsr=nobs
      do 87 no=nobs,1,-1
      if((nnode(no).gt.0).and.(resrdt(no).lt.res4)) goto 87
c      if(nnode(no).eq.0) goto 79
      if(nnode(no).eq.0) then
c        write(16,*) 'nnodej=0'
        goto 79
      endif
      if(resrdt(no).ge.res4) then
c        write(16,*) resrdt(no), res4,'resrdt(no), res4'
      endif
   79 isp=intsprd(no,ne)
      if(isp.eq.0) nordtp=nordtp-1
      if(isp.eq.1) nordts=nordts-1
      nordt=nordt-1
      kobsr1=kobsr
      kobsr=kobsr-1
      if(no.eq.kobsr1) goto 87
      do 85 j=no,kobsr
      intsprd(j,ne)=intsprd(j+1,ne)
      iwrdt(j,ne)=iwrdt(j+1,ne)
      wtrdtobs(j,ne)=wtrdtobs(j+1,ne)
      rdtsec(j,ne)=rdtsec(j+1,ne)
      resrdt(j)=resrdt(j+1)
      ttmrdt(j)=ttmrdt(j+1)
      do 80 i=1,2
        jobsrd(i,j,ne)=jobsrd(i,j+1,ne)
        istord(i,j,ne)=istord(i,j+1,ne)
        nrprdt(i,j,ne)=nrprdt(i,j+1,ne)
   80 continue
   85 continue
   87 continue
      kobsrdt(ne)=kobsr
      if(kobsr.eq.0) return
      nobs=kobsrdt(ne)
   90 continue
c      do 95 noe=1,nobs
c        if(noe.eq.1) write(16,1609) ne
c 1609    format('forwrdtdt resid for rec-pair diff-time: event=',i5,/,
c     *   ' Sta1  ttm1  res1    Sta2  ttm2  res2 : ttmrdt ',
c     *   ' rdtsec   resrdt')
c        j1=jobsrd(1,noe,ne)
c        j2=jobsrd(2,noe,ne)
c        write(16,1611)stn(isto(j1,ne)),ttm(j1),res(j1),stn(isto(j2,ne)),
c     *   ttm(j2),res(j2),ttmrdt(noe),rdtsec(noe,ne),resrdt(noe)
c 1611   format(2(1x,a4,2f6.2),2x,2f6.2,f7.2)
c   95 continue
      return
c***** end of subroutine forwrdrdt *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine ftest(ndf,ndi,ratio)
c
c  interpolate in ftest ratio table
c  ndi is the initial number of degrees of freedom
c  ndf is the final number of degrees of freedom
c
      dimension rattab(4,4),valu(4)
      data rattab/1.69,1.64,1.58,1.51,1.59,1.53,1.47,1.39,
     *            1.50,1.43,1.35,1.25,1.39,1.32,1.22,1.00/
      data valu/0.025,0.016667,0.008333,0.00/
c
      if (ndi.ge.40.and.ndf.ge.40) go to 10
c
c  number of degrees of freedom outside table range
      ratio=2.0
      return
c
c  determine points in table for interploating
   10 continue
      index=ndi/30.
      index1=ndf/30.
      if (index.gt.4) index=4
      if (index1.gt.4) index1=4
c
      go to (11,12,12,14), index
c
   11 continue
      m=1
      m1=2
      go to 20
c
   12 continue
      m=2
      m1=3
      go to 20
c
   14 continue
      m=3
      m1=4
c
   20 continue
c
      go to (21,22,22,24), index1
c
   21 continue
      n=1
      n1=2
      go to 30
c
   22 continue
      n=2
      n1=3
      go to 30
c
   24 continue
      n=3
      n1=4
c
   30 continue
c
c  compute interpolated value from table
      fm=1.0/float(ndi)
      fn=1.0/float(ndf)
c
      vm=valu(m)
      vm1=valu(m1)
c
c  interpolate along m
      vf=rattab(m,n)+(rattab(m1,n)-rattab(m,n))*
     *   (fm-vm)/(vm1-vm)
      vf1=rattab(m,n1)+(rattab(m1,n1)-rattab(m,n1))*
     *    (fm-vm)/(vm1-vm)
c
c
      vn=valu(n)
      vn1=valu(n1)
c  interpolate across n
      ratio=vf+(vf1-vf)*(fn-vn)/(vn1-vn)
c***** end of subroutine ftest *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine getcedt(nc,nobedt)
c
c  for a given cluster,nc, setup earthquake-pair dt for common stations
c       edisedt = maximum distance between edt earthquakes, an input parameter
c       edise12 = distance between 2 specific edt earthquakes
c                 This is used for frayedt and gminedt.
c       frayedt = distance along ray to perturb = facraye * edisedt (edise12)
c       gminedt = minimum slant path length to all edt = 4 * edisedt (edise12)
c       gminedt0 = absolute minimum slant path length for an observation in edt
c       mxiqedt  = max. obs quality for dt observ.
c       nstaepair = maximum number of edt for specific eq pair
c       wtepdt = wt factor for eq-pair dt, like wtsp or wtrdtsht
c
      include 'simul2014_common.inc'
c
c      frayfac= 1.5
      frayfac= facraye
      gminfac= 4.0
c August 2015 changed this to an input parameter
c      gminedt0 = 15.0
c
c max edt could be less for a specific cluster
      mxedtc=mxobsc-nobsc(nc)
      if(mxedtc.gt.maxedt) mxedtc=maxedt
c  Total up edt obs
      if(nc.gt.1) goto 1
      noedt=0
      noedtp=0
      noedts=0
c this is revised to first search for S-P since fewer possible
c allow half S-P
      nstaepair2=nstaepair/2
c
c  Loop through all earthquake pairs, comparing interevent distance and slant 
c  raypath (rdltac(nobs,ne,nc) length to edise12 and gminedt
    1 nobedt=0
c check print to fort.76
c       write(76,7600) nc
 7600 format(' Getedt for cluster:',i5,/,' noed  edtsec ',
     &  'wtedtob isp  dlta1 stn1 ne1  secpc1 secoce1 stn2',
     &  ' ne2 secpc2 secoce2 isp2 ')
      nec=ncev(nc)
c  for event 1, compare stations for event 2
   10 do 400 ne1=1,(nec-1)
      ne2a=ne1+1
      nobs1=kobsce(ne1,nc)
      xe1=cevc(1,ne1,nc)
      ye1=cevc(2,ne1,nc)
      do 350 ne2=ne2a,nec
        xe2=cevc(1,ne2,nc)
        ye2=cevc(2,ne2,nc)
        dx=xe1-xe2
        dy=ye1-ye2
        edise12=sqrt(dx*dx + dy*dy)
        if(edise12.gt.edisedt) goto 350
        frayedt=frayfac*edise12
        gminedt=gminfac*edise12
        nobs2=kobsce(ne2,nc)
c check print
c      print *, 'ne1,ne2,edise12 okay:',ne1,ne2,edise12
c      write(77,7702) 
 7702 format('for ne2: stn intspc iwc secpc rdltac')
c      write(77,7701) ((stn(istoc(j,ne2,nc)),intspc(j,ne2,nc),
c     2  iwc(j,ne2,nc),secpc(j,ne2,nc),rdltac(j,ne2,nc)),j=1,nobs2)
 7701 format(4(1x,a4,2i2,f7.2,f8.2))
        noepair=0
c first loop through for S-P, then for P
        ks=0
        ne1s=ne1
c Loop through obs for event 1
   20   do 200 i=1,nobs1
        isp1=intspc(i,ne1,nc)
        if((ks.eq.0).and.(isp1.eq.0)) goto 200
        if((ks.eq.1).and.(isp1.eq.1)) goto 200
        iq1=iwc(i,ne1,nc)
        if(iq1.gt.mxiqedt) goto 200
        dlta1=rdltac(i,ne1,nc)
        if(dlta1.lt.gminedt) goto 200
        if(dlta1.lt.gminedt0) goto 200
        ista1=istoc(i,ne1,nc)
c        write(76,602)ne1,stn(ista1),iwc(i,ne1,nc),rdltac(i,ne1,nc),
c     &  edise12
  602   format('ne1=',i3,', sta ',a4,',iw=',i4,',dist=',f7.2,
     &  ' edise12=',f7.2)
c
c get number of raypath points for edt partials, 1st eq path
        nrp1=nrp(i)
        nrpedt1=1
        do 50 k=nrp1,2,-1
          k1=k-1
          rx=rp(1,k,i)
          ry=rp(2,k,i)
          rz=rp(3,k,i)
          dx=rp(1,k1,i)-rx
          dy=rp(2,k1,i)-ry
          dz=rp(3,k1,i)-rz
c  compute segment length
          sl=sqrt(dx*dx+dy*dy+dz*dz)
          raylen=raylen+sl
          nrpedt1=nrpedt1+1
          if(raylen.ge.frayedt) goto 55
   50   continue
   55   continue
c
c      write(76,605) ne2,nobs2
  605 format('ne2:',i4,'nobs2=',i4,' istoc(j,ne2,nc),stn')
c      write(76,603) (istoc(j,ne2,nc),stn(istoc(j,ne2,nc)),j=1,nobs2)
  603 format(5(i5,1x,a4))
c
c  loop through obs for eq 2 to match station
        do 150 j=1,nobs2
          ista2=istoc(j,ne2,nc)
          if(ista2.ne.ista1) then
c           if(j.eq.nobs2) write(76,608) j,stn(ista2),stn(ista1)
  608 format(' Last j obs for ne2:',i2,1x,a4,', No match ne1 sta:',a4)
            goto 150
          endif
          if(intspc(j,ne2,nc).ne.isp1) goto 150
          iq2=iwc(j,ne2,nc)
c      write(76,609) j,stn(ista2),iq2,rdltac(j,ne2,nc),isp1
  609 format('Match j=',i4,' :',a4,', iq2=',i2,', dlta2=',f8.2,
     &  ', isp=',i2)
          if(iq2.gt.mxiqedt) goto 155
          dlta2=rdltac(j,ne2,nc)
          if(dlta2.lt.gminedt) goto 155
          if(dlta2.lt.gminedt0) goto 155
c  usable earthquake pair differential time
          nobedt=nobedt+1
          noepair=noepair+1
          noedt=noedt+1
          if(isp1.eq.0) then
            noedtp=noedtp+1
            edtsec(nobedt,nc)=(secpc(i,ne1,nc)-secoce(ne1,nc))
     &      - (secpc(j,ne2,nc)-secoce(ne2,nc))
          else
            noedts=noedts+1
            edtsec(nobedt,nc)=secpc(i,ne1,nc)-secpc(j,ne2,nc)
          endif
          intsped(nobedt,nc)=isp1
          istoed(1,nobedt,nc)=ista1
          istoed(2,nobedt,nc)=ista2
          iwedt(nobedt,nc)=iq1
          if(iq2.gt.iq1) iwedt(nobedt,nc)=iq2
          wtedtobs(nobedt,nc)=1.0/(2.**(iwedt(nobedt,nc)))
          jobsed(1,nobedt,nc)=i
          jobsed(2,nobedt,nc)=j
          jeved(1,nobedt,nc)=ne1
          jeved(2,nobedt,nc)=ne2
c print check
c      write(76,7601)noedt,edtsec(nobedt,nc),wtedtobs(nobedt,nc),
c     & isp1,dlta1,stn(ista1),
c     & ne1,secpc(i,ne1,nc),secoce(ne1,nc),stn(ista2),ne2,
c     & secpc(j,ne2,nc),secoce(ne2,nc),intspc(j,ne2,nc)
 7601 format(i5,f8.3,f6.2,i3,f8.2,2(1x,a4,i5,2f7.3),i3)
c
c get number of raypath points for edt partials, 2nd eq path
          nrp2=nrp(j)
          nrpedt2=1
          do 60 k=nrp2,2,-1
            k1=k-1
            rx=rp(1,k,j)
            ry=rp(2,k,j)
            rz=rp(3,k,j)
            dx=rp(1,k1,j)-rx
            dy=rp(2,k1,j)-ry
            dz=rp(3,k1,j)-rz
c  compute segment length
            sl=sqrt(dx*dx+dy*dy+dz*dz)
            raylen=raylen+sl
            nrpedt2=nrpedt2+1
            if(raylen.ge.frayedt) goto 65
   60     continue
   65     continue
          nrpedt(1,nobedt,nc)=nrpedt1
          nrpedt(2,nobedt,nc)=nrpedt2
c
c  check nobs against maxedt of this cluster
          if(nobedt.eq.mxedtc) goto 410
c  check nobs against max obs for this event pair
          if(noepair.eq.nstaepair) goto 350
c but allow half S-P
          if((ks.eq.0).and.(noepair.eq.nstaepair2)) goto 350
c
        goto 155
  150   continue
  155 continue
  200 continue
c if just did S-P, now loop through P
      if(ks.eq.0) then
        ks=1
        ne1p=ne1
        goto 20
      endif
  350 continue
  400 continue
      goto 415
  410 write(16,1605) nc,nobedt
 1605 format(' QUIT getting eq-pair dt for cluster=',i6,
     * ', nobedt=mxedtc=',i6)
  415 continue
c  total obs
      kobsedt(nc)=nobedt
      write(16,1640) nc,nobedt,ne1p,ne1s,nec
 1640 format(' Get cedt for Cluster:',i4,',  kobsedt=',i6,
     * ', last edt event=',i4,' (P) ',i4,' (S-P), ncev=',i4)
      return
c***** end of subroutine getcedt *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine getgap(azp,dltp,naz,gap,dltmn)
c  variables for gap calculation
      include 'simul2014_common.inc'
      real azp(maxobs),dltp(maxobs),gap
      integer naz,iazst(maxobs)
c
      dltmn=999.0
c  find minimum delta
      do 20 i=1,naz
      if(dltp(i).lt.dltmn) dltmn=dltp(i)
   20 continue
c
c  sort on azimuth
      call sort(azp,iazst,naz)
c      print *,'azp',(azp(j),j=1,naz)
c      print *,'azp(iazst)',(azp(iazst(j)),j=1,naz)
      gap=0.0
      do 100 i=2,naz
      gapi=azp(iazst(i))-azp(iazst(i-1))
      if(gapi.gt.gap) gap=gapi
  100 continue
      az1=azp(iazst(1))+360.0
      gap1=az1-azp(iazst(naz))
      if(gap1.gt.gap) gap=gap1
c      print *,'gap',gap
c
c***** end of subroutine getgap *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine getrdtsht(ne)
c
c  setup receiver dt "shots" for each earthquake
c       rdisrdt = maximum distance between rdt receivers (input)
c       rdisr12 = distance between 2 specific rdt receivers
c                 This is used for gminrdt
c       frayrdt = distance along ray to perturb = 1.5 * rdisrdt
c       frayrdt1 = alternate distance along ray to perturb (use if smaller)
c                = 0.25 * slant distance of individual path
c       gminrdt = minimum slant path length to all rdt = 4 * rdisrdt (rdisr12)
c       gminrdt0 = absolute minimum slant path length for an observation in rdt
c       mxiqrdt  = max. obs quality for dt observ.
c       wtrdtsht = wt for rdtshot, like wtsht
c
      include 'simul2014_common.inc'
c
      gminrdt0= 15.0
c  Loop through all station pairs, comparing interstation distance and slant 
c  raypath (rdlta(nobs,ne)) length to rdisrdt and gminrdt
c  Set up all rdt
c  Total up rdt obs
      if(ne.gt.1) goto 1
      nordt=0
      nordtp=0
      nordts=0
c
    1 nobrdt=0
      nobs=kobs(ne)
c      print *, ' GETRDTSHT, ne=',ne,', nobs=',nobs
      do 200 i=1,nobs-1
        if((iuse2t.ne.0).and.(iclock(isto(i,ne)).ne.0)) goto 200
        iq1=iw(i,ne)
        if(iq1.gt.mxiqrdt) goto 200
        dlta1=rdlta(i,ne)
        if(dlta1.lt.gminrdt0) goto 200
        frayrdt1=0.25*dlta1
        ista1=isto(i,ne)
c          write(6,602)stn(ista1),iw(i,ne),rdlta(i,ne)
  602     format(a4,',iw=',i4,',dist=',f7.2)
        xr1=stc(1,ista1)
        yr1=stc(2,ista1)
        isp1=intsp(i,ne)
c get number of raypath points for rdt partials
        nrp1=nrp(i)
        nrprdt1=1
        do 50 k=nrp1,2,-1
          k1=k-1
          rx=rp(1,k,i)
          ry=rp(2,k,i)
          rz=rp(3,k,i)
          dx=rp(1,k1,i)-rx
          dy=rp(2,k1,i)-ry
          dz=rp(3,k1,i)-rz
c  compute segment length
          sl=sqrt(dx*dx+dy*dy+dz*dz)
          raylen=raylen+sl
          nrprdt1=nrprdt1+1
c segment length to perturb for this obs
          if(raylen.ge.frayrdt) goto 55
          if(raylen.ge.frayrdt1) goto 55
   50   continue
   55   continue
        do 150 j=i+1,nobs
          if(intsp(j,ne).ne.isp1) goto 150
          if((iuse2t.ne.0).and.(iclock(isto(j,ne)).ne.0)) goto 150
          iq2=iw(j,ne)
          if(iq2.gt.mxiqrdt) goto 150
          dlta2=rdlta(j,ne)
          if(dlta2.lt.gminrdt0) goto 150
          frayrdt2=0.25*dlta2
          ista2=isto(j,ne)
          xr2=stc(1,ista2)
          yr2=stc(2,ista2)
          dx=xr1-xr2
          dy=yr1-yr2
          rdisr12=sqrt(dx*dx + dy*dy)
          if(rdisr12.gt.rdisrdt) goto 150
c          PRINT *,'iq2=',iq2,',dlta2=',dlta2,',distr1r2=',rdisr12
c now check for gminrdt=4*rdisr12
          gminrdt=4.0*rdisr12
          if(dlta1.lt.gminrdt) goto 150
          if(dlta2.lt.gminrdt) goto 150
c  usable receiver pair differential time
          nobrdt=nobrdt+1
c          print *, 'usable receiver pair differential time no=',nobrdt
          if(isp1.eq.0) nordtp=nordtp+1
          if(isp1.eq.1) nordts=nordts+1
          nordt=nordt+1
          rdtsec(nobrdt,ne)=secp(i,ne)-secp(j,ne)
          intsprd(nobrdt,ne)=isp1
          istord(1,nobrdt,ne)=ista1
          istord(2,nobrdt,ne)=ista2
          iwrdt(nobrdt,ne)=iq1
          if(iq2.gt.iq1) iwrdt(nobrdt,ne)=iq2
          wtrdtobs(nobrdt,ne)=1.0/(2.**(iwrdt(nobrdt,ne)))
          jobsrd(1,nobrdt,ne)=i
          jobsrd(2,nobrdt,ne)=j
c get number of raypath points for rdt partials
          nrp2=nrp(j)
          nrprdt2=1
          do 60 k=nrp2,2,-1
            k1=k-1
            rx=rp(1,k,j)
            ry=rp(2,k,j)
            rz=rp(3,k,j)
            dx=rp(1,k1,j)-rx
            dy=rp(2,k1,j)-ry
            dz=rp(3,k1,j)-rz
c  compute segment length
            sl=sqrt(dx*dx+dy*dy+dz*dz)
            raylen=raylen+sl
            nrprdt2=nrprdt2+1
            if(raylen.ge.frayrdt) goto 65
            if(raylen.ge.frayrdt2) goto 65
   60     continue
   65     continue
          nrprdt(1,nobrdt,ne)=nrprdt1
          nrprdt(2,nobrdt,ne)=nrprdt2
c  check against maxobs
          if(nobrdt.eq.maxobs) goto 210
  150   continue
  200 continue
      goto 215
  210 write(16,1605) ne,nobrdt
 1605 format(' QUIT getting rec-pair dt for event=',i6,
     * ', nobrdt=maxobs=',i6)
  215 continue
      kobsrdt(ne)=nobrdt
c  total obs
      return
c***** end of subroutine getrdtsht *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine h12(mode,lpivot,l1,m,u,iue,up,c,ice,icv,ncv,
     2               c1,ic1e,ic1v,nc1v)
c  c.l. lawson and r.j. hanson, jpl
c  from "solving least squares problems"
c  modified by c.h. thurber, spring 1980
c  construction and application of a single
c  householder transformation...  q =  i + u*(u**t)/b
      dimension u(m),c(*),c1(*)
      double precision sm,b
      one=1.
c
      if (0.ge.lpivot.or.lpivot.ge.l1.or.l1.gt.m) return
      cl=abs(u(lpivot))
      if (mode.eq.2) go to 60
c
c  construct the transformation
c
      do 10 j=l1,m
         cl=amax1(abs(u(j)),cl)
   10 continue
c      if (cl) 130,130,20
      if(cl.le.0.0) return
   20 clinv=one/cl
      sm=(dble(u(lpivot))*clinv)**2
      do 30 j=l1,m
         sm=sm+(dble(u(j))*clinv)**2
   30 continue
c
c  convert dble prec sm to sngl prec sm1
c
      sm1=sm
      cl=cl*sqrt(sm1)
c      if (u(lpivot)) 50,50,40
      if(u(lpivot).gt.0.0) cl = -cl
   40 cl=-cl
   50 up=u(lpivot)-cl
      u(lpivot)=cl
      go to 70
c
c  apply the transformation i+u*(u**t)/b to matrices c & c1
c
c   60 if (cl) 130,130,70
   60 if(cl.le.0.0) return
   70 if (ncv.le.0) return
      b=dble(up)*u(lpivot)
c
c  b must be nonpositive here.  if b=0 return
c
c      if (b) 80,130,130
      if(b.ge.0.0) return
   80 b=one/b
      i2=1-icv+ice*(lpivot-1)
      incr=ice*(l1-lpivot)
      do 120 j=1,ncv
         i2=i2+icv
         i3=i2+incr
         i4=i3
         sm=c(i2)*dble(up)
         do 90 i=l1,m
            sm=sm+c(i3)*dble(u(i))
            i3=i3+ice
   90    continue
c         if (sm) 100,120,100
         if(sm.ne.0.0) then
  100      sm=sm*b
           c(i2)=c(i2)+sm*dble(up)
           do 110 i=l1,m
             c(i4)=c(i4)+sm*dble(u(i))
             i4=i4+ice
  110      continue
         endif
  120 continue
      i2=1-ic1v+ic1e*(lpivot-1)
      incr=ic1e*(l1-lpivot)
      do 220 j=1,nc1v
         i2=i2+ic1v
         i3=i2+incr
         i4=i3
         sm=c1(i2)*dble(up)
         do 190 i=l1,m
            sm=sm+c1(i3)*dble(u(i))
            i3=i3+ic1e
  190    continue
c         if (sm) 200,220,200
         if(sm.ne.0.0) then
  200      sm=sm*b
           c1(i2)=c1(i2)+sm*dble(up)
           do 210 i=l1,m
             c1(i4)=c1(i4)+sm*dble(u(i))
             i4=i4+ic1e
  210      continue
         endif
  220 continue
  130 return
c***** end of subroutine h12 *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine hyrot  (nmax,x, y, cs, sn, n)
      integer n
      real    x(nmax), y(nmax), cs, sn
c
c
      real    xx
      integer j
c
c
      do 10 j = 1, n
         xx = x(j)
         x(j) = xx*cs + y(j)*sn
         y(j) = y(j)*cs - xx*sn
   10 continue
      return
c***** end of subroutine hyrot *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine input1
c  this routine reads in control parameters and number of eq's
c
c  declaration statements:
      character*24 date24
c
c  common block variables:
      include 'simul2014_common.inc'
c
      call fdate(date24)
      write(16,1605) date24(21:24),date24(4:16)
 1605 format(' Computation began at ',a4,a13)
      write(16,9937)maxpar,mxpari,maxev,maxsta,maxobs
 9937  format(' Program simul2014 (06-Jan-2014) Solves for Vp and ',
     2 'Vp/Vs; Input data is P travel-time and S-P time.',/,
     3 '       Can vary relative weighting of S-P times in ',
     4 'hypocenter location.',/,
     * /,'       OR can solve for Qp using tstar as input',/,/,
     5 '       Allows fixed  and linked nodes (up to',i7,
     *  ' parameters, up to ',i5,' solution parameters);'/,
     6 '       up to ',i4,' events, ',i5,' stations, ',i4,
     *  ' observations per event',
     7 /,'       Psuedo-bending; simul2004 has earth-flattening ',
     8 /,'        transformation in vel3',
     9 '; also improved range of initial arc check for "moho" i3d=3.',
     & /,'simul2010 allows receiver differential times and ',
     & 'earthquake differential times.')
      read(1,*,err=999) neqs,nsht,nbls,wtsht,kout,kout2,kout3
      read(1,*) ntel,wttel,ntobmin,reste1,reste2,reste3
 1033 format(2i3,f4.1)
        if(neqs+nsht+nbls+ntel.le.maxev) goto 20
        write(16,22)
 22     format('0too many events for program arrays.')
        stop
  20    continue
      read(1,*) iuserdt, rdisrdt,facrayr,mxiqrdt,wtrdtsht,res4
      read(1,*) iuseclu, nclu,edisedt,facraye,mxiqedt,nstaepair,wtepdt,
     2  wrmscemx,gminedt0
      if((iuseclu.eq.1).and.(nclu.gt.maxclu)) then
        write(16,23) nclu, maxclu
   23   format(' nclu ',i4,' greater than array, maxclu ',i4,
     2     /,'***** STOP *****')
        write(6,23) nclu,maxclu
        stop
      endif
      read(1,*) nitloc,wtsp,eigtol,rmscut,zmin,dxmax,rderr,ercof
      read(1,*) hitct,dvpmx,dvsmx,cvpmax,cvpmin
      read(1,*) idmp,(vdamp(j),j=1,4),stepl
      read(1,*) ires,i3d,nitmax,snrmct,ihomo,rmstop,ifixl
      read(1,*) delt1,delt2,res1,res2,res3
      read(1,*) ndip,iskip,scale1,scale2
      read(1,*) xfac,tlim,nitpb(1),nitpb(2)
      read(1,*) kttfor,iusep,iuses,invdel,iuse2t
      read(1,*) iuseq, dqmax, qdamp, qrmax, qmin
      if((nitmax.eq.-1).and.((neqs.ne.0).or.(nbls.ne.0))) then
        write(16,1609) neqs,nsht,nbls
 1609   format(/,'******* INPUT ERROR: neqs=',i3,'nsht=',i3,'nbls=',
     &  i3,/,'For synthetic, data must be read in as Shots',
     &  /,'******* STOP *******',/)
        stop
      endif
      if((iuseq.eq.1).and.((neqs.ne.0).or.(nbls.ne.0))) then
        write(16,1611) neqs,nsht,nbls
 1611   format(/,'******* INPUT ERROR: neqs=',i3,'nsht=',i3,'nbls=',
     &  i3,/,'For Q inversion, t* data must be read in as Shots',
     &  /,'******* STOP *******',/)
        stop
      endif
      nebs=neqs+nbls+nsht
      if((iuseq.eq.1).and.(iuses.eq.1)) then
        write(16,1610)
 1610   format(/,'******* INPUT ERROR: Q inversion is only for P',
     2  '*******              Resetting iuses to 0',/)
      endif
c use receiver pair diff time shots
      if(iuserdt.gt.0) then
        frayrdt=facrayr*rdisrdt
      endif
 1051 format(2f7.2,2f5.2)
 1011 format(4i3,f5.3)
 1021 format(i3,3f5.3)
      write(16,1030)
      write(16,1040) kout,kout2,kout3
c      write(16,1031)
c      write(16,1041) neqs,nsht,nbls,wtsht,nitloc,wtsp,zmin,eigtol,
c     2 rmscut,hitct,dvpmx,dvsmx
      write(16,1620)
      write(16,1621) neqs,nsht,nbls,wtsht,nitloc,wtsp,zmin,eigtol,
     2 rmscut,hitct
      write(36,1620)
      write(36,1621) neqs,nsht,nbls,wtsht,nitloc,wtsp,zmin,eigtol,
     2 rmscut,hitct
 1620 format(/,' neqs nsht nbls wtsht',
     2 ' nitloc wtsp   zmin   eigtol  rmscut hitct')
 1621 format(1x,i4,2i5,2x,f4.1,1x,i3,f8.2,f7.2,f8.3,f9.3,f6.0)
      write(16,1629) ntel,wttel,ntobmin,reste1,reste2,reste3
      write(36,1629) ntel,wttel,ntobmin,reste1,reste2,reste3
 1629 format(' ntel wttel ntobmin reste1 reste2 reste3',/,
     2  1x,i4,f6.1,i7,3f7.2)
c
      write(16,1638) kttfor
 1638 format(' kttfor=',i2,' (1=cnsp, 2=phs, 3=6-letter-sta phs)')
c
c changes for 2 origin times
      if(iuse2t.eq.0) then
        nparhy=4
        write(16,1640) iuse2t,nparhy
 1640   format(' iuse2t=',i2,' one origin time, nparhy=',i2)
        write(36,1640) iuse2t,nparhy
      else
        nparhy=5
        write(16,1641) iuse2t,nparhy
        write(36,1641) iuse2t,nparhy
 1641   format(' iuse2t=',i2,' Use 2 origin times (2 master',
     2  ' clocks), nparhy=',i2)
      endif
c
c option to allow receiver-pair differential times
      write(16,1673) iuserdt,rdisrdt,facrayr,mxiqrdt,wtrdtsht,res4
      write(36,1673) iuserdt,rdisrdt,facrayr,mxiqrdt,wtrdtsht,res4
 1673 format('  iuserdt   rdisrdt   facrayr mxiqrdt  wtrdtsht  res4',/,
     2 i7,f11.2,f8.2,i8,f11.2,f8.2)
      if(iuserdt.gt.0) write (16,1674)
 1674 format(' Use receiver-pair differential times ')
c option to allow clusters and earthquake-pair differential times
      write(16,1683) iuseclu, nclu,edisedt,facraye,mxiqedt,nstaepair,
     2  wtepdt,wrmscemx
      write(36,1683) iuseclu, nclu,edisedt,facraye,mxiqedt,nstaepair,
     2  wtepdt,wrmscemx
 1683 format(' iuseclu  nclu   edisedt   facraye mxiqedt nstaepair ',
     2 'wtepdt wrmscemx',/,
     & 2i7,f11.2,f8.2,i6,i10,f9.2,f8.3)
      if(iuseclu.gt.0) write (16,1684)
 1684 format(' Use clusters and earthquake-pair differential times ')
      if((res4.eq.0.0).and.((iuserdt.gt.0).or.(iuseclu.gt.0))) then
        res4=res1
        write(16,1685) res4
 1685   format(' *** ADJUSTING input: res4 must be >0 for rdt or edt',
     &   '   set res4=res1=',f6.3,' ***')
      endif
      write(16,1675) iuseq
 1675 format('  iuseq=',i2)
      if(iuseq.eq.0) then
       if(iuses.gt.0) then
         if(iusep.eq.1) then
           write(36,1676) iusep,iuses
           write(16,1676) iusep,iuses
 1676      format(/,' Both P-velocity and Vp/Vs are included in ',
     2    'inversion  (iusep=',i2,', iuses=',i2,')')
         else
           write(16,1677) iusep,iuses
           write(36,1677) iusep,iuses
 1677      format(/,' Only S velocities are included in ',
     2    'inversion  (iusep=',i2,', iuses=',i2,')')
         endif
       else
         write(16,1678) iusep,iuses
 1678    format(/,' Only P velocities are included in ',
     2    'inversion  (iusep=',i2,', iuses=',i2,')')
c  put station damping in 2nd position if no Vs used
         vdamp(2)=vdamp(3)
       endif
       write(16,1623)
       write(16,1624) dvpmx,dvsmx,cvpmax,cvpmin,idmp,(vdamp(j),j=1,4)
       write(36,1623)
       write(36,1624) dvpmx,dvsmx,cvpmax,cvpmin,idmp,(vdamp(j),j=1,4)
      else
        dvpmx=dqmax
        vdamp(1)=qdamp
        iusep = 1
        iuses = 0
         write(16,1679) iuseq,iusep,iuses
         write(36,1679) iuseq,iusep,iuses
 1679   format(/,' Inversion only for Q ',
     2      ' (iuseq=',i2,' ,iusep=',i2,' ,iuses=',i2,')')
        if(ntel.gt.0) write(16,1689) vdamp(4)
        if(ntel.gt.0) write(36,1689) vdamp(4)
 1689 format('telpad_damp=',f10.3)
        if(invdel.eq.0) then
          write(16,1625)
          write(16,1626) dvpmx,idmp,vdamp(1),qmin,qrmax
          write(36,1625)
          write(36,1626) dvpmx,idmp,vdamp(1),qmin,qrmax
c sta corr for t*
        else
c  put station damping in 2nd position 
          vdamp(2)=vdamp(3)
          write(16,1635)
          write(16,1636) dvpmx,idmp,vdamp(1),vdamp(2)
          write(36,1635)
          write(36,1636) dvpmx,idmp,vdamp(1),vdamp(2)
        endif
      endif
 1623 format(/,'    dvpmx  dvpvsmx corrvpmax corrvpmin ',
     2 'idmp     vpdamp    vpvsdmp    stadamp   telpad_damp')
 1624 format(f9.2,f9.3,2f9.2,i4,f13.3,4f11.3)
 1625 format(/,'     dqmax    idmp        qdamp    qmin    Qin-max')
 1626 format(f9.2,i8,f14.6,f8.2,f10.0)
 1635 format(/,'     dqmax    idmp        qdamp       stadamp')
 1636 format(f9.2,i8,2f14.6)
c
      if(invdel.ne.0) write(16,1080) invdel
 1080 format(/,' Station delays included in inversion',
     2 ' (invdel=',i2,')',/)
c
c      write(16,1032)
c      write(36,1032)
c      write(16,1042) idmp,(vdamp(j),j=1,4),ires,
c     2 nitmax,snrmct,dxmax,rderr,ercof,ihomo,rmstop,ifixl,delt1,delt2,
c     3 res1,res2,res3
c      write(46,1042) idmp,(vdamp(j),j=1,3),ires,
c     2 nitmax,snrmct,dxmax,rderr,ercof,ihomo,rmstop,ifixl,delt1,delt2,
c     3 res1,res2,res3
c
      write(16,1627)
      write(36,1627)
      write(16,1628) ires,nitmax,snrmct,dxmax,rderr,ercof,ihomo,rmstop,
     2 ifixl,delt1,delt2,res1,res2,res3
      write(36,1628) ires,nitmax,snrmct,dxmax,rderr,ercof,ihomo,rmstop,
     2 ifixl,delt1,delt2,res1,res2,res3
 1627 format('  ires nitmax snrmct dxmax rderr ercof')
 1628 format(i3,i7,f9.5,3f6.2,
     * /,' # its. for 1-d vel. model = ',i3,/,' rms for term. = ',f5.3
     * ,/,' fix locations for iterations = ',i3
     * ,/,' distance weighting:',2f7.1,'; residual weighting:',3f7.2)
      write(16,1061) stepl
      write(16,1155)
 1155 format(' parameters for approximate ray tracer',
     * /,6x,'ndip iskip scale1 scale2')
      write(16,1151) ndip,iskip,scale1,scale2
 1151 format(5x,i3,i6,2f7.2)
 1061 format(' step length for integration:',f6.3)
 1030 format(/, ' control parameters',/,' kout kout2 kout3')
 1040 format(1x,i3,i5,i6)
 1031 format(/,' neqs nsht nbls wtsht',
     2 ' nitloc wtsp  zmin eigtol  rmscut hitct',
     3 ' dvpmx dvpvsmx')
 1032 format(/,' idmp    vpdamp   vpvsdmp   stadamp ires ',
     2 'nitmax snrmct dxmax rderr ercof')
 1041 format(1x,3i4,4x,f4.1,1x,i3,f8.2,f6.2,f6.3,f9.3,f6.0,2f6.2)
 1042 format(i3,f12.2,2f10.2,i3,i7,f9.5,3f6.2,
     * /,' # its. for 1-d vel. model = ',i3,/,' rms for term. = ',f5.3
     * ,/,' fix locations for iterations = ',i3
     * ,/,' distance weighting:',2f7.2,'; residual weighting:',3f5.2)
      write(16,1166) i3d,xfac,tlim,nitpb(1),nitpb(2)
 1166 format('  parameters for pseudo-bending:',/,
     * '  i3d xfac   tlim  nitpb(1) nitpb(2)',/,i4,f6.2,f7.4,i6,i10)
c  increment iuses so it can be an index to arrays
      iuses=iuses+1
      ddlt=1.0/(delt2-delt1)
      if(res3.gt.res2) then
        dres12=0.98/(res2-res1)
        dres23=0.02/(res3-res2)
      else
        dres12=1.0/(res2-res1)
        dres23=0.0
      endif
      if(ntel.gt.0) then
        if(reste3.gt.reste2) then
          dreste12=0.98/(reste2-reste1)
          dreste23=0.02/(reste3-reste2)
        else
          dreste12=1.0/(reste2-reste1)
          dreste23=0.0
        endif
      endif
      nevt=neqs+nsht+nbls+ntel
      return
  999 write(16,1699)
 1699 format('error in input1, probably forgot kout3')
c***** end of subroutine input1 *****
      stop
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine input2
c  this routine reads in the station list, sets up the
c  coordinate system, and calculates the stations' cartesian coordinates.
c  (reads from file02 )
c  subroutines required: setorg; disto;
c  common block variables:
      common/shortd/ xltkm,xlnkm,rota,nzco,xlt,xln,snr,csr,
     2 onorth,oeast,rotadg,cmerid,rlatdo,rlondo
      include 'simul2014_common.inc'
c
      character*1 cns, cew
c
      write(16,2000)
 2000 format(/,'  origin :  latitude   longitude   rotation   nzco')
c  read in center of coordinates and rotation angle
c  for converting from lat lon to cartesian (and visa versa).
c  nzco=1 for New Zealand (lat S and long E) NZTM2000,
c    with WGS84 (as in geonet stations DELTA)
c  nzco=2 Alaska (state plane zone4), already had TM and WGS84
c  nzco=3 input central meridian,
c  nzco=0 use origin latitude for central meridian
c  nzco=4 use NZMG as in previous simul versions.
c  nzco=5 use short distance conversion as in original Thurber code
c  note NZ and AK have preset cmerid, which will override cmerid input
c
      read(2,*) ltdo,oltm,lndo,olnm,rota,nzco,cmerid
      if(ltdo.ge.0) then
        rlatdo=float(ltdo)+oltm/60.0
      else
        rlatdo=float(ltdo)-abs(oltm/60.0)
      endif
      if(lndo.ge.0) then
        rlondo=float(lndo)+olnm/60.0
      else
        rlondo=float(lndo)-abs(olnm/60.0)
      endif
      if(nzco.eq.0) cmerid=rlondo
      if(nzco.eq.1) cmerid=173.
c      PRINT *,'input2 rlondo,rlatdo ',rlondo,rlatdo,' cmerid ',cmerid
cdep FH changes assume input origin lat and lon degrees
c       are given as negative for S and E ???
c       So make negative if nzco=1
      if((nzco.eq.1).or.(nzco.eq.4)) then
        ltdo= -1.0*ltdo
        lndo= -1.0*lndo
      endif
c
c  If inverting station corrections, setup output file22 to
c  be used as input file in future runs.
c      if((invdel.eq.0).or.(kout.lt.2)) goto 5
c now write station files always, and have it include TALLY info
      if(kout.lt.2) goto 5
      open(unit=22,file='newstns')
      rewind 22
      write(22,2001) ltdo,oltm,lndo,olnm,rota,nzco
 2001 format(i3,1x,f5.2,i4,1x,f5.2,f8.2,i3)
    5 continue
      write(16,2002) ltdo,oltm,lndo,olnm,rota,nzco
 2002 format(12x,i3,1x,f5.2,2x,i4,1x,f5.2,3x,f7.2,i3)
c  set up short-distance conversion factors, given center of coords
      call setorg(ltdo,oltm,lndo,olnm)
      if(iuseq.eq.0) then
        write(16,2004)
 2004   format('   station     latitude   longitude   elev',
     2   '     x      y      z   pdl s-pdl nfixst iclock')
      else
        write(16,2005)
 2005   format('   station     latitude   longitude   elev',
     2   '     x      y      z   tstar-dl  nfixst iclock')
      endif
c  read in number of stations
      read(2,*) nsts
      if(nsts.le.maxsta)go to 40
      write(16,41)
 41   format('0Too many stations for input arrays.')
      stop
 40   continue
c  read in station list
c  Now read to end of file, skipping sta with lat and lon=0
cfhek  changes made to read east-west, north-south (read in cns, cew and use it)
cfh    - read and format statement changes
      j=1
      nstsi=0
   15 continue
      if(kttfor.eq.3) goto 16
      if(kttfor.eq.4) goto 17
        if(iuseq.eq.0) then
          read(2,2007,end=30) stn(j),ltds(j),cns,sltm(j),lnds(j),cew,
     *    slnm(j),ielev,pdl(j),sdl(j),nfixst(j),iclock(j)
 2007     format(2x,a4,i2,a1,f5.2,i4,a1,f5.2,i5,2f5.2,2i3)
        else
          read(2,2017,err=1017,end=30)  stn(j),ltds(j),cns,sltm(j),
     *    lnds(j),cew,slnm(j),ielev,pdl(j),nfixst(j),iclock(j)
 2017     format(2x,a4,i2,a1,f5.2,i4,a1,f5.2,i5,e11.3,i2,i3)
          goto 1019
 1017     read(2,2018,end=30)  stn(j),ltds(j),cns,sltm(j),
     *    lnds(j),cew,slnm(j),ielev,pdl(j),nfixst(j),iclock(j)
 2018     format(2x,a4,i2,a1,f5.2,i4,a1,f5.2,i5,e12.3,i2,i3)
 1019     continue
        endif
        if((ltds(j).eq.0).and.(lnds(j).eq.0)) goto 15
        z=-ielev*1.0e-3
c need cns and cew in station file
        if((cns.eq.' ').or.(cew.eq.' ')) then
          write(16,1657) stn(j), cns,cew
 1657     format('*** Problem need Lat N or S, Long E or W station: ',
     2     a4,' cns:',a1,' cew;',a1,' Make S E for nzco=1')
          if(nzco.eq.1) cns='S'
          if(nzco.eq.1) cew='E'
        endif
cfhek change made to recognize North South and East West
        if(cns.eq.'S'.or.cns.eq.'s') then
           ltds(j)=-1.*ltds(j)
           sltm(j)=-1.*sltm(j)
        endif
        if(cew.eq.'E'.or.cew.eq.'e') then
           lnds(j)=-1.*lnds(j)
           slnm(j)=-1.*slnm(j)
        endif
        goto 18
c 6 letter station code (Japan) had decimal lat lon
   16   read(2,2008,end=30) stn6(j),rlat,rlon,elevm,pdl(j),sdl(j),
     2    nfixst(j)
 2008   format(a6,f9.4,f10.4,f8.1,2f5.2,2i3)
        ielev=nint(elevm)
        z=-elevm*1.0e-3
c 5-letter station with net (NCEDC,IRIS)
   17   if(kttfor.eq.4) then
          if(iuseq.eq.0) then
            read(2,2003,end=30) stn5(j),net(j),rlat,rlon,ielev,
     2       pdl(j),sdl(j),nfixst(j)
          else
            read(2,2006,end=30) stn5(j),net(j),rlat,rlon,ielev,
     2       pdl(j),nfixst(j)
          endif
 2003     format(a5,3x,a2,1x,f10.4,1x,f10.4,1x,i5,2f5.2,i3)
 2006     format(a5,3x,a2,1x,f10.4,1x,f10.4,1x,i5,e11.3,i2)
          z=-ielev*1.0e-3
        endif
         ltds(j)=ifix(rlat)
         sltm(j)=60.0*(rlat-float(ltds(j)))
         lnds(j)=ifix(rlon)
         slnm(j)=60.0*(rlon-float(lnds(j)))
   18    if(nfixst(j).eq.0) nstsi=nstsi+1
c  calculate cartesian coordinates of station
         call disto(ltds(j),sltm(j),lnds(j),slnm(j),x,y)
         if((kout2.lt.4).or.(j.eq.1)) then
cfhek
           if(ltds(j).ge.0.and.sltm(j).ge.0.) then
                cns='N'
                latdeg=iabs(ltds(j))
                alatmin=abs(sltm(j))
           else
                cns='S'
                latdeg=iabs(ltds(j))
                alatmin=abs(sltm(j))
           endif
           if(lnds(j).ge.0.and.slnm(j).ge.0.) then
                cew='W'
                londeg=iabs(lnds(j))
                alonmin=abs(slnm(j))
           else
                cew='E'
                londeg=iabs(lnds(j))
                alonmin=abs(slnm(j))
           endif
           if(kttfor.ne.3) then
             if(iuseq.eq.0) then
              if(kttfor.ne.4) then
               write(16,2009) j,stn(j),latdeg,cns,alatmin,londeg,cew,
     2         alonmin,ielev,x,y,z,pdl(j),sdl(j),nfixst(j),iclock(j)
              else
               write(16,2013) j,stn5(j),net(j),latdeg,cns,alatmin,
     2         londeg,cew,alonmin,ielev,x,y,z,pdl(j),sdl(j),nfixst(j)
              endif
             else
              if(kttfor.ne.4) then
               write(16,2019) j,stn(j),latdeg,cns,alatmin,londeg,cew,
     2         alonmin,ielev,x,y,z,pdl(j),nfixst(j),iclock(j)
              else
               write(16,2014) j,stn5(j),net(j),latdeg,cns,alatmin,
     2         londeg,cew,alonmin,ielev,x,y,z,pdl(j),nfixst(j)
              endif
            endif
           else
            write(16,2010) j,stn6(j),latdeg,cns,alatmin,londeg,cew,
     2        alonmin,ielev,x,y,z,pdl(j),sdl(j),nfixst(j),iclock(j)
           endif
         endif
 2009    format(1X,i5,3x,a4,1x,i3,a1,f5.2,2x,i4,a1,f5.2,2x,i5,
     *      1x,f7.2,f8.2,f7.2,2f5.2,i3,i7)
 2019    format(1X,i5,3x,a4,1x,i3,a1,f5.2,2x,i4,a1,f5.2,2x,i5,
     *      1x,f7.2,f8.2,f7.2,e11.3,i2,i7)
 2010    format(1X,i5,1x,a6,1x,i3,a1,f5.2,2x,i4,a1,f5.2,2x,i5,
     *      1x,f7.2,f8.2,f7.2,2f5.2,i3,i7)
 2013     format(1x,i5,1x,a5,1x,a2,1x,i3,a1,f5.2,2x,i4,a1,f5.2,2x,i5,
     2      1x,f7.2,f8.2,f7.2,2f5.2,i3)
 2014     format(1x,i5,1x,a5,1x,a2,1x,i3,a1,f5.2,2x,i4,a1,f5.2,2x,i5,
     2      1x,f7.2,f8.2,f7.2,e11.3,i2)
c         call latlon(x,y,latsta,xltsta,lonsta,xlnsta)
c         write(16,1666) latsta,xltsta,lonsta,xlnsta
c 1666    format(2x,'checking latlon ',2(i4,f6.2))
c  store station coordinates
         stc(1,j)=x
         stc(2,j)=y
         stc(3,j)=z
         j=j+1
         goto 15
   30 continue
      nstsrd=j-1
      if(nsts.ne.nstsrd) then
        write(16,2098) nstsrd, nsts,nstsrd
 2098   format('*************************************',/,
     2  'Number station read in',i5,' NOT EQUAL Number of ',
     3  'stations specified',i5,/,
     4  ' Continue with ',i5,/,
     4  '*************************************')
        nsts=nstsrd
      endif
      if(kout2.lt.4) return
c
c  kout2=4, condensed printout
   50 continue
      write(16,1605) nsts
 1605 format('  number of stations read in Input2 =',i5,
     2 ' (print 1st and last only)')
      j=nsts
cfhek
      if(ltds(j).ge.0.and.sltm(j).ge.0.) then
         cns='N'
         latdeg=iabs(ltds(j))
         alatmin=abs(sltm(j))
      else
         cns='S'
         latdeg=iabs(ltds(j))
         alatmin=abs(sltm(j))
      endif
      if(lnds(j).ge.0.and.slnm(j).ge.0.) then
         cew='W'
         londeg=iabs(lnds(j))
         alonmin=abs(slnm(j))
      else
         cew='E'
         londeg=iabs(lnds(j))
         alonmin=abs(slnm(j))
      endif
      if(kttfor.lt.3) then
        write(16,2009) j,stn(j),latdeg,cns,alatmin,londeg,cew,
     2        alonmin,ielev,x,y,z,pdl(j),sdl(j),nfixst(j),iclock(j)
      else
       if(kttfor.eq.3) then
        write(16,2010) j,stn6(j),latdeg,cns,alatmin,londeg,cew,
     2        alonmin,ielev,x,y,z,pdl(j),sdl(j),nfixst(j),iclock(j)
       else
        write(16,2013) j,stn5(j),net(j),latdeg,cns,alatmin,
     2    londeg,cew,alonmin,ielev,x,y,z,pdl(j),sdl(j),nfixst(j)
       endif
      endif
      return
c
c***** end of subroutine input2 *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine input3
c  this routine reads in the initial velocity model in the
c  form of velocity specified on a uniform but not
c  necessarily evenly spaced grid of points
c  (reads from file03 )
c
c  common block variables:
      include 'simul2014_common.inc'
c
c  declaration statements:
      integer ixf(maxpar),iyf(maxpar),izf(maxpar)
      integer ixm(maxpar),iym(maxpar),izm(maxpar)
      integer ixl(maxpar),iyl(maxpar),izl(maxpar)
c      dimension casc(50)
      character*1 casc(50)
      character*110 line
      character*1 vtype(2)
      parameter(zero=0.0,izero=0)
c
      ierror=0
      vtype(1)='P'
      vtype(2)='S'
c
c  put letters for each group of linked nodes
      data  casc/'A','B','C','D','E','F','G','H','I','J','K','L',
     2  'M','N','P','Q','R','S','T','U','V','W','X','Y','Z',
     3  'a','b','c','d','e','f','g','h','i','j','k','l','m','n',
     4  'p','q','r','s','t','u','v','w','x','y','z'/
c      print *,'casc ',casc
c
c  for this version the gridpoints can be unevenly spaced
c  the origin of the coordinate system is at (x,y,z)=(0,0,0)
c  which will not in general correspond to the point
c  xn(1),yn(1),zn(1).
c  xn,yn,zn should be factors of bld (ie: a.0 for bld=1.0 or a.b for bld=0.1)
c
c input the number of gridpoints in x, y and z directions
c  and bld factor (1.0 or 0.1 km) used to set up velocity interpolation grid
c  Now allow comment lines at the beginning of the velocity file
c  comment lines begin with character 'c'
   10 read(3,3111) line
 3111 format(a85)
      if((line(1:1).eq.'c').or.(line(1:1).eq.'C')) goto 10
c      read(line,3002) bld,nx,ny,nz
      read(line,*) bld,nx,ny,nz
 3002 format(f4.1,3i3)
      if((bld.ne.1.0).and.(bld.ne.0.1)) then
        write(16,1625) bld
 1625   format(/, '******** STOP *********, bld must be 1.0 or 0.1,
     2   not ',f6.2)
      endif
        atemp=iuses*(nx-2)*(ny-2)*(nz-2)
        if(atemp.le.maxpar)goto 40
        write(16,42) maxpar,atemp
 42     format('0Too many nodes for program array sizes. maxpar=',i10,
     2    ', atemp=',i10)
        stop
 40     continue
c
      do 123 k=1,atemp
      imerge(k)=0
      jequal(k)=0
 123  continue
c
c  input the x grid, y grid, and z grid
cfh read in free format (makes life easier...)
c        read(3,3004) (xn(i),i=1,nx)
c        read(3,3004) (yn(i),i=1,ny)
c        read(3,3004) (zn(i),i=1,nz)
        read(3,*) (xn(i),i=1,nx)
        read(3,*) (yn(i),i=1,ny)
        read(3,*) (zn(i),i=1,nz)
 3003 format(3i3)
 3004 format(20f6.1)
c
      write(16,3005) bld,nx,ny,nz
 3005 format(//,' velocity grid size:',/,
     * 'bld =',f4.1,5x,' nx =',i4,5x,'ny =',i4,5x,'nz =',i3)
c
cfh give all these numbers the same format
      write(16,3006) (xn(i),i=1,nx)
cfh 3006 format(/,' xgrid',/,3x,12f7.1,8f6.1)
 3006 format(/,' xgrid',/,3x,20f7.1)
      write(16,3007) (yn(i),i=1,ny)
cfh 3007 format(/,' ygrid',/,3x,12f7.1,8f6.1)
 3007 format(/,' ygrid',/,3x,20f7.1)
      write(16,3008) (zn(i),i=1,nz)
cfh 3008 format(/,' zgrid',/,3x,8f6.1,12f7.1/)
 3008 format(/,' zgrid',/,3x,20f7.1/)
c
c  read in which nodes to have fixed velocity
c  end with blank line
      i=1
   50 read(3,3003) ixf(i),iyf(i),izf(i)
      if(ixf(i).le.0) goto 60
      i=i+1
      goto 50
   60 continue
      inf=i-1
c
c  start cht 1998
c  lines moved followed by new code
c  compute total number of gridpoints (nodes)
      nodes=nx*ny*nz
      nxy=nx*ny
      nx2=nx-2             ! number non-edge nodes in row
      nxy2=nx2*(ny-2)      ! number non-edge nodes in layer
      nz2=nz-2
      nodes2=nz2*nxy2
c  peripheral nodes
      nx1=nx-1
      ny1=ny-1
      nz1=nz-1
c
c  read in which nodes have "linked" velocity, "master" node first
c  followed by "linked" nodes - end each group with blank line
c  end with another blank line.  if no linked nodes, just include
c  a blank line
c
      im=1
      il=1
      ican=0
 52   read(3,3003) ixm(im),iym(im),izm(im)
c       print *,ixm(im),iym(im),izm(im)
      if(ixm(im).le.0) goto 62
      izmim=izm(im)
      if (izmim.gt.nz) izmim=izmim-2
      mnode=(izmim-2)*nxy2+(iym(im)-2)*nx2+ixm(im)-1
c
c  link type - constant (1) or linear (2)?
      read(3,*) ltype(mnode)
c
      if(cnode(mnode).ne.'0') then
        write(16,1686) cnode(mnode)
        write(6,1686) cnode(mnode)
 1686   format(' *-*-* ERROR velocity input.  This node has',
     2  ' already been ',/,' *-*-* assigned cnode= ',a1)
        ierror=1
        write(16,1696) ixm(im),iym(im),izm(im),mnode,ltype(mnode)
 1696   format(/,'  master node, ix, iy, iz, #, type: ',3i4,i6,i4)
      endif
      cnode(mnode)='M'
      im=im+1
      ican=ican+1
      if(ican.gt.50) ican=1
      canode(mnode)=casc(ican)
c      print *,'canode(',mnode,')=',canode(mnode),' ican=',ican,
c     2  ' casc(ican)',casc(ican)
c
 54   read(3,3003) ixl(il),iyl(il),izl(il)
c       print *,ixl(il),iyl(il),izl(il)
      if(ixl(il).le.0) goto 52
      izlil=izl(il)
      if (izlil.gt.nz) izlil=izlil-2
      lnode=(izlil-2)*nxy2+(iyl(il)-2)*nx2+ixl(il)-1
      if(cnode(lnode).ne.'0') then
        write(16,1686) cnode(mnode)
        ierror=1
        write(16,1697) ixl(il),iyl(il),izl(il),lnode
 1697   format('  linked node, x, y, z, #: ',3i4,i6)
      endif
      il=il+1
c
      imast=im-1
      ilink=il-1
      imerge(lnode)=1
      jequal(lnode)=mnode
      if(ltype(mnode).eq.1) then
        cnode(lnode)='C'
      else
        cnode(lnode)='L'
c simul2010 allows more than 1 direction
cc Note that the velocity adjustment only does simple linear adjustment in 
cc  x,y OR z.  So check to see if applicable.
c        if((ixl(ilink).eq.ixm(imast)).and.
c     2    (iyl(ilink).eq.iym(imast))) goto 58
c        if((ixl(ilink).eq.ixm(imast)).and.
c     2    (izl(ilink).eq.izm(imast))) goto 58
c        if((iyl(ilink).eq.iym(imast)).and.
c     2    (izl(ilink).eq.izm(imast))) goto 58
c          write(16,1618) ixl(ilink),iyl(ilink),izl(ilink),ixm(imast),
c     2      iym(imast),izm(imast)
c 1618     format(' *** ERROR in INPUT ***',/,' *** Nodes',3i3,' and',
c     2    3i3,' should not be Linearly linked ***',/,
c     3    ' *** Changing to Constant Linking ***')
c          ltype(mnode)=1
c          cnode(lnode)='C'
c  58    continue
      endif
      canode(lnode)=canode(mnode)
c
      goto 54
 62   continue
      write(16,1616) imast,ilink
 1616 format(/,' number of master and linked nodes:',2i7)
c
c  end cht 1998
c
c  now read in the velocity values
   65 write(16,3101)
c     do 38 kv=1,iuses
         kv=1
         do 37 k=1,nz
            k2=k + (kv-1)*nz
            write(16,3015) k,vtype(kv),zn(k)
            do 36 j=1,ny
cfh               read(3,3011) (vel(i,j,k2),i=1,nx)
               read(3,*) (vel(i,j,k2),i=1,nx)
               write(16,3013) (vel(i,j,k2),i=1,nx)
   36       continue
   37    continue
c  38 continue
c CHANGE FOR VP/VS INVERSION
      if((iuses.eq.2).or.(iuseq.eq.1)) then
        do 100 k=1,nz
          if(iuseq .eq. 0) then
              write(16,3016) k,zn(k)
          else
              write(16,3017) k,zn(k)
          endif
           do 99 j=1,ny
             if(iuseq .eq. 0) then
cfh                 read(3,3011) (vpvs(i,j,k),i=1,nx)
                 read(3,*) (vpvs(i,j,k),i=1,nx)
                 write(16,3013) (vpvs(i,j,k),i=1,nx)
             else
               read(3,*) (qval(i,j,k),i=1,nx)
               write(16,3014) (qval(i,j,k),i=1,nx)
             endif
   99     continue
  100  continue
c  compute Vs from Vp and Vp/Vs or compute 1/tstar
        kv=2
       if(iuseq .eq. 0) then
           do 120 k=1,nz
              ks=k+nz
              write(16,3015) k,vtype(kv),zn(k)  
              do 115 j=1,ny
                 do 110 i=1,nx
                    vel(i,j,ks)=vel(i,j,k)/vpvs(i,j,k)
  110           continue
                 write(16,3013) (vel(i,j,ks),i=1,nx)
  115        continue
  120     continue
       else
         do 140 k=1,nz
            ks=k+nz
            write(16,3018) k,zn(k)
            do 135 j=1,ny
               do 130 i=1,nx
                  vel(i,j,ks)=vel(i,j,k)*qval(i,j,k)
  130           continue
                 write(16,3014) (vel(i,j,ks),i=1,nx)
  135        continue
  140     continue
       endif
      endif
c
 3013 format(22f6.2)
 3014 format(20f7.1)
 3015 format(/,' layer',i3,5x,a1,' velocity',10x,'z =',f7.1)
 3016 format(/,' layer',i3,5x,'Vp/Vs',10x,'z =',f7.1)
 3017 format(/,' layer',i3,5x,'Q',10x,'z =',f7.1)
 3018 format(/,' layer',i3,5x,'Q * Vp',10x,'z =',f7.1)
 3011 format(20f5.2)
 3101 format(//,' velocity values on three-dimensional grid')
c  Number of medium parameters to invert for
      npar=nodes2*iuses
      print *,'nodes2 ',nodes2,', iuses ',iuses,', npar',npar,
     2 ', invdel',invdel
      nparv=npar
c      PRINT *, 'nodes2=',nodes2,',npar=',npar,',nparv=',nparv
      if(invdel.ne.0) npar=(npar+nsts*iuses)
      nparvs=npar
      if(ntel.ne.0) npar=npar+2*ntel
c
c  Check to see whether medium parameters fit within array sizes
      if(nparv.gt.maxpar) goto 980
c
cfhdmep
c get number of Vp and Vp/Vs nodes that are free in the inversion
c    nodes reduced by fixednodes
      nparpi=nodes2
      nparsi=nodes2
c  fix specified nodes by setting nfix(k)=1, else=0
      if((inf.eq.0).and.(ilink.eq.0)) goto 496
        do 70 i=1,inf
           iizf=izf(i)-2
c  if s velocity node
cfhdmep           if(izf(i).gt.nz) iizf=izf(i)-4
           if(izf(i).gt.nz) then
              iizf=izf(i)-4
              nparsi=nparsi-1
           else
              nparpi=nparpi-1
           endif
c
           k=iizf*nxy2 + (iyf(i)-2)*nx2 + (ixf(i)-1)
           nfix(k)=1
           if(cnode(k).ne.'0') then
             write(16,1686) cnode(mnode)
             ierror=1
             write(16,1681) mnode,ixf(i),iyf(i),izf(i),k
 1681        format('input3 fixed node',i8,' ixf,iyf,izf:',3i5,
     2       ' node number:',i8)
           endif
           cnode(k)='F'
   70   continue
c
       write(16,1611)
 1610 format(/,' velocity FIXED at the following nodes(1):')
 1611 format(/,' VELOCITY INVERSION GRID    0=free, F=Fixed',/,
     2  '   M=Master, C=Constant Pert. Link, ',
     3  'L=Linear Pert. Link')
  311 do 495 kv=1,iuses
         nz1=nz-1
         ny1=ny-1
         do 320 k=2,nz1
          if(iuseq.eq.0) then
            if(kv.eq.1) write(16,1009) k,vtype(kv),zn(k)
 1009       format(/,' layer',i3,5x,a1,'-velocity nodes',
     2         10x,'z =',f7.1)
            if(kv.eq.2) write(16,3016) k,zn(k)
          else
            write(16,3017) k,zn(k)
          endif
          write(16,1008) (i,i=2,nx1)
 1008     format(' Y X',i2,39i3)
            kk=k+(kv-1)*nz2
            do 310 j=2,ny1
               n1=(kk-2)*nxy2+(j-2)*nx2+1
               n2=n1+nx2-1
               write(16,1006) j,(cnode(i),i=n1,n2)
c               write(16,1005) (nfix(i),i=n1,n2)
  310       continue
  320    continue
c 1005 format('    ',18i6)
 1006 format(' ',i2,40(2x,a1))
  495 continue
c
c  Now print out letters for each set of linked nodes
       write(16,1621)
 1621 format(/,' map of LINKED nodes with each set represented by ',
     2 'letter',/,'  - = nonlinked (free or fixed)')
      do 595 kv=1,iuses
         nz1=nz-1
         ny1=ny-1
         do 520 k=2,nz1
          if(iuseq.eq.0) then
            if(kv.eq.1) write(16,1009) k,vtype(kv),zn(k)
            if(kv.eq.2) write(16,3016) k,zn(k)
          else
            write(16,3017) k,zn(k)
          endif
          write(16,1008) (i,i=2,nx1)
            kk=k+(kv-1)*nz2
            do 510 j=2,ny1
               n1=(kk-2)*nxy2+(j-2)*nx2+1
               n2=n1+nx2-1
               write(16,1006) j,(canode(i),i=n1,n2)
c               write(16,1005) (nfix(i),i=n1,n2)
  510       continue
  520    continue
  595 continue
  496 continue
c
c  ndexfx: index from full nodes to nodes reduced by fixed (invert nodes)
c  mdexfx: index from inversion solution nodes to full velocity nodes
      in=0
c
c  start cht 1998
      infl=inf+ilink
c
      write(16,6161) inf,ilink,infl
 6161 format(/,' number of fixed, linked, fixed+linked nodes: ',3i8)
      write(16,6162) maxpar
 6162 format(/,8x,'maxpar=',i8)
c
      do 80 n=1,nparv
c
c  remove fixed and linked nodes from inversion solution node set
c
         if(nfix(n).eq.1) goto 80
c  start of imerge if-then-else
         if(imerge(n).eq.0) go to 888
c
c  calculate x-y-z indices of linked velocity grid
         k=(n-1)/nxy2+2
         j=2+(n-1+(2-k)*nxy2)/nx2
         i=1+n+nx2*(2-j)+nxy2*(2-k)
         if(k.ge.nz) k=k+2      ! if s velocity node
c
c  calculate x and z indices of "master" velocity grid
         nma=jequal(n)
c
         km=(nma-1)/nxy2+2
         jm=2+(nma-1+(2-km)*nxy2)/nx2
         im=1+nma+nx2*(2-jm)+nxy2*(2-km)
         if(km.ge.nz) km=km+2      ! if s velocity node
c
cDEP It seems better to have "constant" link for velocity derivatives
cDEP  to result in constant perturbation between master and linked,
cDEP  rather than just constant velocity.  So do not change initial
cDEP  model here. Of course a constant velocity initial model could
cDEP  be input if that was desired.
cDEPc  start of constant/linear link if-then-else
cDEP      if (ltype(nma).eq.1) then
cDEP      write(16,2626) i,j,k,n,im,jm,km,nma
cDEP 2626 format('  setting value for constant-type linked node: ',
cDEP     &  4i4,2x,4i4)
cDEPc
cDEP       if(iuseq.eq.0) then
cDEP           vel(i,j,k)=vel(im,jm,km)
cDEP           if (k.gt.nz)
cDEP     &     vpvs(i,j,k-nz)=vpvs(im,jm,km-nz)
cDEPc
cDEP       else
cDEP           vel(i,j,k+nz)=vel(im,jm,km+nz)
cDEP       endif
cDEPc
cDEP      else
cDEPc   continuing constant/linear link if-then-else
cDEP       il=i
cDEP       jl=j
cDEP      if (iuseq.ne.0) then
cDEP       k=k+nz
cDEP       km=km+nz
cDEP      endif
cDEP       kl=k
cDEPc
cDEP       if (i.ne.im) then
cDEP       il=2*i-im
cDEP       vel(i,j,k)=vel(im,jm,km)+(vel(il,jl,kl)-vel(im,jm,km))*
cDEP     & (xn(i)-xn(im))/(xn(il)-xn(im))
cDEP       endif
cDEPc
cDEP       if (j.ne.jm) then
cDEP       jl=2*j-jm
cDEP       vel(i,j,k)=vel(im,jm,km)+(vel(il,jl,kl)-vel(im,jm,km))*
cDEP     & (yn(i)-yn(im))/(yn(il)-yn(im))
cDEP       endif
cDEPc
cDEP       if (k.ne.km) then
cDEP       kl=2*k-km
cDEP       vel(i,j,k)=vel(im,jm,km)+(vel(il,jl,kl)-vel(im,jm,km))*
cDEP     & (zn(i)-zn(im))/(zn(il)-zn(im))
cDEP       endif
cDEPc
cDEP      write(16,2627) i,j,k,n,im,jm,km,nma,il,jl,kl,vel(i,j,k)
cDEP 2627 format('  setting value for linear-type linked nodes: ',/,
cDEP     &  4i4,2x,4i4,2x,3i4,f7.2)
cDEPc
cDEP      if (k.gt.nz) then
cDEP       if (i.ne.im) then
cDEP       il=2*i-im
cDEP       vpvs(i,j,k)=vpvs(im,jm,km)+(vpvs(il,jl,kl)-vpvs(im,jm,km))*
cDEP     & (xn(i)-xn(im))/(xn(il)-xn(im))
cDEP       endif
cDEPc
cDEP       if (j.ne.jm) then
cDEP       jl=2*j-jm
cDEP       vpvs(i,j,k)=vpvs(im,jm,km)+(vpvs(il,jl,kl)-vpvs(im,jm,km))*
cDEP     & (yn(i)-yn(im))/(yn(il)-yn(im))
cDEP       endif
cDEPc
cDEP       if (k.ne.km) then
cDEP       kl=2*k-km
cDEP       vpvs(i,j,k)=vpvs(im,jm,km)+(vpvs(il,jl,kl)-vpvs(im,jm,km))*
cDEP     & (zn(i)-zn(im))/(zn(il)-zn(im))
cDEP       endif
cDEP      endif
cDEPc
cDEP      endif
c
      go to 80
c
c  end of linked node section
c
  888 continue
c  add to index if not linked or fixed
c
         in=in+1
         ndexfx(n)=in
         mdexfx(in)=n
c         write(16,1888) in,n
c 1888    format(' adding inversion node number ',
c     &   i5,' corresponding to gridpoint number ',i5)
c
   80 continue
      inf2=nparv-in
      if(inf2.eq.infl) goto 85
c
c  end cht 1998
c
        write(16,1615) infl,nparv,in,inf2
 1615   format(/,' **** number of fixed and linked nodes input,',i8,
     2  ', does not equal velocity nodes,',i8,', minus invert',
     3  ' nodes,',i8,'.  Continue with inf=',i8,' ****',/)
        infl=inf2
   85 continue
      nparvi=nparv-infl
      npari=nparvi
      if(invdel.eq.0) then
        nparvsi=nparvi
        goto 95
      else
        npari=nparvi+nstsi*iuses
      endif
c  also set indices if station delays are included in inversion
      i1=nparv+1
      do 90 i=i1,nparvs
         is=i-nparv
c s-delay
         if(is.gt.nsts) is=is-nsts
         if(nfixst(is).eq.1) goto 90
         in=in+1
         ndexfx(i)=in
         mdexfx(in)=i
   90 continue
      npari=in
      nparvsi=npari
   95 continue
c tele path delay index, 1st P, then S-P
      if(ntel.eq.0) goto 97
      i1=nparvs+1
      PRINT *,'input3 i1,ntel,npar,in',i1,ntel,npar,in
      PRINT *,'     nparvs ',nparvs
      do 96 i=i1,npar
        in=in+1
        ndexfx(i)=in
        mdexfx(in)=i
C      WRITE(56,*)'   i,ndexfx(i),mdexfx(ndexfx(i)) ',i,ndexfx(i),
C     2 mdexfx(ndexfx(i))
   96 continue
      npari=in
   97 write(16,1620) npar,nparv,npari,nparvi
      WRITE(6,1620) npar,nparv,npari,nparvi
 1620 format(' INPUT3:npar,nparv,npari,nparvi',4i8)
c  Check to see whether medium parameters fit within array sizes
      if(npari.gt.mxpari) goto 990
      if(npar.gt.maxpar) goto 995
c  Stop if there was an error reading in fixed and linked nodes
      if(ierror.ne.0) then
        write(16,1683)
        write(6,1683)
 1683   format(/,'STOP SIMUL2014! , Error in velocity input file ')
        stop
      endif
c
c  Set up an array which is used to point to node indices, for any x,y,z
      call bldmap
c
      return
c
  980 continue
      write(6,1683)
      write(16,1698) nparv,maxpar
 1698 format(/,'  ****** STOP ******',/,i8,' velocity nodes, program',
     2 ' arrays only allow',i6)
      stop
  990 continue
      write(6,1683)
      write(16,1699) npari,mxpari
 1699 format(/,'  ****** STOP ******',/,i8,' parameters to invert for',
     2 ', program arrays only allow',i6)
      stop
  995 continue
      write(6,1683)
      write(16,1695) npar,maxpar
 1695 format(/,'  ****** STOP ******',/,i8,' parameters, arrays are ',
     2 'for ',i8)
      stop
c
c***** end of subroutine input3 *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
c
      subroutine input4(n,nfile)
c  this routine reads in the p travel times for all stations observing
c  the current event.  first trial hypocentral parameters are read, then
c  station names and weights and travel times.  arrivals at stations not
c  in the station list are discarded; normalized weights are calculated.
c  (reads earthquakes from file04, reads shots from file07 )
c  subroutine required: disto;
c
c  declaration statements:
c  local variables:
      real dep,x,y,pwt,tt(6),sta(6)
      character*4 rmki(6)
c      dimension nwav(2),rmki(6),ip(6),tt(6),sta(6),is(6)
c      character*4 sta(6),chk,chksp,cblank,cblnk0
      character*4 chk,chksp,cblank,cblnk0
      character*1 cns,cew,is(6)
      character*6 sta6(6)
      character*5 sta5(6)
      character*3 cmp(6)
      character*2 netin(6)
      integer ip(6),nwav(2)
      character*110 line
c
c  common block variables:
      include 'simul2014_common.inc'
c
      data  blank,iblank,blank0 /4h    ,4h    ,4h0   /
      data  cblank,cblnk0 /4h    ,4h0   /
c
      if (n.gt.1) go to 1
c  print out heading for event cards
      write(16,4000)
 4000 format(//,' trial event locations:',/,
     2 5x,'n    origin time     latitude longitude  depth',
     3 ' mag',5x,'x      y      z',4x,'nob np ns seco2')
      nobt=0
      nobtp=0
      nobts=0
      nordt=0
      nordtp=0
      nordts=0
    1 nobs=0
      nwav(1)=0
      nwav(2)=0
      wsum=0
c  read in event card
c  check for extra travel time card
   10 read(nfile,4111,end=99) line,chk
 4111 format(a85,t1,a4)
c
    5 if((chk.eq.cblank).or.(chk.eq.cblnk0)) goto 10
c
   15 continue
      if(iuse2t.eq.0) then
         read(line,4001,err=5098) iyrmo(n),iday(n),ihr(n),
     2   mino(n),seco(n),ltde(n),cns,eltm(n),lnde(n),cew,elnm(n),
     3   dep,rmag(n)
      else
         read(line,4001,err=5098) iyrmo(n),iday(n),ihr(n),
     2   mino(n),seco(n),ltde(n),cns,eltm(n),lnde(n),cew,elnm(n),
     3   dep,rmag(n),seco2(n)
         if(seco2(n).eq.0.0) then
           seco2(n)=seco(n)
           write(16,1638) seco2(n)
 1638      format(' +++ seco2=0.0, reset to seco ',f7.2,' +++')
         endif
      endif
 4001 format(a4,a2,1x,a2,i2,1x,f5.2,i3,a1,f5.2,1x,i3,a1,f5.2,2f7.2,
     2  t69,f6.2)
      if(cns.eq.'S'.or.cns.eq.'s') then
        ltde(n)=-1.*ltde(n)
        eltm(n)=-1.*eltm(n)
      endif
      if(cew.eq.'E'.or.cew.eq.'e') then
        lnde(n)=-1.*lnde(n)
        elnm(n)=-1.*elnm(n)
      endif
c
      if (iyrmo(n).eq.iblank) go to 10
c  calculate event location in cartesian coordinate system
c      call disto(2,n,x,y)
      call disto(ltde(n),eltm(n),lnde(n),elnm(n),x,y)
c  store event coordinates
      evc(1,n)=x
      evc(2,n)=y
      evc(3,n)=dep
c
c  read in travel times to stations
c  format including phase-remarks
   20 continue
      read(nfile,4007) line
 4007 format(a110)
c
      if(kttfor.eq.1) goto 21
      if(kttfor.eq.3) goto 121
      if(kttfor.eq.4) goto 122
      if(kttfor.ne.2) then
        write(16,1621) kttfor
 1621   format('kttfor=',i2,' SHOULD BE 1=cnsp, 2=phs, 3=phs sta6lt',
     &  /,' Continue with kttfor=1, assume cnsp format ***')
        kttfor=1
        goto 21
      endif
c phs format (SCB), one travel-time per line
c SCB change to read output from "pha2cnsp.pl"
      jend=1
      j=1
      read(line,4992,err=5099) sta(j),rmki(j),tt(j)
      read(line,4993,err=5099) is(j),ip(j)
 4992 format(a4,1x,a4,1x,f9.4)
 4993 format(6x,a1,1x,i1)
      goto 25
c
c phs with 6-letter station (for Japan)
  121 jend=1
      j=1
      read(line,4994,err=5099) sta6(j),rmki(j),tt(j)
      read(line,4995,err=5099) sta(j),is(j),ip(j)
 4994 format(a6,1x,a4,1x,f9.4)
 4995 format(a4,4x,a1,1x,i1)
      goto 25
c
c phs with 5-letter station with net (for NCEDC,IRIS)
c Allow cmp for possible future use
  122 jend=1
      j=1
      read(line,4990,err=5099) sta5(j),cmp(j),netin(j),rmki(j),tt(j)
      read(line,4991,err=5099) sta(j),is(j),ip(j)
C      PRINT *,'sta5(j),cmp(j),netin(j),rmki(j),tt(j)',
C     2 sta5(j),cmp(j),netin(j),rmki(j),tt(j)
C      PRINT *,'sta(j),is(j),ip(j) ',sta(j),' ',is(j),ip(j),
C     2 ' blank= ',blank
 4990 format(a5,1x,a3,1x,a2,1x,a4,1x,f9.4)
 4991 format(a4,10x,a1,1x,i1)
      goto 25
c  
c  read in travel times to stations in groups of six
   21 jend=6
c CHT old format
      if(iuseq.eq.0) then
        read(line,4996,err=921) (sta(j),rmki(j),tt(j),j=1,6),
     2  (is(j),ip(j),j=1,6)
 4996   format(6(a4,a2,f6.2),t1,6(4x,a1,i1,6x))
        goto 25
  921   continue
c cnsp format
        read(line,4006,err=22) (sta(j),rmki(j),tt(j),j=1,6),
     2  (is(j),ip(j),j=1,6)
      else
c tstar data for iuseq =1
        read(line,4029,err=22) (sta(j),rmki(j),tt(j),j=1,6),
     2  (is(j),ip(j),j=1,6)
      endif
 4006 format(6(2a4,f6.2),t1,6(5x,a1,1x,i1,6x))
 4029 format(6(2a4,f10.4),t1,6(5x,a1,1x,i1,10x))
      goto 25
c Alberto Michelini (LBL) format, could be used for tstar
   22 read(line,4008,err=5099) (sta(j),rmki(j),tt(j),j=1,5),
     2  (is(j),ip(j),j=1,5)
 4008 format(5(2a4,f7.4),t1,5(5x,a1,1x,i1,7x))
      jend=5
c
c  loop over the six (or jend) readings
c  terminate if station is a blank (end of observations for event)
c  skip if weight is incorrect in some way
c  search list for current station, and index it if in list
c  also calculate p-arrival time and unnormalized weight
c  and increment number of observations for event
c
   25 continue
c  blank line separates the events
      if ((sta(1).eq.blank).or.(sta(1).eq.blank0)) go to 50
      do 30 j=1,jend
         if(ip(j).ge.4) goto 30
         pwt=1.0/(2.**(ip(j)))
         if (pwt.le.0.0) go to 30
c only use if labeled P or S
         if((is(j).ne.'P').and.(is(j).ne.'p').and.
     2     (is(j).ne.'S').and.(is(j).ne.'s')) goto 30
c search station list
         do 33 k=1,nsts
          if(kttfor.lt.3) then
            if (sta(j).eq.stn(k)) go to 35
          else
            if(kttfor.eq.3) then
              if(sta6(j).eq.stn6(k)) goto 35
            else
              if((sta5(j).eq.stn5(k)).and.(netin(j).eq.net(k))) goto 35
            endif
          endif
   33    continue
c         if(kout2.gt.1) goto 30
         if(sta(j).eq.blank .or. sta(j).eq.blank0) then
c            write(16,1609)
 1609       format(' Blank-station observation ignored')
         else
           if(kttfor.lt.3) then
             write(16,1610) sta(j),iyrmo(n),iday(n),ihr(n),
     *       mino(n),seco(n)
 1610        format(' ** WARNING:  Observed station "',a4,
     *       '" (event ',a4,a2,1x,a2,i2,1x,f5.2,
     *       ') MISSING in station list--observation ignored **')
           else
             if(kttfor.eq.3) then
               write(16,1611)sta6(j),iyrmo(n),iday(n),ihr(n),
     *         mino(n),seco(n)
 1611          format(' ** WARNING:  Observed station "',a6,
     *         '" (event ',a4,a2,1x,a2,i2,1x,f5.2,
     *         ') MISSING in station list--observation ignored **')
            else
              write(16,1631) sta5(j),netin(j),iyrmo(n),iday(n),
     *         ihr(n),mino(n),seco(n)
 1631          format(' ** WARNING:  Observed station "',a5,1x,a2,
     *         '" (event ',a4,a2,1x,a2,i2,1x,f5.2,
     *         ') MISSING in station list--observation ignored **')
            endif
           endif
         endif
         go to 30
c  add to set of observing stations if on list
   35    continue
c  don't use arrivals further away than delt2
         dx=evc(1,n)-stc(1,k)
         dy=evc(2,n)-stc(2,k)
         dz=evc(3,n)-stc(3,k)
c
c  Print out a file of event and receiver coordinates if kout2=0
c  and nitmax=0
      if((kout2.eq.0).and.(nitmax.eq.0)) then
c  Only print out P
      if((is(j).eq.'P').or.(is(j).eq.'p')) then
        eltd=float(ltde(n))+eltm(n)/60.0
        elnd=float(lnde(n))+elnm(n)/60.0
        elnd=-1.0*elnd
        sltd=float(ltds(k))+sltm(k)/60.0
        slnd=float(lnds(k))+slnm(k)/60.0
        slnd=-1.0*slnd
        write(38,3802) eltd,elnd,dep,seco(n),sltd,slnd,stc(3,k),
     2 pdl(k),tt(j)
 3802  format(5x,2f10.4,2f7.2,2x,2f10.4,2f7.2,2x,f8.2)
      endif
      endif
c
      rdelta=sqrt((dx*dx)+(dy*dy)+(dz*dz))
      if((iuseq.eq.1).and.(nitmax.gt.-1)) then
c make quick check of high Q to throw out bad readings
c high Q is an input parameter
         Qobs = rdelta/(6.0*tt(j))
         if(Qobs.ge.qrmax) then
          if(kttfor.lt.3) then
            write(16,1622) stn(k),tt(j),qobs,qrmax
          else
            if(kttfor.eq.3) then
              write(16,1623) stn6(k),tt(j),qobs,qrmax
            else
              write(16,1633) stn5(k),net(k),tt(j),qobs,qrmax
            endif
          endif
 1622      format(' **** NOT USING ',a4,' t* obs=',f7.4,'Qobs=',e12.4,
     2      ' GREATER THAN QRMAX=',f7.0)
 1623      format(' **** NOT USING ',a6,' t* obs=',f7.4,'Qobs=',e12.4,
     2      ' GREATER THAN QRMAX=',f7.0)
 1633      format(' **** NOT USING ',a5,1x,a2,' t* obs=',f7.4,'Qobs=',
     2      e12.4,' GREATER THAN QRMAX=',f7.0)
           goto 30
         endif
       endif
c
c  Use epicentral distance instead of slant distance
         delta=sqrt((dx*dx)+(dy*dy))
         if(delta.lt.delt2) goto 36
         if(kout2.eq.0) write(16,1602) delta,delt2,stn(k),stn6(k),n
 1602    format(' **** delta,delt2=',f15.2,f8.2,2x,a4,a6,
     &     ' not used for event',i4)
         goto 30
c        input iuses=0 (changed to iuses+1 in input1) to disable
c           S-wave processing.
   36    if(iuses.eq.1.and.is(j).eq.'S') goto 30
         if((iusep.eq.0).and.((is(j).eq.'P').or.(is(j).eq.'p')))
     *      goto 30
         nobs=nobs+1
c check for too many observations 
         if(nobs.gt.maxobs) then
            write(6,1603) iyrmo(n),iday(n),ihr(n),mino(n),maxobs
 1603       format('  INPUT ERROR, Too many readings for event:',
     2      a4,a2,1x,a2,i2,/,'    CONTINUE with first',i5,
     3      'observations')
            write(16,1603) iyrmo(n),iday(n),ihr(n),mino(n),maxobs
            nobs=maxobs
            goto 50
         endif
         dlta(nobs,n)=delta
         rdlta(nobs,n)=rdelta
         iw(nobs,n)=ip(j)
         isto(nobs,n)=k
         if((iuse2t.eq.0).or.(iclock(k).eq.0)) then
           secp(nobs,n)=tt(j)+seco(n)
         else
           secp(nobs,n)=tt(j)+seco2(n)
         endif
         rmk(nobs,n)=rmki(j)
c  here are additions by wap to read in s times. If the time is
c  an s time, is(j) will equal 'S'. If this is true, the
c  proper adjustment will be made for the station delay, and
c  a 1 will be put in intsp(nobs,n).
         intsp(nobs,n)=0
         if(rmk(nobs,n).eq.cblank) rmk(nobs,n)=' P  '
         if((is(j).ne.'S').and.(is(j).ne.'s')) goto 37
c  intsp(nobs,n) tells whether the nobs'th observation of the
c  arrival time is a P (0), or an S (1).
         intsp(nobs,n)=1
c** S-P CHANGE
c   For S-P we do not want to add origin time
         secp(nobs,n)=tt(j)
         if(rmk(nobs,n).eq.cblank) rmk(nobs,n)=' S  '
  37     continue
         kv=intsp(nobs,n)+1
         nrd(k,kv)=nrd(k,kv)+1
         nwav(kv)=nwav(kv)+1
         wt(nobs,n)=pwt
         wsum=wsum+pwt
   30 continue
c  continue reading stations and travel times
      go to 20
c
   50 continue
c  end of travel time readings for current event
c  check that there are at least four readings for this event
c  (does not matter for shots)
      if(nobs.eq.0) goto 90
      if((nobs.lt.4).and.(n.le.neqs)) go to 90
      if((nobs.lt.2).and.(n.le.neb)) go to 90
cfhek 
      if(ltde(n).ge.0.and.eltm(n).ge.0.) then
        cns='N'
        latdeg=iabs(ltde(n))
        alatmin=abs(eltm(n))
      else
        cns='S'
        latdeg=iabs(ltde(n))
        alatmin=abs(eltm(n))
      endif
      if(lnde(n).ge.0.and.elnm(n).ge.0.) then
        cew='W'
        londeg=iabs(lnde(n))
        alonmin=abs(elnm(n))
      else
        cew='E'
        londeg=iabs(lnde(n))
        alonmin=abs(elnm(n))
      endif
      if(iuse2t.eq.0) then
        write(16,4004) n,iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
     2  latdeg,cns,alatmin,londeg,cew,alonmin,dep,rmag(n),x,y,dep,
     3  nobs,nwav(1),nwav(2)
C        WRITE(6,4004) n,iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
C     2  latdeg,cns,alatmin,londeg,cew,alonmin,dep,rmag(n),x,y,dep,
C     3  nobs,nwav(1),nwav(2)
      else
        write(16,4004) n,iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
     2  latdeg,cns,alatmin,londeg,cew,alonmin,dep,rmag(n),x,y,dep,
     3  nobs,nwav(1),nwav(2),seco2(n)
      endif
cek      write(16,4004) n,iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
cek     2 ltde(n),cns,eltm(n),lnde(n),cew,elnm(n),dep,rmag(n),x,y,dep,
cek     3 nobs,nwav(1),nwav(2)
 4004 format(2h *,i5,1x,a4,a2,1x,a2,i2,f6.2,i3,a1,f5.2,i4,a1,f5.2,
     2 f7.2,f5.2,1x,3f7.2,1x,3i3,f6.2)
c
c  normalize reading weights
c   change in normalizing factor 19-jul-1983
      wfac=nobs/wsum
c     wfac=sqrt(nobs/wsum)
c  check for too many readings
      if(nobs.le.maxobs) goto 55
        write(16,1698) n,nobs,maxobs,maxobs
 1698   format(' **** ERROR **** in event ',i4,', TOO MANY ',
     2  'OBSERVATIONS, nobs=',i6,' array size allows',i6,/,
     2  '   continuing with',i6,' observations',/)
        nobs=maxobs
   55 kobs(n)=nobs
      kobps(1,n)=nwav(1)
      kobps(2,n)=nwav(2)
      nobt=nobt+nobs
      nobtp=nobtp+nwav(1)
      nobts=nobts+nwav(2)
      if(n.gt.netemp) then
        nobtex=nobtex+nobs
        nobt=nobt+nint((wtsht-1.0)*(float(nobs)))
      else
        nobteq=nobteq+nobs
      endif
      do 60 j=1,nobs
      wt(j,n)=wt(j,n)*wfac
   60 continue
c
c  FIX TO CONVERT S TO S-P DATA, CHT
c  For input S-P data characters 2,3 of phase remark are 'Sp'
c     ie: "polarity" field is 'p'
c
      do 987 j=1,nobs
      if (intsp(j,n).eq.1) then
      write(chksp,1620) rmk(j,n)
 1620 format(a4)
      if(chksp(2:3).eq.'Sp') goto 987
      if(iuse2t.eq.0) then
        secp(j,n)=secp(j,n)+seco(n)
      else
        secp(j,n)=secp(j,n)+seco2(n)
      endif
c
c  find matching P observation
      imatch=0
      do 986 jj=1,nobs
      if ((isto(jj,n).eq.isto(j,n)).and.
     & (intsp(jj,n).eq.0)) then
        imatch=imatch+1
        if(imatch.eq.1) then
          secs=secp(j,n)
          secp(j,n)=secp(j,n)-secp(jj,n)
          rmk(j,n)(2:3)='Sp'
c         WRITE(16,1619) stn(isto(j,n)),rmk(j,n)
 1619 format('input4 matching P for S-P ',a4,'  rmk ',a4)
c         WRITE(06,1619) stn(isto(j,n)),rmk(j,n)
        endif
        if(imatch.gt.1) then
          secsp2=secs-secp(jj,n)
          dsp=abs(secsp2-secp(j,n))
          if(dsp.gt.0.55) then
            write(16,9876) stn(isto(jj,n)),stn(isto(j,n)),
     &      secp(jj,n),secp(j,n),secs-secp(jj,n),dsp
 9876       format('  WARNING - second S-P match! ',a4,1x,a4,1x,3f8.3,
     2      ' diff',f5.2,' >0.55s  SKIP')
            wt(j,n)=0.
          endif
        endif
      endif
 986  continue
      endif
 987  continue
c
      return
c
c  not enough readings for event n
   90 continue
      write(16,4009)  iyrmo(n),iday(n),ihr(n),mino(n)
 4009 format(' *** event ',a4,a2,1x,a2,i2,' should be discarded - ',
     2 'too few observations ***',/,' *** SKIPPING TO NEXT EVENT ***')
      if(nfile.eq.7) then       ! shot datafile
        nsht=nsht-1
      else
        if(nfile.eq.4) then     !earthquake datafile
          if(ifixl.le.0) then
            neqs=neqs-1
          else
            nbls=nbls-1
          endif
          netemp=netemp-1
        else                    ! blast datafile
          nbls=nbls-1
          nbtemp=nbtemp-1
        endif
        neb=neb-1
      endif
      nevt=nevt-1
      goto 1
c  ran out of event cards - error]
   99 write(16,4099) nfile
 4099 format('**** end of data file',i2,' - neqs too large ****')
      if(nfile.eq.4) then
c  assume input error if neqs was less than 10% off, and continue with smaller neqs
c  assume other error,possibly wrong file, if more than 10% off, then stop
        if((abs(netemp-n)).gt.(0.10*netemp)) goto 199
        netemp=n-1
        neb=nbtemp+netemp
        nevt=neb+nsht
          if(ifixl.le.0) then
            neqs=n-1
          else
            nbls=neb
          endif
        write(16,1699) netemp
 1699   format(' continue with neqs= ',i4)
        return
      else
c blasts
        if(nfile.eq.8) then
          if(abs(nbtemp-(n-neqs)).gt.(0.10*nbtemp)) goto 199
          nbtemp=n-1-neqs
          if(ifixl.le.0) then
            nbls=n-1-neqs
          else
            nbls=n-1
          endif
          write(16,1696) nbtemp
 1696     format(' continue with nbls= ',i4)
          return
c shots
        else
          nshtin=n-1-neqs-nbls
          if((abs(nsht-nshtin)).gt.(0.10*nsht)) goto 199
          nsht=nshtin
          write(16,1697) nsht
 1697     format(' continue with nsht= ',i4)
          return
        endif
      endif
  199 continue
      write(16,1690)
 1690 format(' stop since neqs more than 10% off, or problem '
     2 ,'was nbls or nsht')
      stop
 5098 write(16,5050) n
      write(16,4012) nchar,line,chk
 4012 format(2x,i5,' characters in line',/,2x,a85,/,2x,'first 4 ',
     2 'char = ',a4)
cfhek
      if(ltde(n).ge.0.and.eltm(n).ge.0.) then
        cns='N'
        latdeg=iabs(ltde(n))
        alatmin=abs(eltm(n))
      else
        cns='S'
        latdeg=iabs(ltde(n))
        alatmin=abs(eltm(n))
      endif
      if(lnde(n).ge.0.and.elnm(n).ge.0.) then
        cew='W'
        londeg=iabs(lnde(n))
        alonmin=abs(elnm(n))
      else
        cew='E'
        londeg=iabs(lnde(n))
        alonmin=abs(elnm(n))
      endif
cek
      if(iuse2t.eq.0) then
        write(16,4011) iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
     2  latdeg,cns,alatmin,londeg,cew,alonmin,dep,rmag(n)
      else
        write(16,4011) iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
     2  latdeg,cns,alatmin,londeg,cew,alonmin,dep,rmag(n),seco2(n)
      endif
cek      write(16,4011) iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
cek     2 ltde(n),eltm(n),lnde(n),elnm(n),dep,rmag(n)
 4011 format(1x,a4,a2,a3,i2,1x,f5.2,1x,i3,a1,f5.2,1x,i3,a1,f5.2,2f7.2,
     2  f6.2)
c
      goto 10
c     stop
 5099   write(16,5050) n,line
      write(16,4015) (sta(j),is(j),ip(j),tt(j),j=1,6)
 4015 format(1x,6(a4,a1,i1,f6.2))
 5050   format('0error in data in input4. event= ',i8,/,'** ',a110)
        stop
c****** end of subroutine input4 ******
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
c
      subroutine input9(nc,nfile)
c this routine reads in a cluster of earthquakes from fort.9 (nfile)
c cluster starts with 'BEGIN CLUSTER', ends with 'END CLUSTER'
c cluster earthquake data is read in same format as other 
c event data (like input4), but different variable names.
c  reads in the p travel times for all stations observing
c  the current event.  first trial hypocentral parameters are read, then
c  station names and weights and travel times.  arrivals at stations not
c  in the station list are discarded; normalized weights are calculated.
c  subroutine required: disto;
c
c  declaration statements:
c  local variables:
      real dep,x,y,pwt,tt(6),sta(6)
      character*4 chk,chksp,cblank,cblnk0,rmki(6)
      character*1 cns,cew,is(6)
      integer ip(6),nwav(2)
      character*5 sta5(6)
      character*6 sta6(6)
      character*2 netin(6)
      character*3 cmp(6)
      character*110 line
c
c  common block variables:
      include 'simul2014_common.inc'
c
      data  blank,iblank,blank0 /4h    ,4h    ,4h0   /
      data  cblank,cblnk0 /4h    ,4h0   /
c
      if((nc.gt.1).or.(nevt.gt.0)) goto 401
      nobt=0
      nobtp=0
      nobts=0
      nordt=0
      nordtp=0
      nordts=0
c  Start cluster read
  401 read(nfile,4111,end=499) line,chk
      if(line(1:13).eq.'BEGIN CLUSTER') goto 402
      if((chk.eq.cblank).or.(chk.eq.cblnk0)) goto 401
      write(16,1641)nc,line(1:10)
 1641 format(' *Input9 for Cluster:',i3,'SKIPPING Line:',a10)
      goto 401
  499 nclu=nc-1
      write(16,1642) nc,nclu
 1642 format('*** EOF:Input9 for Cluster:',i3,
     & ' CORRECTED nclu=',i3,' ***')
      if(nclu.lt.1) stop
      goto 500
c      return
c
  402 write(16,1643) nc
cc PRINT
      if(kout2.eq.2) write(6,1643) nc
 1643 format(/,' Cluster:',i3,' Reading events')
      n=0
      nobsc(nc)=0
c  print out heading for event cards
      write(16,4000)
cc PRINT
c      write(6,4000)
 4000 format(//,' trial event locations:',/,
     2 5x,'n    origin time     latitude longitude  depth',
     3 ' mag',5x,'x      y      z',4x,'nob np ns seco2')
  410 n=n+1
    1 nobs=0
      nwav(1)=0
      nwav(2)=0
      wsum=0
c  read in event card
c  check for extra travel time card or End Cluster
   10 read(nfile,4111,end=99) line,chk
 4111 format(a85,t1,a4)
      if(line(1:11).ne.'END CLUSTER') goto 5
      if(n.le.2) then
        write(16,1645) nc
 1645   format(' ** TOO FEW events in cluster:',i3,', Skipping to ',
     &  'next Cluster')
        goto 401
      endif
      ncev(nc)=n-1
      goto 500
c      return
    5 if((chk.eq.cblank).or.(chk.eq.cblnk0)) goto 10
c
   15 continue
      read(line,4001,err=5098) iyrmoce(n,nc),idayce(n,nc),ihrce(n,nc),
     2 minoce(n,nc),secoce(n,nc),ltdce(n,nc),cns,celtm(n,nc),
     3 lndce(n,nc),cew,celnm(n,nc),dep,rmagce(n,nc)
 4001 format(a4,a2,1x,a2,i2,1x,f5.2,i3,a1,f5.2,1x,i3,a1,f5.2,2f7.2,
     2  t69,f6.2)
      if(cns.eq.'S'.or.cns.eq.'s') then
        ltdce(n,nc)=-1.*ltdce(n,nc)
        celtm(n,nc)=-1.*celtm(n,nc)
      endif
      if(cew.eq.'E'.or.cew.eq.'e') then
        lndce(n,nc)=-1.*lndce(n,nc)
        celnm(n,nc)=-1.*celnm(n,nc)
      endif
c
      if(iyrmoce(n,nc).eq.iblank) goto 10
c  calculate event location in cartesian coordinate system
      call disto(ltdce(n,nc),celtm(n,nc),lndce(n,nc),celnm(n,nc),x,y)
c  store event coordinates
      cevc(1,n,nc)=x
      cevc(2,n,nc)=y
      cevc(3,n,nc)=dep
c
c  read in travel times to stations
c  format including phase-remarks
   20 continue
      read(nfile,4007) line
 4007 format(a110)
c
      if(iuseq.eq.0) then
      if(kttfor.eq.1) goto 21
      if(kttfor.eq.3) goto 121
      if(kttfor.eq.4) goto 122
      if(kttfor.ne.2) then
        write(16,1621) kttfor
 1621   format('kttfor=',i2,' SHOULD BE 1=cnsp, 2=phs, 3=phs sta6lt',
     &  /,' Continue with kttfor=1, assume cnsp format ***')
        kttfor=1
        goto 21
      endif
c phs format (SCB), one travel-time per line
c SCB change to read output from "pha2cnsp.pl"
      jend=1
      j=1
      read(line,4992,err=5099) sta(j),rmki(j),tt(j)
      read(line,4993,err=5099) is(j),ip(j)
 4992 format(a4,1x,a4,1x,f6.3)
 4993 format(6x,a1,1x,i1)
      goto 25
c
c phs with 6-letter station (for Japan)
  121 jend=1
      j=1
      read(line,4994,err=5099) sta6(j),rmki(j),tt(j)
      read(line,4995,err=5099) sta(j),is(j),ip(j)
 4994 format(a6,1x,a4,1x,f9.4)
 4995 format(a4,4x,a1,1x,i1)
      goto 25
c
c phs with 5-letter station with net (for NCEDC,IRIS)
c Allow cmp for possible future use
  122 jend=1
      j=1
      read(line,4990,err=5099) sta5(j),cmp(j),netin(j),rmki(j),tt(j)
      read(line,4991,err=5099) sta(j),is(j),ip(j)
 4990 format(a5,1x,a3,1x,a2,1x,a4,1x,f9.4)
 4991 format(a4,10x,a1,1x,i1)
      goto 25
c  
c  read in travel times to stations in groups of six
   21 jend=6
c CHT old format
      read(line,4996,err=921) (sta(j),rmki(j),tt(j),j=1,6),
     2  (is(j),ip(j),j=1,6)
 4996 format(6(a4,a2,f6.2),t1,6(4x,a1,i1,6x))
      goto 25
 921  continue
c cnsp format
        read(line,4006,err=22) (sta(j),rmki(j),tt(j),j=1,6),
     2  (is(j),ip(j),j=1,6)
c tstar data for iuseq =1
      else
        if(kttfor.eq.4) then
          j=1
          jend=1
          read(line,4990,err=5099) sta5(j),cmp(j),netin(j),rmki(j),tt(j)
C      PRINT *,'sta5(j),cmp(j),netin(j),rmki(j),tt(j)',
C     2 sta5(j),cmp(j),netin(j),rmki(j),tt(j)
          read(line,4991,err=5099) sta(j),is(j),ip(j)
C      PRINT *,'sta(j),is(j),ip(j) ',sta(j),' ',is(j),ip(j),
C     2 ' blank= ',blank
        else
          jend=6
          read(line,4029,err=22) (sta(j),rmki(j),tt(j),j=1,6),
     2    (is(j),ip(j),j=1,6)
        endif
      endif
 4006 format(6(2a4,f6.2),t1,6(5x,a1,1x,i1,6x))
 4029 format(6(2a4,f10.4),t1,6(5x,a1,1x,i1,10x))
      goto 25
c Alberto Michelini (LBL) format, could be used for tstar
   22 read(line,4008,err=5099) (sta(j),rmki(j),tt(j),j=1,5),
     2  (is(j),ip(j),j=1,5)
 4008 format(5(2a4,f7.4),t1,5(5x,a1,1x,i1,7x))
      jend=5
c
c  loop over the six (or jend) readings
c  terminate if station is a blank (end of observations for event)
c  skip if weight is incorrect in some way
c  search list for current station, and index it if in list
c  also calculate p-arrival time and unnormalized weight
c  and increment number of observations for event
c
   25 continue
c  blank line separates the events
      if ((sta(1).eq.blank).or.(sta(1).eq.blank0)) goto 50
      do 30 j=1,jend
         if(ip(j).ge.4) goto 30
         pwt=1.0/(2.**(ip(j)))
         if (pwt.le.0.0) goto 30
c only use if labeled P or S
         if((is(j).ne.'P').and.(is(j).ne.'p').and.
     2     (is(j).ne.'S').and.(is(j).ne.'s')) goto 30
c search station list
         do 33 k=1,nsts
          if(kttfor.lt.3) then
            if (sta(j).eq.stn(k)) go to 35
          else
            if(kttfor.eq.3) then
              if(sta6(j).eq.stn6(k)) goto 35
            else
              if((sta5(j).eq.stn5(k)).and.(netin(j).eq.net(k))) goto 35
            endif
          endif
   33    continue
c         if(kout2.gt.1) goto 30
         if(sta(j).eq.blank .or. sta(j).eq.blank0) then
c            write(16,1609)
 1609       format(' Blank-station observation ignored')
         else
           if(kttfor.lt.3) then
             write(16,1610)sta(j),iyrmoce(n,nc),idayce(n,nc),
     *       ihrce(n,nc),minoce(n,nc),secoce(n,nc)
 1610        format(' ** WARNING:  Observed station "',a4,
     *       '" (event ',a4,a2,1x,a2,i2,1x,f5.2,
     *       ') MISSING in station list--observation ignored **')
           else
             if(kttfor.eq.3) then
               write(16,1611)sta6(j),iyrmoce(n,nc),idayce(n,nc),
     *         ihrce(n,nc),minoce(n,nc),secoce(n,nc)
 1611          format(' ** WARNING:  Observed station "',a6,
     *         '" (event ',a4,a2,1x,a2,i2,1x,f5.2,
     *         ') MISSING in station list--observation ignored **')
            else
              write(16,1631)sta5(j),netin(j),iyrmoce(n,nc),idayce(n,nc),
     *         ihrce(n,nc),minoce(n,nc),secoce(n,nc)
 1631          format(' ** WARNING:  Observed station "',a5,1x,a2,
     *         '" (event ',a4,a2,1x,a2,i2,1x,f5.2,
     *         ') MISSING in station list--observation ignored **')
             endif
           endif
         endif
         goto 30
c  add to set of observing stations if on list
   35    continue
c Don't use 2nd clock stations for cluster events
      if(iuse2t.eq.0) goto 135
      if(iclock(k).eq.0) goto 135
      write(16,1612) stn(k)
 1612 format('*** SKIPPING 2nd clock station: ',a4,' for ',
     & 'Cluster event ***')
      goto 30
c  don't use arrivals further away than delt2
  135    dx=cevc(1,n,nc)-stc(1,k)
         dy=cevc(2,n,nc)-stc(2,k)
         dz=cevc(3,n,nc)-stc(3,k)
c
c  Print out a file of event and receiver coordinates if kout2=0
c  and nitmax=0
      if((kout2.eq.0).and.(nitmax.eq.0)) then
c  Only print out P
      if((is(j).eq.'P').or.(is(j).eq.'p')) then
        eltd=float(ltdce(n,nc))+celtm(n,nc)/60.0
        elnd=float(lndce(n,nc))+celnm(n,nc)/60.0
        elnd=-1.0*elnd
        sltd=float(ltds(k))+sltm(k)/60.0
        slnd=float(lnds(k))+slnm(k)/60.0
        slnd=-1.0*slnd
        write(38,3802) eltd,elnd,dep,secoce(n,nc),sltd,slnd,stc(3,k),
     2 pdl(k),tt(j)
 3802  format(5x,2f10.4,2f7.2,2x,2f10.4,2f7.2,2x,f8.2)
      endif
      endif
c
      rdelta=sqrt((dx*dx)+(dy*dy)+(dz*dz))
      if(iuseq.eq.1) then
c make quick check of high Q to throw out bad readings
c high Q is an input parameter
         Qobs = rdelta/(6.0*tt(j))
         if(Qobs.ge.qrmax) then
          if(kttfor.lt.3) then
           write(16,1622) stn(k),tt(j),qobs,qrmax
          else
            if(kttfor.eq.3) then
              write(16,1623) stn6(k),tt(j),qobs,qrmax
            else
              write(16,1633) stn5(k),net(k),tt(j),qobs,qrmax
            endif
          endif
 1622      format(' **** NOT USING ',a4,' t* obs=',f7.4,'Qobs=',e12.4,
     2      ' GREATER THAN QRMAX=',f7.0)
 1623      format(' **** NOT USING ',a6,' t* obs=',f7.4,'Qobs=',e12.4,
     2      ' GREATER THAN QRMAX=',f7.0)
 1633      format(' **** NOT USING ',a5,1x,a2,' t* obs=',f7.4,'Qobs=',
     2      e12.4,' GREATER THAN QRMAX=',f7.0)
          goto 30
         endif
       endif
c
c  Use epicentral distance instead of slant distance
         delta=sqrt((dx*dx)+(dy*dy))
         if(delta.lt.delt2) goto 36
         if(kout2.eq.0) write(16,1602) delta,delt2,stn(k),stn6(k),n
 1602    format(' **** delta,delt2=',f15.2,f8.2,2x,a4,a6,
     &     ' not used for event',i4)
         goto 30
c        input iuses=0 (changed to iuses+1 in input1) to disable
c           S-wave processing.
   36    if(iuses.eq.1.and.is(j).eq.'S') goto 30
         if((iusep.eq.0).and.((is(j).eq.'P').or.(is(j).eq.'p')))
     *      goto 30
         nobs=nobs+1
c check for too many observations 
         if(nobs.gt.maxobs) then
            write(6,1603) iyrmoce(n,nc),idayce(n,nc),ihrce(n,nc),
     2      minoce(n,nc),maxobs
            write(16,1603) iyrmoce(n,nc),idayce(n,nc),ihrce(n,nc),
     2      minoce(n,nc),maxobs
 1603       format('  INPUT ERROR, Too many readings for event:',
     2      a4,a2,1x,a2,i2,/,'    CONTINUE with first',i5,
     3      'observations')
            nobs=maxobs
            goto 50
         endif
         dltac(nobs,n,nc)=delta
         rdltac(nobs,n,nc)=rdelta
         iwc(nobs,n,nc)=ip(j)
         istoc(nobs,n,nc)=k
         secpc(nobs,n,nc)=tt(j)+secoce(n,nc)
         rmkc(nobs,n,nc)=rmki(j)
c  here are additions by wap to read in s times. If the time is
c  an s time, is(j) will equal 'S'. If this is true, the
c  proper adjustment will be made for the station delay, and
c  a 1 will be put in intspc(nobs,n,nc).
         intspc(nobs,n,nc)=0
         if(rmkc(nobs,n,nc).eq.cblank) rmkc(nobs,n,nc)=' P  '
         if((is(j).ne.'S').and.(is(j).ne.'s')) goto 37
c  intspc(nobs,n,nc) tells whether the nobs'th observation of the
c  arrival time is a P (0), or an S (1).
         intspc(nobs,n,nc)=1
c** S-P CHANGE
c   For S-P we do not want to add origin time
         secpc(nobs,n,nc)=tt(j)
         if(rmkc(nobs,n,nc).eq.cblank) rmkc(nobs,n,nc)=' S  '
  37     continue
         kv=intspc(nobs,n,nc)+1
         nrd(k,kv)=nrd(k,kv)+1
         nwav(kv)=nwav(kv)+1
         wtc(nobs,n,nc)=pwt
         wsum=wsum+pwt
   30 continue
c  continue reading stations and travel times
      goto 20
c
   50 continue
c  end of travel time readings for current event
c  check that there are at least four readings for this event
      if(nobs.lt.4) goto 90
cfhek 
      if(ltdce(n,nc).ge.0.and.celtm(n,nc).ge.0.) then
        cns='N'
        latdeg=iabs(ltdce(n,nc))
        alatmin=abs(celtm(n,nc))
      else
        cns='S'
        latdeg=iabs(ltdce(n,nc))
        alatmin=abs(celtm(n,nc))
      endif
      if(lndce(n,nc).ge.0.and.celnm(n,nc).ge.0.) then
        cew='W'
        londeg=iabs(lndce(n,nc))
        alonmin=abs(celnm(n,nc))
      else
        cew='E'
        londeg=iabs(lndce(n,nc))
        alonmin=abs(celnm(n,nc))
      endif
        write(16,4004) n,iyrmoce(n,nc),idayce(n,nc),ihrce(n,nc),
     2  minoce(n,nc),secoce(n,nc),
     3  latdeg,cns,alatmin,londeg,cew,alonmin,dep,rmagce(n,nc),x,y,dep,
     4  nobs,nwav(1),nwav(2)
cc PRINT if kout=2
        if(kout2.eq.2) write(6,4004) n,iyrmoce(n,nc),idayce(n,nc),
     2  ihrce(n,nc),minoce(n,nc),secoce(n,nc),
     3  latdeg,cns,alatmin,londeg,cew,alonmin,dep,rmagce(n,nc),x,y,dep,
     4  nobs,nwav(1),nwav(2)
 4004 format(2h *,i5,1x,a4,a2,1x,a2,i2,f6.2,i3,a1,f5.2,i4,a1,f5.2,
     2 f7.2,f5.2,1x,3f7.2,1x,3i3,f6.2)
c
c  check for too many readings
      if(nobs.le.maxobs) goto 55
        write(16,1698) n,nobs,maxobs,maxobs
 1698   format(' **** ERROR **** in event ',i4,', TOO MANY ',
     2  'OBSERVATIONS, nobs=',i6,' array size allows',i6,/,
     2  '   continuing with',i6,' observations',/)
        nobs=maxobs
c  normalize reading weights
   55 wfac=nobs/wsum
      kobsce(n,nc)=nobs
      nobsc(nc)=nobsc(nc)+nobs
      nobt=nobt+nobs
      nobtp=nobtp+nwav(1)
      nobts=nobts+nwav(2)
      nobtce=nobtce+nobs
      nobteq=nobteq+nobs
      do 60 j=1,nobs
      wtc(j,n,nc)=wtc(j,n,nc)*wfac
   60 continue
c
c  FIX TO CONVERT S TO S-P DATA, CHT
c  For input S-P data characters 2,3 of phase remark are 'Sp'
c     ie: "polarity" field is 'p'
c
      do 987 j=1,nobs
      if (intspc(j,n,nc).eq.1) then
      write(chksp,1620) rmkc(j,n,nc)
 1620 format(a4)
c check print
c      write(77,7708) stn(istoc(j,n,nc)),rmkc(j,n,nc),chksp, chksp(2:3),
c     2  secpc(j,n,nc)
 7708 format('sta(j),rmkc:',a4,1x,a4,' chksp:',a4,', chksp(2:3):',a2,
     2  'secpc=',f8.2)
      if(chksp(2:3).eq.'Sp') goto 987
c      write(77,7709) stn(istoc(j,n,nc)), n
 7709 format('For sta ',a4,' cluevent',i3,' rmk shows s intead of s-p?')
      secpc(j,n,nc)=secpc(j,n,nc)+secoce(n,nc)
c
c  find matching P observation
      imatch=0
      do 986 jj=1,nobs
      if (stn(istoc(jj,n,nc)).eq.stn(istoc(j,n,nc)).and.
     & intspc(jj,n,nc).eq.0) then
        imatch=imatch+1
        if(imatch.eq.1) then
          secpc(j,n,nc)=secpc(j,n,nc)-secpc(jj,n,nc)
          rmkc(j,n,nc)(2:3)='Sp'
        endif
        if(imatch.gt.1) write(16,9876) stn(istoc(jj,n,nc)),
     &  stn(istoc(j,n,nc)),secpc(jj,n,nc),secpc(j,n,nc),
     &  secpc(j,n,nc)-secpc(jj,n,nc)
 9876 format('  WARNING - second S-P match! ',a4,1x,a4,1x,3f8.3,' SKIP')
        if(imatch.gt.1) wtc(j,n,nc)=0.
      endif
 986  continue
      endif
  987 continue
c
c check print
c      write(77,4004) n,iyrmoce(n,nc),idayce(n,nc),ihrce(n,nc),
c     2  minoce(n,nc),secoce(n,nc),
c     3  latdeg,cns,alatmin,londeg,cew,alonmin,dep,rmagce(n,nc),x,y,dep,
c     4  nobs,nwav(1),nwav(2)
c      write(77,7601) ((stn(istoc(j,n,nc)),intspc(j,n,nc),iwc(j,n,nc),
c     2  secpc(j,n,nc),rdltac(j,n,nc)),j=1,nobs)
 7601 format(4(1x,a4,2i2,f7.2,f8.2))
c
c
c Get next event in this Cluster
      goto 410
c
c  not enough readings for event n
   90 continue
      write(16,4009) iyrmoce(n,nc),idayce(n,nc),ihrce(n,nc),minoce(n,nc)
 4009 format(' *** event ',a4,a2,1x,a2,i2,' should be discarded - ',
     2 'too few observations ***',/,' *** SKIPPING TO NEXT EVENT ***')
      goto 1
c
c  ran out of event cards before 'END CLUSTER' - error]
   99 write(16,4099) nfile,nc
 4099 format('**** Input9:end of data file',i2,', reading cluster:',i3,
     & ' ****')
      if((nc.lt.nclu).or.(n.le.20)) then
        write(16,1644) nlcu
 1644   format(' **** nclu=',i3,' STOP ****')
        stop
      endif
      ncev(nc)=n-1
      nclu=nc
      write(16,1642) nc,nclu
c
c compute centroid of cluster
  500 continue
      ne=ncev(nc)
      do 515 i=1,3
      csum=0.0
        do 510 j=1,ne
        csum=csum+cevc(i,j,nc)
  510   continue
      ccenc(i,nc)=csum/float(ne)
  515 continue
      return
c
 5098 write(16,5050) n
      write(16,4012) nchar,line,chk
 4012 format(2x,i5,' characters in line',/,2x,a85,/,2x,'first 4 ',
     2 'char = ',a4)
cfhek
      if(ltdce(n,nc).ge.0.and.celtm(n,nc).ge.0.) then
        cns='N'
        latdeg=iabs(ltdce(n,nc))
        alatmin=abs(celtm(n,nc))
      else
        cns='S'
        latdeg=iabs(ltdce(n,nc))
        alatmin=abs(celtm(n,nc))
      endif
      if(lndce(n,nc).ge.0.and.celnm(n,nc).ge.0.) then
        cew='W'
        londeg=iabs(lndce(n,nc))
        alonmin=abs(celnm(n,nc))
      else
        cew='E'
        londeg=iabs(lndce(n,nc))
        alonmin=abs(celnm(n,nc))
      endif
cek
      write(16,4011) iyrmoce(n,nc),idayce(n,nc),ihrce(n,nc),
     & minoce(n,nc),secoce(n,nc),
     2  latdeg,cns,alatmin,londeg,cew,alonmin,dep,rmagce(n,nc)
 4011 format(1x,a4,a2,a3,i2,1x,f5.2,1x,i3,a1,f5.2,1x,i3,a1,f5.2,2f7.2,
     2  f6.2)
      goto 10
c
 5099   write(16,5051) n,kttfor,line
      write(16,4015) (sta(j),is(j),ip(j),tt(j),j=1,6)
 4015 format(1x,6(a4,a1,i1,f6.2))
 5050   format('0error in hyp data input9 event= ',i8,/,'** ',a110)
 5051   format('0error in phase data input9 event= ',i8,' kttfor=',i2,
     2    /,'line:',a110)
        stop
c****** end of subroutine input9 ******
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine input10(n,nfile)
c  this routine reads in the p travel times for all stations observing
c  reads in the teleseismic travel-time residuals for all stations
c  observing the current tele event (file 10).
c  first teleseismic hypocenter and tele path delay term(s) are read.
c  then station names, weights, travel-time residuals, and piercing
c  points (at base of inversion grid).  The piercing points are input
c  and not calculated, DMEP used Kohler's code to compute the input
c  file.  The tele path delay is also updated in the velocity inversion.
c  plan for both S and P tele path delays (telpad,telpads)
c  No distance wt'g for tele.
c
c  subroutine required: disto;
c
c  declaration statements:
c  local variables:
      real dep,x,y,pwt,tt(6),sta(6),xc(3)
      character*4 rmki(6)
c      dimension nwav(2),rmki(6),ip(6),tt(6),sta(6),is(6)
c      character*4 sta(6),chk,chksp,cblank,cblnk0
      character*4 chk,chksp,cblank,cblnk0
      character*1 cns,cew,is(6)
      character*6 sta6(6)
      character*5 sta5(6)
      character*2 netin(6)
      integer ip(6),nwav(2)
      character*110 line
c
c  common block variables:
      include 'simul2014_common.inc'
c
      data  blank,iblank,blank0 /4h    ,4h    ,4h0   /
      data  cblank,cblnk0 /4h    ,4h0   /
c
      nt=n-nebs
      if (n.gt.1) go to 1
c  print out heading for event cards
      write(16,4000)
 4000 format(//,' trial event locations:',/,
     2 5x,'n    origin time     latitude longitude  depth',
     3 ' mag',5x,'x      y      z',4x,'nob np ns seco2')
      nobt=0
      nobtp=0
      nobts=0
      nordt=0
      nordtp=0
      nordts=0
    1 nobs=0
      nwav(1)=0
      nwav(2)=0
      wsum=0
c  read in event card
c  check for extra travel time card
   10 read(nfile,4111,end=99) line,chk
 4111 format(a85,t1,a4)
c
    5 if((chk.eq.cblank).or.(chk.eq.cblnk0)) goto 10
c
   15 continue
      read(line,4003,err=16) iyrmo(n),iday(n),ihr(n),
     2   mino(n),seco(n),elat(n),elon(n),
     3   dep,rmag(n),telpad(nt),telpads(nt)
        call ltlndm(elat(n),elon(n),ltde(n),eltm(n),ecns(n),lnde(n),
     2   elnm(n),ecew(n),1)
      goto 17
   16 read(line,4002,err=5098) iyrmo(n),iday(n),ihr(n),
     2   mino(n),seco(n),ltde(n),ecns(n),eltm(n),lnde(n),cew,elnm(n),
     3   dep,rmag(n),telpad(nt),telpads(nt)
 4002 format(a4,a2,1x,a2,i2,1x,f5.2,i3,a1,f5.2,1x,i3,a1,f5.2,2f7.2,
     2  t57,2f9.2)
 4003 format(a4,a2,1x,a2,i2,1x,f5.2,f9.4,f10.4,2f7.2,
     2  t57,2f9.2)
        call ltlndm(elat(n),elon(n),ltde(n),eltm(n),ecns(n),lnde(n),
     2   elnm(n),ecew(n),2)
c
   17 if(iyrmo(n).eq.iblank) go to 10
c do not put tele hyp into local coords, save depth there
      evc(3,n)=dep
c
c  read in travel times to stations
c  format including phase-remarks
   20 continue
      read(nfile,4007) line
 4007 format(a110)
c  blank line separates the events
      read(line,4997) sta(1)
      if ((sta(1).eq.blank).or.(sta(1).eq.blank0)) go to 50
c
      nobs=nobs+1
c check for too many observations 
         if(nobs.gt.maxobs) then
            write(6,1603) iyrmo(n),iday(n),ihr(n),mino(n),maxobs
 1603       format('  INPUT ERROR, Too many readings for event:',
     2      a4,a2,1x,a2,i2,/,'    CONTINUE with first',i5,
     3      'observations')
            write(16,1603) iyrmo(n),iday(n),ihr(n),mino(n),maxobs
            nobs=maxobs
            goto 50
         endif
c use m for nobs in arrays for easier editing
      m=nobs
c phs format for tele, one travel-time residual per line, with piercing point
      jend=1
      j=1
c  decimal degrees
      if(kttfor.eq.4) then
        read(line,4993,err=5099) is(j),ip(j)
        read(line,4992,err=21) sta5(j),netin(j),rmki(j),secte(m,nt),
     2  teplat(m,nt),teplon(m,nt),tepp(3,m,nt)
 4992   format(a5,1x,a2,1x,a4,f10.4,2x,2f10.4,f8.2)
 4993   format(10x,a1,1x,i1)
        goto 22
c if error on reading decimal degrees, try deg minute format
   21   read(line,4991,err=5099)sta5(j),netin(j),rmki(j),secte(m,nt),
     2  ltdtep(m,nt),tecns(m,nt),tepltm(m,nt),lndtep(m,nt),
     3  tecew(m,nt),teplnm(m,nt),tepp(3,m,nt)
 4991   format(a5,1x,a2,1x,a4,f10.4,2x,2(i4,a1,f5.2),f8.2)
        goto 14
      endif
      read(line,4997,err=23) sta(j),rmki(j),secte(m,nt),teplat(m,nt),
     2  teplon(m,nt),tepp(3,m,nt)
c      WRITE(86,699) sta(j),rmki(j),
c     2  secte(m,nt),teplat(m,nt),teplon(m,nt),tepp(3,m,nt)
  699 FORMAT('sta rmk secte tep lt ln z ',a4,3x,a4,
     2 f8.2,2f10.4,f8.2)
      read(line,4998,err=5099) is(j),ip(j)
   22 call ltlndm(teplat(m,nt),teplon(m,nt),ltdtep(m,nt),tepltm(m,nt),
     2  tecns(m,nt),lndtep(m,nt),teplnm(m,nt),tecew(m,nt),1)
 4997 format(a4,3x,a4,f10.4,2x,2f10.4,f8.2)
 4998 format(8x,a1,1x,i1)
      goto 24
c if error on reading decimal degrees, try deg minute format
   23 read(line,4999,err=5099) sta(j),rmki(j),secte(m,nt),ltdtep(m,nt),
     2 tecns(m,nt),tepltm(m,nt),lndtep(m,nt),tecew(m,nt),teplnm(m,nt),
     3 tepp(3,m,nt)
      PRINT *,'21:sta rmk secte tep lt ln z',sta(j),rmki(j),
     2  secte(m,nt),teplat(m,nt),teplon(m,nt),tepp(3,m,nt)
c      WRITE(86,*) '21:sta rmk secte tep lt ln z',sta(j),rmki(j),
c     2  secte(m,nt),teplat(m,nt),teplon(m,nt),tepp(3,m,nt)
      read(line,4998,err=5099) is(j),ip(j)
 4999 format(a4,3x,a4,f10.4,2x,2(i4,a1,f5.2),f8.2)
   14 call ltlndm(teplat(m,nt),teplon(m,nt),ltdtep(m,nt),tepltm(m,nt),
     2  tecns(m,nt),lndtep(m,nt),teplnm(m,nt),tecew(m,nt),2)
c phs with 6-letter station (for Japan)
   24 if(kttfor.eq.3) read(line,4994) sta6(j)
 4994 format(a6,1x,a4,1x,f9.4)
      call disto(ltdtep(m,nt),tepltm(m,nt),lndtep(m,nt),teplnm(m,nt),
     2  x,y)
      tepp(1,m,nt)=x
      tepp(2,m,nt)=y
c
c  skip if weight is incorrect in some way
c  search list for current station, and index it if in list
c  also calculate p-arrival time and unnormalized weight
c
   25 continue
c only 1 obs per line for tele so this is all in the one read loop
c      do 30 j=1,jend
         if(ip(j).ge.4) goto 30
         pwt=1.0/(2.**(ip(j)))
         if (pwt.le.0.0) go to 30
c only use if labeled P or S
         if((is(j).ne.'P').and.(is(j).ne.'p').and.
     2     (is(j).ne.'S').and.(is(j).ne.'s')) goto 30
c search station list
         do 33 k=1,nsts
          if(kttfor.ne.3) then
            if (sta(j).eq.stn(k)) go to 35
          else
            if(sta6(j).eq.stn6(k)) goto 35
          endif
   33    continue
c         if(kout2.gt.1) goto 30
         if(sta(j).eq.blank .or. sta(j).eq.blank0) then
c            write(16,1609)
 1609       format(' Blank-station observation ignored')
         else
           if(kttfor.lt.3) then
             write(16,1610) sta(j),iyrmo(n),iday(n),ihr(n),
     *       mino(n),seco(n)
 1610        format(' ** WARNING:  Observed station "',a4,
     *       '" (event ',a4,a2,1x,a2,i2,1x,f5.2,
     *      ') MISSING in station list--observation ignored **')
           else
             if(kttfor.eq.3) then
               write(16,1611)sta6(j),iyrmo(n),iday(n),ihr(n),
     *         mino(n),seco(n)
 1611          format(' ** WARNING:  Observed station "',a6,
     *         '" (event ',a4,a2,1x,a2,i2,1x,f5.2,
     *         ') MISSING in station list--observation ignored **')
             else
               write(16,1631)sta5(j),netin(j),iyrmo(n),iday(n),
     *          ihr(n),mino(n),seco(n)
 1631          format(' ** WARNING:  Observed station "',a5,1x,a2,
     *         '" (event ',a4,a2,1x,a2,i2,1x,f5.2,
     *         ') MISSING in station list--observation ignored **')
             endif
           endif
         endif
         nobs=nobs-1
         go to 30
c  add to set of observing stations if on list
   35    continue
         dx=tepp(1,m,nt)-stc(1,k)
         dy=tepp(2,m,nt)-stc(2,k)
         dz=tepp(3,m,nt)-stc(3,k)
c
c  Print out a file of piercing pt and receiver coordinates if kout2=0
c  and nitmax=0
      if((kout2.eq.0).and.(nitmax.eq.0)) then
c  Only print out P
        if((is(j).eq.'P').or.(is(j).eq.'p')) then
          write(38,3802) teplat(m,nt),teplon(m,nt),tepp(3,m,nt),telpad,
     2    slat(k),slon(k),stc(3,k),pdl(k),secte(m,nt)
 3802     format(5x,2f10.4,2f7.2,2x,2f10.4,2f7.2,2x,f8.2)
        endif
      endif
c
      rdelta=sqrt((dx*dx)+(dy*dy)+(dz*dz))
      if(iuseq.eq.1) then
c make quick check of high Q to throw out bad readings
c high Q is an input parameter
c for tele Q, only check very far off, since also telpad
c use 2*qrmax
         Qobs = rdelta/(6.0*tt(j))
         qrmax2=2.0*qrmax
         if(Qobs.ge.qrmax2) then
          if(kttfor.lt.3) then
            write(16,1622) stn(k),tt(j),qobs,qrmax2
          else
            if(kttfor.eq.3) then
              write(16,1623) stn6(k),tt(j),qobs,qrmax2
            else
              write(16,1633) stn5(k),net(k),tt(j),qobs,qrmax2
            endif
          endif
 1622   format(' **** TELE NOT USING ',a4,' t* obs=',f7.4,'Qobs=',e12.4,
     2      ' GREATER THAN QRMAX2=',f7.0,'  TELE')
 1623   format(' **** TELE NOT USING ',a6,' t* obs=',f7.4,'Qobs=',e12.4,
     2      ' GREATER THAN QRMAX2=',f7.0,'  TELE')
 1633      format(' **** NOT USING ',a5,1x,a2,' t* obs=',f7.4,'Qobs=',
     2      e12.4,' GREATER THAN QRMAX2=',f7.0)
           nobs=nobs-1
           goto 30
         endif
       endif
c
c  Use epicentral distance instead of slant distance
         delta=sqrt((dx*dx)+(dy*dy))
c        input iuses=0 (changed to iuses+1 in input1) to disable
c           S-wave processing.
   36    if(iuses.eq.1.and.is(j).eq.'S') goto 30
         if((iusep.eq.0).and.((is(j).eq.'P').or.(is(j).eq.'p')))
     *      goto 30
         dlta(nobs,n)=delta
         rdlta(nobs,n)=rdelta
         iw(nobs,n)=ip(j)
         isto(nobs,n)=k
         rmk(nobs,n)=rmki(j)
c  here are additions by wap to read in s times. If the time is
c  an s time, is(j) will equal 'S'. If this is true, the
c  proper adjustment will be made for the station delay, and
c  a 1 will be put in intsp(nobs,n).
         intsp(nobs,n)=0
         if(rmk(nobs,n).eq.cblank) rmk(nobs,n)=' P  '
         if((is(j).ne.'S').and.(is(j).ne.'s')) goto 37
c  intsp(nobs,n) tells whether the nobs'th observation of the
c  arrival time is a P (0), or an S (1).
         intsp(nobs,n)=1
c** S-P CHANGE
         if(rmk(nobs,n).eq.cblank) rmk(nobs,n)=' S  '
  37     continue
         kv=intsp(nobs,n)+1
         nrd(k,kv)=nrd(k,kv)+1
         nwav(kv)=nwav(kv)+1
         wt(nobs,n)=pwt
         wsum=wsum+pwt
   30 continue
c  continue reading stations and travel times
      go to 20
c
   50 continue
c  end of travel time readings for current event
c  check for minimum number of tele obs
      if(nobs.lt.ntobmin) goto 90
c get centroid
      do 53 i=1,3
        xc(i)=0.0
        do 52 j=1,nobs
   52   xc(i)=xc(i)+tepp(i,j,nt)
        tecen(i,nt)=xc(i)/float(nobs)
   53 continue
      call latlon(tecen(1,nt),tecen(2,nt),lat,xlat,lon,xlon)
      call ltlndm(tecnlt(nt),tecnln(nt),lat,xlat,tecns(1,nt),
     2  lon,xlon,tecew(1,nt),2)
cfhek 
      if(ltde(n).ge.0.and.eltm(n).ge.0.) then
        cns='N'
        latdeg=iabs(ltde(n))
        alatmin=abs(eltm(n))
      else
        cns='S'
        latdeg=iabs(ltde(n))
        alatmin=abs(eltm(n))
      endif
      if(lnde(n).ge.0.and.elnm(n).ge.0.) then
        cew='W'
        londeg=iabs(lnde(n))
        alonmin=abs(elnm(n))
      else
        cew='E'
        londeg=iabs(lnde(n))
        alonmin=abs(elnm(n))
      endif
        write(16,4004) n,iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
     2  latdeg,cns,alatmin,londeg,cew,alonmin,dep,rmag(n),
     3  (tecen(i,nt),i=1,3),nobs,nwav(1),nwav(2),telpad(nt)
 4004 format(2h *,i5,1x,a4,a2,1x,a2,i2,f6.2,i3,a1,f5.2,i4,a1,f5.2,
     2 f7.2,f5.2,1x,3f7.2,1x,3i3,2f9.2)
c
c  normalize reading weights
c   change in normalizing factor 19-jul-1983
      wfac=nobs/wsum
c     wfac=sqrt(nobs/wsum)
c  check for too many readings
      if(nobs.le.maxobs) goto 55
        write(16,1698) n,nobs,maxobs,maxobs
 1698   format(' **** ERROR **** in event ',i4,', TOO MANY ',
     2  'OBSERVATIONS, nobs=',i6,' array size allows',i6,/,
     2  '   continuing with',i6,' observations',/)
        nobs=maxobs
   55 kobs(n)=nobs
      kobps(1,n)=nwav(1)
      kobps(2,n)=nwav(2)
      nobt=nobt+nobs
      nobtp=nobtp+nwav(1)
      nobts=nobts+nwav(2)
c this is tele input
      nobtte=nobtte+nobs
      nobt=nobt+nint((wtsht-1.0)*(float(nobs)))
      do 60 j=1,nobs
      wt(j,n)=wt(j,n)*wfac
   60 continue
c
c  FIX TO CONVERT S TO S-P DATA, CHT
c  For input S-P data characters 2,3 of phase remark are 'Sp'
c     ie: "polarity" field is 'p'
c
      do 987 j=1,nobs
      if (intsp(j,n).eq.1) then
      write(chksp,1620) rmk(j,n)
 1620 format(a4)
      if(chksp(2:3).eq.'Sp') goto 987
c
c  find matching P observation
      imatch=0
      do 986 jj=1,nobs
      if ((isto(jj,n).eq.isto(j,n)).and.
     & (intsp(jj,n).eq.0)) then
        imatch=imatch+1
        if(imatch.eq.1) then
          secte(j,nt)=secte(j,nt)-secte(jj,nt)
          rmk(j,n)(2:3)='Sp'
        endif
        if(imatch.gt.1) write(16,9876) stn(isto(jj,n)),stn(isto(j,n)),
     &  secte(jj,nt),secte(j,nt),secte(j,nt)-secte(jj,nt)
 9876 format('  WARNING - second S-P match! ',a4,1x,a4,1x,3f8.3,' SKIP')
        if(imatch.gt.1) wt(j,n)=0.
      endif
 986  continue
      endif
 987  continue
c
      return
c
c  not enough readings for event n
   90 continue
      write(16,4009)  iyrmo(n),iday(n),ihr(n),mino(n)
 4009 format(' *** event ',a4,a2,1x,a2,i2,' should be discarded - ',
     2 'too few observations ***',/,' *** SKIPPING TO NEXT EVENT ***')
      if(nfile.eq.10) then       ! tele datafile
        ntel=ntel-1
      endif
      nevt=nevt-1
      goto 1
c  ran out of event cards - error]
   99 write(16,4099) nfile
 4099 format('**** end of data file',i2,' - ntel too large ****')
      if(nfile.eq.10) then
          ntelin=n-1-neqs-nbls-nsht
          if((abs(ntel-ntelin)).gt.(0.10*ntel)) goto 199
          ntel=ntelin
          write(16,1697) ntel
 1697     format(' continue with ntel= ',i4)
          return
      endif
  199 continue
      write(16,1690)
 1690 format(' stop since neqs more than 10% off, or problem '
     2 ,'was nbls or nsht or ntel')
      stop
 5098 write(16,5050) n
      write(16,4012) nchar,line,chk
 4012 format(2x,i5,' characters in line',/,2x,a85,/,2x,'first 4 ',
     2 'char = ',a4)
cfhek
      if(ltde(n).ge.0.and.eltm(n).ge.0.) then
        cns='N'
        latdeg=iabs(ltde(n))
        alatmin=abs(eltm(n))
      else
        cns='S'
        latdeg=iabs(ltde(n))
        alatmin=abs(eltm(n))
      endif
      if(lnde(n).ge.0.and.elnm(n).ge.0.) then
        cew='W'
        londeg=iabs(lnde(n))
        alonmin=abs(elnm(n))
      else
        cew='E'
        londeg=iabs(lnde(n))
        alonmin=abs(elnm(n))
      endif
cek
      if(iuse2t.eq.0) then
        write(16,4011) iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
     2  latdeg,cns,alatmin,londeg,cew,alonmin,dep,rmag(n)
      else
        write(16,4011) iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
     2  latdeg,cns,alatmin,londeg,cew,alonmin,dep,rmag(n),seco2(n)
      endif
cek      write(16,4011) iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
cek     2 ltde(n),eltm(n),lnde(n),elnm(n),dep,rmag(n)
 4011 format(1x,a4,a2,a3,i2,1x,f5.2,1x,i3,a1,f5.2,1x,i3,a1,f5.2,2f7.2,
     2  f6.2)
c
      goto 10
c     stop
 5099   write(16,5050) n,line
      write(16,4015) (sta(j),is(j),ip(j),tt(j),j=1,6)
 4015 format(1x,6(a4,a1,i1,f6.2))
 5050   format('0error in data in input10. event= ',i8,/,'** ',a110)
        stop
c****** end of subroutine input10 ******
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine intmap(x,y,z,ip,jp,kp)
c  Modified by W. Prothero so a single call can get the indices
c  common block variables:
      include 'simul2014_common.inc'
c
      ip=int((x+xl)/bld)
      ip=ixloc(ip)
      jp=int((yl+y)/bld)
      jp=iyloc(jp)
      kp=int((z+zl)/bld)
      kp=izloc(kp)
c  If an array element=0, the position is off the map.
      return
c***** end of subroutine intmap *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
       subroutine latlon(x,y,lat,xlat,lon,xlon)
c  Subroutine to convert from Cartesian coordinates back to
c  latitude and longitude.
c
      real x,y,xlat,xlon
      integer lat,lon
      common/shortd/ xltkm,xlnkm,rota,nzco,xlt,xln,snr,csr,
     2 onorth,oeast,rotadg,cmerid,rlatdo,rlondo
c
      rad=1.7453292e-2
      rlt=9.9330647e-1
c
      fy=csr*y-snr*x
      fx=snr*y+csr*x
c
      if(nzco.lt.5) goto 100
c
      fy=fy/xltkm
      plt=xlt+fy
c
      xlt1=atan(rlt*tan(rad*(plt+xlt)/120.))
      fx=fx/(xlnkm*cos(xlt1))
      pln=xln+fx
c
      lat=plt/60.
      xlat=plt-lat*60.
      lon=pln/60.
      xlon=pln-lon*60.
c
      return
c
c convert with transverse mercator, NZTM2000, or input cmerid,
c  or TM State Plane Coordinates for Alaska
c  or nzco=4 use old New Zealand map grid coordinates
  100 continue
c get easting and northing
c since SIMUL has W=positive, must subtract point from origin for easting.
      xeast=oeast-fx
      ynorth=fy+onorth
      if(nzco.eq.4) call nzmgll(ynorth,xeast,lat,xlat,lon,xlon)
      if(nzco.eq.2) call spc2ll(ynorth,xeast,lat,xlat,lon,xlon)
      if((nzco.le.1).or.(nzco.eq.3)) call tmllxy(lat,xlat,lon,
     2    xlon,ynorth,xeast,2)
c NZ south latitude so negative in simul prog
      if((nzco.eq.1).or.(rlatdo.lt.0.0)) then
        lat= -1 * lat
        xlat = -1.0 * xlat
      endif
c Note that simul assumes W is positive (for USA) so long also neg
      if(nzco.eq.1) then
        lon= -1 * lon
        xlon = -1.0 * xlon
      endif
c      WRITE(6,605) x,y,lat,xlat,lon,xlon
  605 format('latlon x,y',2f8.2,' lat lon',2(i5,f8.4))
c
      return
c
c***** end of subroutine latlon ******
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine ltlndm(plt,pln,ltdeg,pltmin,chns,lndeg,plnmin,chew,
     2   kopt)
c  Simul coords have N and W positive for California, which is fine with
c  input deg minutes and cns,cew
c  Make the decimal degrees in arrays have E positive
c       kopt=1, input is decimal lat,lon
c       kopt=2, input is deg,minute
      character*1 chns,chew
      if(kopt.eq.1) then
        ltdeg=ifix(plt)
        pltmin=60.0*(plt-float(ltdeg))
        chns='N'
        if(plt.lt.0.0) chns='S'
        lndeg=ifix(pln)
        plnmin=60.0*(pln-float(lndeg))
        if(pln.ge.0.0) then
          chew='E'
          lndeg=-1*abs(lndeg)
          plnmin=-1.*abs(plnmin)
        else
          chew='W'
          lndeg=abs(lndeg)
          plnmin=abs(plnmin)
        endif
      else if(kopt.eq.2) then
        plt=float(ltdeg)+pltmin/60.0
        if(chns.eq.'S'.or.chns.eq.'s') then
          plt=-1.0*(abs(plt))
          ltdeg=-1*ltdeg
          pltmin=-1.*pltmin
        endif
        pln=float(lndeg)+plnmin/60.0
        if(chew.eq.'W'.or.chew.eq.'w') pln=-1.0*(abs(pln))
        if(chew.eq.'E'.or.chew.eq.'e') then
          lndeg=-1*lndeg
          plnmin=-1.*plnmin
        endif
      else
        PRINT *,'in ltlndm, kopt must be 1 or 2'
        PRINT *, 'plt,pln,ltdeg,pltmin,chns,lndeg,plnmin,chew'
        PRINT *,plt,pln,ltdeg,pltmin,chns,lndeg,plnmin,chew
        WRITE(16,*) 'in ltlndm, kopt must be 1 or 2'
        WRITE(16,*) 'plt,pln,ltdeg,pltmin,chns,lndeg,plnmin,chew'
        WRITE(16,*) plt,pln,ltdeg,pltmin,chns,lndeg,plnmin,chew
        stop
      endif
      return
c****** end of subroutine ltlndm ******
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine locclu(nc,niti,nwrce)
c  this routine locates the events in a cluster using both the 
c  travel-times and earthquake differential times.
c  this routine locates the given event in the initial velocity
c  structure using approximate ray tracing
c
c  common block variables:
      common/shortd/ xltkm,xlnkm,rota,nzco,xlt,xln,snr,csr,
     2 onorth,oeast,rotadg,cmerid,rlatdo,rlondo
      include 'simul2014_common.inc'
c
c  declaration statements:
      real dr5(3),dr5i(3),s(mxhparc),v(mxhparc,mxhparc),
     *  adj(mxhparc),ster(mxhparc),cov(mxhparc)
      real adja(5)
      real wc(mxobsc),xmovc(maxcev)
      double precision ahc(mxobsc,mxobsc)
cc  variables for error ellipse calculations
c      double precision ac(mxhparc,mxhparc)
c      integer iaz(3),idip(3)
c      real covar(mxhparc,mxhparc),v3(mxhparc,mxhparc),serr(mxhparc)
c      parameter (drad=1.7453292d-02)
c  variables for gap calculation
      real azp(maxobs),dltp(maxobs),gapp
      integer naz,nwav(2),nwavce(2,maxcev),numbad(maxcev)
c  for cluster will be earthquake
       nparhy=4
c
c  Define small change in hypocenter, depending on bld factor
      small=0.5
      if(bld.eq.0.1) small=0.1
c
   11 nec=ncev(nc)
      jahy=nparhy*nec
      do 200 ne=1,nec
      xmovc(ne)=900.0
  200 continue
c  nitl=counter for location iterations
      nitl=0
c  niti=counter for inversion iterations, niti=0 for initial locations
      if(niti.eq.0) then
        do 205 ne=1,nec
        jflc(ne)=0
  205   continue
      else
        do 208 ne=1,nec
        jflc(ne)=2
        if(ihomo.eq.niti) jflc(ne)=1
  208   continue
      end if
      if (nitloc.eq.0) go to 9999
      rmswt=0.0
      rmslst=10.0
      rmsbfl=20.0
c
c  location iteration loop
    1 continue
c
c zero out ahc
      do 35 i=1,mxobsc
      do 35 j=1,mxobsc
      ahc(i,j)=0.0
   35 continue
c
c Loop over events in the cluster
      nobs1=0
      nwrce=0
      ssqrwce=0.0
      wnormce=0.0
      do 250 ne=1,nec
      jflag=jflc(ne)
      nobs=kobsce(ne,nc)
      if(ne.gt.1) nobs1=nobs1+kobsce(ne-1,nc)
c  event coordinates
      xe=cevc(1,ne,nc)
      ye=cevc(2,ne,nc)
      ze=cevc(3,ne,nc)
c
      naz=0
c  loop over all observations of this event
      nwavce(1,ne)=0
      nwavce(2,ne)=0
      do 10 no=1,nobs
c  check for P or S reading 
      isp=intspc(no,ne,nc)
      nwavce(isp+1,ne)=nwavce(isp+1,ne) +1
      ns=istoc(no,ne,nc)
      xr=stc(1,ns)
      yr=stc(2,ns)
      zr=stc(3,ns)
c  If small change in hypocenter location, start with
c  old pb 3D path for additional hypocenter only iterations
      if(i3d.lt.2) goto 120
      if(plce(no,ne).eq.0.0) plce(no,ne)=rdltac(no,ne,nc)
c  small change is defined at beginning of loceqk
c      if(no.eq.1) write(78,*) 'ne,nc',ne,nc,' xmovc(ne) ',
c     & xmovc(ne),' small ',small
      if(xmovc(ne).gt.small) goto 120
c  Move end of raypath
      dr5(1)=xe-rpce(1,1,no,ne)
      dr5(2)=ye-rpce(2,1,no,ne)
      dr5(3)=ze-rpce(3,1,no,ne)
      do 105 k=1,3
        dr5i(k)=dr5(k)/5.0
  105 continue
      do 115 j=1,5
        dj=float(5-(j-1))
        do 110 k=1,3
          rpce(k,j,no,ne)=rpce(k,j,no,ne)+dr5i(k)*dj
  110   continue
  115 continue
      ttime=ttcce(no,ne)
      goto 125
c  determine 3-d path using circular raypaths.
c  determine approximate path
  120 ncrold=ncoc(no,ne,nc)
      ndpold=ndoc(no,ne,nc)
c      PRINT *,'BEFORE Rayweb'
      call rayweb(1,ne,no,isp,xe,ye,ze,xr,yr,zr,
     * fstime,jflag,ncrold,ndpold)
      ncoc(no,ne,nc)=ncrold
      ndoc(no,ne,nc)=ndpold
c
c  calculate residual and hypocentral partial derivatives
      ttime=fstime
cc check 
c      write(76,7601) nc,ne,no,stn(ns),isp,ttime,xe,ye,ze,xr,yr,zr
 7601 format(' Locclu for nc,ne,no',3i4,' stn:',a4,' isp',i2,' ttime '
     & ,'after rayweb=',f8.3,/,
     & '  xe,ye,ze=',3f8.2,', xr,yr,zr=',3f8.2)
c
      if (i3d.lt.2) go to 12
  125 continue
c  number of pb iter depends on distance
      nitpbu=nitpb(1)
      if(rdltac(no,ne,nc).gt.delt1) nitpbu=nitpb(2)
c      PRINT *,'BEFORE Minima'
      call minima(1,ne,no,isp,ttime,nitpbu,jpb)
      if(jpb.lt.nitpbu) goto 130
        write(26,2601) no,stn(ns),jpb,nitpbu,stn6(ns)
 2601   format(' Minima: no=',i4,2x,a4,', used maximum ',
     2  'number PB iter.: j=',i3,', nitpb=',i3,2x,a6)
  130 continue
   12 continue
c
      call ttmderce(nc,ne,no,0,ttime,niti,nnodej)
c        write(16,305)ns,ttime
c305    format(' for station ',i8,' the ttime= ',f7.3)
      ttcce(no,ne)=ttime
      if(niti.lt.nitmax) goto 10
        x1=rpce(1,1,no,ne)
        x2=rpce(1,2,no,ne)
        y1=rpce(2,1,no,ne)
        y2=rpce(2,2,no,ne)
        z1=rpce(3,1,no,ne)
        z2=rpce(3,2,no,ne)
        call aztoa(x1,x2,y1,y2,z1,z2,xr,yr,azim,tkofan)
        az(no)=azim
        toa(no)=tkofan
        if(intsp(no,ne).eq.1) goto 10
        naz=naz+1
        azp(naz)=az(no)
        dltp(naz)=dltac(no,ne,nc)
   10 continue
c
c      if(niti.ge.nitmax) call getgap(azp,dltp,naz,gapp,dltmn)
      if(niti.ge.nitmax) then
        call getgap(azp,dltp,naz,gapp,dltmn)
        igapce(ne,nc)=nint(gapp)
        idltmnce(ne,nc)=nint(dltmn)
      endif
c
      jflc(ne)=2
      jflag=2
c  calculate weights and apply to hypocenter matrix;
c  calculate rms residual
c  place derivs + resids into a-matrix for inversion
      call wthypc(nc,ne,nobs1,nobs,ahc,rmswte,nwr,ssqrwe,wnorme,wc)
      ssqrwce=ssqrwce+ssqrwe
      wnormce=wnormce+wnorme
      wrmsrce(ne,nc)=rmswte
      if((kout2.eq.0).or.(kout2.eq.2).or.(nwr.lt.4.).or. 
     2 ((rmswte.gt.wrmscemx).and.(nitl.gt.0)))  
     3 write(16,1001)ne,nitl,xe,ye,ze,secoce(ne,nc),rmswte
 1001 format(i3,' location, o.t., rmswt on iteration',i3,' are -',3f7.2,
     2 ',',f7.2,',',f6.3)
      if((kout2.eq.0).or.(kout2.eq.2).or.(nwr.lt.4.).or.
     2 ((rmswte.gt.wrmscemx).and.(nitl.gt.0)))  
     3 write(6,1011)nc,ne,nitl,xe,ye,ze,secoce(ne,nc),rmswte
 1011 format(i4,i3,' location, o.t., rmswt on iteration',i3,
     2 ' are -',3f7.2,',',f7.2,',',f6.3)

c  check for too few readings
      if(nwr.lt.4) then
        write(16,1664) nwr
 1664   format(' * * * * STOP * * * * Too few obs: nwr=',i3)
        stop
      endif
      nwrce=nwrce+nwr
c end of cevents loop
  250 continue
      rmswt=sqrt(ssqrwce/wnormce)
c
c If this is initial locclu loop, then set up edt obs
      if((niti.eq.0).and.(nitl.eq.0)) then
        call getcedt(nc,nobedt)
        nobsct(nc)=nobsc(nc)+nobedt
      endif
c for t* eq are shots and do not need location
      if(iuseq.eq.1) return
c nobsct is total obs tt+edt
      nobst=nobsct(nc)
c
c Next set up residuals and partials for eq dt
c  wthypedt section to get residuals and partials for edt obs
      nobsedt=kobsedt(nc)
      ssqrwcd=0.0
      wnormcd=0.0
c print check
c      write(76,7602) nc
 7602 format('edt obs for cluster:',i4,/,
     & ' noed   wc  ttmedt  resedt jev1 ttmcev1 jev2 ttmcev2')
      do 300 noe=1,nobsedt
      k=nobsc(nc)+noe
      job1=jobsed(1,noe,nc)
      job2=jobsed(2,noe,nc)
      jev1=jeved(1,noe,nc)
      jev2=jeved(2,noe,nc)
c revise edtsec for p ttimes
      if(intsped(noe,nc).eq.0) then
        edtsec(noe,nc)=(secpc(job1,jev1,nc)-secoce(jev1,nc))
     &     - (secpc(job2,jev2,nc)-secoce(jev2,nc))
      endif
      ttmedt(noe)=ttmc(job1,jev1)-ttmc(job2,jev2)
      resedt(noe)=edtsec(noe,nc)-ttmedt(noe)
c  make wt be average of the 2 wts
      wc(k)=(wtcombc(job1,jev1,nc)+wtcombc(job2,jev2,nc))*0.5
c apply wtepdt to edt obs
      wc(k)=wc(k)*wtepdt
c print out for checking
c      write(76,7603) noe,wc(k),ttmedt(noe),resedt(noe),
c     &  jev1,ttmc(job1,jev1),jev2,ttmc(job2,jev2)
 7603 format(i5,f5.1,2f8.3,i5,f8.3,i5,f8.3)
c      write(76,7604) stn(istoc(job1,jev1,nc)),secpc(job1,jev1,nc),
c     & ttcce(job1,jev1),stn(istoc(job2,jev2,nc)),
c     & secpc(job2,jev2,nc),ttcce(job2,jev2)
 7604 format(' jev1:stn ',a4,' secp=',f7.2,', ttcce=',f7.2,
     & ', jev2:stn ',a4,' secp=',f7.2,', ttcce=',f7.2)
      ja1=nparhy*(jev1-1)
      ja2=nparhy*(jev2-1)
      wnormcd=wnormcd+wc(k)*wc(k)
      wres=resedt(noe)*wc(k)
      ssqrwcd=ssqrwcd+wres*wres
      do 285 j=1,nparhy
        ahc(k,ja1+j)=dthc(job1,j,jev1)*wc(k)
        ahc(k,ja2+j)= -1.0*dthc(job2,j,jev2)*wc(k)
  285 continue
      ahc(k,jahy+1)=wres
  300 continue
      rmswtcd=sqrt(ssqrwcd/wnormcd)
      wnormclu=wnormce+wnormcd
      rmswtclu=sqrt((ssqrwce+ssqrwcd)/wnormclu)
      wrmsclu(nc)=rmswtclu
c  write out rms of cluster
      if((kout2.eq.0).or.(kout2.le.3)) write(16,1698) nc,nitl,rmswt,
c      if((kout2.eq.0).or.(kout2.eq.2)) write(16,1698) nc,nitl,rmswt,
     2 wnormce,rmswtcd,wnormcd,rmswtclu,wnormclu
 1698 format('Cluster:',i3,' Iter=',i2,' tt obs: wrms=',f6.3,
     2 ' wnormce=',f9.1,'; eq dt obs: wrms=',f6.3,' wnormcd=',f9.1,
     3 '; For all obs: wrms=',f6.3,' wnormclu=',f9.1)
c  terminate or continue iterating?
c  check for worsened location, over 2 iterations
      if(rmswtclu.gt.rmslst) then
        if(rmswtclu.gt.rmsbfl) then
          write(16,1633) rmsbfl, rmslst,rmswtclu
 1633     format(' Rms Increasing:',3f7.3,'; stop iterating')
          goto 999
        endif
      endif
      if (nitl.ge.nitloc) go to 999
c ************* wap *****************************
c  if rapid convergence, continue despite rmscut
c     write(16,1666) rmswtclu,rmslst,rmscut
c 1666 format(' rmswtclu,rmslst,rmscut ',3f6.3)
      if((rmslst-rmswtclu).ge.0.02) goto 15
c  quit if 2 iterations better than rmscut
      if (rmswtclu.le.rmscut) then
        if(rmslst.le.rmscut) then
          if(rmsbfl.le.rmscut) goto 999
        endif
      endif
c  WP added this to end iteration when convergence in reached.
c  terminate locclu when minor improvement in rms, since then
c  can get poorer hyps
      if((rmslst-rmswtclu).lt.0.007) then
          write(16,1634) rmsbfl, rmslst,rmswtclu
 1634     format(' Rms Not Improving:',3f7.3,'; stop iterating')
          goto 999
      endif
c  find new hypocenters for cluster
   15 continue
      rmsbfl=rmslst
      rmslst=rmswtclu
      nitl=nitl+1
c      PRINT *,'Locclu nitl=',nitl
c  determine singular value decomposition of hypocenter matrix
c      call fksvd(ah,s,v,maxobs,5,nobs,nparhy,1,.false.,.true.)
c      write(78,*) 'mxobsc,mxhparc,nobst,jahy',mxobsc,mxhparc,nobst,jahy
      call fksvd(ahc,s,v,mxobsc,mxhparc,nobst,jahy,1,.false.,.true.)
c  calculate adjustments and apply to hypocentral parameters
c  determine number of non-zero singular values
      nfre=0
      do 20 i=1,jahy
      if (s(i).gt.eigtol) nfre=nfre+1
   20 continue
      if(nfre.ne.jahy) write(78,*) 'nc=',nc,' nfre=',nfre,' jahy=',jahy
c  calculate adjustments for earthquakes
c put in 4 rather than nparhy
      do 330 ne=1,nec
        do 30 i=1,4
          ia=(ne-1)*4+i
          adj(ia)=0.0
          adja(i)=0.0
          do 28 j=1,nfre
            adja(i)=adja(i)+v(ia,j)*ahc(j,jahy+1)/s(j)
   28     continue
          adj(ia)=adja(i)
   30   continue
c  check for poor depth constraint
        i4=ne*4
        adj(i4)=adj(i4)*0.75    !Cliff uses 0.75
        if(s(i4).le.eigtol) write(16,9876) adj(i4)
 9876 format(' *** poorly constrained depth ** Do not adjust',
     & ' (adj=',f6.2,') ***')
        if (s(i4).le.eigtol) adj(i4)=0.0
c  apply adjustments to hypocentral parameters
        i1=(ne-1)*4+1
        secoce(ne,nc)=secoce(ne,nc)+adj(i1)
        dot=adj(i1)
        xmov=dot*dot
        do 40 i=1,3
          ia=(ne-1)*4+1+i
c  check to see that hyp. adjustment is less than dxmax
c          if(adj(ia)) 32,38,34
c   32      if(adj(ia).lt.(-dxmax)) adj(ia)= -dxmax
c           goto 38
c   34      if(adj(ia).gt.dxmax) adj(ia)=dxmax
          if(adj(ia).lt.0.0) then
            if(adj(ia).lt.(-dxmax)) adj(ia)= -dxmax
          else if(adj(ia).gt.0.0) then
            if(adj(ia).gt.dxmax) adj(ia)=dxmax
          endif
   38     cevc(i,ne,nc)=cevc(i,ne,nc)+adj(ia)
          xmov=xmov+adj(ia)*adj(ia)
   40   continue
c  check for depth less than zmin
      if(cevc(3,ne,nc).gt.zmin) goto 45
        write(16,1610) zmin,cevc(3,ne,nc)
 1610 format(5x,'***** Depth Fixed at ZMIN *****',' zmin, zadj= ',2f7.2)
        cevc(3,ne,nc)=zmin
   45 continue
      xmov=sqrt(xmov)
      xmovc(ne)=xmov
      if(xmovc(ne).gt.5.0) write(78,*) 'ne,nc',ne,nc,' Big xmovc(ne) ',
     & xmovc(ne)
  330 continue
      go to 1
  999  continue
ccc  if location only, calculate ssqr here instead of parsep
cc      if(nitmax.gt.0) goto 48
ccc  compute contribution to ssqr
cc      sqwtsh=sqrt(wtsht)
cc      totrms(ne)=0.0
cc      do 18 i=1,nobs
cc         resi=res(i)
cc         resi2=resi*resi
cc         totrms(ne)=totrms(ne)+resi2
cc         rrw=resi2*w(i)
cc         ssqrw=ssqrw+rrw
cc         wnobt=wnobt+w(i)
cc         if(intsp(i,ne).eq.0) then
ccc  P arrival
cc            ssqrwp=ssqrwp+rrw
cc            wnobtp=wnobtp+w(i)
cc         else
ccc  S arrival
cc            ssqrws=ssqrws+rrw
cc            wnobts=wnobts+w(i)
cc         endif
cc         if (ne.gt.netemp) resi=resi*sqwtsh
cc         ssqr=ssqr+resi*resi
cc   18 continue
cc      totrms(ne)=sqrt(totrms(ne)/nobs)
c
c  Covariance calculations for hypocenter parameters
c  Compute diagonal elements of covariance matrix
c  set reading error (sec).
c 16feb2010 compute like loceqk, but later could evaluate for 
c           cluster earthquake hypocenter
   48 continue
      sig=rderr
      sig2=sig*sig
      ibad=0
      nbad=0
      do 350 ne=1,nec
c  initially set st.er. to zero (so correct for blasts)
      do 50 i=1,nparhy
         sum=0.0
         ia=(ne-1)*4+i
         ster(ia)=0.0
         do 60 j=1,nparhy
            js=(ne-1)*4+j
            sss=s(js)
c fix s so won't get arithmetic fault
            if(sss.eq.0.0) sss=0.00000001
            sum=sum+v(ia,js)*v(ia,js)/(sss*sss)
   60    continue
         cov(i)=sig2*sum
         ster(i)=sqrt(cov(i))
   50 continue
        write(16,1099)ne,nitl,wrmsrce(ne,nc),secoce(ne,nc),
     2  (cevc(j,ne,nc),j=1,3),(ster(i),i=1,4)
cc PRINT
c        write(6,1099)ne,nitl,wrmsrce(ne,nc),secoce(ne,nc),
c     2  (cevc(j,ne,nc),j=1,3),(ster(i),i=1,4)
 1099   format(1x,'Cl-Event ',i4,';nit=',i2,'; rmswte=',f6.3,
     2  ';ot,x,y,z =',f6.2,2f8.2,f7.2,
     3  '; St.Er.=',f8.4,'(ot)',f7.4,'(x)',f7.4,'(y)',
     4  f9.4,'(z)')
c check for too high wrms of this event
c delete bad events and then restart locclu
      if(wrmsrce(ne,nc).gt.wrmscemx) then
        if(niti.eq.0) then
          ibad=1
          nbad=nbad+1
          numbad(nbad)=ne
        else
          write(16,1657) ne,nc,iyrmoce(ne,nc),idayce(ne,nc),
     *    ihrce(ne,nc),minoce(ne,nc),secoce(ne,nc)
 1657     format(' ** WARNING SHOULD DELETE:  Cluster:',i5,' DELETING ',
     2   'EVENT',i5,', with bad wrms: ',a4,a2,1x,a2,i2,1x,f5.2)
        endif
      endif
  350 continue
      if(ibad.eq.1) then
        do 355 i=1,nbad
          neb=numbad(i)
          nwav(1)=nwavce(1,neb)
          nwav(2)=nwavce(2,neb)
          call delcev(neb,nc,nwav)
  355   continue
        write(16,1642) nc,nbad,wrmscemx
 1642   format(' * * *',/,' Restarting LOCCLU for cluster',i5,' after',
     2    12x,'deleted',i4,' bad (high wrms) events. wrmsce_max=',
     3    f7.3,/,' * * *')
        goto 11
      endif
c  write out rms of cluster
      write(16,1699) nc,rmswt,wnormce,rmswtcd,wnormcd,rmswtclu,wnormclu
 1699 format('Cluster:',i3,' tt obs: wrms=',f6.3,' wnormce=',f9.1,
     2 '; eq dt obs: wrms=',f6.3,' wnormcd=',f9.1,
     3 '; For all obs: wrms=',f6.3,' wnormclu=',f9.1)
c  (Fred Klein's stuff deleted for cluster)
 9999 continue
c  update dlta
      do 85 ne=1,nec
        nobs=kobsce(ne,nc)
        do 80 no=1,nobs
        ns=istoc(no,ne,nc)
        dx=cevc(1,ne,nc)-stc(1,ns)
        dy=cevc(2,ne,nc)-stc(2,ns)
        dz=cevc(3,ne,nc)-stc(3,ns)
        dltac(no,ne,nc)=sqrt(dx*dx+dy*dy)
        rdltac(no,ne,nc)=sqrt(dx*dx+dy*dy+dz*dz)
   80   continue
   85 continue
      return
c***** end of subroutine locclu *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine loceqk(ne,niti,nwr)
c  this routine locates the given event in the initial velocity
c  structure using approximate ray tracing
c
c  common block variables:
      common/shortd/ xltkm,xlnkm,rota,nzco,xlt,xln,snr,csr,
     2 onorth,oeast,rotadg,cmerid,rlatdo,rlondo
      include 'simul2014_common.inc'
c
c  declaration statements:
      real dr5(3),dr5i(3),s(5),v(5,5),adj(5),ster(5),cov(5)
      real w(maxobs)
      double precision ah(maxobs,maxobs)
c  variables for error ellipse calculations
      double precision ac(3,3)
      integer iaz(3),idip(3)
      real covar(5,5),v3(5,5),serr(3)
      parameter (drad=1.7453292d-02)
c  variables for gap calculation
      real azp(maxobs),dltp(maxobs),gapp
      integer naz
c
c  Define small change in hypocenter, depending on bld factor
      small=0.5
      if(bld.eq.0.1) small=0.1
c
      xmov=900.0
c  nitl=counter for location iterations
      nitl=0
c  niti=counter for inversion iterations, niti=0 for initial locations
      if(niti.eq.0) then
        jfl=0
      else
        jfl=2
        if(ihomo.eq.niti) jfl=1
      end if
      jflag=jfl
      if (nitloc.eq.0) go to 9999
      rmswt=0.0
      rmslst=99.0
      rmsbfl=20.0
      nobs=kobs(ne)
c  location iteration loop
    1 continue
c  event coordinates
      xe=evc(1,ne)
      ye=evc(2,ne)
      ze=evc(3,ne)
c
      naz=0
c  loop over all observations of this event
      do 10 no=1,nobs
c  check for P or S reading 
      isp=intsp(no,ne)
      ns=isto(no,ne)
      xr=stc(1,ns)
      yr=stc(2,ns)
      zr=stc(3,ns)
c  If small change in hypocenter location, start with
c  old pb 3D path for additional hypocenter only iterations
      if(i3d.lt.2) goto 120
      if(pl(no).eq.0.0) pl(no)=rdlta(no,ne)
c     if(xmov.gt.(0.05*pl(no))) goto 120
c  small change is defined at beginning of loceqk
      if(xmov.gt.small) goto 120
c  Move end of raypath
      dr5(1)=xe-rp(1,1,no)
      dr5(2)=ye-rp(2,1,no)
      dr5(3)=ze-rp(3,1,no)
      do 105 k=1,3
        dr5i(k)=dr5(k)/5.0
  105 continue
      do 115 j=1,5
        dj=float(5-(j-1))
        do 110 k=1,3
          rp(k,j,no)=rp(k,j,no)+dr5i(k)*dj
  110   continue
  115 continue
      ttime=ttc(no)
      goto 125
c  determine 3-d path using circular raypaths.
c  determine approximate path
  120 ncrold=nco(no,ne)
      ndpold=ndo(no,ne)
      call rayweb(0,ne,no,isp,xe,ye,ze,xr,yr,zr,
     * fstime,jflag,ncrold,ndpold)
      nco(no,ne)=ncrold
      ndo(no,ne)=ndpold
c
c  calculate residual and hypocentral partial derivatives
      ttime=fstime
c
      if (i3d.lt.2) go to 12
  125 continue
c  number of pb iter depends on distance
      nitpbu=nitpb(1)
      if(rdlta(no,ne).gt.delt1) nitpbu=nitpb(2)
      call minima(0,ne,no,isp,ttime,nitpbu,jpb)
      if(jpb.lt.nitpbu) goto 130
        write(26,2601) no,stn(ns),jpb,nitpbu,stn6(ns)
 2601   format(' Minima: no=',i4,2x,a4,', used maximum ',
     2  'number PB iter.: j=',i3,', nitpb=',i3,2x,a6)
  130 continue
   12 continue
c
      call ttmder(ne,no,0,ttime,niti,nnodej)
c        write(16,305)ns,ttime
c305    format(' for station ',i8,' the ttime= ',f7.3)
      ttc(no)=ttime
      if(niti.lt.nitmax) goto 10
        x1=rp(1,1,no)
        x2=rp(1,2,no)
        y1=rp(2,1,no)
        y2=rp(2,2,no)
        z1=rp(3,1,no)
        z2=rp(3,2,no)
        call aztoa(x1,x2,y1,y2,z1,z2,xr,yr,azim,tkofan)
        az(no)=azim
        toa(no)=tkofan
        if(intsp(no,ne).eq.1) goto 10
        naz=naz+1
        azp(naz)=az(no)
        dltp(naz)=dlta(no,ne)
   10 continue
c
      if(niti.ge.nitmax) call getgap(azp,dltp,naz,gapp,dltmn)
      igap(ne)=nint(gapp)
      idltmn(ne)=nint(dltmn)
c
      jfl=2
      jflag=2
c  calculate weights and apply to hypocenter matrix;
c  calculate rms residual
c  place derivs + resids into a-matrix for inversion
      call wthyp(ne,nobs,ah,rmswt,nwr,w)
      if((kout2.eq.0).or.(kout2.eq.2)) write(16,1001) nitl,xe,ye,ze,
     2 seco(ne),rmswt,ne,nobs
      if((kout2.eq.0).or.(kout2.eq.2))
     2 WRITE(6,1001) nitl,xe,ye,ze,seco(ne),rmswt,ne,nobs
 1001 format('  location, o.t., rmswt on iteration',i3,' are -',3f7.2,
     2 ',',f7.2,',',f6.3,' ne',i6,' nobs',i6)
c  change for nwr.lt.2 for blasts
      if(ne.le.neqs) then
        if(nwr.lt.4) return
      else
        if(nwr.lt.2) return
      endif
c  terminate or continue iterating?
c  check for worsened location, over 2 iterations
      if(rmswt.gt.rmslst) then
        if(rmswt.gt.rmsbfl) then
          evc(1,ne)=xbflst
          evc(2,ne)=ybflst
          evc(3,ne)=zbflst
          seco(ne)=otbflt
          seco2(ne)=ot2bflt
          rmswt=rmsbfl
          goto 999
        endif
      endif
      if (nitl.ge.nitloc) go to 999
c ************* wap *****************************
c  if rapid convergence, continue despite rmscut
c     write(16,1666) rmswt,rmslst,rmscut
 1666 format(' rmswt,rmslst,rmscut ',3f6.3)
      if((rmslst-rmswt).ge.0.02) goto 15
c  quit if 2 iterations better than rmscut
      if (rmswt.le.rmscut) then
        if(rmslst.le.rmscut) then
          if(rmsbfl.le.rmscut) goto 999
        endif
      endif
c  WP added this to end iteration when convergence in reached.
c  GE terminate when no improvement of rms in 2 iterations
      if((rmslst-rmswt).lt.0.0002) then
        if((rmsbfl-rmslst).lt.0.0002) goto 999
      endif
      if(ne.gt.neqs) goto 15
      if(xmov.lt.0.05) goto 999
c  find new hypocenter
   15 continue
      xbflst=xlst
      ybflst=ylst
      zbflst=zlst
      otbflt=otlst
      ot2bflt=ot2lst
      rmsbfl=rmslst
      xlst=evc(1,ne)
      ylst=evc(2,ne)
      zlst=evc(3,ne)
      otlst=seco(ne)
      ot2lst=seco2(ne)
      rmslst=rmswt
      nitl=nitl+1
c  determine singular value decomposition of hypocenter matrix
      call fksvd(ah,s,v,maxobs,5,nobs,nparhy,1,.false.,.true.)
c  calculate adjustments and apply to hypocentral parameters
c  determine number of non-zero singular values
      nfre=0
      do 20 i=1,nparhy
      if (s(i).gt.eigtol) nfre=nfre+1
   20 continue
c      print *,'nfre=',nfre,'nparhy=',nparhy
c  calculate adjustments
c   only adjust origin time for blast
      if(ne.gt.neqs) then
c    calculate adjustment
        adj(1)=0.0
        do 25 j=1,nfre
        adj(1)=adj(1)+v(1,j)*ah(j,nparhy+1)/s(j)
        if(iuse2t.gt.0) adj(5)=adj(5)+v(5,j)*ah(j,nparhy+1)/s(j)
   25   continue
c    apply adjustment
        seco(ne)=seco(ne)+adj(1)
        if(iuse2t.gt.0) seco2(ne)=seco2(ne)+adj(5)
c  calculate adjustments for earthquakes
      else
        do 30 i=1,nparhy
           adj(i)=0.0
           do 29 j=1,nfre
              adj(i)=adj(i)+v(i,j)*ah(j,nparhy+1)/s(j)
   29      continue
   30   continue
c  check for poor depth constraint
        adj(4)=adj(4)*0.75    !Cliff uses 0.75
        if(s(4).le.eigtol) write(16,9876) adj(4)
 9876 format(' *** poorly constrained depth ** Do not adjust',
     & ' (adj=',f6.2,') ***')
        if (s(4).le.eigtol) adj(4)=0.0
c  apply adjustments to hypocentral parameters
        seco(ne)=seco(ne)+adj(1)
        if(nparhy.eq.5) seco2(ne)=seco2(ne)+adj(5)
        do 40 i=1,3
           i1=i+1
c  check to see that hyp. adjustment is less than dxmax
           if(adj(i1)) 32,38,34
   32      if(adj(i1).lt.(-dxmax)) adj(i1)= -dxmax
           goto 38
   34      if(adj(i1).gt.dxmax) adj(i1)=dxmax
   38      evc(i,ne)=evc(i,ne)+adj(i1)
   40   continue
c  check for depth less than zmin
      if(evc(3,ne).gt.zmin) goto 45
        write(16,1610) zmin,evc(3,ne)
 1610 format(5x,'***** Depth Fixed at ZMIN *****',' zmin, zadj= ',2f7.2)
        evc(3,ne)=zmin
        goto 999
   45 end if
      dx=xlst-evc(1,ne)
      dy=ylst-evc(2,ne)
      dz=zlst-evc(3,ne)
      dot=otlst-seco(ne)
      xmov=sqrt(dx*dx + dy*dy + dz*dz + dot*dot)
      go to 1
  999  continue
c  residual reduction achieved - proceed to parameter separation
c  if location only, calculate ssqr here instead of parsep
      if(nitmax.gt.0) goto 48
c  compute contribution to ssqr
      sqwtsh=sqrt(wtsht)
      totrms(ne)=0.0
      do 18 i=1,nobs
         resi=res(i)
         resi2=resi*resi
         totrms(ne)=totrms(ne)+resi2
         rrw=resi2*w(i)
         ssqrw=ssqrw+rrw
         wnobt=wnobt+w(i)
         if(intsp(i,ne).eq.0) then
c  P arrival
            ssqrwp=ssqrwp+rrw
            wnobtp=wnobtp+w(i)
         else
c  S arrival
            ssqrws=ssqrws+rrw
            wnobts=wnobts+w(i)
         endif
         if (ne.gt.netemp) resi=resi*sqwtsh
         ssqr=ssqr+resi*resi
   18 continue
      totrms(ne)=sqrt(totrms(ne)/float(nobs))
c
c  Covariance calculations for hypocenter parameters
c  Compute diagonal elements of covariance matrix
c  set reading error (sec).
   48 sig=rderr
      sig2=sig*sig
c  initially set st.er. to zero (so correct for blasts)
      do 47 i=1,nparhy
         ster(i)=0.0
   47 continue
      ni=nparhy
      if (ne.gt.neqs) ni=1
      do 50 i=1,ni
         sum=0.0
         do 60 j=1,nparhy
            sss=s(j)
c fix s so won't get arithmetic fault
            if(sss.eq.0.0) sss=0.00000001
            sum=sum+v(i,j)*v(i,j)/(sss*sss)
   60    continue
         cov(i)=sig2*sum
         ster(i)=sqrt(cov(i))
   50 continue
      wrmsr(ne)=rmswt
      if(iuse2t.eq.0) then
c        WRITE(6,1099)ne,nitl,rmswt,seco(ne),(evc(j,ne),j=1,3),
c     2  (ster(i),i=1,4)
        write(16,1099)ne,nitl,rmswt,seco(ne),(evc(j,ne),j=1,3),
     2  (ster(i),i=1,4)
 1099   format(1x,'event ',i5,';nit=',i2,';wtd rms res=',f6.3,
     2  ';ot,x,y,z =',f6.2,2f8.2,f7.2,
     3  '; St.Er.=',f8.4,'(ot)',f7.3,'(x)',f7.3,'(y)',
     4  f9.3,'(z)')
      else
        write(16,1097)ne,nitl,rmswt,seco(ne),(evc(j,ne),j=1,3),
     2  seco2(ne),ster
 1097   format(1x,'event ',i4,';nit=',i3,';wrms=',f6.3,
     2  ';ot,x,y,z,ot2=',f6.2,2f8.2,f7.2,f6.2,
     3  ';St.Er.=',f6.3,'(ot)',f7.3,'(x)',f7.3,'(y)',
     4  f9.3,'(z)',f6.3,'(ot2)')
      endif
c  this is from Fred Klein's program HYPOINVERSE (6mar1989)
c  set Fred's "n" to 4
      n=4
C--SKIP THE REMAINING CALCS IF THEY ARE NOT NEEDED FOR FINAL OUTPUT
      if((kout.lt.5).or.(ne.gt.neqs)) goto 9999
C--CALCULATE COVARIANCE MATRIX AS SIGMA**2 * V * EIGVAL**-2 * VT
C--ESTIMATED ARRIVAL TIME ERROR
        if(ercof.le.0.0) then
          write(6,1612) ercof
          write(16,1612) ercof
 1612     format('*** STOP! ercof=',f5.2,'  For kout=5, computing ',/,
     5    '   hypoellipse-style errors, ERCOF should be >0   ***')
          stop
        endif
	SIGSQ=RDERR*rderr+ERCOF*rmswt*rmswt
	TEMP=EIGTOL**2
	DO 252 I=1,4
	DO 250 J=1,I
	COVAR(I,J)=0.
	IF (I.GT.N .OR. J.GT.N) THEN
	  IF (I.EQ.J) COVAR(I,J)=999.
	  GO TO 250
	END IF
	DO 245 L=1,N
c 245 COVAR(I,J)=COVAR(I,J)+V(I,L)*V(J,L)/(s(L)**2+TEMP)
245	COVAR(I,J)=COVAR(I,J)+V(I,L)*V(J,L)/(s(L)*s(L))
	COVAR(I,J)=SIGSQ*COVAR(I,J)
250	COVAR(J,I)=COVAR(I,J)
252	CONTINUE

C--EVALUATE THE HYPOCENTER ERROR ELLIPSE BY DIAGONALIZING
C--THE SPATIAL PART OF THE COVARIANCE MATRIX
C--USE ac AS TEMPORARY STORAGE
	DO 257 I=1,3
	DO 255 J=1,3
255	ac(I,J)=COVAR(I+1,J+1)
257	CONTINUE
        mmax=5
	CALL fksvd (ac,SERR,V3,MMAX,3,3,3,0,.FALSE.,.TRUE.)
	DO 260 I=1,3
	SERR(I)=SQRT(SERR(I))
	IF (SERR(I).GT.99.) SERR(I)=99.
260	CONTINUE

C--COMPUTE ERH AND ERZ AS THE LARGEST OF THE HORIZ AND VERTICAL
C--PROJECTIONS OF THE PRINCIPAL STANDARD ERRORS
	ERH=0.
	ERZ=0.
	DO 265 I=1,3
	TEMP=SERR(I)*SQRT(V3(1,I)**2+V3(2,I)**2)
	IF (TEMP.GT.ERH) ERH=TEMP
	TEMP=SERR(I)*ABS(V3(3,I))
	IF (TEMP.GT.ERZ) ERZ=TEMP
265	CONTINUE
	IF (ERZ.GT.99.) ERZ=99.
	IF (ERH.GT.99.) ERH=99.
c     write(16,1098)ne,nitl,rmswt,seco(ne),(evc(j,ne),j=1,3),
c    2 erh,erz
 1098 format(2x,'event ',i4,'; nit=',i3,'; wtd rms res=',f6.3,
     2 '; ot,x,y,z =',4f6.2,
     3 '; Hypinv.code St.Er.=',f10.3,' ERH',6x,f9.3,'ERZ')

C--NOW CALC THE ORIENTATIONS OF THE PRINCIPAL STD ERRORS
	DO 290 J=1,3
	IAZ(J)=0
	IDIP(J)=90
	TEMP=SQRT(V3(1,J)**2+V3(2,J)**2)
	IF (TEMP.EQ.0.) GO TO 290
c  Change indices because different coordinates. Hypoinverse array has 1=lat, 2=long(west=positive)
c     IAZ(J)=ATAN2(-V3(2,J),V3(1,J))/drad
c   Simul has evc(1)=west, evc(2)=north
c   and allows rotation
      iaz(j)=(atan2(-v3(1,j),v3(2,j))-rota)/drad
	IDIP(J)=ATAN2(V3(3,J),TEMP)/drad
	IF (IDIP(J).LT.0) THEN
	  IDIP(J)=-IDIP(J)
	  IAZ(J)=IAZ(J)+180
	END IF
	IF (IAZ(J).LT.0) IAZ(J)=IAZ(J)+360
        if(iaz(j).gt.360) iaz(j)=iaz(j)-360
290	CONTINUE
      call out(ne,niti,nwr,serr,iaz,idip,erh,erz)
 9999 continue
c  update dlta
      do 300 no=1,nobs
      ns=isto(no,ne)
      dx=evc(1,ne)-stc(1,ns)
      dy=evc(2,ne)-stc(2,ns)
      dlta(no,ne)=sqrt(dx*dx+dy*dy)
      dz=evc(3,ne)-stc(3,ns)
      rdlta(no,ne)=sqrt(dx*dx+dy*dy+dz*dz)
  300 continue
      return
c***** end of subroutine loceqk *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine ludecp(a,ul,n,ier)
c  routine to perform cholesky decomposition
c
c  common block variables:
      include 'simul2014_common.inc'
c
      dimension a(mgsol),ul(mgsol)
      data zero,one,four,sixtn,sixth/0.0,1.,4.,16.,.0625/
      d1=one
      d2=zero
      rn=one/(n*sixtn)
      ip=1
      ier=0
      do 45 i=1,n
         iq=ip
         ir=1
            do 40 j=1,i
            x=a(ip)
            if(j.eq.1) go to 10
               do 5 k=iq,ip1
               x=x-ul(k)*ul(ir)
               ir=ir+1
    5       continue
   10       if (i.ne.j) go to 30
            d1=d1*x
            if (a(ip)+x*rn.le.a(ip)) go to 50
   15       if (abs(d1).le.one) go to 20
            d1=d1*sixth
            d2=d2+four
            go to 15
   20       if (abs(d1).ge.sixth) go to 25
            d1=d1*sixtn
            d2=d2-four
            go to 20
   25       ul(ip)=one/sqrt(x)
            go to 35
   30       ul(ip)=x*ul(ir)
   35       ip1=ip
            ip=ip+1
            ir=ir+1
   40    continue
   45 continue
      go to 9005
   50 ier=129
 9000 continue
 9005 return
c* * * end of subroutine ludecp * * *
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine luelmp(a,b,n,x)
c  routine to perform elimination part of solution of ax=b
c
c  common block variables:
      include 'simul2014_common.inc'
c
      dimension a(mgsol),b(mxpri1),x(mxpri1)
      data zero/0./
c  solution of ly=b
      ip=1
      kw=0
      do 15 i=1,n
         t=b(i)
         im1=i-1
         if (kw.eq.0) go to 9
         ip=ip+kw-1
         do 5 k=kw,im1
            t=t-a(ip)*x(k)
            ip=ip+1
    5    continue
         go to 10
    9    if (t.ne.zero) kw=i
         ip=ip+im1
   10    x(i)=t*a(ip)
         ip=ip+1
   15 continue
c  solution of ux=y
      n1=n+1
      do 30 i=1,n
         ii=n1-i
         ip=ip-1
         is=ip
         iq=ii+1
         t=x(ii)
         if (n.lt.iq) go to 25
         kk=n
         do 20 k=iq,n
            t=t-a(is)*x(kk)
            kk=kk-1
            is=is-kk
   20    continue
   25    x(ii)=t*a(is)
   30 continue
      return
c***** end of subroutine luelmp *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine medder(nc,ne,kopt,nwr)
c routine to add medium derivatives from shot data
c  to the medium matrix dtmp for inversion
c  kopt =1 for shot
c  kopt =2 for receiver-pair differential time "shot"
c  kopt =3 for cluster event (treated as shot in velocity inversion)
c  kopt =4 for cluster earthquake-pair differential times
c  kopt =5 for teleseismic event
c
c  common block variables:
      include 'simul2014_common.inc'
      common/wtpars/w(mxobsa)
c
      if(kopt.eq.1.or.kopt.eq.5) nobs=kobs(ne)
      if(kopt.eq.2) then
        nobs=kobsrdt(ne)
c   Put rdt residuals into res array before medder and parsep
        do 95 no=1,nobs
        res(no)=resrdt(no)
   95   continue
      endif
      if(kopt.eq.3) then
        if(ne.eq.1) write(16,1632)
 1632   format('Getting medium derivatives, and parameter separation:')
        nobs=kobsce(ne,nc)
        do 100 no=1,nobs
        res(no)=resc(no,ne)
  100   continue
      endif
      if(kopt.eq.4) then
        nobs=kobsedt(nc)
        do 105 no=1,nobs
        res(no)=resedt(no)
  105   continue
      endif
      rms=0.0
c  compute weighting
c     as in parsep
c
c PRINT
c      WRITE(56,601) nobs,(res(j),w(j),dlta(j,ne),j=1,nobs)
  601 FORMAT('nobs,1 medder',i4,' res w delta',4(f8.2,f7.3,f10.2))
      nwr=0
      wnorm=0.0
      do 13 j=1,nobs
c  reading weight
         if(kopt.eq.1.or.kopt.eq.5) w(j)=wt(j,ne)
         if(kopt.eq.2) w(j)=wtrdtobs(j,ne)
         if(kopt.eq.3) w(j)=wtc(j,ne,nc)
         if(kopt.eq.4) w(j)=wtedtobs(j,nc)
c  residual weighting
c  downweighting(linear) 0 to 98% res1 to res2, 98 to 100% res2 to res3
         ares=abs(res(j))
         if(kopt.eq.5) goto 8
         if(ares.le.res2) then
            wr=1.0-(ares-res1)*dres12
            if (wr.gt.1.0) wr=1.0
         else
            if(res3.gt.res2) then
               wr=0.02-(ares-res2)*dres23
               if (wr.lt.0.0) wr=0.0
            else
               wr=0.0
            endif
         endif
         goto 9
    8    continue
         if(ares.le.reste2) then
            wr=1.0-(ares-reste1)*dreste12
            if (wr.gt.1.0) wr=1.0
         else
            if(reste3.gt.reste2) then
               wr=0.02-(ares-reste2)*dreste23
               if (wr.lt.0.0) wr=0.0
            else
               wr=0.0
            endif
         endif
         wd=1.0
         goto 11
c  distance weighting
    9    if(kopt.eq.1) dltaob=dlta(j,ne)
c  for receiver-pair dt, stations are close so okay to use 1st for dlta
         if(kopt.eq.2) dltaob=dlta(jobsrd(1,j,ne),ne)
         if(kopt.eq.3) dltaob=dltac(j,ne,nc)
c  For eq-pair, sources are close so okay to use 1st for dlta
         if(kopt.eq.4) then
           job1=jobsed(1,j,nc)
           jev1=jeved(1,j,nc)
           dltaob=dltac(job1,jev1,nc)
         endif
         wd=1.0-(dltaob-delt1)*ddlt
         if (wd.gt.1.0) wd=1.0
         if (wd.lt.0.0) wd=0.0
c  unnormalized weight
   11    w(j)=w(j)*wr*wd
         wnorm=wnorm+w(j)
         if (w(j).gt.0.0) nwr=nwr+1
  13  continue
C      PRINT *,'nobs,2 medder w(j) ',nobs, (w(j),j=1,nobs)
c PRINT
      if(kopt.eq.4) then
CPRINT        write(91,9101) nc,nobs
 9101 format('Medder cluster',i5,', nobs=',i6,', Printout res(j),wt(j)',
     2  ' wtedtobs(j,nc)')
CPRINT        write(91,9102)(res(j),w(j),wtedtobs(j,nc),j=1,nobs)
 9102   format(5(f9.3,2f7.3))
      endif
c  check to be sure 1 or more readings with nonzero weight
      if(nwr.ge.1) goto 12
      write(16,1612)
 1612 format(' ***** SKIPPING EVENT !!! No readings with nonzero ',
     &  'weights *****')
c  normalize weights
   12 wfac=nwr/wnorm
      if((kopt.eq.1).and.(ne.gt.netemp)) wfac=wfac*wtsht
      if(kopt.eq.2) wfac=wfac*wtrdtsht
      if(kopt.eq.4) wfac=wfac*wtepdt
      if(kopt.eq.5) wfac=wfac*wttel
      nswrt=nswrt+nwr
      do 14 j=1,nobs
         w(j)=w(j)*wfac
   14 continue
c  store residuals and compute rms residual
c  compute weighted rms residual as in wthyp
      wnorm=0.0
      rmswt=0.0
      do 10 ii=1,nobs
         resp(ii)= res(ii)*w(ii)
         w2=w(ii)*w(ii)
         rmswt= rmswt + w2*res(ii)* res(ii)
         wnorm=wnorm+w2
   10 continue
      rmswt=sqrt(rmswt/wnorm)
      if(kopt.eq.1) then
        wrmsr(ne)=rmswt
        write(16,1001) ne,rmswt,nwr
 1001 format(' explosion ',i4,'; weighted rms res =',f7.3,'; nwr=',i6)
      endif
      if(kopt.eq.2) then
        wrmsrdt(ne)=rmswt
        write(16,1002) ne,rmswt,nwr
 1002 format(' receiver-pair diff-time shot ',i4,'; weighted rms res =',
     *  f7.3,'; nwr=',i6)
      endif
      if(kopt.eq.3) then
        wrmsrce(ne,nc)=rmswt
        write(16,1003) nc,ne,rmswt,nwr
 1003 format(' Cluster:',i4,' Cl_Earthquake:',i4,'; weighted rms res =',
     *  f7.3,'; nwr=',i6)
      endif
      if(kopt.eq.4) then
        wrmsedt(nc)=rmswt
        write(16,1004) nc,rmswt,nwr
 1004 format(' Cluster:',i4,' Eq-diff-times; weighted rms res =',
     *  f7.3,'; nwr=',i6)
      endif
      if(kopt.eq.5) then
        wrmsr(ne)=rmswt
        write(16,1005) ne,rmswt,nwr
 1005 format(' tele event ',i4,'; weighted rms res =',f7.3,'; nwr=',i6)
      endif
c  store partial derivatives
      nono=npari*nobs
      do 20 j=1,nono
         i=(j-1)/npari+1
         dtmp(j)=dtm(j)*w(i)
c check station corrections
        im=j-(i-1)*npari
        is=mdexfx(im)-nparv
C        if((is.le.0).or.(dtm(j).ne.1.0000)) goto 20
        if(dtm(j).eq.0.00) goto 20
c        WRITE(16,1620) is,j,dtm(j),dtmp(j)
c1620    format(' MEDDER:is,j,dtm,dtmp',i5,i10,2f14.6)
   20 continue
c  zero out dtm matrix
      do 45 m=1,nono
         dtm(m)=0.0
   45 continue
cDEP Save weight
      do 60 j=1,nobs
      if(kopt.eq.1.or.kopt.eq.5) wtcomb(j,ne)=w(j)
      if(kopt.eq.2) wtcombrd(j,ne)=w(j)
      if(kopt.eq.3) wtcombc(j,ne,nc)=w(j)
      if(kopt.eq.4) wtcombed(j,nc)=w(j)
   60 continue
      return
c***** end of subroutine medder *****
      end
c
c-------------------------------------------------------
c
      subroutine minima(kopt,ne,no,isp,fstime,npbmax,jpb)
c
c*****this routine finds the minimum path using pseudo-bending
c if kopt=0, then standard earthquake
c    kopt=1, then cluster earthquake
c
c  common block variables:
      common/pathm/x(260),y(260),z(260),v(260),vq(260),tra,qtra,n,nn
      common/temp/xtemp(260),ytemp(260),ztemp(260),rtemp(260),ttemp(260)
      common/pat/xt(260),yt(260),zt(260),rt(260),rm,rs,rtm,rts,nt,j
      common/upb/ upbtot,nupb
      include 'simul2014_common.inc'
c
      tra=fstime
c
      if(kopt.eq.0) then
        n=nrp(no)
        do 15 i=1,n
         x(i)=rp(1,i,no)
         y(i)=rp(2,i,no)
         z(i)=rp(3,i,no)
   15   continue
      else
        n=nrpce(no,ne)
        do 20 i=1,n
         x(i)=rpce(1,i,no,ne)
         y(i)=rpce(2,i,no,ne)
         z(i)=rpce(3,i,no,ne)
   20   continue
      endif
c
      do 21 i=1,n
         xp=x(i)
         yp=y(i)
         zp=z(i)
         xtemp(i)=xp
         ytemp(i)=yp
         ztemp(i)=zp
         call vel3eft(isp,xp,yp,zp,vp)
         v(i)=vp
         call vel3eft(1,xp,yp,zp,vpp)
         vq(i)=vpp
   21 continue
      call travel
c
      nn=n-1
      nupb=nupb+1
c     write(6,6000) nupb
 6000 format(' nubp=',i8)
      do 100 j=1,npbmax
         ta=tra
         call bend(isp,xfac)
         call travel
         deltat=ta-tra
c
         if (deltat.lt.0.0) go to 102
         do 22 i=1,n
            x(i)=xtemp(i)
            y(i)=ytemp(i)
            z(i)=ztemp(i)
   22    continue
         if(deltat.le.tlim) go to 102
100   continue
c
  102 continue
      if(iuseq .eq. 1) then
        call qtravel
      endif
c
      if(j.gt.npbmax) j=npbmax
      upbtot=upbtot + float(j)
      jpb=j
      if(kopt.eq.0) then
        do 300 i=1,n
         rp(1,i,no)=x(i)
         rp(2,i,no)=y(i)
         rp(3,i,no)=z(i)
  300   continue
      else
        do 305 i=1,n
         rpce(1,i,no,ne)=x(i)
         rpce(2,i,no,ne)=y(i)
         rpce(3,i,no,ne)=z(i)
  305   continue
      endif
      if(iuseq .eq. 0) then
        fstime=tra
      else
      fstime=qtra
      endif
c
      return
c ***** end of subroutine minima *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine out(ne,nit,nwr,serr,iaz,idip,erh,erz)
c
c  declaration statements:
      integer iaz(3),idip(3)
      real serr(3)
	DIMENSION KSIG(3)
	CHARACTER KSFL*1
      character*1 is,ie
      character*24 date24
c
c  common block variables:
      include 'simul2014_common.inc'
c
      call fdate(date24)
c
c  Write new hypocenters to file12 in hypoinverse format
      if(ne.eq.1) write(12,1201) nit,date24(21:24),date24(4:16)
 1201 format(' iteration step=',i3,'  computed',1x,a4,a13,
     1 ' hyp error ellipse',
     2 /,' Date HrMn Sec   Lat    Long DepthMgNwr',7x,'Rms',
     3 'Az1D1Ser1Az2D2Ser2Mg',3x,'Ser3',
     4 5x,'Erh Erz')
      i=ne
      call latlon(evc(1,i),evc(2,i),lat,xlat,lon,xlon)
c       convert origin time into gmt
        sec= seco(i)
        nin=mino(i)
 210    if(sec.lt.0) goto 230
 220    if(sec.lt.60) goto 240
        sec=sec-60
        nin=nin+1
        goto 220
  230   sec=sec+60
        nin=nin-1
        goto 210
  240 continue
c       write(12,260) iyrmo(i),iday(i),ihr(i),nin,sec,lat,xlat,lon,
c    * xlon,evc(3,i),rmag(i)
c260  format(a4,a2,1x,a2,i2,1x,f5.2,i3,1x,f5.2,1x,i3,1x,f5.2,2f7.2)
c
c  use Fred Klein's HYPOINVERSE format from his routine HYSUM
C--CALLED BY HYPOINVERSE TO OUTPUT SUMMARY DATA
c  set latitude to be north and longitude to be west
      is='n'
      if(lat.lt.0) then
        is='s'
        lat=-1*lat
        xlat=-1.0*xlat
      endif
      ie='w'
      if(lon.lt.0) then
        ie='e'
        lon=-1*lon
        xlon=-1.0*xlon
      endif
c  use hypoinverse format
c  NOTE: I've left all the variables in the HYPOINVERSE output
c  for possible later use.
c  However those that are not calculated will be zero.
      ih71s=0
      iunit=12
c  Find minimum source-receiver distance
      do 5001 j=1,kobs(i)
      if(dlta(j,i).lt.dmin) dmin=dlta(j,i)
 5001 continue
C--CONVERT SOME DATA TO INTEGER FOR OUTPUT
	KLTM=xlat*100.+.5
	KLNM=xlon*100.+.5
	KQ=sec*100.+.5
	KZ=evc(3,i)*100.+.5
	LFMAG=rmag(i)*10.+.5
	LXMAG=rmag(i)*10.+.5
	KDMIN=DMIN+.5
	IF (KDMIN.GT.999) KDMIN=999
	KRMS=wrmsr(i)*100.+.5
	KERH=ERH*100.+.5
	KERZ=ERZ*100.+.5
c  nfm not used
c	NFM=NFRM
c	IF (NFM.GT.99) NFM=99
c use rmag(i) as input from for004 for both mxmag and mfmag
        mxmag=nint(10.0*rmag(i))
        mfmag=nint(10.0*rmag(i))
	IXMAG=.1*MXMAG+.5
	IF (IXMAG.GT.999) IXMAG=999
	IFMAG=.1*MFMAG+.5
	IF (IFMAG.GT.999) IFMAG=999
	DO 10 k=1,3
10	KSIG(k)=SERR(k)*100.
C--WRITE A SUMMARY RECORD
C--HYPO71 FORMAT
c  maxgap is not used although it could be calculated from az 
c  array (in AZTOA and LOCEQK)
	IF (IH71S.EQ.2 .AND. IUNIT.EQ.12) THEN
	  KSFL=' '
c	  IF (NWS.GT.0) KSFL='S'
	  WRITE (12,1002) iyrmo(i),iday(i),ihr(i),nin,sec,LAT,IS,
c	2 xlat,LON,IE,xlon,evc(3,i),rmag(i),nwr,MAXGAP,DMIN,
     2 xlat,LON,IE,xlon,evc(3,i),rmag(i),nwr,DMIN,
     3 wrmsr(i),erh,erz,ksfl
 1002 FORMAT (a4,a2,1x,a2,I2,F6.2,I3,A1,F5.2,I4,A1,F5.2,2F7.2,
c     2 I3,I4,F5.1,F5.2,2F5.1,1X,A1)
     2 I3,4x,F5.1,F5.2,2F5.1,1X,A1)
      ELSE

C--HYPOINVERSE FORMAT
        WRITE (IUNIT,1001) iyrmo(i),iday(i),ihr(i),nin, KQ,LAT,IS,
c    2 KLTM,LON,IE,KLNM,KZ, LXMAG,nwr,MAXGAP, KDMIN,KRMS,
     2 KLTM,LON,IE,KLNM,KZ, LXMAG,nwr,KDMIN,KRMS,
c    3 (IAZ(k),IDIP(k),KSIG(k),k=1,2), LFMAG,REMK,KSIG(3),
     3 (IAZ(k),IDIP(k),KSIG(k),k=1,2), LFMAG,KSIG(3),
c    4 RMK1,RMK2,NWS, KERH,KERZ,NFM, IXMAG,IFMAG,KXMMAD,KFMMAD,
     4 KERH,KERZ, IXMAG,IFMAG
1001   FORMAT (a4,a2,a2,i2, I4,I2,A1,
c     2 I4,I3,A1,I4,I5, I2,3I3,I4,
     2 I4,I3,A1,I4,I5, I2,i3,3x,i3,I4,
c    3 2(I3,I2,I4), I2,A3,I4,
     3 2(I3,I2,I4), I2,3x,I4,
c    4 2A1,I2, 2I4,I2, 4I3,
     4 2x,2x,2I4,2x, 2i3)
c    5 A3,A1, 3A1, I1,I3)
      END IF
c***** end of subroutine out *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine outadj(nit,istop,istop1)
c
c  common block variables:
      include 'simul2014_common.inc'
c
c  declaration statements:
      real fillin(maxnx)
      real x(15000)
      parameter(zero=0.0,izero=0)
      character*5 vtype(3)
c
      vtype(1)='P-Vel'
      vtype(2)='Vp/Vs'
      vtype(3)='Q    '
c
      if(nit.eq.0) goto 45
      write(16,1000) ssqr,nobtt,nobtp,nobts,nobteq,nobtex,wtsht,
     2  nobtte,wttel,nwrt,nswrt,
     3  nrdtsh,nordt,nordtp,nordts,wtrdtsht
       write(36,1000) ssqr,nobtt,nobtp,nobts,nobteq,nobtex,wtsht,
     2  nobtte,wttel,nwrt,nswrt,
     3  nrdtsh,nordt,nordtp,nordts,wtrdtsht
 1000 format(/,' sum of squared residuals =',f12.3,
     2 '; total number of observations =',i9,
     3 ', P obs =',i6,', S obs =',i6,
     4 /,20x,'earthquake obs.=',i9,', explosion obs.=',i9,
     5 ' (wtsht=',f5.2,'), tele obs=',i9,' (wttel=',f5.2,')',
     6 /,20x,'number obs with non-zero wt:  (eq-parsep) nwrt=',i9,
     7 ',   (sh-medder) nswrt=',i9,
     8 /,8x,'nev with rec-pair dt (nrdtsh)=',i6,
     9 /,14x,'receiver-pair dt: nobs=',i9,', P obs=',i9,', S obs=',i9,
     * ', wtrdtsht=', f5.2)
      if(iuseclu.gt.0) then
        write(16,1020) nclu,nobtce,noedt,noedtp,noedts,wtepdt
        write(36,1020) nclu,nobtce,noedt,noedtp,noedts,wtepdt
 1020   format(8x,'For Cluster: nclu=',i4,',  clus_event tt obs=',i9,
     &  /,14x,'earthquake-pair dt: nobs=',i9,', P obs=',i9,', S obs=',
     &  i9,', wtepdt=',f5.2)
      endif
      rmsres=sqrt(ssqr/float(nobtt))
      rmssep=0.0
      do 88 ne=1,nevt
         rmssep=rmssep+seprms(ne)
   88 continue
      rmssep=sqrt(rmssep/float(nobtt))
      write(16,1999) rmsres,rmssep
 1999 format(' rms residual =',f8.5,'   model rms =',f8.5)
      rmsw=sqrt(ssqrw/wnobt)
      write(16,1600) ssqrw,wnobt,rmsw
 1600 format(/,' weighted sum of squared residuals =',f12.3,
     * '; weighted total number of obs =',f10.1,
     * '; weighted rms res = ',f10.5)
c  store ssqr and number of velocity parameters for f-test
c  unless last step was backup
      if(istop1.ge.2) goto  5
      ssqrw1=ssqrw
      ssqr1=ssqr
      mbl1=mbl
      var1=var
      varw1=varw
      ndof1=ndof
      wndof1=wndof
    5 if(invdel.eq.0)go to 123
c=================================break================================
c  Output station corrections
      write(16,111)
  111 format(/,6x,'STATION CORRECTIONS AND ADJUSTMENTS',
     2 /,' station  Pdelay Adjusted S-Pdelay Adjusted')
      do 40 n=1,nsts
         nv=n+nparv
         nv1=nv+nsts
         if(hit(nv).lt.hitct.and.hit(nv1).lt.hitct) go to 40
         if(hit(nv).lt.hitct) vadj(nv)=0.
         if(hit(nv1).lt.hitct) vadj(nv1)=0.
c        write(16,120)stn(n),vadj(nv),vadj(nv1)
        if(kttfor.lt.3) then
          write(16,130) stn(n),pdl(n),vadj(nv),sdl(n),vadj(nv1)
        else
          if(kttfor.eq.3) then
            write(16,131) stn6(n),pdl(n),vadj(nv),sdl(n),vadj(nv1)
          else
            write(16,132)stn5(n),net(n),pdl(n),vadj(nv),sdl(n),vadj(nv1)
          endif
        endif
  130    format(1x,a4,2(f10.7,f10.7))
  131    format(1x,a6,2(f10.7,f9.7))
  132    format(1x,a5,1x,a2,2(f12.7,f11.7))
  120    format(3(1x,a4,1x,2f10.3))
   40 continue
c     write(16,110)
  110 format(/,' station corrections')
c     write(16,120)(stn(i),pdl(i),sdl(i),i=1,nsts)
 123  continue
c
      if(nparvi.eq.0) goto 999
c
      do 919 i=1,nx
         fillin(i)=0.0
  919 continue
      nx1=nx-1
      ny1=ny-1
      nz1=nz-1
c
   45 do 25 kv=1,iuses
         ivtype=kv
         if(iuseq.eq.1) ivtype=3
         if((kv.eq.1).and.(iusep.eq.0)) goto 25
c
c        write(16,1688) nit
c1688    format(' outadj, nit=',i5)
         if(nit.eq.0) goto 24
         if(nit.gt.1) goto 488
         write(16,1633)
 1633    format(/,' DERIVATIVE WEIGHT SUM ')
         do 486 k=2,nz1
            k2=k + (kv-1)*nz
c          if(iuseq.eq.0) then
c              write(16,1009) k,vtype(kv),zn(k)
c          else
c              write(16,1009) k,vtype(3),zn(k)
c          endif
            write(16,1009) k,vtype(ivtype),zn(k)
            kk=k+(kv-1)*nz2
            do 485 j=2,ny1
               n1=(kk-2)*nxy2+(j-2)*nx2+1
               n2=n1+nx2-1
               write(16,1635)(hit(i),i=n1,n2),zero
 1635          format(' 0.',20f7.0)
               do 484 i=n1,n2
                  sumhit(k2)=sumhit(k2)+hit(i)
  484          continue
  485       continue
  486    continue
  488    continue
c
         write(16,1001)
 1001    format(/,'  velocity model changes')
         do 12 k=2,nz1
            k2=k + (kv-1)*nz
            if(sumhit(k2).eq.0.0) goto 12
c          if(iuseq.eq.0) then
c              write(16,1009) k,vtype(kv),zn(k)
c          else
c              write(16,1009) k,vtype(3),zn(k)
c          endif
            write(16,1009) k,vtype(ivtype),zn(k)
            write(16,6002) (fillin(i),i=1,nx)
            kk=k + (kv-1)*nz2
            do 10 j=2,ny1
               j1=(kk-2)*nxy2+(j-2)*nx2+1
               j2=j1+nx2-1
               if(iuseq.eq.0) then
                 if(kv.eq.2) then
                   write(16,1004) (vadj(i),i=j1,j2),fillin(1)
                 else
                   write(16,1003) (vadj(i),i=j1,j2),fillin(1)
                 endif
               else
                 write(16,1003) (qadj(i),i=j1,j2),fillin(1)
               endif
   10       continue
            write(16,6002) (fillin(i),i=1,nx)
   12    continue
 6002    format(20f7.2)
 6003    format(20f7.1)
 1003    format('   0.00',20f7.2)
 1004    format('  0.000',20f7.3)
         if(iuseq.eq.0) then
           write(16,2001)
         else
           write(16,2002)
         endif
 2001    format(/,' corrected velocity model',/)
 2002    format(/,' corrected Q model',/)
         do 23 k=2,nz1
            k2=k + (kv-1)*nz
            if(sumhit(k2).eq.0.0) goto 23
c          if(iuseq.eq.0) then
c              write(16,1009) k,vtype(kv),zn(k)
c          else
c              write(16,1009) k,vtype(3),zn(k)
c          endif
            write(16,1009) k,vtype(ivtype),zn(k)
 1009       format(/,' layer',i3,5x,a5,' nodes',10x,'z =',f7.1)
            do 22 j=1,ny
               if(kv.eq.1) then
               if(iuseq.eq.0) then
                 write(16,6002) (vel(i,j,k2),i=1,nx)
               else
                 write(16,6003) (qval(i,j,k),i=1,nx)
               endif
             endif
               if(kv.eq.2) write(16,6002) (vpvs(i,j,k),i=1,nx)
  22        continue
  23     continue
c  Find average and variance (section taken from COMPVEL.FOR)
   24    write(16,3030)
 3030    format(/,'  MODEL STATISTICS, by Depth, excluding edge nodes')
         if((kv.eq.1).and.(nit.eq.0)) write(36,3034) nit
 3034    format(/,'  MODEL STATISTICS on nit=',i2)
         iv=0
         varv=0.0
         do 200 k=2,nz1
c          if(iuseq.eq.0) then
c              write(16,3015) zn(k),vtype(kv)
c          else
c              write(16,3015) zn(k),vtype(3)
c          endif
            write(16,3015) zn(k),vtype(ivtype)
            k2=k+(kv-1)*nz
            ix=0
            do 190 j=2,ny1
               do 180 i=2,nx1
                  ix=ix+1
                  if(kv.eq.1) then
                   if(iuseq.eq.0) then
                       x(ix)=vel(i,j,k2)
                   else
                       x(ix)=qval(i,j,k)
                   endif
                  else
                     x(ix)=vpvs(i,j,k)
                  endif
  180          continue
  190       continue
            call avsd(0.00,x,ix,sd,av,devtot)
            iv=iv+ix 
            varv=varv+devtot
            write(16,3017) ix,av,sd
  200    continue
         varv=varv/iv
         sdv=sqrt(varv)
       if(iuseq.eq.0) then
           write(16,3031)
       else
           write(16,3032)
       endif
 3031    format(/,' VARIANCE OF THIS VELOCITY MODEL')
 3032    format(/,' VARIANCE OF THIS Q MODEL')
         write(36,3018) iv,vtype(ivtype),varv,sdv
         write(16,3018) iv,vtype(ivtype),varv,sdv
 3015    format(/,' grid z=',f8.2,5x,a5)
 3017    format('  For ', i5,' gridpts, av=',f7.2,', sd=',f7.3)
 3018    format('  For all',i7,1x,a5,' gridpts, variance=',f14.7,
     2    ', sd=',f10.3)
c
   25 continue
  999 continue
c***** end of subroutine outadj *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine outend
c  Routine for final output.  See description of output files
c  and "kout" output control parameter at the beginning of
c  the program.
c
c  common block variables:
      include 'simul2014_common.inc'
      common/machin/ eta,tol
      common/shortd/ xltkm,xlnkm,rota,nzco,xlt,xln,snr,csr,
     2 onorth,oeast,rotadg,cmerid,rlatdo,rlondo
      common/upb/ upbtot,nupb
      common/nerr/errnml(1000),kerr,rderrq(5)
c
c  declaration statements:
      integer lat,lon
      real xlat,xlon
      real zeroi(maxnx)
      real velav(maxnz*2)           !for average layer velocities
c      real type(3),ttime(maxobs),sekms(maxnx)
      real ttime(maxobs),sekms(maxnx)
      character*1 type(3)
      data type(1),type(2),type(3) /1hP,1hS,1hQ/
      parameter(zero=0.0,izero=0)
c
c  cht 1998:
      real resout(maxnx,maxny,maxnz*2)
c
      character*50 formvl,formrs,line,frmrs1,frmrs2,frmvl1,frmvl2
      character*2 ch2nx2
      character*1 ch2nx3,ch2nx4,cns,cew,cn(maxev),ce(maxev)
      character*5 vtype(3)
      character*24 date24
c
      vtype(1)='P-Vel'
      vtype(2)='Vp/Vs'
      vtype(3)='Q'
      call fdate(date24)
      formrs='(1x,3hres,5x,10(1h(,f5.2,1h),5x),3hres)'
      frmrs1='(1x,3hres,5x,10(1h(,f5.2,1h),5x))'
cfh format problems
c      frmrs2='(9x,2(12x),8(1h(,f5.2,1h),5x),3hres)'
      frmrs2='(9x,02(12x),08(1h(,f5.2,1h),5x),3hres)'
      if(iuseq.eq.0) then
        formvl='(1h0,f4.2,10(f6.2,1h+,f5.4),f6.2)'
        frmvl1='(1h0,f4.2,10(f6.2,1h+,f5.4))'
cfh format problems
c        frmvl2='(1h0,4x,2(12x),8(f6.2,1h+,f5.4),f6.2)'
      frmvl2='(1h0,4x,02(12x),08(f6.2,1h+,f5.4),f6.2)'
      else
        formvl='(    f5.0,10(f6.0,1h+,f5.1),f6.0)'
        frmvl1='(    f5.0,10(f6.0,1h+,f5.0))'
        frmvl2='(1h0,4x,2(12x),8(f6.0,1h+,f5.0),f6.0)'
      endif
      do 2 i=1,maxnx
         zeroi(i)=zero
    2 continue
c
      if(nitmax.gt.-1) goto 3
c get normal distribution error if synthetic (not for Q)
      iseed=123456789
      do i=1,1000
        errnml(i)=rnormal(iseed)
      enddo
c kerr is counter
      kerr=1
c assign rderr according to reading weight
c this gave set of terrsy with sd=rderr for test data
      rderrq(1)=rderr*0.3
      rderrq(2)=rderr*0.7
      rderrq(3)=rderr*1.2
      rderrq(4)=rderr*2.2
      rderrq(5)=5.0*rderr
    3 write(16,100) ssqr,nobtt,nobtp,nobts,nobteq,nobtex,wtsht,
     2  nobtte,wttel,nordt,nordtp,nordts,wtrdtsht,
     3 nwrt,nswrt
  100 format(/,' sum of squared residuals =',f10.3,
     2 '; total number of observations =',i9,
     3 ', P obs =',i8,', S obs =',i8,
     4 /,53x,'earthquake obs.=',i9,', explosion obs.=',i9,
     5 ' (wtsht=',f5.2,'), tele obs=',i9,' (wttel=',f5.2,')',
     6 /,53x,'receiver-pair dt: nobs=',i9,', P obs=',i9,', S obs=',i9,
     8 ', wtrdtsht=', f5.2,
     9 /,20x,'with non-zero wt:',t65,'nwrt=',i9,',',10x,
     * 'nswrt=',i9)
c** change for creating synthetic data, nitmax= -1
      if(nitmax.eq.-1) goto 76
      rmsres=sqrt(ssqr/float(nobtt))
      write(16,1999) rmsres
 1999 format(' rms residual =',f8.5)
      rmsw=sqrt(ssqrw/wnobt)
      write(16,1600) ssqrw,wnobt,rmsw
 1600 format(/,' weighted sum of squared residuals =',f12.3,
     * '; weighted total number of obs =',f10.1,
     * '; weighted rms res = ',f10.5)
c
c  Write new hypocenters to output(file16)
   76 open(unit=33,file='hyp_prt.out')
      if(iuse2t.eq.0) then
        write(16,1650)
        write(33,1650)
      else
        write(16,1651)
        write(33,1651)
      endif
 1650 format(//,15x,'FINAL LOCATIONS',//,' EVT YRMODY HRMN  SEC',
     2 2x,'LATITUDE LONGITUDE  DEPTH     MAG  NO RMSRES',4x,'x',6x,
     3 'y',6x,'z',4x,'GAP DMIN RZDM  NP NS')
 1651 format(//,15x,'FINAL LOCATIONS',//,' EVT YRMODY HRMN  SEC',
     2 2x,'LATITUDE LONGITUDE  DEPTH     MAG  NO RMSRES',4x,'x',6x,
     3 'y',6x,'z',4x,'GAP DMIN RZDM  NP NS   SEC2    DSEC')
c
c  cht 1998:
      open(32,file='hypo.gmt')
c
c  Map 0,1,2,3,4,5 to -1,0,1,2,2,3
      kkout=nint(kout*0.85)-1
c      if(kkout) 301,295,290
c  Write new travel time file24 to :q be used as input in future runs
      if(kkout.lt.0) goto 301
      if(kkout.eq.0) goto 295
  290 if(neqs.gt.0) open(unit=24,file='tteq.out')
      if(nbls.gt.0) open(unit=28,file='ttbls.out')
      if((nitmax.ge.0).and.(nsht.gt.0)) open(unit=27,file='ttsht.out')
      if(nitmax.eq.-1) open(unit=27,file='ttsyn.out')
      if(ntel.gt.0) then
        open(unit=40,file='tttel.out')
        write(40,4008)
 4008   format(50x,'Tele hyps and ttobs at base of inversion layers',/)
      endif
      if(ntel.gt.0) open(unit=41,file='telpad.out')
c  Write new hypocenters to file13 in hypo71 format
  295 open(unit=13,file='finalsmpout')
      rewind(13)
      write(13,1301) neqs,date24(21:24),date24(4:16)
 1301 format(i5,75x,'computed',2x,a4,a13)
  301 neqbl=neqs+nbls
      iunit=24
       rmsall = 0.0
c      do 5 i=1,nevt
c write out tele later
      do 5 i=1,nebs
         if(i.gt.neqs) iunit=28
         if(i.gt.neqbl) iunit=27
c  convert origin time into gmt
         sec= seco(i)
         nin=mino(i)
 210     if(sec.lt.0) goto 230
 220     if(sec.lt.60) goto 240
         sec=sec-60
         nin=nin+1
         goto 220
 230     sec=sec+60
         nin=nin-1
         goto 210
  240    continue
         mina(i)=nin
         seca(i)=sec
         sec2=0.0
         dsec=0.0
         if((iuse2t.eq.1).and.(seco2(i).ne.0.0)) then
           dmin=mino(i)-nin
           sec2=seco2(i)+dmin*60.0
           dsec=seco(i)-seco2(i)
         endif
         call latlon(evc(1,i),evc(2,i),lat,xlat,lon,xlon)
cfhek
              if(((lat.ge.0).and.(xlat.ge.0.)).or.(nzco.eq.2)) then
                   cns='N'
                   latdeg=iabs(lat)
                   alatmin=abs(xlat)
              else
                   cns='S'
                   latdeg=iabs(lat)
                   alatmin=abs(xlat)
              endif
              if(((lon.ge.0).and.(xlon.ge.0.)).or.(nzco.eq.2)) then
                   cew='W'
                   londeg=iabs(lon)
                   alonmin=abs(xlon)
              else
                   cew='E'
                   londeg=iabs(lon)
                   alonmin=abs(xlon)
              endif
c If nzco.eq.1, assume S and E, even if left blank in the input
         if(nzco.eq.1) cns='S'
         if(nzco.eq.1) cew='E'
         ltde(i)=latdeg
         eltm(i)=alatmin
         cn(i)=cns
         lnde(i)=londeg
         elnm(i)=alonmin
         ce(i)=cew
c
c ratio depth to deltmin
      if(idltmn(i).gt.0) rzdm=evc(3,i)/float(idltmn(i))
c
c  start cht 1998
c  output gmt-format epicenters
c
      xla=lat*1.0+xlat/60.
      xlo=-1.0*lon-xlon/60.
      write(32,8765) xlo,xla,evc(3,i)
 8765 format(f12.5,1x,f12.5,1x,f7.2)
c
c  end cht 1998
c
c         if(kkout) 304,303,302
         if(kkout.lt.0) goto 304
         if(kkout.eq.0) goto 303
  302    continue
c
         if(i.eq.1) then
           write(iunit,2404) iyrmo(i),iday(i),ihr(i),
     2      nin,sec,latdeg,cns,alatmin,londeg,cew,alonmin,evc(3,i),
     3       rmag(i),igap(i),idltmn(i),sec2,date24(21:24),date24(4:16)
         else
           if(iuse2t.eq.0) then
             write(iunit,2401) iyrmo(i),iday(i),ihr(i),
     2       nin,sec,latdeg,cns,alatmin,londeg,cew,alonmin,evc(3,i),
     3       rmag(i),igap(i),idltmn(i)
           else
             write(iunit,2401) iyrmo(i),iday(i),ihr(i),
     2       nin,sec,latdeg,cns,alatmin,londeg,cew,alonmin,evc(3,i),
     3       rmag(i),igap(i),idltmn(i),sec2
           endif
 2401      format(a4,a2,1x,a2,i2,1x,f5.2,i3,a1,f5.2,1x,i3,a1,f5.2,
     2     2f7.2,2i4,t69,f6.2)
 2404      format(a4,a2,1x,a2,i2,1x,f5.2,i3,a1,f5.2,1x,i3,a1,f5.2,
     2     2f7.2,2i4,t69,f6.2,10x,'computed',2x,a4,a13)
         end if
c  write travel time observations for each event
         do 30 jobs=1,kobs(i)
            ttime(jobs)=secp(jobs,i)
           if((iuse2t.eq.0).or.(iclock(isto(jobs,i)).eq.0)) then
             if(intsp(jobs,i).eq.0) ttime(jobs)=secp(jobs,i)-seco(i)
           else
             if(intsp(jobs,i).eq.0) ttime(jobs)=secp(jobs,i)-seco2(i)
           endif
   30    continue
c
c add error to synthetic travel times (not t*)
        if((nitmax.eq.-1).and.(iuseq.eq.0)) then
        IF(I.EQ.1) OPEN(78,file='temp.synerr')
      WRITE(78,7801) i
 7801 FORMAT(' synthetic error for event ',i5,' iqual, intsp, tterr')
 7802 FORMAT(i4,i3,f7.3)
        do 32 jobs=1,kobs(i)
          irq=iw(jobs,i)+1
          terrsy= rderrq(irq)*errnml(kerr)
          if(intsp(jobs,i).eq.1) then
            kerr=kerr+1
            terrsy=terrsy-rderrq(1)*errnml(kerr)
          endif
      WRITE(78,7802) iw(jobs,i),intsp(jobs,i),terrsy
          ttime(jobs)=ttime(jobs)+terrsy
          kerr=kerr+1
          if(kerr.gt.1000) kerr=1
   32   continue
        endif
c  sometimes have commented this out to get old cnsp from simphs input data
ctemp      if(kttfor.eq.2) then
ctemp        do 35 j=1,kobs(i)
ctemp   35   write(iunit,4992) stn(isto(j,i)),rmk(j,i),ttime(j)
ctemp        write(iunit,2406)
ctemp        goto 303
ctemp      endif
 4992 format(a4,1x,a4,1x,f9.4)
      if(kttfor.eq.3) then
        do 37 j=1,kobs(i)
   37   write(iunit,4994) stn6(isto(j,i)),rmk(j,i),ttime(j)
        write(iunit,2406)
        goto 303
      endif
 4994 format(a6,1x,a4,1x,f9.4)
      if(kttfor.eq.4) then
        do 137 j=1,kobs(i)
  137   write(iunit,4990)stn5(isto(j,i)),net(isto(j,i)),
     2  rmk(j,i),ttime(j)
        write(iunit,2406)
        goto 303
      endif
 4990 format(a5,5x,a2,1x,a4,1x,f9.4)
       if(iuseq.eq.0) then
         write(iunit,2402) (stn(isto(j,i)),rmk(j,i),
     2      ttime(j),j=1,kobs(i))
         write(iunit,2403)
       else
         write(iunit,2405) (stn(isto(j,i)),rmk(j,i),
     2        ttime(j),j=1,kobs(i))
         write(iunit,2406)
       endif
 303  continue
      if(iuse2t.eq.0) then
         write(13,1302) iyrmo(i),iday(i),ihr(i),nin,sec,
     2    latdeg,cns,alatmin,londeg,cew,alonmin,evc(3,i),
     3    rmag(i),kobs(i),wrmsr(i),igap(i),idltmn(i)
      else
         write(13,1302) iyrmo(i),iday(i),ihr(i),nin,sec,
     2    latdeg,cns,alatmin,londeg,cew,alonmin,evc(3,i),
     3    rmag(i),kobs(i),wrmsr(i),igap(i),idltmn(i),sec2
      endif
c
cek Now write it on file 16
  304 continue
      if(iuse2t.eq.0) then
           write(16,1602) i,iyrmo(i),iday(i),ihr(i),nin,sec,
     2      latdeg,cns,alatmin,londeg,cew,alonmin,evc(3,i),rmag(i),
     3      kobs(i),wrmsr(i),(evc(j,i),j=1,3),igap(i),idltmn(i),rzdm,
     4      kobps(1,i),kobps(2,i)
           write(33,1602) i,iyrmo(i),iday(i),ihr(i),nin,sec,
     2      latdeg,cns,alatmin,londeg,cew,alonmin,evc(3,i),rmag(i),
     3      kobs(i),wrmsr(i),(evc(j,i),j=1,3),igap(i),idltmn(i),rzdm,
     4      kobps(1,i),kobps(2,i)
      else
           write(16,1612) i,iyrmo(i),iday(i),ihr(i),nin,sec,
     2      latdeg,cns,alatmin,londeg,cew,alonmin,evc(3,i),rmag(i),
     3      kobs(i),wrmsr(i),(evc(j,i),j=1,3),igap(i),idltmn(i),rzdm,
     4      kobps(1,i),kobps(2,i),sec2,dsec
           write(33,1612) i,iyrmo(i),iday(i),ihr(i),nin,sec,
     2      latdeg,cns,alatmin,londeg,cew,alonmin,evc(3,i),rmag(i),
     3      kobs(i),wrmsr(i),(evc(j,i),j=1,3),igap(i),idltmn(i),rzdm,
     4      kobps(1,i),kobps(2,i),sec2,dsec
      endif
cek  304    write(16,1602) i,iyrmo(i),iday(i),ihr(i),nin,sec,
cek     2      lat,xlat,lon,xlon,evc(3,i),rmag(i),kobs(i),wrmsr(i),
cek     3      (evc(j,i),j=1,3)
c
      rmsall = rmsall+wrmsr(i)
    5 continue
      rmsall = rmsall/float(nevt)
      write(16,2407) rmsall
 2407 format('RMSALL:  ',f8.6)
 2402 format(6(a4,a4,f6.2))
 2403 format(4h    )
 2405 format(6(a4,a4,f10.4))
 2406 format(6h     0)
 2408 format(6(a6,a4,f6.2))
 2409 format(6(a6,a4,f10.4))
 1302 format(a4,a2,1x,a2,i2,f6.2,1x,i2,a1,f5.2,1x,i3,a1,f5.2,f7.2,
     2 2x,f5.2,i3,9x,f5.2,2i4,t79,f6.2)
cekfh
cek 1602 format(1x,i3,1x,a4,a2,1x,a2,i2,f6.2,1x,i2,'n',f5.2,1x,i3,'w',
cek     2 f5.2,f7.2,2x,f5.2,i4,f6.2,3f7.2)
 1602 format(i4,1x,a4,a2,1x,a2,i2,f6.2,i3,a1,f6.2,i4,a1,
     2 f6.2,f7.2,2x,f5.2,i4,f6.2,3f7.2,2i5,f5.1,i4,i3)
 1612 format(i4,1x,a4,a2,1x,a2,i2,f6.2,i3,a1,f6.2,i4,a1,
     2 f6.2,f7.2,2x,f5.2,i4,f6.2,3f7.2,2i5,f5.1,i4,i3,2f7.2)
c
c    8 if(kkout) 310,310,305
    8 if(kkout.le.0) goto 310
  305 close(unit=24)
  310 continue
c
c output receiver-pair shots
      if(nrdtsh.eq.0) goto 390
      write(16,1652) 
 1652 format(//,15x,'RECEIVER-PAIR SHOTS FINAL LOCATIONS',/,
     1 t64,'RDT',/,' RSH YRMODY HRMN  SEC',
     2 2x,'LATITUDE LONGITUDE  DEPTH     MAG  NO RMSRES',4x,'x',6x,
     3 'y',6x,'z',4x,'GAP DMIN')
      do 370 i=1,nevt
        if(kobsrdt(i).eq.0) goto 370
        write(16,1602) i,iyrmo(i),iday(i),ihr(i),mina(i),seca(i),
     2   ltde(i),cn(i),eltm(i),lnde(i),ce(i),elnm(i),evc(3,i),rmag(i),
     3   kobsrdt(i),wrmsrdt(i),(evc(j,i),j=1,3),igap(i),idltmn(i)
        rmsall=rmsall+wrmsrdt(i)
  370 continue
      rmsall=rmsall/float(nrdtsh)
      write(16,2407) rmsall
  390 continue
c
c output cluster events
      if((iuseclu.eq.0).or.(nclu.eq.0)) goto 395
      write(16,1655)
 1655 format(//,15x,'CLUSTER_EVENT FINAL LOCATIONS', /,
     1 ' EVT YRMODY HRMN  SEC',
     2 2x,'LATITUDE LONGITUDE  DEPTH     MAG  NO RMSRES',4x,'x',6x,
     3 'y',6x,'z',4x,'GAP DMIN RZDM  NCLU RMSCLU')
       if(kout.ge.2) open(unit=29,file='ttcev.out')
       iunit=29
c
      do 157 nc=1,nclu
       if(kout.ge.2) write(iunit,2901) nc
 2901  format('BEGIN CLUSTER',10x,i6)
        do 155 i=1,ncev(nc) 
c  convert origin time into gmt
         sec= secoce(i,nc)
         nin=minoce(i,nc)
 1210    if(sec.lt.0) goto 1230
 1220    if(sec.lt.60) goto 1240
         sec=sec-60
         nin=nin+1
         goto 1220
 1230    sec=sec+60
         nin=nin-1
         goto 1210
 1240    continue
         minace(i,nc)=nin
         secace(i,nc)=sec
         sec2=0.0
         dsec=0.0
         call latlon(cevc(1,i,nc),cevc(2,i,nc),lat,xlat,lon,xlon)
              if(((lat.ge.0).and.(xlat.ge.0.)).or.(nzco.eq.2)) then
                   cns='N'
                   latdeg=iabs(lat)
                   alatmin=abs(xlat)
              else
                   cns='S'
                   latdeg=iabs(lat)
                   alatmin=abs(xlat)
              endif
              if(((lon.ge.0).and.(xlon.ge.0.)).or.(nzco.eq.2)) then
                   cew='W'
                   londeg=iabs(lon)
                   alonmin=abs(xlon)
              else
                   cew='E'
                   londeg=iabs(lon)
                   alonmin=abs(xlon)
              endif
c If nzco.eq.1, assume S and E, even if left blank in the input
         if(nzco.eq.1) cns='S'
         if(nzco.eq.1) cew='E'
         ltdce(i,nc)=latdeg
         celtm(i,nc)=alatmin
         cn(i)=cns
         lndce(i,nc)=londeg
         celnm(i,nc)=alonmin
         ce(i)=cew
c
c ratio depth to deltmin
      if(idltmnce(i,nc).gt.0) rzdm=cevc(3,i,nc)/float(idltmnce(i,nc))
c
c  output gmt-format epicenters
      xla=lat*1.0+xlat/60.
      xlo=-1.0*lon-xlon/60.
      write(32,8765) xlo,xla,cevc(3,i,nc)
c
      if(kout.lt.2) goto 151
         if(i.eq.1) then
           write(iunit,2404) iyrmoce(i,nc),idayce(i,nc),ihrce(i,nc),
     2      nin,sec,latdeg,cns,alatmin,londeg,cew,alonmin,cevc(3,i,nc),
     3       rmagce(i,nc),igapce(i,nc),idltmnce(i,nc),sec2,
     4       date24(21:24),date24(4:16)
         else
           write(iunit,2401) iyrmoce(i,nc),idayce(i,nc),ihrce(i,nc),
     2     nin,sec,latdeg,cns,alatmin,londeg,cew,alonmin,cevc(3,i,nc),
     3     rmagce(i,nc),igapce(i,nc),idltmnce(i,nc)
         endif
c  write travel time observations for each event
         do 150 jobs=1,kobsce(i,nc)
           ttime(jobs)=secpc(jobs,i,nc)
           if(intspc(jobs,i,nc).eq.0) ttime(jobs)=secpc(jobs,i,nc)
     &       -secoce(i,nc)
  150    continue
c add error to synthetic travel times (not t*)
         if((nitmax.eq.-1).and.(iuseq.eq.0)) then
           do 153 jobs=1,kobsce(i,nc)
             irq=iwc(jobs,i,nc)+1
             terrsy= rderrq(irq)*errnml(kerr)
             if(intspc(jobs,i,nc).eq.1) then
               kerr=kerr+1
               terrsy=terrsy-rderrq(1)*errnml(kerr)
             endif
             ttime(jobs)=ttime(jobs)+terrsy
             kerr=kerr+1
             if(kerr.gt.1000) kerr=1
  153     continue
        endif
c
      if(kttfor.eq.1) then
        if(iuseq.eq.0) then
         write(iunit,2402) (stn(istoc(j,i,nc)),rmkc(j,i,nc),
     2     ttime(j),j=1,kobsce(i,nc))
         write(iunit,2403)
        else
         write(iunit,2405) (stn(istoc(j,i,nc)),rmkc(j,i,nc),
     2     ttime(j),j=1,kobsce(i,nc))
         write(iunit,2406)
        endif
        goto 158
      endif
      if(kttfor.eq.2) then
         write(iunit,4992) (stn(istoc(j,i,nc)),rmkc(j,i,nc),
     2     ttime(j),j=1,kobsce(i,nc))
         goto 156
      endif
      if(kttfor.eq.3) then
         write(iunit,4994) (stn6(istoc(j,i,nc)),rmkc(j,i,nc),
     2     ttime(j),j=1,kobsce(i,nc))
         goto 156
      endif
      if(kttfor.eq.4) then
         write(iunit,4990) (stn5(istoc(j,i,nc)),net(istoc(j,i,nc)),
     2     rmkc(j,i,nc),ttime(j),j=1,kobsce(i,nc))
         goto 156
      endif
  156  write(iunit,2406)
  158  if(i.eq.ncev(nc)) write(iunit,2902) nc
 2902  format('END CLUSTER',10x,i6)
  151  continue
        write(13,1303) iyrmoce(i,nc),idayce(i,nc),ihrce(i,nc),nin,sec,
     2    latdeg,cns,alatmin,londeg,cew,alonmin,cevc(3,i,nc),
     3    rmagce(i,nc),kobsce(i,nc),wrmsrce(i,nc),igapce(i,nc),
     4    idltmnce(i,nc),nc
         write(16,1603) i,iyrmoce(i,nc),idayce(i,nc),ihrce(i,nc),
     2      nin,sec,latdeg,cns,alatmin,londeg,cew,alonmin,
     3      cevc(3,i,nc),rmagce(i,nc),kobsce(i,nc),wrmsrce(i,nc),
     4      (cevc(j,i,nc),j=1,3),igapce(i,nc),idltmnce(i,nc),rzdm,
     5      nc,wrmsclu(nc)
         write(33,1603) i,iyrmoce(i,nc),idayce(i,nc),ihrce(i,nc),
     2      nin,sec,latdeg,cns,alatmin,londeg,cew,alonmin,
     3      cevc(3,i,nc),rmagce(i,nc),kobsce(i,nc),wrmsrce(i,nc),
     4      (cevc(j,i,nc),j=1,3),igapce(i,nc),idltmnce(i,nc),rzdm,
     5      nc,wrmsclu(nc)
  155   continue
  157 continue
      if(kout.ge.2) close(29)
 1303 format(a4,a2,1x,a2,i2,f6.2,1x,i2,a1,f5.2,1x,i3,a1,f5.2,f7.2,
     2 2x,f5.2,i3,9x,f5.2,2i4,t79,i6)
 1603 format(i4,1x,a4,a2,1x,a2,i2,f6.2,i3,a1,f6.2,i4,a1,
     2 f6.2,f7.2,2x,f5.2,i4,f6.2,3f7.2,2i5,f5.1,i6,f6.2)
  395 continue
c
      if(ntel.eq.0) goto 399
c output tele events, updated tel path delay(s);
c    centroid piercing points with tel path delays
      write(41,4041)
 4041 format('avg_pierce_pt  Lat  Long  Depth  P_tele_path_del S-P_pad',
     2  '   x      y       z')
      do 397 nt=1,ntel
        n=nebs+nt
        write(16,1642) n,iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
     2  ltde(n),ecns(n),eltm(n),lnde(n),ecew(n),elnm(n),evc(3,n),
     3  rmag(n),kobs(n),wrmsr(n),telpad(nt),telpads(nt)
 1642   format(i4,1x,a4,a2,1x,a2,i2,f6.2,i3,a1,f6.2,i4,a1,
     2   f6.2,f7.2,2x,f5.2,i4,f6.2,t121,2f7.2)
        write(40,4003) iyrmo(n),iday(n),ihr(n),
     2   mino(n),seco(n),elat(n),elon(n),
     3   evc(3,n),rmag(n),telpad(nt),telpads(nt)
c 4003 format('      ',/,a4,a2,1x,a2,i2,1x,f5.2,f9.4,f10.4,f7.2,f7.2,
 4003 format(a4,a2,1x,a2,i2,1x,f5.2,f9.4,f10.4,f7.2,f7.2,
     2  t57,2f9.2)
         call latlon(tecen(1,nt),tecen(2,nt),lat,xlat,lon,xlon)
              if(((lat.ge.0).and.(xlat.ge.0.)).or.(nzco.eq.2)) then
                   cns='N'
                   latdeg=iabs(lat)
                   alatmin=abs(xlat)
              else
                   cns='S'
                   latdeg=iabs(lat)
                   alatmin=abs(xlat)
              endif
              if(((lon.ge.0).and.(xlon.ge.0.)).or.(nzco.eq.2)) then
                   cew='W'
                   londeg=iabs(lon)
                   alonmin=abs(xlon)
              else
                   cew='E'
                   londeg=iabs(lon)
                   alonmin=abs(xlon)
              endif
        write(16,1643) n,iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
     2   latdeg,cns,alatmin,londeg,cew,alonmin,tecen(3,nt),rmag(n),
     3   kobs(n),wrmsr(n),(tecen(j,nt),j=1,3),telpad(nt),telpads(nt)
        write(33,1643) n,iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
     2   latdeg,cns,alatmin,londeg,cew,alonmin,tecen(3,nt),rmag(n),
     3   kobs(n),wrmsr(n),(tecen(j,nt),j=1,3),telpad(nt),telpads(nt)
        write(13,1302) iyrmo(n),iday(n),ihr(n),mino(n),seco(n),
     2    latdeg,cns,alatmin,londeg,cew,alonmin,tecen(3,nt),
     4    rmag(n),kobs(n),wrmsr(n)
 1643   format(i4,1x,a4,a2,1x,a2,i2,f6.2,i3,a1,f6.2,i4,a1,
     2   f6.2,f7.2,2x,f5.2,i4,f6.2,3f7.2,t121,2f7.2)
        write(41,4165) tecnlt(nt),tecnln(nt),tecen(3,nt),telpad(nt),
     2   telpads(nt),(tecen(i,nt),i=1,3)
 4165 format(f12.5,1x,f12.5,1x,f7.2,1x,2f9.2,2x,2f10.2,f7.2)
        do 396 j=1,kobs(n)
          write(40,4997) stn(isto(j,n)),rmk(j,n),secte(j,nt),
     2     teplat(j,nt),teplon(j,nt),tepp(3,j,nt)
 4997     format(a4,3x,a4,f10.4,2x,2f10.4,f8.2)
  396   continue
        write(40,2406)
  397 continue
      close(40)
      close(41)
c
c
  399 close(unit=33)
      close(unit=13)
c  output number of observations by station
      if(kttfor.lt.3) then
        write(16,1640)
 1640   format(//,15x,'TALLY OF OBSERVATIONS',//,1x,
     2   10('station obs  '),/,10(7x,'P  S-P'),/,
     3   10(' ==== === ==='))
        write(16,1645) (stn(i),nrd(i,1),nrd(i,2),i=1,nsts)
 1645   format(10(1x,a4,2i4))
      else
        if(kttfor.eq.3) then
          write(16,1641)
 1641     format(//,15x,'TALLY OF OBSERVATIONS',//,1x,
     2    10('station   obs  '),/,10(7x,'P  S-P'),/,
     3    10(' ====== === ==='))
          write(16,1644) (stn6(i),nrd(i,1),nrd(i,2),i=1,nsts)
 1644     format(10(1x,a6,2i4))
        else
          write(16,1646)
 1646     format(//,15x,'TALLY OF OBSERVATIONS',//,1x,
     2    8(' sta net obs '),/,8(10x,'P  S-P'),/,
     3    8(' ===== == === ==='))
          write(16,1647) (stn5(i),net(i),nrd(i,1),nrd(i,2),i=1,nsts)
 1647     format(8(1x,a5,1x,a2,2i4))
        endif
      endif
c
c
c** Change for synthetic data, nitmax= -1
   95 if(nitmax.le.0) goto 900
      if(invdel.eq.0)go to 96
c       output final station corrections
        write(16,110)
 110  format(//,' FINAL STATION CORRECTIONS',/,' stn  typ Delay  StErr',
     2 ' Resol   nobs  stn  stn grdpt  inv soln',/,29x,
     3         'nrd hit name num   num arry arry',/,51x,
     4                           '  num  num',/,
     5 '   _____________ _____ _____ ___ ___',5(' ____'))
        rat=ssqrw/ssqrw1
c  also write new station corrections to file22 to be used as input
c  in future runs.
   96 if(kout.ge.2) write(22,2201) nsts,date24(21:24),date24(4:16)
c 2201 format(i4,24x,'computed',2x,a4,a13)
 2201 format(i5,1x,'latitude longitude elev',
     2 ' pdl s-pdl fix iclock      x      y',
     3 8x,'Pnobs  S-Pnobs,',3x,'Computed ',a4,a13)
      do 120 i=1,nsts
         ielev= stc(3,i)*(-1.0e+03)
cfhek
              if(ltds(i).ge.0.and.sltm(i).ge.0.) then
                   cns='N'
                   latdeg=iabs(ltds(i))
                   alatmin=abs(sltm(i))
              else
                   cns='S'
                   latdeg=iabs(ltds(i))
                   alatmin=abs(sltm(i))
              endif
              if(lnds(i).ge.0.and.slnm(i).ge.0.) then
                   cew='W'
                   londeg=iabs(lnds(i))
                   alonmin=abs(slnm(i))
              else
                   cew='E'
                   londeg=iabs(lnds(i))
                   alonmin=abs(slnm(i))
              endif
       if(kttfor.lt.3) then
         if(iuseq.eq.0) then
           write(22,2207) stn(i),latdeg,cns,alatmin,londeg,cew,
     2      alonmin,ielev,pdl(i),sdl(i),nfixst(i),iclock(i),
     3      stc(1,i),stc(2,i),nrd(i,1),nrd(i,2)
         else
           write(22,2208) stn(i),latdeg,cns,alatmin,londeg,cew,
     2      alonmin,ielev,pdl(i),nfixst(i),
     3      stc(1,i),stc(2,i),nrd(i,1),nrd(i,2)
         endif
 2207    format(2x,a4,i2,a1,f5.2,i4,a1,f5.2,i5,2f5.2,2i3,4x,2f9.3,2i9)
 2208    format(2x,a4,i2,a1,f5.2,i4,a1,f5.2,i5,e11.3,i2,7x,2f9.3,2i9)
       else
         call ltlndm(rlat,rlon,latdeg,alatmin,cns,londeg,alonmin,
     2       cew,2)
         if(kttfor.eq.3) then
           elevm=float(ielev)
           write(22,2008) stn6(i),rlat,rlon,elevm,pdl(i),sdl(i),
     2       nfixst(i),stc(1,i),stc(2,i),nrd(i,1),nrd(i,2)
 2008       format(a6,f9.4,f10.4,f8.1,2f5.2,4x,2f9.3,2i9)
         else
           if(iuseq.eq.0) then
             write(22,2003) stn5(i),net(i),rlat,rlon,ielev,pdl(i),
     2        sdl(i),nfixst(i),stc(1,i),stc(2,i),nrd(i,1),nrd(i,2)
           else
             write(22,2006) stn5(i),net(i),rlat,rlon,ielev,pdl(i),
     2        nfixst(i),stc(1,i),stc(2,i),nrd(i,1),nrd(i,2)
          endif
 2003     format(a5,3x,a2,1x,f10.4,1x,f10.4,1x,i5,2f5.2,i3,
     2         4x,2f9.3,2i9)
 2006     format(a5,3x,a2,1x,f10.4,1x,f10.4,1x,i5,e11.3,i2,
     2         7x,2f9.3,2i9)
        endif
      endif
c
   97 if(invdel.eq.0) goto  120
         if(nfixst(i).eq.1) goto 118
         nvp=nparv+i
         if(hit(nvp).eq.zero) goto 115
         stderr(nvp)=rat*stderr(nvp)
         nvpi=ndexfx(nvp)
         nvpa=index(nvpi)
         ihit=nint(hit(nvp))
         if(kttfor.ne.3) then
           write(16,125) stn(i),pdl(i),stderr(nvp),drm(nvp),
     2      nrd(i,1),ihit,stn(i),i,nvp,nvpi,nvpa
  125      format(1x,a4,' P  ',f7.3,2f6.3,2i4,1x,a4,2i5,2i5)
         else
           write(16,126)stn6(i),pdl(i),stderr(nvp),drm(nvp),
     2      nrd(i,1),ihit,stn(i),i,nvp,nvpi,nvpa
  126      format(1x,a6,' P  ',f7.3,2f6.3,2i4,1x,a4,2i5,2i5)
         endif
  115    continue
         if(iuses.eq.2) then
            nvs=nvp+nsts
            if(khit(nvs).eq.izero) goto 120
            stderr(nvs)=rat*stderr(nvs)
            nvsi=ndexfx(nvs)
            nvsa=index(nvsi)
            ihit=nint(hit(nvs))
            if(kttfor.ne.3) then
              write(16,195) sdl(i),stderr(nvs),drm(nvs),
     2         nrd(i,2),ihit,stn(i),i,nvs,nvsi,nvsa
  195         format(5x,' S-P',f7.3,2f6.3,2i4,1x,a4,2i5,2i5)
            else
              write(16,196) sdl(i),stderr(nvs),drm(nvs),
     2         nrd(i,2),ihit,stn6(i),i,nvs,nvsi,nvsa
  196         format(5x,' S-P',f7.3,2f6.3,2i4,1x,a6,2i5,2i5)
            endif
         endif
  118    if(kout.lt.2) goto 120
  120 continue
      if(kout.ge.2) close(22)
c
 353  continue
      write(16,1692) nparvi,iuses
 1692 format('in outend after 353, nparvi=',i8,', iuses=',i3)
c  skip velocity output if all fixed
      if(nparvi.eq.0) goto 900
c  compute vp/vs
c     do 360 k=1,nz
c        ks=k+nz
c        do 358 j=1,ny
c           do 356 i=1,nx
c              vpvs(i,j,k)=vel(i,j,k)/vel(i,j,ks)
c 356       continue
c 358    continue
  360 continue
c
c  start cht 1998
c
      close(32)
c
c  output x-y-z-v file!!
c
      write(16,9111)
 9111 format(/,' velocity X-Y-Z grid in file vlxyzltln.out')
      open(53,file='vlxyzltln.out')
      xltd=xlt/60.0
      xlnd=xln/60.0
      if(nzco.eq.1) xlnd=-1.0*xlnd
      write(53,9001) xltd,xlnd,rotadg,nzco,cmerid
 9001 format(' velocity X-Y-Z grid, Origin=',2f10.4,' rota=',f8.2,
     2  ' nzco=',i3,' cmerid=',f10.4)
      if(iuseq.eq.0) then
        write(53,9002)
 9002   format(' latitude longitude      X_km      Y_km    Z_km   ',
     2   '    Vp    Vp/Vs     Vs    Vp-DWS  Vp-DRE  vpvs-DWS vpvs-DRE')
      else
        write(53,9003)
 9003   format(' latitude longitude      X_km      Y_km    Z_km   ',
     2  '    Q        V       Q-DWS    Q-DRE')
      endif
      nz2=nz-2
      ny2=ny-2
      nx2=nx-2
      do 911 k=1,nz2
      k1=k+1
      do 911 j1=1,ny
        j=j1-1
      do 911 i1=1,nx
        i=i1-1
c index for hit,drm
      n1=(k-1)*nxy2+(j-1)*nx2+i
      n2=(nz2+k-1)*nxy2+(j-1)*nx2+i
c
      call latlon(xn(i1),yn(j1),latg,xlatg,long,xlong)
c      print *, 'latlon',xn(i1),yn(j1),latg,xlatg,long,xlong
      rlatg=float(latg)+xlatg/60.0
      if(nzco.eq.1) rlatg= -1.0* rlatg
      rlong=float(long)+xlong/60.0
      if((i1.gt.1).and.(i1.lt.nx).and.(j1.gt.1).and.(j1.lt.ny)) then
        if(iuseq.eq.0) then
         write(53,9116) rlatg,rlong,xn(i1),yn(j1),zn(k1),vel(i1,j1,k1),
     2   vpvs(i1,j1,k1),vel(i1,j1,k1+nz),hit(n1),drm(n1),hit(n2),drm(n2)
        else
         write(53,9115) rlatg,rlong,xn(i1),yn(j1),zn(k1),qval(i1,j1,k1),
     2   vel(i1,j1,k1),hit(n1),drm(n1)
        endif
      else
        if(iuseq.eq.0) then
         write(53,9116) rlatg,rlong,xn(i1),yn(j1),zn(k1),vel(i1,j1,k1),
     2   vpvs(i1,j1,k1),vel(i1,j1,k1+nz),zero,zero,zero,zero
        else
         write(53,9115) rlatg,rlong,xn(i1),yn(j1),zn(k1),qval(i1,j1,k1),
     2   vel(i1,j1,k1),zero,zero
        endif
      endif
 911  continue
 9116 format(2f10.4,2f10.2,f8.2,f8.4,f8.2,2(f11.0,f9.5))
 9115 format(2f10.4,2f10.2,f8.2,f10.1,f8.2,f11.0,f9.5)
      if(iuses.eq.1) goto 364
C-ABOVEc  repeat for vp/vs
C-ABOVE      do 919 k=1,nz2
C-ABOVE      k1=k+1
C-ABOVE      do 919 j=1,ny2
C-ABOVE      j1=j+1
C-ABOVE      do 919 i=1,nx2
C-ABOVE      i1=i+1
C-ABOVEc index for hit,drm
C-ABOVE      n1=(nz2+k-1)*nxy2+(j-1)*nx2+i
C-ABOVEc
C-ABOVE      call latlon(xn(i1),yn(j1),latg,xlatg,long,xlong)
C-ABOVEc      print *, 'latlon',xn(i1),yn(j1),latg,xlatg,long,xlong
C-ABOVE      rlatg=float(latg)+xlatg/60.0
C-ABOVE      rlong=float(long)+xlong/60.0
C-ABOVE      write(53,9116) rlatg,rlong,xn(i1),yn(j1),zn(k1),vpvs(i1,j1,k1),
C-ABOVE     2  hit(n1),drm(n1)
C-ABOVE 919  continue
c
  364 close(53)
c
      if(ires.eq.0) goto 365
c  output diagonal resolution element file, in velocity file format
c
      write(16,9112)
 9112 format(/,' resolution diagonal element grid in file resDRE.out')
      open(59,file='resDRE.out')
c
c  vp first
c
      do 925 k=1,nz2
        k1=k+1
      do 925 j=1,ny2
        j1=j+1
      do 925 i=1,nx2
        i1=i+1
c
c  try to find right part of drm matrix!
      n1=(k-1)*nxy2+(j-1)*nx2+i
      resout(i1,j1,k1)=drm(n1)
  925 continue
      do 926 j1=1,ny
      do 926 i1=1,nx
        resout(i1,j1,1)=0.0
        resout(i1,j1,nz)=0.0
  926 continue
c
      write(59,5902) bld,nx,ny,nz,iuses,date24(21:24),date24(4:16)
 5902 format(f4.1,4i3,7x,'DRE computed',2x,a4,a13)
      if(bld.eq.1.0) then
        write(59,2303) (xn(i),i=1,nx)
        write(59,2303) (yn(i),i=1,ny)
        write(59,2303) (zn(i),i=1,nz)
      else
        write(59,2304) (xn(i),i=1,nx)
        write(59,2304) (yn(i),i=1,ny)
        write(59,2304) (zn(i),i=1,nz)
      endif
      write(59,2305) izero,izero,izero
      do 912 k=1,nz
      do 912 j=1,ny
        write(59,5911) (resout(i,j,k),i=1,nx)
 5911   format(14f9.5)
 912  continue
c
c  now do vp/vs resolution output
      if((iuses.eq.1).or.(iuseq.eq.0)) goto 930
c
      do 929 k=1,nz2
      do 929 j=1,ny2
      do 929 i=1,nx2
c
c  try to find right part of drm matrix!
      n1=(k-1)*nxy2+(j-1)*nx2+i+nx2*ny2*nz2
      resout(i,j,k)=drm(n1)
 929  continue
c
      do 921 k=1,nz
      do 921 j=1,ny
        write(59,5911) (resout(i,j,k),i=1,nx)
  921 continue
c
  930 close(59)
c
c  end cht 1998
c
  365 if(kout.lt.2) goto 311
c  Write the final velocity model to file23 in a format that
c  can be used as input for another run.                     
c  and to file25 in station format for plotting.
      open(unit=23,file='velomod.out')
      if(kout.gt.3)
     2 open(unit=25,file='vl_layer.out')
      write(23,2302) bld,nx,ny,nz,iuses,date24(21:24),date24(4:16)
 2302 format(f4.1,4i3,7x,'computed',2x,a4,a13)
      if(bld.eq.1.0) then
        write(23,2303) (xn(i),i=1,nx)
        write(23,2303) (yn(i),i=1,ny)
        write(23,2303) (zn(i),i=1,nz)
      else
        write(23,2304) (xn(i),i=1,nx)
        write(23,2304) (yn(i),i=1,ny)
        write(23,2304) (zn(i),i=1,nz)
      endif
 2303 format(20f6.0)
 2304 format(20f6.1)
      write(23,2305) izero,izero,izero
 2305 format(3i3)
c      do 38 kv=1,iuses
c Now write vs after vp/vs
      kv=1
         do 36 k=1,nz
            k2=k+(kv-1)*nz
            if((kout.gt.3).and.(k.gt.1).and.(k.lt.nz)) then     
cfhek - introduce average layer velocity for output to file 25
cek average layer velocity:  velav
              velav(k2)=0.0
              do ii1=2,nx-1
                 do kk1=2,ny-1
                    velav(k2)=velav(k2)+vel(ii1,kk1,k2)
                 enddo
              enddo
              velav(k2)=velav(k2)/((nx-2)*(ny-2))
              write(25,2501) k,zn(k),velav(k2)
cfhek 2501       format(1x,'LAYR',1x,i2,8x,f7.2,'k')
 2501         format(1x,'LAYER',1x,i2,8x,f7.2,'km',3x,f7.2,1x,'kmsec-1')
              write(25,4455)
 4455      format(1x,' long, lat, percent off avg-vel, abs. velocity')
cek average layer velocity end
           endif
c
          if(iuseq.eq.0) then
            if(vel(1,1,k2).lt.8.4) then
              write(23,2311) (vel(i,1,k2),i=1,nx)
            else
              write(23,2315) (vel(i,1,k2),i=1,nx)
            endif
          else
              write(23,2511) (qval(i,1,k),i=1,nx)
          endif
            do 34 j=2,ny-1
             if(iuseq.eq.0) then
               if(vel(1,j,k2).lt.8.4) then
                 write(23,2311) (vel(i,j,k2),i=1,nx)
               else
                 write(23,2315) (vel(i,j,k2),i=1,nx)
               endif
             else
                 write(23,2511) (qval(i,j,k),i=1,nx)
             endif
               if((k.eq.1).or.(k.eq.nz)) goto 34
               if(kout.le.3) goto 34
               do 33 i=2,nx-1
                  call latlon(xn(i),yn(j),lat,xlat,lon,xlon)
cfhek
                if(lat.ge.0.and.xlat.ge.0.) then
                   cns='N'
                   lat=iabs(lat)
                   xlat=abs(xlat)
                else
                   cns='S'
                   lat=iabs(lat)
                   xlat=abs(xlat)
                endif
                if(lon.ge.0.and.xlon.ge.0.) then
                   cew='W'
                   lon=iabs(lon)
                   xlon=abs(xlon)
                else
                   cew='E'
                   lon=iabs(lon)
                   xlon=abs(xlon)
                endif
cek
                if(iuseq.eq.0) then
cek                  write(25,2502) vel(i,j,k2),lat,cns,xlat,lon,cew,xlon,
cek     2            xn(i),yn(j),zn(k)
cek 2502             format(2x,f4.2,i2,a1,f5.2,1x,i3,a1,f5.2,10x,3f7.2)
cek  velodiff= total percent velocity change relative to average layer velo
                  velodiff=100.*(vel(i,j,k2)-velav(k2))/velav(k2)
                  alonb=abs(lon+xlon/60.)
                  if(cew.eq.'E') then
                    alona=alonb
                  else
                    alona=-1.*alonb
                  endif
                  alatb=abs(lat+xlat/60.)
                  if(cns.eq.'N') then
                    alata=alatb
                  else
                    alata=-1.*alatb
                  endif
                  write(25,2503) alona,alata,velodiff,vel(i,j,k2)
 2503             format(1x,f9.4,1x,f9.4,1x,f6.2,2x,f6.2)
cfh q-output not changed here for coordinates...
                else
                  write(25,2502) qval(i,j,k),lat,cns,xlat,lon,cew,xlon,
     2              xn(i),yn(j),zn(k)
                endif
 2502           format(2x,f4.2,i2,1x,a1,f5.2,1x,i3,1x,a1,f5.2,10x,3f7.2)
   33          continue
   34       continue
          if(iuseq.eq.0) then
            if(vel(1,ny,k2).lt.8.4) then
              write(23,2311) (vel(i,ny,k2),i=1,nx)
            else
              write(23,2315) (vel(i,ny,k2),i=1,nx)
            endif
          else
              write(23,2511) (qval(i,ny,k),i=1,nx)
          endif
   36    continue
   38 continue
c  write Vp/Vs at end of new velocity file
      if(iuses.eq.1) goto 50
      do 42 k=1,nz
         write(23,2314) (vpvs(i,1,k),i=1,nx)
         do 134 j=2,ny-1
            write(23,2314) (vpvs(i,j,k),i=1,nx)
            if((k.eq.1).or.(k.eq.nz)) goto 134
            if(kout.le.3) goto 134
            do 133 i=2,nx-1
               call latlon(xn(i),yn(j),lat,xlat,lon,xlon)
cekfh
              if(lat.ge.0.and.xlat.ge.0.) then
                   cns='N'
                   lat=iabs(lat)
                   xlat=abs(xlat)
              else
                   cns='S'
                   lat=iabs(lat)
                   xlat=abs(xlat)
              endif
              if(lon.ge.0.and.xlon.ge.0.) then
                   cew='W'
                   lon=iabs(lon)
                   xlon=abs(xlon)
              else
                   cew='E'
                   lon=iabs(lon)
                   xlon=abs(xlon)
              endif
              write(25,2502) vpvs(i,j,k),lat,cns,xlat,lon,cew,xlon,
     2             xn(i),yn(j),zn(k)
cek
  133       continue
  134    continue
         write(23,2314) (vpvs(i,ny,k),i=1,nx)
   42 continue
c write out Vs now, after vpvs
      kv=2
      do 46 k=1,nz
        k2=k+(kv-1)*nz
        write(23,2311) (vel(i,1,k2),i=1,nx)
        do 44 j=2,ny-1
          write(23,2311) (vel(i,j,k2),i=1,nx)
   44   continue
        write(23,2311) (vel(i,ny,k2),i=1,nx)
   46 continue
c
 2311 format(22f5.2)
 2312 format(22f5.1)
 2314 format(22f6.3)
 2315 format(22f6.2)
 2511 format(20f7.1)
   50 close(unit=23)
      if(kout.gt.3) close(unit=25)
c
  311 do 495 kv=1,iuses
         if((kv.eq.1).and.(iusep.eq.0)) goto 495
  312    if(iuseq.eq.0) then
           write(16,1000) type(kv)
       else
           write(16,1501) type(3)
       endif
 1000    format(//,' FINAL ',a1,'-VELOCITY MODEL')
 1501    format(//,' FINAL ',a1,' MODEL')
         do 10 k=1,nz
            k2=k+(kv-1)*nz
            if(sumhit(k2).eq.0.0) goto 10
c write out vp/vs
            if(kv.eq.2) then
               write(16,1009) k,vtype(kv),zn(k)
               do 9 j=1,ny
                  write(16,1001) (vpvs(i,j,k),i=1,nx)
    9          continue
            endif
          if(iuseq.eq.0) then
              write(16,1010) k,type(kv),zn(k)
          else
              write(16,1511) k,type(3),zn(k)
          endif
 1009       format(/,' layer',i3,5x,a5,' nodes',10x,'z =',f7.1)
 1010       format(/,' layer',i3,5x,a1,'-velocity nodes',10x,
     2         'z =',f7.1)
 1511       format(/,' layer',i3,5x,a1,' nodes',10x,
     2         'z =',f7.1)
            do 11 j=1,ny
             if(iuseq.eq.0) then
                 write(16,1001) (vel(i,j,k2),i=1,nx)
             else
                 write(16,1521) (qval(i,j,k),i=1,nx)
             endif
   11       continue
   10    continue
c
         write(16,1003) 
         nz1=nz-1
         ny1=ny-1
         do 22 k=2,nz1
            k2=k+(kv-1)*nz
            if(sumhit(k2).eq.0.0) goto 22
          if(iuseq.eq.0) then
              write(16,1009) k,vtype(kv),zn(k)
          else
              write(16,1009) k,vtype(3),zn(k)
          endif
            kk=k+(kv-1)*nz2
            do 20 j=2,ny1
               n1=(kk-2)*nxy2+(j-2)*nx2+1
               n2=n1+nx2-1
               write(16,1005) (khit(i),i=n1,n2),izero
   20       continue
   22    continue
c
 1003    format(/,' OBSERVATION MATRIX - KHIT - (will be 0 for',
     2      ' fixed nodes)')
 1005    format('     0',20i6)
 1001    format(22f6.2)
 1521    format(20f7.1)
c
         write(16,1633)
 1633    format(/,' DERIVATIVE WEIGHT SUM -- only print non-zero ',
     2    'rows -- ')
         do 487 k=2,nz1
            k2=k+(kv-1)*nz
            if(sumhit(k2).eq.0.0) goto 487
          if(iuseq.eq.0) then
              write(16,1009) k,vtype(kv),zn(k)
          else
              write(16,1009) k,vtype(3),zn(k)
          endif
            kk=k+(kv-1)*nz2
            do 486 j=2,ny1
              n1=(kk-2)*nxy2+(j-2)*nx2+1
              n2=n1+nx2-1
              ih0=0
              do 485 i=n1,n2
                if(hit(i).gt.0.0) ih0=1
  485         continue
              if(ih0.eq.0) goto 486
              n12=n2-n1+1
c      PRINT *,'j',j,' n1,n2,n12',n1,n2,n12
              if(n12.le.18) then
                write(16,1635)  j,(hit(i),i=n1,n2),zero
 1635           format('j',i3,' 0.',19f7.0)
              else
                n118=n1+18
                n119=n1+19
                write(16,1635)  j,(hit(i),i=n1,n118)
c      PRINT *,'n118,n119',n118,n119
                write(16,1639)  (hit(i),i=n119,n2),zero
 1639           format(20f7.0)
              endif
  486       continue
  487    continue
c
         if(ires.eq.0) goto 495
c
c  assign -1.0 to diagonal resolution for fixed or linked gridpoint
         do 490 n=1,nparv
            if(nfix(n).eq.1.or.imerge(n).eq.1) drm(n)=-1.0
  490    continue
c
       if(iuseq.eq.0) then
           write(16,1634)
       else
           write(16,1534)
       endif
 1634    format(/,' RESOLUTION : GRIDOINT NUMBER, DIAGONAL RESOLUTION',
     2      ' ELEMENT (-1.0 indicates fixed velocity gridpoint)')
 1534    format(/,' RESOLUTION : GRIDOINT NUMBER, DIAGONAL RESOLUTION',
     2      ' ELEMENT (-1.0 indicates fixed Q gridpoint)')
       if(ires.eq.1) then
            write(47,1637) nx-2,ny-2,nz-2,date24(21:24),date24(4:16)
 1637       format('1 ',3i3,6x,'Resolution Matrix, Computed',2x,a4,a13)
            write(47,2304) (xn(i),i=2,nx-1)
            write(47,2304) (yn(i),i=2,ny-1)
            write(47,2304) (zn(i),i=2,nz-1)
            write(47,2305) izero,izero,izero
       endif

         do 489 k=2,nz1
            k2=k+(kv-1)*nz
c
c          to make the Q output quite similar to the velmodel
c
c
          if(iuseq.eq.0) then
              if(sumhit(k2).eq.0.0) goto 489
              write(16,1009) k,vtype(kv),zn(k)
          else
              write(16,1009) k,vtype(3),zn(k)
          endif
            kk=k+(kv-1)*nz2
            do 488 j=2,ny1
               n1=(kk-2)*nxy2+(j-2)*nx2+1
               n2=n1+nx2-1
               write(16,1636) (i,drm(i),i=n1,n2)
             if(ires.eq.1) then
                 write(47,1638) (drm(i),i=n1,n2)
             endif
c 1636          format(1x,17(i5,':',f7.4))
 1636          format(1x,13(i5,':',f7.4))
 1638          format(20f8.4)
  488       continue
  489    continue
  495 continue
c
      rat=ssqrw/ssqrw1
      if(ires.eq.0) goto 900
c
c-ALWAYS WRITE OUT NODES FILEc  write out pointer from full nodes to inversion nodes
c-ALWAYS WRITE OUT NODES FILE      if(inf.le.0) goto 545
cfhdmep - unit 18 now opened in main; npar... written in rescov 
cfh      open(unit=18,carriagecontrol='list',status='new')
cfh      write(18,1800) npar,nparv,npari,nparvi
cfh 1800 format(4i7,'   0.0001')
      do 540 l=1,npar
         write(18,1801) l,ndexfx(l)
 1801    format(2i7)
  540 continue
      close(18)
c
  545 do 550 l=1,npar
         if(khit(l).eq.0) goto 550
         if ((hit(l).lt.hitct).or.(nfix(l).eq.1).or.
     & (imerge(l).eq.1)) go to 550
         stderr(l)=rat*stderr(l)
c        write(16,2004) l,drm(l),stderr(l)
 2004    format('  parameter number ',i3,': resolution:',f8.4,
     *      '  error (%):',f8.4)
c
  550 continue
c  Also write resolution in  more readable format
c   problems with this printout format
      if(iusesq.eq.0) goto 80
      if(nx.gt.20) goto 80
      if(iuseq.eq.0) then
        write(16,1026)
      else
        write(16,1526)
      endif
 1026 format(/,10x,'VELOCITY WITH STANDARD ERROR(km/s) AND RESOLUTION'
     2 /,'  where error and res are both 0.00, the node was not',
     3 ' inverted for',
     4 ', where resol.=-1.00, gridpoint velocity was fixed')
 1526 format(/,10x,'Q WITH STANDARD ERROR(km/s) AND RESOLUTION'
     2 /,'  where error and res are both 0.00, the node was not',
     3 ' inverted for',
     4 ', where resol.=-1.00, gridpoint Q was fixed')
      nx2=nx-2
      nxy2=nx2*(ny-2)
cfh problems if nx.ge.20!!!
c  create proper sized format for nx nodes
      if(nx2.le.10) then
        write(line,27) nx2
        read(line,28) ch2nx2
   27   format(1x,i2)
   28   format(1x,a2)
        formvl(11:12)=ch2nx2
        formrs(14:15)=ch2nx2
      else
        nx3=nx2-10
        write(line,25) nx3
        read(line,26) ch2nx3
cfh next 2 lines
c   25   format(1x,i1)
c   26   format(1x,a1)
   25   format(1x,i2.2)
   26   format(1x,a2)
        nx4=10-nx3
        write(line,25) nx4
        read(line,26) ch2nx4
        frmvl2(9:9)=ch2nx4
        frmvl2(16:16)=ch2nx3
        frmrs2(5:5)=ch2nx4
        frmrs2(12:12)=ch2nx3
      endif
      do 75 kv=1,iuses
         do 70 k=2,nz-1
            k2=k+(kv-1)*nz 
            if(sumhit(k2).eq.0.0) goto 70
            kk=k+(kv-1)*nz2
            if(iuseq.eq.0) then
              write(16,1009) k,vtype(kv),zn(k)
            else
              write(16,1009) k,vtype(3),zn(k)
            endif
            if(kv.eq.1) then
             if(iuseq.eq.0) then
                 write(16,1021) (vel(i,1,k2),i=1,nx)
             else
                 write(16,1022) (qval(i,1,k),i=1,nx)
             endif
            else
               write(16,1021) (vpvs(i,1,k),i=1,nx)
            endif
            kp=(kk-2)*nxy2
            do 60 j=2,ny-1
               jp=(j-2)*nx2
               kjp=kp+jp
c  convert standard error from % to km/s
               do 55 i=1,nx2
                  v=vel(i+1,j,k2)
                  if(kv.eq.2) v=vpvs(i+1,j,k)
                  if(iuseq.eq.1) v=qval(i+1,j,k)
                  sekms(i)=stderr(kjp+i)*v
                  if(stderr(kjp+i).gt.tol) then
                     covdi(kjp+i)=1.0/(sekms(i)*sekms(i))
                  endif
   55          continue
               if(nx2.le.10) then
                  if(kv.eq.1) then
                   if(iuseq.eq.0) then
                       write(16,formvl) vel(1,j,k2),(vel(i+1,j,k2),
     2                    sekms(i),i=1,nx2),vel(nx,j,k2)
                   else
                       write(16,formvl) qval(1,j,k),(qval(i+1,j,k),
     2                    sekms(i),i=1,nx2),qval(nx,j,k)
                   endif
                  else
                     write(16,formvl) vpvs(1,j,k),(vpvs(i+1,j,k),
     2                  sekms(i),i=1,nx2),vpvs(nx,j,k)
                  endif
                  write(16,formrs) (drm(kjp+i),i=1,nx2)
               else
                  if(kv.eq.1) then
                   if(iuseq.eq.0) then
                       write(16,frmvl1) vel(1,j,k2),(vel(i+1,j,k2),
     2                    sekms(i),i=1,10),vel(nx,j,k2)
                   else
                       write(16,frmvl1) qval(1,j,k),(qval(i+1,j,k),
     2                    sekms(i),i=1,10),qval(nx,j,k)
                   endif
                  else
                     write(16,frmvl1) vpvs(1,j,k),(vpvs(i+1,j,k),
     2                  sekms(i),i=1,10),vpvs(nx,j,k)
                  endif
                  write(16,frmrs1) (drm(kjp+i),i=1,10)
                  if(kv.eq.1) then
                   if(iuseq.eq.0) then
                       write(16,frmvl2,err=60) (vel(i+1,j,k2),sekms(i),
     2                    i=11,nx2),vel(nx,j,k2)
                   else
                       write(16,frmvl2) (qval(i+1,j,k),sekms(i),
     2                    i=11,nx2),qval(nx,j,k)
                   endif
                  else
                     write(16,frmvl2,err=60) (vpvs(i+1,j,k),sekms(i),
     2                  i=11,nx2),vpvs(nx,j,k)
                  endif
                  write(16,frmrs2) (drm(kjp+i),i=11,nx2)
               endif
   60       continue
            if(kv.eq.1) then
             if(iuseq.eq.0) then
                 write(16,1024) (vel(i,ny,k2),i=1,nx)
             else
                 write(16,1025) (qval(i,ny,k),i=1,nx)
             endif
            else
               write(16,1024) (vpvs(i,ny,k),i=1,nx)
            endif
   70    continue
   75 continue
   16 format('0',f4.2,12(f7.2,'+',f4.2))
 1021 format(f5.2,f7.2,10f12.2)
 1022 format(2f6.0,10f12.2)
 1023 format(1x,'res',6x,10('(',f4.2,')',6x))
 1024 format('0',f4.2,f7.2,10f12.2)
 1025 format(2f6.0,10f12.2)
c
   80 continue
c  write out 1/covariance(diag) file in format like velocity file
c  so can be plotted
      if((kout2.eq.5).and.(ires.gt.0)) then
        open(unit=45,file='covarid.out')
        write(45,4502) nx,ny,nz,date24(21:24),date24(4:16)
 4502   format(3i3,6x,'1/Diag. Covariance, Computed',2x,a4,a13)
        write(45,2304) (xn(i),i=1,nx)
        write(45,2304) (yn(i),i=1,ny)
        write(45,2304) (zn(i),i=1,nz)
        write(45,2305) izero,izero,izero
      write(16,1739)
 1739 format(/,' 1/(Covariance Diag. )')
      do 746 kv=1,iuses
c  write peripheral grids to file45 also
         do 725 j=1,ny
            write(45,1737) (zeroi(i),i=1,nx)
  725    continue
         do 736 k=2,nz1
            k2=k+(kv-1)*nz
            if(sumhit(k2).gt.0.0) then
            if(iuseq.eq.0) then
                 write(16,1009) k,vtype(kv),zn(k)
            else
                 write(16,1009) k,vtype(3),zn(k)
            endif
          endif
            kk=k+(kv-1)*nz2
            write(45,1737) (zeroi(i),i=1,nx)
            do 735 j=2,ny1
               n1=(kk-2)*nxy2+(j-2)*nx2+1
               n2=n1+nx2-1
               if(sumhit(k2).gt.0.0)
     2            write(16,1737) zero,(covdi(i),i=n1,n2),zero
 1737          format(14e9.2)
               write(45,1737) zero,(covdi(i),i=n1,n2),zero
  735       continue
  736    continue
         write(45,1737) (zeroi(i),i=1,nx)
         do 745 j=1,ny
            write(45,1737) (zeroi(i),i=1,nx)
  745    continue
  746 continue
      close(45)
      endif
c
  900 continue
      if(i3d.le.0) goto 950
      upbavg=upbtot/(float(nupb))
      write(16,1665) nupb, upbavg
 1665 format(/,' For',i10,' calls to MINIMA, average number',
     2 ' pseudo-bending iter. =',f7.2)
c  950 call time(tm)
  950 call fdate(date24)
      write(16,1660) date24(21:24),date24(4:16)
 1660 format(/,' Computation finished at ',a4,a13)
c***** end of subroutine outend *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine outhit(nit)
c write out DWS for all nodes (including fixed and linked) when nitmax=1
c
c write to file 43 like velocity file format
c write to file 44 with layers, ysec and xsec printed out
c
c  common block variables:
      include 'simul2014_common.inc'
c
      real zeroi(maxnx),hitx(maxny)
      character*24 date24
      character*5 vtype(3)
      character*1 chi,chj,chk
      parameter(zero=0.0,izero=0)
      chi='i'
      chj='j'
      chk='k'
      vtype(1)='P-Vel'
      vtype(2)='Vp/Vs'
      vtype(3)='Q'
      call fdate(date24)
      do 2 i=1,maxnx
         zeroi(i)=zero
    2 continue
      write(16,1623)
 1623 format(/,' DWS for ALL nodes. Use for planning fixed, linked.')
      if(nit.eq.0) write(16,1624)
 1624 format('   Note: for iter=0 only, DWS uses wtsht for estimated ',
     2 'combined wt for shots.')
c Write the DWS for ALL nodes to a file in format like velocity file.
c This can then be used to plot or calculate linked nodes.
      open(unit=43,file='DWSALL')
      write(43,4302) bld,nx,ny,nz,iuses,date24(21:24),date24(4:16)
c      open(unit=44,file='DWSALL.GRD')
c      write(44,4302) nx,ny,nz,date24(21:24),date24(4:16)
 4302 format(f4.1,4i3,6x,'DWS - All nodes, Computed',2x,a4,a13)
      if(bld.eq.1.0) then
        write(43,2303) (xn(i),i=1,nx)
        write(43,2303) (yn(i),i=1,ny)
        write(43,2303) (zn(i),i=1,nz)
      else
        write(43,2304) (xn(i),i=1,nx)
        write(43,2304) (yn(i),i=1,ny)
        write(43,2304) (zn(i),i=1,nz)
      endif
      write(43,2305) izero,izero,izero
c
 2303 format(20f6.0)
 2304 format(20f6.1)
 2305 format(3i3)
      do 446 kv=1,iuses
c  write peripheral grids to file43 also
         do 425 j=1,ny
            write(43,4301) (zeroi(i),i=1,nx)
  425    continue
         do 436 k=2,nz1
            if(iuseq.eq.0) then
                 write(16,1009) k,vtype(kv),zn(k)
c                 write(44,1009) k,vtype(kv),zn(k)
            else
                 write(16,1009) k,vtype(3),zn(k)
c                 write(44,1009) k,vtype(3),zn(k)
            endif
            kk=k+(kv-1)*nz2
            write(43,4301) (zeroi(i),i=1,nx)
c            write(44,4401) chj,chi,(i,i=2,nx1)
            do 435 j=2,ny1
               n1=(kk-2)*nxy2+(j-2)*nx2+1
               n2=n1+nx2-1
               write(16,1635) (hitall(i),i=n1,n2)
c               write(44,4402) j,(hitall(i),i=n1,n2)
               write(43,4301) zero,(hitall(i),i=n1,n2),zero
  435       continue
            write(43,4301) (zeroi(i),i=1,nx)
  436    continue
c         write(43,4301) (zeroi(i),i=1,nx)
         do 445 j=1,ny
            write(43,4301) (zeroi(i),i=1,nx)
  445    continue
c
c  now print out ysec dws to file 44
         do 456 j=2,ny1
          if(iuseq.eq.0) then
c               write(44,1010) j,vtype(kv),yn(j)
          else
c               write(44,1010) j,vtype(3),yn(j)
          endif
c          write(44,4401) chk,chi,(i,i=2,nx1)
          do 450 k=2,nz1
            kk=k+(kv-1)*nz2
            n1=(kk-2)*nxy2+(j-2)*nx2+1
            n2=n1+nx2-1
c            write(44,4402) k,(hitall(i),i=n1,n2)
  450       continue
  456    continue
c
c
c  now print out xsec dws to file 44
         do 466 i=2,nx1
          if(iuseq.eq.0) then
c               write(44,1011) i,vtype(kv),xn(i)
          else
c               write(44,1011) i,vtype(3),xn(i)
          endif
c          write(44,4401) chk,chj,(j,j=2,ny1)
          do 462 k=2,nz1
            kk=k+(kv-1)*nz2
            n1=(kk-2)*nxy2+(j-2)*nx2+1
            do 460 j=2,ny1
              n=(kk-2)*nxy2+(j-2)*nx2+i
              hitx(j)=hitall(n)
  460       continue
c            write(44,4402) k,(hitx(j),j=2,ny1)
  462       continue
  466    continue
  446 continue
c
 1009       format(/,' layer',i3,5x,a5,' nodes',10x,'z =',f7.1)
 1010       format(/,' y-grid',i3,5x,a5,' nodes',10x,'y =',f7.1)
 1011       format(/,' x-grid',i3,5x,a5,' nodes',10x,'x =',f7.1)
 1635          format(20f7.0,f6.0)
 4301 format(20f9.0)
 4401 format(4x,a1,3x,a1,i3,20i10)
 4402 format(i5,2x,20f10.0)
      close(43)
c      close(44)
c
c***** end of subroutine outhit *****
      return
      end
c
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine outlis(ne)
c  write out each event in a format similar to HYPO71 listing format
c  so that it can be used as input to FPFIT fault-plane solution
c  program
c
c  declaration statements:
      parameter (one=1.00)
      character*1 phs(2)
      character*1 cew,cns
      character*24 date24
c
c  common block variables:
      include 'simul2014_common.inc'
      common/shortd/ xltkm,xlnkm,rota,nzco,xlt,xln,snr,csr,
     2 onorth,oeast,rotadg,cmerid,rlatdo,rlondo
      common/wtpars/w(mxobsa)
c
      data phs/'P','S'/
cfhek distiguish North-South and East-West by N,S,E,W
      cew='W'
      cns='N'
c
      if (ne.gt.1) goto 5
c  print header line
      call fdate(date24)
      write(34,3407) date24(21:24),date24(4:16)
 3407 format(/,' AZIM and TOA calculated with 3D velocity model ',
     2 '(SIMUL2014) ',a4,a13,/)
    5 write(34,3401)
 3401 format(2x,'DATE',4x,'ORIGIN',3X,'LATITUDE LONGITUDE  DEPTH',4X,
     2 'MAG NO',11X,'RMS')
      call latlon(evc(1,ne),evc(2,ne),lat,xlat,lon,xlon)
cek
              if(lat.ge.0.and.xlat.ge.0.) then
                   cns='N'
                   lat=iabs(lat)
                   xlat=abs(xlat)
              else
                   cns='S'
                   lat=iabs(lat)
                   xlat=abs(xlat)
              endif
              if(lon.ge.0.and.xlon.ge.0.) then
                   cew='W'
                   lon=iabs(lon)
                   xlon=abs(xlon)
              else
                   cew='E'
                   lon=iabs(lon)
                   xlon=abs(xlon)
              endif
cek
cek  
      if((nzco.eq.1).or.(nzco.eq.4)) then
        cew='E'
        cns='S'
      endif
      if(iuse2t.eq.0) then
        write(34,3402) iyrmo(ne),iday(ne),ihr(ne),mino(ne),seco(ne),
     2  lat,cns,xlat,lon,cew,xlon,evc(3,ne),rmag(ne),kobs(ne),wrmsr(ne)
      else
        write(34,3402) iyrmo(ne),iday(ne),ihr(ne),mino(ne),seco(ne),
     2  lat,cns,xlat,lon,cew,xlon,evc(3,ne),rmag(ne),kobs(ne),wrmsr(ne),
     3  seco2(ne)
      endif
 3402 format(1x,a4,a2,1x,a2,i2,f6.2,1x,i2,a1,f5.2,1x,i3,a1,f5.2,f7.2,
     2 2x,f5.2,i3,9x,f5.2,t69,f6.2)
cek  303 write(34,3402) iyrmo(ne),iday(ne),ihr(ne),mino(ne),seco(ne),
cek     2 lat,xlat,lon,xlon,evc(3,ne),rmag(ne),kobs(ne),wrmsr(ne)
cek 3402 format(1x,a4,a2,1x,a2,i2,f6.2,1x,i2,'n',f5.2,1x,i3,'w',f5.2,f7.2,
cek     2 2x,f5.2,i3,9x,f5.2)
c Actual hypo71 list format is 3403, but need 6-char station
c
c      write(34,3403)
 3403 format(/,2x,'STN  DIST  AZ TOA PRMK HRMN  PSEC TPOBS',
     2 14X,'PRES  PWT')
      write(34,3404)
 3404 format(/,2x,'STN    DIST  AZ TOA PRMK HRMN   PSEC  TPOBS',
     2  10x,'PRES  PWT')
      nobs=kobs(ne)
      do 100 i=1,nobs
c  don't use S arrivals for fault-plane solution
C         if(intsp(i,ne).eq.1)  goto 100
        iaz=nint(az(i))
        itoa=nint(toa(i))
        if(intsp(i,ne).eq.0) then
         if(kttfor.lt.3) then
           write(34,3408) stn(isto(i,ne)),dlta(i,ne),iaz,
     2      itoa,rmk(i,ne),ihr(ne),mino(ne),secp(i,ne),
     3      (secp(i,ne)-seco(ne)),res(i),wtcomb(i,ne)
         else
           write(34,3409) stn6(isto(i,ne)),dlta(i,ne),iaz,
     2      itoa,rmk(i,ne),ihr(ne),mino(ne),secp(i,ne),
     3      (secp(i,ne)-seco(ne)),res(i),wtcomb(i,ne)
         endif
        else
c write out S-P in list format
         if(kttfor.lt.3) then
           write(34,3410) stn(isto(i,ne)),dlta(i,ne),iaz,
     2      rmk(i,ne),ihr(ne),mino(ne),secp(i,ne),
     3      res(i),wtcomb(i,ne)
         else
           write(34,3411) stn6(isto(i,ne)),dlta(i,ne),iaz,
     2      rmk(i,ne),ihr(ne),mino(ne),secp(i,ne),
     3      res(i),wtcomb(i,ne)
         endif
        endif
  100 continue
c hypo71
 3405 format(1x,a4,1x,f5.1,1x,i3,1x,i3,1x,a4,
     2 1x,a2,i2,f6.2,f6.2,
     3 13x,f5.2,1x,f5.2)
c new 4-character
 3408 format(1x,a4,2x,f6.1,1x,i3,1x,i3,1x,a4,
     2 1x,a2,i2,f7.2,f7.2,
     3 9x,f5.2,1x,f5.2)
 3410 format(1x,a4,2x,f6.1,1x,i3,5x,a4,
     2 1x,a2,i2,7x,f7.2,
     3 9x,f5.2,1x,f5.2)
c 6-character
 3409 format(1x,a6,f6.1,1x,i3,1x,i3,1x,a4,
     2 1x,a2,i2,f7.2,f7.2,
     3 9x,f5.2,1x,f5.2)
 3411 format(1x,a6,f6.1,1x,i3,5x,a4,
     2 1x,a2,i2,7x,f7.2,
     3 9x,f5.2,1x,f5.2)
      write(34,3406)
 3406 format(/)
      return
c ****** end of subroutine outlis ******
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine outres(ne)
c  write out station residuals and weights ( from parsep) for each event
c
c  common block variables:
      include 'simul2014_common.inc'
      common/wtpars/w(mxobsa)
c
c  declaration statements:
      real ttobs(maxobs)
      character*1 phs(2)
      parameter (one=1.00)
c
      data phs/'P','S'/
c
      do 10 j=1,kobs(ne)
        if(ne.gt.nebs) then
          nt=ne-nebs
          ttobs(j)=secte(j,nt)-telpad(nt)
          if(intsp(j,ne).eq.1) ttobs(j)=secte(j,nt)-telpads(nt)
        else
          ttobs(j)=secp(j,ne)
          if(intsp(j,ne).eq.0) then
            if((iuse2t.eq.0).or.(iclock(isto(j,ne)).eq.0)) then
              ttobs(j)=secp(j,ne)-seco(ne)
            else
              ttobs(j)=secp(j,ne)-seco2(ne)
            endif
          endif
        endif
   10 continue
      if(ne.gt.(neqs+nbls)) goto 200
      write(16,3) ne,iyrmo(ne),iday(ne),ihr(ne),mino(ne),seco(ne),
     2 rmag(ne)
      write(20,3) ne,iyrmo(ne),iday(ne),ihr(ne),mino(ne),seco(ne),
     2 rmag(ne)
 3    format(1h0,' residuals and weights(parsep) for event=',i4, 
     2 3x,a4,a2,1x,a2,i2,1x,f5.2,'  mag ',f4.2)
      if(kttfor.ne.3) then
        write(16,51)
51      format(1x,4('sta  ph  wt  res:O-C ttobs delta',1x))
        write(16,53) (stn(isto(j,ne)),rmk(j,ne),w(j),res(j),
     2  ttobs(j),dlta(j,ne),j=1,kobs(ne))
53      format(4(1x,a4,a4,f5.2,f7.3,2f6.2))
      else
        write(16,41)
   41   format(1x,4('sta6   ph  wt  res:O-C ttobs delta',1x))
        write(16,43) (stn6(isto(j,ne)),rmk(j,ne),w(j),res(j),
     2  ttobs(j),dlta(j,ne),j=1,kobs(ne))
   43   format(4(1x,a6,a4,f5.2,f7.3,2f6.2))
      endif
      write(20,55)
   55 format(1x,'sta  ph  wt  res:O-C ttobs ttcal  delta',
     2 4x,'x ev',4x,'y ev',4x,'z ev',4x,'x st',4x,'y st',
     3 4x,'z st')
      do 20 j=1,kobs(ne)
        ttcal=ttobs(j)-res(j)
        if(kttfor.ne.3) then
          write(20,54) stn(isto(j,ne)),rmk(j,ne),w(j),res(j),
     2      ttobs(j),ttcal,dlta(j,ne),
     3      (evc(ie,ne),ie=1,3),(stc(is,isto(j,ne)),is=1,3)
        else
          write(20,44)stn6(isto(j,ne)),rmk(j,ne),w(j),res(j),
     2      ttobs(j),ttcal,dlta(j,ne),
     3      (evc(ie,ne),ie=1,3),(stc(is,isto(j,ne)),is=1,3)
        endif
   20 continue
   44 format(a6,a4,f5.2,f7.3,2f7.3,f7.2,6f8.2)
   54 format(1x,a4,a4,f5.2,f8.3,2f8.3,f8.2,6f8.2)
      write(16,1601)
 1601 format(/)
      return
  200 continue
      write(16,33) ne,iyrmo(ne),iday(ne),ihr(ne),mino(ne),seco(ne),
     2 rmag(ne)
      write(20,33) ne,iyrmo(ne),iday(ne),ihr(ne),mino(ne),seco(ne),
     2 rmag(ne)
 33   format(1h0,' residuals and weights(medder) for event=',i4, 
     2 3x,a4,a2,1x,a2,i2,1x,f5.2,'  mag ',f4.2)
      if(kttfor.ne.3) then
        write(16,51)
        write(16,53) (stn(isto(j,ne)),rmk(j,ne),w(j),res(j),
     2  ttobs(j),dlta(j,ne),j=1,kobs(ne))
      else
        write(16,41)
        write(16,43)(stn6(isto(j,ne)),rmk(j,ne),w(j),res(j),
     2  ttobs(j),dlta(j,ne),j=1,kobs(ne))
      endif
      write(20,55)
      do 30 j=1,kobs(ne)
        ttcal=ttobs(j)-res(j)
        if(kttfor.ne.3) then
          if(ne.le.nebs) then
            write(20,54) stn(isto(j,ne)),rmk(j,ne),w(j),res(j),
     2      ttobs(j),ttcal,dlta(j,ne),
     3      (evc(ie,ne),ie=1,3),(stc(is,isto(j,ne)),is=1,3)
          else
            write(20,54) stn(isto(j,ne)),rmk(j,ne),w(j),res(j),
     2      ttobs(j),ttcal,dlta(j,ne),
     3      (tepp(ie,j,nt),ie=1,3),(stc(is,isto(j,ne)),is=1,3)
          endif
        else
          if(ne.le.nebs) then
            write(20,44)stn6(isto(j,ne)),rmk(j,ne),w(j),res(j),
     2      ttobs(j),ttcal,dlta(j,ne),
     3      (evc(ie,ne),ie=1,3),(stc(is,isto(j,ne)),is=1,3)
          else
            write(20,44)stn6(isto(j,ne)),rmk(j,ne),w(j),res(j),
     2      ttobs(j),ttcal,dlta(j,ne),
     3      (tepp(ie,j,nt),ie=1,3),(stc(is,isto(j,ne)),is=1,3)
          endif
        endif
   30 continue
      write(16,1601)
      return
c ****** end of subroutine outres ******
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine outresce(nc,ne)
c  write out station residuals and weights ( from parsep) for each 
c  cluster_event
c
c  common block variables:
      include 'simul2014_common.inc'
      common/wtpars/w(mxobsa)
c
c  declaration statements:
      real ttobs(maxobs)
      character*1 phs(2)
      parameter (one=1.00)
c
      data phs/'P','S'/
c
      do 10 j=1,kobsce(ne,nc)
         ttobs(j)=secpc(j,ne,nc)
         if(intspc(j,ne,nc).eq.0) then
           ttobs(j)=secpc(j,ne,nc)-secoce(ne,nc)
         endif
   10 continue
      write(16,33) ne,iyrmoce(ne,nc),idayce(ne,nc),ihrce(ne,nc),
     2 minoce(ne,nc),secoce(ne,nc),rmagce(ne,nc)
      write(30,33) ne,iyrmoce(ne,nc),idayce(ne,nc),ihrce(ne,nc),
     2 minoce(ne,nc),secoce(ne,nc),rmagce(ne,nc)
      if(kttfor.ne.3) then
        write(16,51)
        write(16,53) (stn(istoc(j,ne,nc)),rmkc(j,ne,nc),w(j),res(j),
     2  ttobs(j),dltac(j,ne,nc),j=1,kobsce(ne,nc))
      else
        write(16,41)
        write(16,43) (stn6(istoc(j,ne,nc)),rmkc(j,ne,nc),w(j),res(j),
     2  ttobs(j),dltac(j,ne,nc),j=1,kobsce(ne,nc))
      endif
      write(30,55)
      do 30 j=1,kobsce(ne,nc)
        ttcal=ttobs(j)-res(j)
        if(kttfor.ne.3) then
          write(30,54) stn(istoc(j,ne,nc)),rmkc(j,ne,nc),w(j),res(j),
     2    ttobs(j),ttmc(j,ne),dltac(j,ne,nc),
     3    (cevc(ie,ne,nc),ie=1,3),(stc(is,istoc(j,ne,nc)),is=1,3)
        else
          write(30,44)stn6(istoc(j,ne,nc)),rmkc(j,ne,nc),w(j),res(j),
     2    ttobs(j),ttmc(j,ne),dltac(j,ne,nc),
     3    (cevc(ie,ne,nc),ie=1,3),(stc(is,istoc(j,ne,nc)),is=1,3)
        endif
   30 continue
   33 format(1h0,' residuals and weights(medder) cluster_event=',i4, 
     2 3x,a4,a2,1x,a2,i2,1x,f5.2,'  mag ',f4.2)
   41 format(1x,4('sta6   ph  wt  res:O-C ttobs delta',1x))
   51 format(1x,4('sta  ph  wt  res:O-C ttobs delta',1x))
   55 format(1x,'sta  ph  wt  res:O-C ttobs ttcal  delta',
     2 4x,'x ev',4x,'y ev',4x,'z ev',4x,'x st',4x,'y st',
     3 4x,'z st')
   43 format(4(1x,a6,a4,f5.2,f7.3,2f6.2))
   53 format(4(1x,a4,a4,f5.2,f7.3,2f6.2))
   44 format(a6,a4,f5.2,f7.3,2f7.3,f7.2,6f8.2)
   54 format(1x,a4,a4,f5.2,f7.3,2f7.3,f7.2,6f8.2)
      write(16,1601)
 1601 format(/)
      return
c ****** end of subroutine outresce ******
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine outresedt(nc)
c  write out station residuals and weights (from medder) for each event
c  only write first 60 to print output, all to File 31
c
c  common block variables:
      include 'simul2014_common.inc'
      common/wtpars/w(mxobsa)
c
      write(16,51) nc,(ccenc(i,nc),i=1,3),wrmsedt(nc)
      write(31,51) nc,(ccenc(i,nc),i=1,3),wrmsedt(nc)
   51 format(1h0,'residuals for eq-pair diff-time: cluster=',i5,
     & ' centroid x,y,z=',3f8.3,' wrmsedt=',f7.3)
      write(16,53)
   53 format(1x,3('sta ph  eq1  eq2  wt  resedt edtobs',2x))
      write(31,55)
   55 format(1x,'sta  ph  eq1  eq2  wt  ttmedt edtobs resedt  delta',
     2 '   x_cev   y_cev   z_cev')
      if(kttfor.ne.3) then
        write(16,57) (stn(istoed(1,j,nc)),
     &  rmkc(jobsed(1,j,nc),jeved(1,j,nc),nc),jeved(1,j,nc),
     2  jeved(2,j,nc),w(j),resedt(j),edtsec(j,nc),
     3  j=1,60)
      else
        write(16,47) (stn6(istoed(1,j,nc)),
     &  rmkc(jobsed(1,j,nc),jeved(1,j,nc),nc),jeved(1,j,nc),
     2  jeved(2,j,nc),w(j),resedt(j),edtsec(j,nc),
     3  j=1,60)
      endif
      write(16,56) kobsedt(nc)
   56 format(' ------ Only 1st 60 edtobs printed here.  For all',
     2 i5,' obs, see File 31 (edtres.out) -----')
   57 format(3(1x,a4,a4,2i3,f5.2,f7.3,f7.3,1x))
   47 format(3(1x,a6,a4,2i3,f5.2,f7.3,f7.3,1x))
   58 format(1x,a4,1x,a4,2i3,f5.2,3f7.2,f9.2,3f8.2)
   48 format(a6,1x,a4,2i3,f5.2,3f7.2,f9.2,3f8.2)
      write(16,1601)
 1601 format(/)
      do 30 j=1,kobsedt(nc)
        j1=jobsed(1,j,nc)
        j2=jobsed(2,j,nc)
        je1=jeved(1,j,nc)
        je2=jeved(2,j,nc)
        if(kttfor.ne.3) then
          write(31,58) stn(istoc(j1,je1,nc)),rmkc(j1,je1,nc),je1,je2,
     2    w(j),ttmedt(j),edtsec(j,nc),resedt(j),dltac(j1,je1,nc),
     3   (cevc(ie,je1,nc),ie=1,3)
        else
          write(31,48) stn(istoc(j1,je1,nc)),rmkc(j1,je1,nc),je1,je2,
     2    w(j),ttmedt(j),edtsec(j,nc),resedt(j),dltac(j1,je1,nc),
     3    (cevc(ie,je1,nc),ie=1,3)
        endif
   30 continue
      return
c ****** end of subroutine outresedt ******
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine outresrdt(ne)
c  write out station residuals and weights (from medder) for each event
c
c  common block variables:
      include 'simul2014_common.inc'
      common/wtpars/w(mxobsa)
c
      write(16,3) ne,iyrmo(ne),iday(ne),ihr(ne),mino(ne),seco(ne),
     2 rmag(ne)
      write(21,3) ne,iyrmo(ne),iday(ne),ihr(ne),mino(ne),seco(ne),
     2 rmag(ne)
    3 format(1h0,'residuals for rec-pair diff-time: event=',i5,
     2 3x,a4,a2,1x,a2,i2,1x,f5.2,'  mag ',f4.2)
      write(16,52)
   52 format(1x,3('sta1 ph  sta2  wt  resrdt rdtobs',2x))
      write(21,56)
   56 format(1x,'sta1  ph  sta2  wt  ttmrdt rdtobs resrdt  delta',
     2 4x,'z ev')
      if(kttfor.ne.3) then
      write(16,53) (stn(isto(jobsrd(1,j,ne),ne)),rmk(jobsrd(1,j,ne),ne),
     2 stn(isto(jobsrd(2,j,ne),ne)),w(j),resrdt(j),rdtsec(j,ne),
     3 j=1,kobsrdt(ne))
      else
      write(16,43)(stn6(isto(jobsrd(1,j,ne),ne)),rmk(jobsrd(1,j,ne),ne),
     2 stn(isto(jobsrd(2,j,ne),ne)),w(j),resrdt(j),rdtsec(j,ne),
     3 j=1,kobsrdt(ne))
      endif
   53 format(3(1x,a4,a4,1x,a4,f5.2,f7.3,f7.2,1x))
   43 format(3(1x,a6,a4,1x,a4,f5.2,f7.3,f7.2,1x))
   54 format(1x,a4,1x,a4,1x,a4,f5.2,3f7.2,2f8.2)
   44 format(1x,a6,1x,a4,1x,a4,f5.2,3f7.2,2f8.2)
      write(16,1601)
 1601 format(/)
      do 30 j=1,kobsrdt(ne)
        j1=jobsrd(1,j,ne)
        j2=jobsrd(2,j,ne)
        if(kttfor.ne.3) then
         write(21,54) stn(isto(j1,ne)),rmk(j1,ne),stn(isto(j2,ne)),w(j),
     2     ttmrdt(j),rdtsec(j,ne),resrdt(j),dlta(j1,ne),evc(3,ne)
        else
         write(21,44)stn6(isto(j1,ne)),rmk(j1,ne),stn(isto(j2,ne)),w(j),
     2     ttmrdt(j),rdtsec(j,ne),resrdt(j),dlta(j1,ne),evc(3,ne)
        endif
   30 continue
      return
c ****** end of subroutine outresrdt ******
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine parsep(nc,ne,kopt,nwr)
c  routine to perform parameter separation using qr
c  kopt =1 for event
c  kopt =2 for receiver-pair differential time "shot"
c  kopt =3 for cluster event (treated as shot in velocity inversion)
c  kopt =4 for cluster earthquake-pair differential times
c  common block variables:
      include 'simul2014_common.inc'
      common/wtpars/w(mxobsa)
c
c  declaration statements:
      real c1(mxobsa,mxobsa),u(mxobsa),up
      parameter(zero=0.0,one=1.0)
c
      if(kopt.eq.1) nobs=kobs(ne)
      if(kopt.eq.2) nobs=kobsrdt(ne)
      if(kopt.eq.3) nobs=kobsce(ne,nc)
      if(kopt.eq.4) nobs=kobsedt(nc)
c
c  for shots, don't compute dtmp
      if (ne.gt.(neqs+nbls)) go to 100
      if(kopt.gt.1) goto 100
c
c  compute weighting
c
      nwr=0
      wnorm=0.0
      do 13 j=1,nobs
c  reading weight
         w(j)=wt(j,ne)
c  residual weighting
c  downweighting(linear) 0 to 98% res1 to res2, 98 to 100% res2 to res3
         ares=abs(res(j))
         if(ares.le.res2) then
            wr=1.0-(ares-res1)*dres12
            if (wr.gt.1.0) wr=1.0
         else
            if(res3.gt.res2) then
               wr=0.02-(ares-res2)*dres23
               if (wr.lt.0.0) wr=0.0
            else
               wr=0.0
            endif
         endif
c  distance weighting
         wd=1.0-(dlta(j,ne)-delt1)*ddlt
         if (wd.gt.1.0) wd=1.0
         if (wd.lt.0.0) wd=0.0
c  unnormalized weight
         w(j)=w(j)*wr*wd
c  19-jul-1983 Change normalizing factor
         wnorm=wnorm+w(j)
c        wnorm=wnorm+w(j)*w(j)
         if (w(j).gt.0.0) nwr=nwr+1
   13 continue
c  check to be sure 4 or more readings with nonzero weight
c only need 2 readings for blast
      if(ne.le.neqs) then
        if(nwr.ge.4) goto 12
        write(16,1612)
 1612   format(' ****** less than 4 readings with nonzero weights',
     2   /,20x,'SKIP TO NEXT EVENT *********')
      else
        if(nwr.ge.2) goto 12
        write(16,1613)
 1613   format(' ****** less than 2 readings with nonzero weights',
     2   /,20x,'SKIP TO NEXT EVENT *********')
      endif
      call outres(ne)
      return
c  normalize weights
c    19-jul-1983 Change in normalizing factor
c     wfac=sqrt(nwr/wnorm)
   12 wfac=nwr/wnorm
      if(ne.gt.netemp) wfac=wfac*wtsht
      resmax=zero
      do 14 j=1,nobs
         w(j)=w(j)*wfac
         ares=abs(res(j))
         if(ares.gt.resmax) resmax=ares
   14 continue
cDEP Save weight
      do 60 j=1,nobs
      wtcomb(j,ne)=w(j)
   60 continue
      if(resmax.gt.(3.0*res2)) call outres(ne)
      nwrt=nwrt+nwr
      if(nitmax.le.0) return
c
      do 10 i=1,nobs
         do 15 j=1,nobs
            c1(j,i)=zero
   15    continue
         c1(i,i)=w(i)
   10 continue
c  apply weighting to hypocenter matrix (1=o.t., 2,3,4=hypocenter)
      nsep=nparhy
c  for blast, only do origin time
      if(ne.gt.neqs) nsep=1
      do 17 i=1,nsep
         do 16 j=1,nobs
            dth(j,i)=dth(j,i)*w(j)
   16    continue
   17 continue
c
c  perform qr decomposition
c  prepare parameters for input to h12 routine
c  loop over the four columns of dth for decomposition
      do 20 i=1,nsep
         irow=i
         irow1=i+1
c  put column of dth into u-vector
         do 22 j=1,nobs
            u(j)=dth(j,i)
   22    continue
c  compute single transformation
         call h12(1,irow,irow1,nobs,u,1,up,dth,1,maxobs,nparhy,
     2      c1,1,maxobs,nobs)
   20 continue
c
c
c  qr decomposition complete - multiply residuals and medium matrix
c  by u0 matrix (q matrix minus its first four columns)
      nobs4=nobs-nsep
      do 35 j=1,nobs
         resj=res(j)
         do 34 i=1,nobs4
c  28 jan84 Changed from 4 to nsep to match CT. dmep
            i4=i+nsep
            resp(i)=resp(i)+c1(i4,j)*resj
   34    continue
   35 continue
c  compute separated medium matrix and store
      do 30 n=1,nobs
         nn=(n-1)*npari
         do 37 i=1,nobs4
            ii=(i-1)*npari
            i4=i+nsep
            c4n=c1(i4,n)
            do 36 m=1,npari
               im=ii+m
               nm=nn+m
               dtmp(im)=dtmp(im)+c4n*dtm(nm)
   36       continue
   37    continue
   30 continue
c
c  add to elements of g matrix
c
  100 continue
c
c  compute contribution to ssqr
      sqwtsh=sqrt(wtsht)
      sqwtte=sqrt(wttel)
      if(kopt.eq.2)sqwtrdt=sqrt(wtrdtsht)
      if(kopt.eq.4) sqwtedt=sqrt(wtepdt)
      totrms(ne)=0.0
      do 18 i=1,nobs
         resi=res(i)
         resi2=resi*resi
         totrms(ne)=totrms(ne)+resi2
         rrw=resi2*w(i)
         ssqrw=ssqrw+rrw
         wnobt=wnobt+w(i)
c**testprint*          write(16,1677) resi,w(i),ssqrw, wnobt
c**testprint* 1677 format( 'resi=',e15.3,',w(i)=',f5.2,',ssqrw=',e15.3,',wnobt=',
c**testprint*     2 e15.3)
         if(kopt.eq.1) isp=intsp(i,ne)
         if(kopt.eq.2)isp=intsp(jobsrd(1,i,ne),ne)
         if(kopt.eq.3) isp=intspc(i,ne,nc)
         if(kopt.eq.4) isp=intsped(i,nc)
         if(isp.eq.0) then
c  P arrival
            ssqrwp=ssqrwp+rrw
            wnobtp=wnobtp+w(i)
         else
c  S arrival
            ssqrws=ssqrws+rrw
            wnobts=wnobts+w(i)
         endif
         if ((kopt.eq.1).and.(ne.gt.netemp)) resi=resi*sqwtsh
         if ((kopt.eq.1).and.(ne.gt.nebs)) resi=resi*sqwtte
         if(kopt.eq.2) resi=resi*sqwtrdt
         if(kopt.eq.4) resi=resi*sqwtedt
         ssqr=ssqr+resi*resi
   18 continue
      totrms(ne)=sqrt(totrms(ne)/float(nobs))
c
c  loop over all observations (in separated form)
      n1=1
      n2=nobs4
      if (ne.gt.(neqs+nbls)) n2=nobs
      if(kopt.gt.1) n2=nobs
c     write(16,359)n1,n2,ne,neqs
c359  format(' parsep: n1,n2,ne,neqs= ',4i8)
c  loop over observations
      do 140 i=n1,n2
c        write(16,364)i,n1,n2
c364     format(' parsep: i,n1,n2= ',3i8)
         ii=(i-1)*npari
c  for a given observation, loop over all nodes
         do 130 j=1,npari
            jj=j+ii
            if (dtmp(jj).eq.zero) go to 130
c  add unknown for first observation only
            if (khit(mdexfx(j)).gt.0) go to 125
            mbl=mbl+1
            index(j)=mbl
C      WRITE(56,5601) j,mbl,mdexfx(j),index(j)
 5601 format('parsep j,mbl,mdexfx(j),index(j) ',4i6)
c  zero next row of g matrix
            mbtot=nbtot+1
            nbtot=nbtot+mbl
            do 120 l=mbtot,nbtot
               g(l)=0.0
  120       continue
  125       khit(mdexfx(j))=khit(mdexfx(j))+1
  126       k=index(j)
            kk2=((k-1)*k)/2
c  build rhs
            rhs(k)=rhs(k)+resp(i)*dtmp(jj)
c  build g matrix
            do 128 l1=1,j
               lj=l1+ii
               if (dtmp(lj).eq.zero) go to 128
               l=index(l1)
               m=l+kk2
               if (k.lt.l) m=k+((l-1)*l)/2
               g(m)=g(m)+dtmp(jj)*dtmp(lj)
c              if(g(m).lt.0.0) write(16,1620) g(m),m,dtmp(jj),
c    2            jj,dtmp(lj),lj,k,l,l1,i,ii,j
c1620          format(' g(m)=',f6.3,'m=',i6,',dtmp(jj)=',f6.3,
c    2            'jj=',i6,',dtmp(lj)=',f6.3,'lj=',i6,',k=',i4,
c    3            ',l=',i4,',l1=',i4,',i=',i4,',ii=',i6,',j=',i4)
  128       continue
  130    continue
  140 continue
c  zero out dtm and dtmp matrices for next event
      nono=nobs*npari
      seprms(ne)=zero
      do 45 m=1,nono
         dtm(m)=zero
         dtmp(m)=zero
   45 continue
      do 945 i=1,nobs
         seprms(ne)=seprms(ne)+resp(i)*resp(i)
         resp(i)=zero
  945 continue
      rmssep=sqrt(seprms(ne)/float(nobs))
      if(ne.gt.(neqs+nbls)) then
         write(16,1698) totrms(ne),rmssep
 1698    format('  ** total rms =',f8.5,
     2      '  **  model rms =',f8.5,' **')
      else
         write(16,1699) nwr,totrms(ne),rmssep
 1699    format(5x,'nwr=',i3,', ','** total rms =',f8.5,
     2      '  **  model rms =',f8.5,' **')
      endif
      if(totrms(ne).ge.(2.0*res1)) then
        if(kopt.eq.1) call outres(ne)
        if(kopt.eq.2) call outresrdt(ne)
        if(kopt.eq.3) call outresce(nc,ne)
        if(kopt.eq.4) call outresedt(nc)
      endif
      return
c***** end of subroutine parsep *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine path(ne,no,xe,ye,ze,ttime)
c  this routine determines the minimum-time ray path
c  in two steps:  first, an approximate path is
c  determined using approximate ray tracing
c  then the approximate path is used as a starting point in
c  shooting ray tracing routine to determine the path in 3-d
c  ***  note - current version does not contain full 3-d ray
c  ***  tracing - routines are under development
c
c  declaration statements:
c
c  common block variables:
      include 'simul2014_common.inc'
c
c  set up parameters to call ray tracing routine
c
c  S or P reading
      isp=intsp(no,ne)
c  receiver coordinates
      ns=isto(no,ne)
      xr=stc(1,ns)
      yr=stc(2,ns)
      zr=stc(3,ns)
c  determine 3-d path using circular raypaths.
c  determine approximate path
  120 jflag=jfl
      ncrold=nco(no,ne)
      ndpold=ndo(no,ne)
c
      call rayweb(0,ne,no,isp,xe,ye,ze,xr,yr,zr,
     * fstime,jflag,ncrold,ndpold)
      ttime=fstime
      ttc(no)=ttime
      nco(no,ne)=ncrold
      ndo(no,ne)=ndpold
c
c  do pseudo-bending if i3d>0
      if (i3d.le.0) return
c
  125 continue
c  number of pb iter depends on distance
      nitpbu=nitpb(1)
      if(rdlta(no,ne).gt.delt1) nitpbu=nitpb(2)
      call minima(0,ne,no,isp,ttime,nitpbu,jpb)
      if(jpb.lt.nitpbu) goto 130
        write(26,2601) no,stn(ns),jpb,nitpbu,stn6(ns)
 2601   format(' Minima: no=',i4,2x,a4,', used maximum ',
     2  'number PB iter.: j=',i3,', nitpb=',i3,1x,a6)
  130 continue
      ttc(no)=ttime
      if(kout3.eq.0) return
c  write out travel-time differences to file 19
      tdif=fstime-ttime
      if(kttfor.ne.3) then
        write(19,1900) ne,stn(isto(no,ne)),dlta(no,ne),fstime,ttime,tdif
      else
        write(19,1901)ne,stn6(isto(no,ne)),dlta(no,ne),fstime,ttime,tdif
      endif
 1900 format(i4,1x,a4,f7.2,3f8.4)
 1901 format(i4,1x,a6,f7.2,3f8.4)
c
      return
c***** end of subroutine path *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine pathce(nc,ne,no,xe,ye,ze,ttime)
c  this routine determines the minimum-time ray path
c  in two steps:  first, an approximate path is
c  determined using approximate ray tracing
c  then the approximate path is used as a starting point in
c  shooting ray tracing routine to determine the path in 3-d
c  ***  note - current version does not contain full 3-d ray
c  ***  tracing - routines are under development
c
c  declaration statements:
c
c  common block variables:
      include 'simul2014_common.inc'
c
c  set up parameters to call ray tracing routine
c
c  S or P reading
      isp=intspc(no,ne,nc)
c  receiver coordinates
      ns=istoc(no,ne,nc)
      xr=stc(1,ns)
      yr=stc(2,ns)
      zr=stc(3,ns)
c  determine 3-d path using circular raypaths.
c  determine approximate path
  120 jflag=jflc(ne)
      ncrold=ncoc(no,ne,nc)
      ndpold=ndoc(no,ne,nc)
c
      call rayweb(1,ne,no,isp,xe,ye,ze,xr,yr,zr,
     * fstime,jflag,ncrold,ndpold)
      ttime=fstime
      ttcce(no,ne)=ttime
      ncoc(no,ne,nc)=ncrold
      ndoc(no,ne,nc)=ndpold
c
c  do pseudo-bending if i3d>0
      if (i3d.le.0) return
c
  125 continue
c  number of pb iter depends on distance
      nitpbu=nitpb(1)
      if(rdltac(no,ne,nc).gt.delt1) nitpbu=nitpb(2)
      call minima(1,ne,no,isp,ttime,nitpbu,jpb)
      if(jpb.lt.nitpbu) goto 130
        write(26,2601) no,stn(ns),jpb,nitpbu,stn6(ns)
 2601   format(' Minima: no=',i4,2x,a4,', used maximum ',
     2  'number PB iter.: j=',i3,', nitpb=',i3,1x,a6)
  130 continue
      ttcce(no,ne)=ttime
      if(kout3.eq.0) return
c  write out travel-time differences to file 19
      tdif=fstime-ttime
      if(kttfor.ne.3) then
        write(19,1900) nc,ne,stn(istoc(no,ne,nc)),rdltac(no,ne,nc),
     2   fstime,ttime,tdif
 1900 format(i4,i5,1x,a4,f7.2,3f8.4)
      else
        write(19,1901) nc,ne,stn6(istoc(no,ne,nc)),rdltac(no,ne,nc),
     2   fstime,ttime,tdif
      endif
 1901 format(i4,i5,1x,a6,f7.2,3f8.4)
c
      return
c***** end of subroutine pathce *****
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine qtravel
c
      common/temp/xtemp(260),ytemp(260),ztemp(260),rtemp(260),ttemp(260)
c
      common/pathm/x(260),y(260),z(260),v(260),vq(260),tra,qtra,n,nn
c
      qtra=0
      do 60 i=2,n
         i1=i-1
         xd=xtemp(i)-xtemp(i1)
         yd=ytemp(i)-ytemp(i1)
         zd=ztemp(i)-ztemp(i1)
         ds=sqrt(xd*xd+yd*yd+zd*zd)
         tv=ds*(1.0/vq(i)+1.0/vq(i1))
         qtra=qtra+tv
  60  continue
      qtra=0.5*qtra
c
      return
c ***** end of subroutine qtravel *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine rayweb(kopt,ne,no,isp,xe,ye,ze,xr,yr,zr,
     *            fstime,jflag,ncrold,ndpold)
c  approximate ray tracing package art2
c    with fast ray tracing code
c     by Cliff Thurber (from his simul3l version)
c
c if kopt=0, then standard earthquake
c    kopt=1, then cluster earthquake
c  common block variables:
      include 'simul2014_common.inc'
c
c  declaration statements:
c  parameters
      real xe,ye,ze,xr,yr,zr,fstime
c  local variables
      real delx,dely,delz,sep,delsep,pthsep(260),strpth(780),
     * fstpth(780),dipvec(3,9),disvec(780,9),trpath(780,9),
     * trtime(9),tmin,tt,trpth1(780),disvcm(780,9)
      integer nd,ns,npt,ncr,i,ic,ndp,nc,n1,n2,n3,nn,np,ndpfst
c
c      print *,'Rayweb: no,xe,ye,ze,xr,yr,zr:',no,xe,ye,ze,xr,yr,zr
c  compute source-receiver separation
      delx=xr-xe
      dely=yr-ye
      delz=zr-ze
      sep=sqrt(delx*delx+dely*dely+delz*delz)
c  determine integer parameters for set of curves to be constructed
      call setup(sep,scale1,scale2,nd,ns,npt,ncr,n2exp)
c
c  set up pthsep array for straight-line travel time calculation
      sn1=1.0/float(ns)
      delsep=sep*sn1
      do 20 i=1,ns
         pthsep(i)=delsep
   20 continue
c
c  determine points along straight-line path
      xstep=delx*sn1
      ystep=dely*sn1
      zstep=delz*sn1
c
      ic=0
      ns1=ns+1
      do 25 ii=1,ns1
c
         i=ii-1
c
         ic=ic+1
         strpth(ic)=xe+xstep*i
         fstpth(ic)=strpth(ic)
         ic=ic+1
         strpth(ic)=ye+ystep*i
         fstpth(ic)=strpth(ic)
         ic=ic+1
         strpth(ic)=ze+zstep*i
         fstpth(ic)=strpth(ic)
   25 continue
c
c  compute travel time along straight-line path
      call ttime(isp,ns,npt,strpth,pthsep,fstime)
c
c only straight rays for i3d = -1
      if(i3d.eq.-1) goto 65
c
      if (ncr.eq.1) go to 65
c
c  compute the dip vectors of length scale2
      call cmpdpv(xe,ye,ze,xr,yr,zr,scale2,ndip,dipvec)
c
c  compute the basic set of displacement vectors
      call cmpdsv(ndip,iskip,ns,dipvec,disvec)
      call cmdsvm(ndip,iskip,ns,dipvec,disvcm)
c
c  set first and last points of all trial paths to source and receiver
      n1=3*npt-2
      n2=n1+1
      n3=n1+2
      ndip1=1+iskip
      ndip2=ndip-iskip
c
c  fast ray tracing code
c
      nz1=nz-1
      ncr0=1
      ncr1=ncr-1
      if (jflag.eq.0) go to 28
      ncr0=ncrold-1
      if (ncr0.lt.1) ncr0=1
      ncr1=ncrold+1
      if (ncr1.gt.ncr) ncr1=ncr
c  ndip was changed (in main) for ihomo iterations., make ndpold=vertical plane.
c     if(jfl.eq.1) ndpold=(ndip+1)/2     ! 21-feb-86, this is done in strt, so unnecessary here
      ndip1=ndpold-1
      if (ndip1.lt.1+iskip) ndip1=1+iskip
      ndip2=ndpold+1
      if (ndip2.gt.ndip-iskip) ndip2=ndip-iskip
   28 continue
c  set "old" values to straight line
      ncrold=0
      ndpold=(ndip+1)/2
c
c
      do 30 ndp=ndip1,ndip2
         trpath(1,ndp)=xe
         trpath(2,ndp)=ye
         trpath(3,ndp)=ze
         trpath(n1,ndp)=xr
         trpath(n2,ndp)=yr
         trpath(n3,ndp)=zr
   30 continue
      trpth1(1)=xe
      trpth1(2)=ye
      trpth1(3)=ze
      trpth1(n1)=xr
      trpth1(n2)=yr
      trpth1(n3)=zr
c
c  loop over the curve sets
      do 40 nc=ncr0,ncr1
         iz0=0
c
c  loop over different dips for one set of curves
         do 42 ndp=ndip1,ndip2
c
            npt2=npt-2
c  loop to determine points along one path
            do 44 np=1,npt2
               n1=3*np+1
               n3=n1+2
               do 43 nn=n1,n3
                  trpath(nn,ndp)=nc*disvec(nn,ndp)+strpth(nn)
                  trpth1(nn)=trpath(nn,ndp)
   43          continue
   44       continue
c
c  set up pthsep array for travel time calculations
            if (ndp.eq.ndip1) call cmpsep(trpth1,pthsep,ns)
c  compute travel time along one path
            call ttime(isp,ns,npt,trpth1,pthsep,tt)
            trtime(ndp)=tt
   42    continue
c
c  sort through trtime to find fastest path from current set
         tmin=1.0e15
         do 50 ndp=ndip1,ndip2
            if (trtime(ndp).gt.tmin) go to 50
            tmin=trtime(ndp)
            ndpfst=ndp
   50    continue
c
c  compare fastest trtime to current value of fstime
c  replace fstime and fstpth if needed
         if (tmin.ge.fstime) go to 40
c         print *,'fstime=',fstime
         fstime=tmin
c  reset "old" values
         ncrold=nc
         ndpold=ndpfst
c
         npt3=3*npt
         do 52 np=1,npt3
            fstpth(np)=trpath(np,ndpfst)
   52    continue
c
   40 continue
c
      if(i3d.lt.3) goto 146
c      print *,'Loop with the more squashed arc (third-root)'
c  Do another loop with the more squashed arc (third-root)
c  that may be better initial path for long ray paths
c  loop over the curve sets
      do 140 nc=ncr0,ncr1
         iz0=0
c
c  loop over different dips for one set of curves
         do 142 ndp=ndip1,ndip2
c
            npt2=npt-2
c  loop to determine points along one path
            do 144 np=1,npt2
               n1=3*np+1
               n3=n1+2
               do 143 nn=n1,n3
                  trpath(nn,ndp)=nc*disvcm(nn,ndp)+strpth(nn)
                  trpth1(nn)=trpath(nn,ndp)
  143          continue
  144       continue
c
c  set up pthsep array for travel time calculations
            if (ndp.eq.ndip1) call cmpsep(trpth1,pthsep,ns)
c  compute travel time along one path
            call ttime(isp,ns,npt,trpth1,pthsep,tt)
            trtime(ndp)=tt
  142    continue
c
c  sort through trtime to find fastest path from current set
         tmin=1.0e15
         do 150 ndp=ndip1,ndip2
            if (trtime(ndp).gt.tmin) go to 150
            tmin=trtime(ndp)
            ndpfst=ndp
  150    continue
c
c  compare fastest trtime to current value of fstime
c         print *,'tmin from curve set=',tmin
c  replace fstime and fstpth if needed
         if (tmin.ge.fstime) go to 140
         fstime=tmin
cc_i3d=3      print *,'Rayweb: no,xe,ye,ze,xr,yr,zr:',no,xe,ye,ze,xr,yr,zr
cc_i3d=3         print *,'squashed arc fstime=',fstime
c  reset "old" values
         ncrold=nc
         ndpold=ndpfst
c
         npt3=3*npt
         do 152 np=1,npt3
            fstpth(np)=trpath(np,ndpfst)
  152    continue
c
  140 continue
  146 continue
c
c  put fstpth into rp array
   65 continue
      if(kopt.eq.0) then
      do 60 np=1,npt
         n3=3*np
         n1=n3-2
         n2=n1+1
         rp(1,np,no)=fstpth(n1)
         rp(2,np,no)=fstpth(n2)
         rp(3,np,no)=fstpth(n3)
   60 continue
      nrp(no)=npt
      else
      do 70 np=1,npt
        n3=3*np
        n1=n3-2
        n2=n1+1
        rpce(1,np,no,ne)=fstpth(n1)
        rpce(2,np,no,ne)=fstpth(n2)
        rpce(3,np,no,ne)=fstpth(n3)
   70 continue
      nrpce(no,ne)=npt
      endif
c
c***** end of subroutine rayweb *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine rescov
c  compute resolution and covariance
c  common block variables:
      include 'simul2014_common.inc'
c
      varnce=ssqrw/float(nobt-mbl-4*neqs)
c     write(16,1001) varnce
 1001 format(/,' RESCOV: data variance is ',f11.6)
c  resolution and error calculations
      if(ires.eq.2) rewind (17)
c  avoid printing zeroth element of rhs array
c  index point from nodes(input) to solution arrays(rhs)
      do 5 kc=1,npari
         if(khit(mdexfx(kc)).eq.0) index(kc)=npari+1
         if(hit(mdexfx(kc)).lt.hitct) index(kc)=npari+1
         if(index(kc).eq.0) index(kc)=npari+1
    5 continue
c  zero out drm and stderr elements
      do 300 l=1,npar
         drm(l)=0.0
         stderr(l)=0.0
  300 continue
c
   10 do 500 l=1,npari
         l1=l
         l1n=mdexfx(l1)
         if(khit(l1n).eq.0) goto 500
         if(hit(l1n).lt.hitct) goto 500
         k=index(l1)
         jj=0
         ii=0
c put appropriate elements of g1 into rhs
         do 435 i=1,mbl
            do 430 j=1,i
               jj=jj+1
               if (i.eq.k.or.j.eq.k) go to 420
               go to 430
  420          ii=ii+1
               rhs(ii)=g1(jj)
c CHECK
c              if(g1(jj).lt.-5.0) write(16,1620)
c    2            g1(jj),jj,ii,i,j,mbl,l,k,npar
c1620          format(' g1(jj)=',f9.3,'jj=',i6,',ii=',i4,',i=',i4,
c    2            ',j=',i4,',mbl=',i4,',l=',i4,',k=',i4,',npar=',i4)
  430       continue
  435    continue
c  compute vector of resolution matrix for this parameter
         call luelmp(g,rhs,mbl,rhs)
c print out full resolution matrix if ires=2 or 3
         if(ires.lt.2) goto 440
cfdmep write out number of free nodes and inverted nodes to unit 18 for
c later use with resolution matrix
         rewind (18)
         write(18,1800) npar,nparv,npari,nparvi,nparpi,nparsi,nrowp,
     2     nrows,nrowst
 1800    format(9i7,'   0.0001')
c
         write(17,1720) l1,mdexfx(l1)
 1720    format(/,'row=',i5,', gridpoint no.=',i5,/)
         write(17,1721) (rhs(index(kc)),kc=1,npari)
 1721    format(20f7.4)
c  store diagonal element
  440    drm(l1n)=rhs(k)
c  standard error and covariance calculation
         do 450 i=1,mbl
            rhs(i)=rhs(i)*varnce
  450    continue
c  compute vector of covariance matrix for this parameter
         call luelmp(g,rhs,mbl,rhs)
c  store standard error of slowness perturbation
         rhs(k)=abs(rhs(k))
         stderr(l1n)=sqrt(rhs(k))
  500 continue
  900 return
c***** end of subroutine rescov *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      function rnormal ( seed )
c*********************************************************************
c obtained from:http://people.sc.fsu.edu/~jburkardt/f77_src/normal/normal.html
c    9-Dec-2014  DMEP
c      function r4_normal_01 ( seed )
c
cc R4_NORMAL_01 returns a unit pseudonormal real R4.
c  Discussion:
c    The standard normal probability distribution function (PDF) has
c    mean 0 and standard deviation 1.
c
c  Licensing:
c    This code is distributed under the GNU LGPL license.
c  Modified:
c    06 August 2013
c  Author:
c    John Burkardt
c  Parameters:
c    Input/output, integer SEED, a seed for the random number generator.
c    Output, real R4_NORMAL_01, a sample of the standard normal PDF.
c
      implicit none

      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real r1
      real r2
c      real r4_normal_01
      real rnormal
      real r4_uniform_01
      integer seed
      real x

      r1 = r4_uniform_01 ( seed )
      r2 = r4_uniform_01 ( seed )
      x = sqrt ( -2.0E+00 * log ( r1 ) ) * cos ( 2.0E+00 * pi * r2 )

c      r4_normal_01 = x
      rnormal = x

      return
c***** end of function rnormal *****
      end

c*********************************************************************
c
      function r4_uniform_01 ( seed )
c
cc R4_UNIFORM_01 returns a unit pseudorandom R4.
c  Discussion:
c    This routine implements the recursion
c      seed = 16807 * seed mod ( 2^31 - 1 )
c      r4_uniform_01 = seed / ( 2^31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R4_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c    This code is distributed under the GNU LGPL license.
c  Modified:
c    17 July 2006
c  Author:
c    John Burkardt
c
c  Reference:
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R4_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer k
      integer seed
      real r4_uniform_01

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r4_uniform_01 = real ( dble ( seed ) * 4.656612875D-10 )

      return
c***** end of function r4_uniform_01 *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine setorg(ltdo,oltm,lndo,olnm)
c  this routine establishes the short distance conversion factors
c  given the origin of coordinates
c  the rotation angle is converted to radians also
c  common block variables:
      common/shortd/ xltkm,xlnkm,rota,nzco,xlt,xln,snr,csr,
     2 onorth,oeast,rotadg,cmerid,rlatdo,rlondo
c  local variables:
      double precision dlt1,dlt2,dxlt,drad,drlt
      parameter ( re=6378.163, ell=298.26)
      parameter (drad=1.7453292d-2)
      parameter (drlt=9.9330647d-1)
      rotadg=rota
      rota=rota*drad       !convert from degrees to radians
c
cfhek
cek correction for S and E: minutes must also be taken as negative!
cek
cek   EK  12.6.96
cek
cek        dxlt=dble(60.*ltdo+oltm)
cek        xln=60.*lndo+olnm
      dxlt=dble(60.*abs(ltdo)+oltm)
      xln=60.*abs(lndo)+olnm
      if(ltdo.lt.0) dxlt=-1.*dxlt
      if(lndo.lt.0) xln=-1.*xln
cek
      xlt=sngl(dxlt)
c  conversion factor for latitude
      dlt1=datan(drlt*dtan(dxlt*drad/60.d0))
      dlt2=datan(drlt*dtan((dxlt+1.d0)*drad/60.d0))
      del=sngl(dlt2-dlt1)
      r=re*(1.0-sngl(dsin(dlt1)**2)/ell)
      xltkm=r*del
c  conversion factor for longitude
      del=sngl(dacos(1.0d0-(1.0d0-dcos(drad/60.d0))*dcos(dlt1)**2))
      bc=r*del
      xlnkm=bc/sngl(dcos(dlt1))
      write(16,3001) xltkm,bc
 3001 format(5x,'short distance conversion factors',/,
     2       10x,'one min lat  ',f7.4,' km',/,
     3       10x,'one min lon  ',f7.4,' km',/)
c
c  convert coordinates with rotation cosines
      snr=sin(rota)
      csr=cos(rota)
c
c  convert to NZ Map Grid if nzco=4
      if(nzco.eq.0) return
c  convert to NZ Map Grid if nzco=4
      if(nzco.eq.4) then
        call llnzmg(ltdo,oltm,lndo,olnm,onorth,oeast)
        write(16,1601) nzco,onorth,oeast
 1601   format(' nzco=',i1,', use NZ Map Grid distance conversions',
     2  /,'    Origin: Easting(km)',f12.3,', Northing(km)',f12.3)
        return
      endif
c  convert to AK state Plane coords if nzco=2
      if(nzco.eq.2) then
        call ll2spc(ltdo,oltm,lndo,olnm,onorth,oeast)
        write(16,1602)nzco,onorth,oeast
 1602   format(' nzco=',i1,', use Alaska StatePlane Coords for ',
     2  'distance converstions',
     3  /,'    Origin: Easting(km)',f12.3,', Northing(km)',f12.3)
        return
      endif
c  use transverse mercator if nzco=0,1,3
c MAKE NZ negative
      if(nzco.eq.1) then
        ltdo=-1*abs(ltdo)
        lndo=-1*abs(lndo)
      endif
c make sure minutes are negative if degrees are
      if(ltdo.lt.0) oltm=-1.0*abs(oltm)
      if(lndo.lt.0) olnm=-1.0*abs(olnm)
      latdp=ltdo
      platm=oltm
      londp=lndo
      plonm=olnm
        if((nzco.le.1).or.(nzco.eq.3)) call tmllxy(latdp,platm,londp,
     2    plonm,ynorth,xeast,1)
      onorth=ynorth
      oeast=xeast
      WRITE(16,*) 'nzco=',nzco,' cmerid= ',cmerid
c      WRITE(6,*) 'setorg: ltdo,oltm ',ltdo,oltm,' lndo,olnm ',
c     2 lndo,olnm,' oeast,onorth ',oeast,onorth
        if(nzco.eq.0) write(16,1604)nzco,cmerid,onorth,oeast
        if(nzco.eq.1) write(16,1603)nzco,onorth,oeast
        if(nzco.eq.3) write(16,1604)nzco,cmerid,onorth,oeast
 1603   format(' nzco=',i1,', use NZTM2000 distance conversions',
     2  /,'    Origin: Easting(km)',g16.4,', Northing(km)',g16.4)
 1604  format(' nzco=',i1,', use TM, central meridian=',f13.4,
     2  /,'    Origin: Easting(km)',g16.4,', Northing(km)',g16.4)
      return
c***** end of subroutine setorg *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine tmllxy(latd,xlat,lond,xlon,ynorkm,xeaskm,iway)
c setup parameters for transverse mercator to use in
c utm_geo code
c    iway=1 to go from lat-lon to x-y
c        =2 to go from x-y to lat-lon
c not using UTM zones, but that could be done
c because simul code has west positive, the working longitudes
c are negative for east, and so central meridian is negative
c
      double precision dbcmer,dblon,dblat,dbxm,dbym
      common/shortd/ xltkm,xlnkm,rota,nzco,xlt,xln,snr,csr,
     2 onorth,oeast,rotadg,cmerid,rlatdo,rlondo
      iutmzn=0
      inorth=0
      if((nzco.eq.1).or.(rlatdo.lt.0.0)) inorth=1
      dbcmer=dble(-1.0*cmerid)
      if(iway.eq.1) then
        dblat=dble(latd)+dble(xlat)/60.D0
        dblon=dble(lond)+dble(xlon)/60.D0
        call utm_geo(dblon,dblat,dbxm,dbym,dbcmer,
     2   iutmzn,iway,inorth)
        ynorkm=dbym/1000.0
c seems that need this multiplied by -1.0 for simul coords
        xeaskm= -1.0*(dbxm/1000.0)
      endif
      if(iway.eq.2) then
        dbxm=dble(-1000.0*xeaskm)
        dbym=dble(1000.0*ynorkm)
        call utm_geo(dblon,dblat,dbxm,dbym,dbcmer,
     2   iutmzn,iway,inorth)
        latd=int(dblat)
        xlat=60.0*(dblat-float(latd))
        lond=int(dblon)
        xlon=60.0*(dblon-float(lond))
      endif
      return
c***** end of subroutine tmllxy *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine setup(sep,scale1,scale2,nd,ns,npt,ncr,n2exp)
c
c  parameters
      real sep,scale1,scale2
c
      integer nd,ns,npt,ncr
c
c  determine the number of path divisions - use scale1
      nd=1+nint(3.32193*log10(sep/scale1))
      if (nd.gt.n2exp) nd=n2exp
c  number of segments along path
      ns=2**nd
c  number of points on path
      npt=ns+1
c
c  determine the number of curves - use scale2
      ncr=1+0.5*sep/scale2
c
      if (sep.gt.scale1) return
      nd=0
      ns=1
      npt=2
      ncr=1
c
c***** end of subroutine setup *****
      return
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine sort(x,ixst,n)
c  indirect sort routine from Meissner&Organick p352
      dimension ixst(1),x(1)
      do 5 i=1,n
      ixst(i)=i
    5 continue
      n1=n-1
      do 10 j=1,n1
        next=ixst(j+1)
        do 20 i=j,1,-1
c       write(6,400) j,i,ixst(i)
  400   format(' j=',i4,', i=',i4,', ixst(i)=',i4)
          if(x(next).gt.x(ixst(i))) goto 9
          ixst(i+1)=ixst(i)
   20     continue
c***** end of subroutine sort *****
    9   ixst(i+1)=next
   10   continue
      return
      end
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine strt(nit)
c
c  declaration statements:
      parameter(zero=0.0)
c
c  common block variables:
      common/machin/ eta,tol
      include 'simul2014_common.inc'
c
      iarr=maxpar
      iarri=mxpari
c
      nwrt=zero
      nswrt=zero
      ssqr=zero
      ssqrw=zero
      ssqrwp=zero
      ssqrws=zero
      wnobt=zero
      wnobtp=zero
      wnobts=zero
      do 75 n=1,iarr
        cnode(n)='0'
        canode(n)='-'
   75 continue
      do 77 n=1,iarri
         index(n)=0
         jndex(n)=0
         rhs(n)=zero
   77 continue
c
c
c  reset fast art arrays?
      if (nit.gt.ihomo) go to 80
      do 790 m=1,maxev
         do 79 n=1,maxobs
            ndo(n,m)=(ndip+1)/2
   79    continue
  790 continue
c
      if (nit.ge.1) go to 80
      do 780 m=1,maxev
         do 78 n=1,maxobs
            nco(n,m)=1
   78    continue
  780 continue
c
c
   80 mbl=0
      nbtot=0
c  zero hit counter arrays
      do 25 j=1,iarr
         hit(j)=zero
         hitall(j)=zero
         khit(j)=0
   25 continue
c
c  zero out partial derivative arrays
      mxobs1=maxobs-1
      do 30 j=1,iarri
         jj2=maxobs*j
         jj1=jj2-mxobs1
         do 35 i=jj1,jj2
            dtm(i)=zero
            dtmp(i)=zero
   35    continue
   30 continue
      do 39 i=1,maxobs
         resp(i)=zero
   39 continue
c
      if (nit.gt.0) return
c
c  machine constants
      halfu=0.5
   50 temp1=1.0+halfu
      if (temp1.le.1.0) go to 100
      halfu=0.5*halfu
      go to 50
  100 eta=2.0*halfu
c
      temp2=0.5
  150 continue
      if (temp2.le.0.0) go to 200
      tol=temp2
      temp2=0.5*temp2
      go to 150
  200 continue
c
c     tol=tol/eta
c for unix make tol bigger
      tol=tol/(eta*eta)
c
      write(16,1055) eta,tol
 1055 format(/,'  computed machine constants eta and tol:',2e16.6)
c
      return
c***** end of subroutine strt *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine travel
c
      common/temp/xtemp(260),ytemp(260),ztemp(260),rtemp(260),ttemp(260)
c
      common/pathm/x(260),y(260),z(260),v(260),vq(260),tra,qtra,n,nn
c
      tra=0
      do 60 i=2,n
         i1=i-1
         xd=xtemp(i)-xtemp(i1)
         yd=ytemp(i)-ytemp(i1)
         zd=ztemp(i)-ztemp(i1)
         ds=sqrt(xd*xd+yd*yd+zd*zd)
         tv=ds*(1.0/v(i)+1.0/v(i1))
         tra=tra+tv
  60  continue
      tra=0.5*tra
c
      return
c ***** end of subroutine travel *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine ttime(isp,ns,npt,pathr,pthsep,tt)
c
c  travel time along path via trapezoidal rule integration
c
c  parameters
      real pathr(780),pthsep(260),tt
c
      integer ns,npt
c  local variables
      real vpt(260),x,y,z,v
c
      integer ip,np
c
c  loop over points along path to determine velocities
      ip=0
      do 10 np=1,npt
         ip=ip+1
         x=pathr(ip)
         ip=ip+1
         y=pathr(ip)
         ip=ip+1
         z=pathr(ip)
         call vel3eft(isp,x,y,z,v)
   10 vpt(np)=v
c
      tt=0.0
c  sum up travel time - use trapezoidal rule
      do 20 np=1,ns
         np1=np+1
c  Check for value outside defined area
         if((vpt(np).le.0.0).or.(vpt(np1).le.0.0)) goto 99
         tt=tt+pthsep(np)/(vpt(np)+vpt(np1))
   20 continue
      tt=2.0*tt
c
      return
   99 write(16,100) np,vpt(np),np1,vpt(np1)
  100 format(' **** ERROR IN TTIME ****, velocity le 0.0',
     2 /,'    np   vpt(np)   np1   vpt(np1)',/,1x,2(i5,f10.2))
      stop
c***** end of subroutine ttime *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine ttmder(ne,noe,kopt,ttime,nit,nnodej)
c
c  declaration statements:
      dimension in(8)
c
c  common block variables:
      common/weight/ wv(8),ip,jp,kp,kpg
      include 'simul2014_common.inc'
c
c   Modified 4-Apr-2014 to have S-P residuals get partials for 
c   both vp and vpvs.
c   This subroutine has been modified for Q inversion 13-March-1998.
c   Following Andreas Rietbrock's version.
c   Then (when iuseq=1), the second part of the model array
c   is used to store Vp*Qp (not Vp/Vs).
c   The solution is only for Q, not for Vp though.  
c   Thus the solution arrays that would otherwise be related
c   to the first part of the model array (Vp) are used.
c   Therefore in ttmder, vel3 is called with "1" to retrieve
c   the V*Q value; and vel3 is called with "0" to get the 
c   indices for the solution arrays.
c
c    kopt=0 only earthquake location
c    kopt=1 getting velocity partials
c    kopt=2 earthquake pair
c
c nnodej= counts nodes (no-fixed) sampled by this observation
      nnodej=0
c teleseismic event counter
      nt=ne-nebs
c
c  for receiver-pair dt, 2 paths to consider
      npath=1
      if(kopt.eq.2) npath=2
      no=noe
      do 65 ipath=1,npath
      if(kopt.eq.2) no=jobsrd(ipath,noe,ne)
c  resegment ray path
c  calculate travel time derivatives with respect to hypocentral
c  parameters - use coords of first two points on raypath
c  determine slowness at source
      xe=rp(1,1,no)
      ye=rp(2,1,no)
      ze=rp(3,1,no)
c  check for P or S reading
      isp=intsp(no,ne)
      if(kopt.eq.2) goto 40
      if(iuseq.eq.0) then
              call vel3eft(isp,xe,ye,ze,v)
      else
              call vel3eft(1,xe,ye,ze,v)
      endif
       
      us=1.0/v
c**  The program has been modified to include S-P times for the inversion.
c**  The corresponding travel-time derivatives are calculated.(06/29/92)
      if (isp.eq.1) then
              isp=0
              call vel3eft(isp,xe,ye,ze,vp)
              up=1.0/vp
              us=us-up
        isp=1
      endif
c**
c  determine cartesian derivatives from direction cosines
      dx=rp(1,2,no)-rp(1,1,no)
      dy=rp(2,2,no)-rp(2,1,no)
      dz=rp(3,2,no)-rp(3,1,no)
      ds=sqrt(dx*dx+dy*dy+dz*dz)
      uds=-us/ds
c  hypocentral derivatives
      dth(no,2)=uds*dx
      dth(no,3)=uds*dy
      dth(no,4)=uds*dz
c  origin time derivative
      if(iuse2t.eq.0) then
        dth(no,1)=1.0
      else
        nclock=iclock(isto(no,ne))
        if(nclock.eq.0) then
          dth(no,1)=1.0
          dth(no,5)=0.0
        else
          dth(no,1)=0.0
          dth(no,5)=1.0
        endif
      endif
c** S-P CHANGE
c   zero the origin time derivative for S-P
      if (isp.eq.1) then
        dth(no,1)=0.0
        dth(no,5)=0.0
      endif
c**
c  zero cartesian derivatives for quarry blasts
      if((ne.le.neqs).or.(ne.gt.(neqs+nbls))) goto 20
      dth(no,2)=0.0
      dth(no,3)=0.0
      dth(no,4)=0.0
   20 continue
c
c  skip next section if all velocity nodes are fixed (nparvi=0)
      if(nparvi.eq.0) goto 88
c  skip next section if only doing earthquake location
        if(kopt.eq.0) go to 88
c
   40 continue
c  travel time and velocity partial derivatives
      tt=0.0
      half=0.5
c  loop over segments comprising the ray path
      nrp1=nrp(no)-1
      pl(no)=0.0
      ibegrp=1
      if(kopt.eq.2) ibegrp=nrp(no)-nrprdt(ipath,noe,ne)+1
      do 60 i=ibegrp,nrp1
         i1=i+1
         rx=rp(1,i,no)
         ry=rp(2,i,no)
         rz=rp(3,i,no)
         dx=rp(1,i1,no)-rx
         dy=rp(2,i1,no)-ry
         dz=rp(3,i1,no)-rz
c  compute segment length
         sl=sqrt(dx*dx+dy*dy+dz*dz)
         pl(no)=pl(no)+sl
c  decide on number of subsegments and compute length
         nseg=nint(sl/stepl)+1
         fnsegi=1.0/float(nseg)
         ssl=sl*fnsegi
         dxs=dx*fnsegi
         dys=dy*fnsegi
         dzs=dz*fnsegi
c
         xp=rx-half*dxs
         yp=ry-half*dys
         zp=rz-half*dzs
c  loop over subsegments
         do 55 is=1,nseg
            xp=xp+dxs
            yp=yp+dys
            zp=zp+dzs
c
          if(iuseq.eq.0) then
               call vel3eft(isp,xp,yp,zp,v)
          else
c  v=V*Q, vv=V
               call vel3eft(1,xp,yp,zp,v)
               call vel3eft(0,xp,yp,zp,vv)
          endif
          dt=ssl/v
c
c For S-P, loop over vpvs and vp
          isp2=isp+1
          do 50 jder=1,isp2
c  The next section is a change from 'block' to 'linear'
c   partial derivatives, by C.Thurber,may10,1983.
c  Nodes with non-zero weight
c For jder=2, make dvp relate to S-P
            in(1)=ip-1+nx2*(jp-2)+nxy2*(kp-2)-nxy2*(2*(isp))
            if(jder.eq.2) in(1)=ip-1+nx2*(jp-2)+nxy2*(kpg-2)
            in(2)=in(1)+1
            in(3)=in(1)+nx2
            in(4)=in(3)+1
            in(5)=in(1)+nxy2
            in(6)=in(5)+1
            in(7)=in(5)+nx2
            in(8)=in(7)+1
C PRINT        IF((NO.EQ.1).OR.(ISP.EQ.1)) THEN
C PRINT          WRITE(56,5601) NO,isp,jder,IN,ip,jp,kp,kpg
C PRINT   5601   FORMAT('no=',i4,' isp=',i1,' jder=',i1,' in =',8i8,
C PRINT       2  ' ip jp kp kpg=',4i6)
C PRINT          WRITE(56,5602) i,nrp1,is,nseg,xp,yp,zp
C PRINT   5602   FORMAT('i,nrp1=',2i4,' is,nseg=',2i4,' xp,yp,zp=',3f8.2)
C PRINT        ENDIF
c
c  Assign zero weight to boundary nodes (these nodes are not 
c  included in the inversion, but are in the velocity array,
c  thus we want to avoid writing to negative or incorrect 
c  elements of the partial derivative matrix)
c
            if(ip.eq.1) then
c              write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
               wv(1)=0.0
               wv(3)=0.0
               wv(5)=0.0
               wv(7)=0.0
            else
               if(ip.eq.nx1) then
c                 write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
                  wv(2)=0.0
                  wv(4)=0.0
                  wv(6)=0.0
                  wv(8)=0.0
               end if
            endif
c
            if(jp.eq.1) then
c              write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
               wv(1)=0.0
               wv(2)=0.0
               wv(5)=0.0
               wv(6)=0.0
            else
               if(jp.eq.ny1) then
c                 write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
                  wv(3)=0.0
                  wv(4)=0.0
                  wv(7)=0.0
                  wv(8)=0.0
               endif
            endif
c
            if((kpg.eq.1).or.(kpg.eq.(nz1+1))) then
c              write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
               do 30 izg=1,4
                  wv(izg)=0.0
   30          continue
            else
               if((kpg.eq.nz1).or.(kpg.eq.(2*nz1))) then
c                 write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
                  do 35 izg=5,8
                     wv(izg)=0.0
   35             continue
               endif
            endif
 1610       format(' ASSIGNING ZERO WEIGHTS IN TTMDER',
     2         ' no=',i3,',xp=',f7.2,',yp=',f7.2,',zp=',
     3         f7.2,',v=',f5.3,',ip=',i2,',jp=',i2,',kpg=',i2)
c
c  Accumulate model partial derivatives
            do 48 kk=1,2
               kk1=kk-1
               do 47 jj=1,2
                  jj1=jj-1
                  do 46 ii=1,2
                     ii1=ii-1
                     ijk=ii+2*jj1+4*kk1
c skip boundary nodes
                     if(wv(ijk).lt.0.05) goto 46
c DEP
c write out DWS for all nodes (including fixed and linked) when nitmax=1
c (useful for planning fixed and linked)
c Include weight factor like for inversion
c Note that for shots on nit=0, combined weight, including residual
c wt'g is not yet calculated so use wtsht.
             if(nitmax.eq.1) then
               if((ne.gt.(neqs+nbls)).and.(nit.eq.0)) then
                 if(kopt.eq.2) then
                 hitall(in(ijk))=hitall(in(ijk))+wv(ijk)*wtsht*wtrdtsht
                 else if(ne.gt.nebs) then
                 hitall(in(ijk))=hitall(in(ijk))+wv(ijk)*wttel
                 else
                 hitall(in(ijk))=hitall(in(ijk))+wv(ijk)*wtsht
                 endif
               else
                 if(kopt.eq.2) then
                 hitall(in(ijk))=hitall(in(ijk))+wv(ijk)*wtcombrd(no,ne)
                 else
                 hitall(in(ijk))=hitall(in(ijk))+wv(ijk)*wtcomb(no,ne)
                 endif
               endif
             endif
c  skip fixed nodes
                     if(nfix(in(ijk)).eq.1) goto 46
                     ini=ndexfx(in(ijk))
                     nnodej=nnodej+1
c
c  start cht 1998
c  also save index of master for accumulating DWS(hit)
      inmas=in(ijk)
      if (imerge(in(ijk)).eq.1) then
        ini=ndexfx(jequal(in(ijk)))
        inmas=jequal(in(ijk))
      endif
c  end cht 1998
c
c check for writing to an array element that is outside of the inversion
c  solution array
                     if((ini.lt.1).or.(ini.gt.nparvi)) then
                        write(16,1606) ini,ijk,npari,nparvi
 1606                   format(' *** Error in TTMDER, accessing',
     2                     ' gridpoint outside of velocity inversion',
     3                     ' gridpoints, ini=',i5,', ijk=',i5,/,
     4                     22x,'Probably boundary gridpoints are',
     5                     ' too close (TTMDER tries to write DTM',
     6                     ' elements with wv >= 0.05',/,' npari=',
     7                     i7,', nparvi=',i7)
                        write(16,1603) ne,no,xp,yp,zp,v,ip,jp,kp,kpg
 1603                   format(' ne=',i5,', no=',i5,', xp=',f8.2,
     2                     ', yp=',f8.2,', zp=',f8.2,', v=',f8.3,/,
     3                     21x,'ip=',i6,',   jp=',i6,',   kp=',i6,
     4                     ',   kpg=',i6)
                        write(16,1607) (j,in(j),j,wv(j),j=1,8)
 1607                   format(' in(',i1,')=',i6,' wv(',i1,')=',e15.5)
                        write(16,1608)
 1608                   format(' * * * * STOP * * * * (to avoid',
     2                     ' writing outside of defined DTM array)')
                        stop
                     end if
                     inp=ini+(no-1)*npari
cDEP Now include weight factor for hit(DWS)
                     if((ne.gt.(neqs+nbls)).and.(nit.eq.0)) then
                       if(kopt.eq.2) then
                       hit(inmas)=hit(inmas)+wv(ijk)*wtsht*wtrdtsht
                       else if(ne.gt.nebs) then
                       hit(inmas)=hit(inmas)+wv(ijk)*wttel
                       else
                       hit(inmas)=hit(inmas)+wv(ijk)*wtsht
                       endif
                     else
                       if(kopt.eq.2) then
                       hit(inmas)=hit(inmas)+wv(ijk)*wtcombrd(no,ne)
                       else
                       hit(inmas)=hit(inmas)+wv(ijk)*wtcomb(no,ne)
                       endif
                     endif
c** S-P CHANGE
c   set up equation to solve for delta-Vp/Vs
                   if (isp.eq.0) then
                     if(iuseq.eq.0) then
                       dtmip = dt*wv(ijk)*
     &                     vel(ip+ii1,jp+jj1,kp+kk1)/v
                     else
                       call vel3eft(1,xp,yp,zp,vv)
                       dtmip =dt*wv(ijk)*
     &                     vel(ip+ii1,jp+jj1,kp+kk1)/v
                       call vel3eft(0,xp,yp,zp,vv)
                     endif
                   endif
                   if (isp.eq.1) then
                     isp=0
                     call vel3eft(isp,xp,yp,zp,vp)
                     if(jder.eq.2) then
                       vpvsi=vp/v
                       dtmip = dt*wv(ijk)*(vpvsi-1.0)*
     &                     vel(ip+ii1,jp+jj1,kp+kk1)/vp
                     else
                       dtmip =ssl*wv(ijk)/vp
                     endif
                     isp=1
C      PRINT *,'ip,jp,kp ',ip,jp,kp,' v,vv ',v,vv,' dtmip ',dtmip
                   endif
c  for receiver-pair differential time, res1-res2 
c    so obs2 has its dtmp subtracted
                   if(ipath.eq.2) dtmip= -1.0*dtmip
                   dtm(inp)=dtm(inp)+dtmip
c***
   46             continue
   47          continue
   48       continue
   50     continue
c
   55    continue
   60 continue
c End of path partial derivative loop
   65 continue
c for receiver-pair dt, indiv. sta. residuals already calculated
      if(kopt.eq.2) goto 91
c  arrival time residual
  88    continue
      is=isto(no,ne)
      if(nt.le.0) goto 89
c get tele residual
      if(isp.eq.0) then
        res(no)=secte(no,nt)-telpad(nt)-ttime
      else
        intsp(no,ne)=0
        call path(ne,no,xe,ye,ze,ptime)
        intsp(no,ne)=1
        smptime=ttime-ptime
        res(no)=secte(no,nt)-telpads(nt)-smptime
      endif
      goto 90
   89 if((iuse2t.eq.0).or.(nclock.eq.0)) then
        res(no)=secp(no,ne)-seco(ne)-ttime-pdl(is)
      else
        res(no)=secp(no,ne)-seco2(ne)-ttime-pdl(is)
      endif
      ttm(no)=ttime
c** S-P CHANGE
c**     Estimate the calculated P and S times and then the S-P time
c       Then estimate the residuals between the observed and calculated.
      if (isp.eq.1) then
        intsp(no,ne)=0
        call path(ne,no,xe,ye,ze,ptime)
        intsp(no,ne)=1
        smptime=ttime-ptime
        res(no)=secp(no,ne)-smptime - sdl(is)
        ttm(no)=smptime
      endif
c**
   90 if(kopt.le.1) goto 93
c  Compute receiver-pair dt residual
   91 continue
      ttmrdt(noe)=ttm(jobsrd(1,noe,ne))-ttm(jobsrd(2,noe,ne))
      resrdt(noe)=rdtsec(noe,ne)-ttmrdt(noe)
c      if((kout2.lt.2).and.(nnodej.gt.0)) then
c        if(noe.eq.1) write(16,1609) ne
c 1609    format('ttmder resid for rec-pair diff-time: event=',i5,/,
c     *   ' Sta1  ttm1  res1    Sta2  ttm2  res2 : ttmrdt ',
c     *   ' rdtsec   resrdt')
c        j1=jobsrd(1,noe,ne)
c        j2=jobsrd(2,noe,ne)
c        write(16,1611)stn(isto(j1,ne)),ttm(j1),res(j1),stn(isto(j2,ne)),
c     *   ttm(j2),res(j2),ttmrdt(noe),rdtsec(noe,ne),resrdt(noe)
c 1611   format(2(1x,a4,2f6.2),2x,2f6.2,f7.2)
c      endif
c
   93  continue
c compute tele path delay derivatives
      if((ntel.eq.0).or.(ne.le.nebs)) goto 96
c full index
      ncor=nparv+2*nsts*invdel+nt
      if(intsp(no,ne).ne.0) ncor=ncor+ntel
      ini=ndexfx(ncor)
      inp=(no-1)*npari+ini
      dtm(inp)=1.0
c      hit(ncor)=hit(ncor)+1
      if(nit.eq.0) then
        hit(ncor)=hit(ncor)+wttel
      else
        hit(ncor)=hit(ncor)+wtcomb(no,ne)
      endif
c      PRINT *,'ttmder ne,nt,ntel,nebs',ne,nt,ntel,nebs
c      PRINT *,'       ncor,ini,hit(ncor)',ncor,ini,hit(ncor),
c     2 ' inp,dtm(inp)',inp,dtm(inp)
c
   96 if(invdel.eq.0)return
c  Calculate station delay partial derivatives (unless doing location only)
      if(kopt.eq.0) return
      do 100 i=1,npath
        no=noe
        if(kopt.eq.2) no=jobsrd(i,noe,ne)
        is=isto(no,ne)
        if(nfixst(is).eq.1) return
        ncor=nparv+is
c  If the observation is an S-wave, put the derivative nsts further up.
        if(intsp(no,ne).ne.0)ncor=ncor+nsts
        ini=ndexfx(ncor)
c   The derivative array, dtm(), is a long, linear array with all
c   of the observations and their derivatives end to end.
c   hit is an array which is used for this observation only. Khit
c   is set in subroutine parsep.
c   Each observation occupies npar locations, which = 1+(no-1)*npar to
c   no*npar
c   npar=nodes2, or nodes2+2*nsts when stn delays are being inverted, too.
        inp=(no-1)*npari+ini
        dtm(inp)=1.0
        if(i.eq.2) dtm(inp) = -1.0
        hit(ncor)=hit(ncor)+1.0
  100 continue
      return
c***** end of subroutine ttmder *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine ttmderce(nc,ne,noe,kopt,ttime,nit,nnodej)
c
c  declaration statements:
      dimension in(8)
c
c  common block variables:
      common/weight/ wv(8),ip,jp,kp,kpg
      include 'simul2014_common.inc'
c
c  This version is for cluster events, different variable names
c  but code is not really changed much
c    kopt=0 only earthquake location
c    kopt=1 getting velocity partials
c    kopt=2 earthquake pair
c
c   This subroutine has been modified for Q inversion 13-March-1998.
c   Following Andreas Rietbrock's version.
c   Then (when iuseq=1), the second part of the model array
c   is used to store Vp*Qp (not Vp/Vs).
c   The solution is only for Q, not for Vp though.  
c   Thus the solution arrays that would otherwise be related
c   to the first part of the model array (Vp) are used.
c   Therefore in ttmder, vel3 is called with "1" to retrieve
c   the V*Q value; and vel3 is called with "0" to get the 
c   indices for the solution arrays.
c
c
c nnodej= counts nodes (no-fixed) sampled by this observation
      nnodej=0
c
c  for earthquake-pair dt, 2 paths to consider
      npath=1
      if(kopt.eq.2) npath=2
      no=noe
      do 65 ipath=1,npath
      if(kopt.eq.2) no=jobsed(ipath,noe,nc)
c  resegment ray path
c  calculate travel time derivatives with respect to hypocentral
c  parameters - use coords of first two points on raypath
c  determine slowness at source
      xe=rpce(1,1,no,ne)
      ye=rpce(2,1,no,ne)
      ze=rpce(3,1,no,ne)
c  check for P or S reading
      isp=intspc(no,ne,nc)
      if(kopt.eq.2) goto 40
      if(iuseq.eq.0) then
              call vel3eft(isp,xe,ye,ze,v)
      else
              call vel3eft(1,xe,ye,ze,v)
      endif
       
      us=1.0/v
c**  The program has been modified to include S-P times for the inversion.
c**  The corresponding travel-time derivatives are calculated.(06/29/92)
      if (isp.eq.1) then
              isp=0
              call vel3eft(isp,xe,ye,ze,vp)
              up=1.0/vp
              us=us-up
        isp=1
      endif
c**
c  determine cartesian derivatives from direction cosines
      dx=rpce(1,2,no,ne)-rpce(1,1,no,ne)
      dy=rpce(2,2,no,ne)-rpce(2,1,no,ne)
      dz=rpce(3,2,no,ne)-rpce(3,1,no,ne)
      ds=sqrt(dx*dx+dy*dy+dz*dz)
      uds=-us/ds
c  hypocentral derivatives
      dthc(no,2,ne)=uds*dx
      dthc(no,3,ne)=uds*dy
      dthc(no,4,ne)=uds*dz
c  origin time derivative
c only use true clock stations for cluster events
c  and thus dthc array has 4 hyp param, not 5
c      if(iuse2t.eq.0) then
        dthc(no,1,ne)=1.0
c      else
c        nclock=iclock(isto(no,ne))
c        if(nclock.eq.0) then
c          dthc(no,1,ne)=1.0
c          dthc(no,5,ne)=0.0
c        else
c          dthc(no,1,ne)=0.0
c          dthc(no,5,ne)=1.0
c        endif
c      endif
c** S-P CHANGE
c   zero the origin time derivative for S-P
      if (isp.eq.1) then
        dthc(no,1,ne)=0.0
c        dthc(no,5,ne)=0.0
      endif
c
c don't need this for cluster earthquakes
cc  zero cartesian derivatives for quarry blasts
c      if((ne.le.neqs).or.(ne.gt.(neqs+nbls))) goto 20
c      dthc(no,2,ne)=0.0
c      dthc(no,3,ne)=0.0
c      dthc(no,4,ne)=0.0
c   20 continue
c
c  skip next section if all velocity nodes are fixed (nparvi=0)
      if(nparvi.eq.0) goto 88
c  skip next section if only doing earthquake location
        if(kopt.eq.0) go to 88
        if(nitmax.le.0) goto 88
c
   40 continue
c  travel time and velocity partial derivatives
      tt=0.0
      half=0.5
c  loop over segments comprising the ray path
      nrp1=nrpce(no,ne)-1
      plce(no,ne)=0.0
      ibegrp=1
      if(kopt.eq.2) nrp1=nrpedt(ipath,noe,nc)-1
      do 60 i=ibegrp,nrp1
         i1=i+1
         rx=rpce(1,i,no,ne)
         ry=rpce(2,i,no,ne)
         rz=rpce(3,i,no,ne)
         dx=rpce(1,i1,no,ne)-rx
         dy=rpce(2,i1,no,ne)-ry
         dz=rpce(3,i1,no,ne)-rz
c  compute segment length
         sl=sqrt(dx*dx+dy*dy+dz*dz)
         plce(no,ne)=plce(no,ne)+sl
c  decide on number of subsegments and compute length
         nseg=nint(sl/stepl)+1
         fnsegi=1.0/float(nseg)
         ssl=sl*fnsegi
         dxs=dx*fnsegi
         dys=dy*fnsegi
         dzs=dz*fnsegi
c
         xp=rx-half*dxs
         yp=ry-half*dys
         zp=rz-half*dzs
c  loop over subsegments
         do 55 is=1,nseg
            xp=xp+dxs
            yp=yp+dys
            zp=zp+dzs
c
          if(iuseq.eq.0) then
               call vel3eft(isp,xp,yp,zp,v)
          else
c  v=V*Q, vv=V
               call vel3eft(1,xp,yp,zp,v)
               call vel3eft(0,xp,yp,zp,vv)
          endif
          dt=ssl/v
c
c For S-P, loop over vpvs and vp
          isp2=isp+1
          do 50 jder=1,isp2
c  The next section is a change from 'block' to 'linear'
c   partial derivatives, by C.Thurber,may10,1983.
c  Nodes with non-zero weight
c For jder=2, make dvp relate to S-P
            in(1)=ip-1+nx2*(jp-2)+nxy2*(kp-2)-nxy2*(2*isp)
            if(jder.eq.2) in(1)=ip-1+nx2*(jp-2)+nxy2*(kpg-2)
            in(2)=in(1)+1
            in(3)=in(1)+nx2
            in(4)=in(3)+1
            in(5)=in(1)+nxy2
            in(6)=in(5)+1
            in(7)=in(5)+nx2
            in(8)=in(7)+1
c
c  Assign zero weight to boundary nodes (these nodes are not 
c  included in the inversion, but are in the velocity array,
c  thus we want to avoid writing to negative or incorrect 
c  elements of the partial derivative matrix)
c
            if(ip.eq.1) then
c              write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
               wv(1)=0.0
               wv(3)=0.0
               wv(5)=0.0
               wv(7)=0.0
            else
               if(ip.eq.nx1) then
c                 write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
                  wv(2)=0.0
                  wv(4)=0.0
                  wv(6)=0.0
                  wv(8)=0.0
               end if
            endif
c
            if(jp.eq.1) then
c              write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
               wv(1)=0.0
               wv(2)=0.0
               wv(5)=0.0
               wv(6)=0.0
            else
               if(jp.eq.ny1) then
c                 write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
                  wv(3)=0.0
                  wv(4)=0.0
                  wv(7)=0.0
                  wv(8)=0.0
               endif
            endif
c
            if((kpg.eq.1).or.(kpg.eq.(nz1+1))) then
c              write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
               do 30 izg=1,4
                  wv(izg)=0.0
   30          continue
            else
               if((kpg.eq.nz1).or.(kpg.eq.(2*nz1))) then
c                 write(16,1610) no,xp,yp,zp,v,ip,jp,kpg
                  do 35 izg=5,8
                     wv(izg)=0.0
   35             continue
               endif
            endif
 1610       format(' ASSIGNING ZERO WEIGHTS IN TTMDER',
     2         ' no=',i3,',xp=',f7.2,',yp=',f7.2,',zp=',
     3         f7.2,',v=',f5.3,',ip=',i2,',jp=',i2,',kpg=',i2)
c
c  Accumulate model partial derivatives
            do 48 kk=1,2
               kk1=kk-1
               do 47 jj=1,2
                  jj1=jj-1
                  do 46 ii=1,2
                     ii1=ii-1
                     ijk=ii+2*jj1+4*kk1
c skip boundary nodes
                     if(wv(ijk).lt.0.05) goto 46
c DEP
c write out DWS for all nodes (including fixed and linked) when nitmax=1
c (useful for planning fixed and linked)
c Include weight factor like for inversion
           if(nitmax.eq.1) then
             if(kopt.eq.2) then
               hitall(in(ijk))=hitall(in(ijk))+wv(ijk)*wtcombed(noe,nc)
             else
               hitall(in(ijk))=hitall(in(ijk))+wv(ijk)*wtcombc(no,ne,nc)
             endif
           endif
c  skip fixed nodes
                     if(nfix(in(ijk)).eq.1) goto 46
                     ini=ndexfx(in(ijk))
                     nnodej=nnodej+1
c
c  start cht 1998
c  also save index of master for accumulating DWS(hit)
      inmas=in(ijk)
      if (imerge(in(ijk)).eq.1) then
        ini=ndexfx(jequal(in(ijk)))
        inmas=jequal(in(ijk))
      endif
c  end cht 1998
c
c check for writing to an array element that is outside of the inversion
c  solution array
                     if((ini.lt.1).or.(ini.gt.nparvi)) then
                        write(16,1606) ini,ijk,npari,nparvi
 1606                   format(' *** Error in TTMDER, accessing',
     2                     ' gridpoint outside of velocity inversion',
     3                     ' gridpoints, ini=',i5,', ijk=',i5,/,
     4                     22x,'Probably boundary gridpoints are',
     5                     ' too close (TTMDER tries to write DTM',
     6                     ' elements with wv >= 0.05',/,' npari=',
     7                     i7,', nparvi=',i7)
                        write(16,1603) ne,no,xp,yp,zp,v,ip,jp,kp,kpg
 1603                   format(' ne=',i5,', no=',i5,', xp=',f8.2,
     2                     ', yp=',f8.2,', zp=',f8.2,', v=',f8.3,/,
     3                     21x,'ip=',i6,',   jp=',i6,',   kp=',i6,
     4                     ',   kpg=',i6)
                        write(16,1607) (j,in(j),j,wv(j),j=1,8)
 1607                   format(' in(',i1,')=',i6,' wv(',i1,')=',e15.5)
                        write(16,1608)
 1608                   format(' * * * * STOP * * * * (to avoid',
     2                     ' writing outside of defined DTM array)')
                        stop
                     end if
                     inp=ini+(no-1)*npari
cDEP Now include weight factor for hit(DWS)
                     if((ne.gt.(neqs+nbls)).and.(nit.eq.0)) then
                       if(kopt.eq.2) then
                       hit(inmas)=hit(inmas)+wv(ijk)*wtsht*wtrdtsht
                       else
                       hit(inmas)=hit(inmas)+wv(ijk)*wtsht
                       endif
                     else
                       if(kopt.eq.2) then
                       hit(inmas)=hit(inmas)+wv(ijk)*wtrdtsht
                       else
                       hit(inmas)=hit(inmas)+wv(ijk)*wtcomb(no,ne)
                       endif
                     endif
c** S-P CHANGE
c   set up equation to solve for delta-Vp/Vs
                   if (isp.eq.0) then
                     if(iuseq.eq.0) then
                       dtmip = dt*wv(ijk)*
     &                     vel(ip+ii1,jp+jj1,kp+kk1)/v
                     else
                       dtmip =dt*wv(ijk)*
     &                     vel(ip+ii1,jp+jj1,kp+kk1)/v
                     endif
                   endif
                   if (isp.eq.1) then
                     isp=0
                     call vel3eft(isp,xp,yp,zp,vp)
                     if(jder.eq.2) then
                       vpvsi=vp/v
                       dtmip = dt*wv(ijk)*(vpvsi-1.0)*
     &                     vel(ip+ii1,jp+jj1,kp+kk1)/vp
                     else
                       dtmip =ssl*wv(ijk)/vp
                     endif
                     isp=1
                   endif
c  for earthquake-pair differential time, res1-res2 
c    so obs2 has its dtmp subtracted
                   if(ipath.eq.2) dtmip= -1.0*dtmip
                   dtm(inp)=dtm(inp)+dtmip
c***
   46             continue
   47          continue
   48       continue
   50     continue
c
   55   continue
   60 continue
c End of path partial derivative loop
   65 continue
c for earthquake-pair dt, indiv. sta. residuals already calculated
      if(kopt.eq.2) goto 90
c  arrival time residual
  88    continue
      is=istoc(no,ne,nc)
      resc(no,ne)=secpc(no,ne,nc)-secoce(ne,nc)-ttime-pdl(is)
      ttmc(no,ne)=ttime
c** S-P CHANGE
c**     Estimate the calculated P and S times and then the S-P time
c       Then estimate the residuals between the observed and calculated.
      if (isp.eq.1) then
        intspc(no,ne,nc)=0
        call pathce(nc,ne,no,xe,ye,ze,ptime)
        intspc(no,ne,nc)=1
        smptime=ttime-ptime
        resc(no,ne)=secpc(no,ne,nc)-smptime - sdl(is)
        ttmc(no,ne)=smptime
      endif
c**
      if(kopt.le.1) goto 93
c  Compute earthquake-pair dt residual
   90 continue
      job1=jobsed(1,noe,nc)
      job2=jobsed(2,noe,nc)
      jev1=jeved(1,noe,nc)
      jev2=jeved(2,noe,nc)
      if(intsped(noe,nc).eq.0) 
     &  edtsec(noe,nc)=(secpc(job1,jev1,nc)-secoce(jev1,nc))
     &  -(secpc(job2,jev2,nc)-secoce(jev2,nc))
      ttmedt(noe)=ttmc(job1,jev1)-ttmc(job2,jev2)
      resedt(noe)=edtsec(noe,nc)-ttmedt(noe)
c Print Check
c      if((kout2.lt.2).and.(nnodej.gt.0)) then
c        if(noe.eq.1) write(16,1609) nc
c 1609    format('ttmder resid for eq-pair diff-time: cluster=',i5,/,
c     *   ' Sta1  ttm1  res1    Sta2  ttm2  res2 : ttmedt ',
c     *   ' edtsec   resedt')
c        write(16,1611)stn(istoed(1,noe,nc)),ttmc(job1,jev1),
c     *   resc(job1,jev1),stn(istoed(2,noe,nc)),ttmc(job2,jev2),
c     *   resc(job2,jev2),ttmedt(noe),edtsec(noe,nc),resedt(noe)
c 1611   format(2(1x,a4,2f6.2),2x,2f6.2,f7.2)
c      endif
c
  93  continue
      if(invdel.eq.0)return
c  Calculate station delay partial derivatives (unless doing location only)
c  Eq-pair diff times (to common station) do not contribute to sta delays
      if((kopt.ne.1).or.(nitmax.le.0)) return
        if(nfixst(is).eq.1) return
        ncor=nparv+is
c  If the observation is an S-wave, put the derivative nsts further up.
        if(intspc(no,ne,nc).ne.0) ncor=ncor+nsts
        ini=ndexfx(ncor)
c   The derivative array, dtm(), is a long, linear array with all
c   of the observations and their derivatives end to end.
c   hit is an array which is used for this observation only. Khit
c   is set in subroutine parsep.
c   Each observation occupies npar locations, which = 1+(no-1)*npar to
c   no*npar
c   npar=nodes2, or nodes2+2*nsts when stn delays are being inverted, too.
        inp=(no-1)*npari+ini
        dtm(inp)=1.0
        if(i.eq.2) dtm(inp) = -1.0
        hit(ncor)=hit(ncor)+1.0
  100 continue
      return
c***** end of subroutine ttmderce *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine vel3(isp,x,y,z,v)
c  This routine is Cliff Thurber's
c  common block variables:
      common/weight/ wv(8),ip,jp,kp,kpg
      include 'simul2014_common.inc'
c
c  use Prothero's intmap here
      call intmap(x,y,z,ip,jp,kp)
c
      ip1=ip+1
      jp1=jp+1
      kp1=kp+1
c
c      write(16,100)x,xl,ip,y,yl,jp,z,zl,kp
c100 format(3(2f7.3,i3))
      xf=(x-xn(ip))/(xn(ip1)-xn(ip))
      yf=(y-yn(jp))/(yn(jp1)-yn(jp))
      zf=(z-zn(kp))/(zn(kp1)-zn(kp))
      xf1=1.0-xf
      yf1=1.0-yf
      zf1=1.0-zf
c
      wv(1)=xf1*yf1*zf1
      wv(2)=xf*yf1*zf1
      wv(3)=xf1*yf*zf1
      wv(4)=xf*yf*zf1
      wv(5)=xf1*yf1*zf
      wv(6)=xf*yf1*zf
      wv(7)=xf1*yf*zf
      wv(8)=xf*yf*zf
c  calculate velocity
c  S-velocity is stored after P-velocity
c  (or V*Q if iuseq=1)
      kpg=kp
      if(isp.eq.1) kp=kp+nz
      kp1=kp+1
      v=wv(1)*vel(ip,jp,kp)+wv(2)*vel(ip1,jp,kp)
     2 +wv(3)*vel(ip,jp1,kp)+wv(4)*vel(ip1,jp1,kp)
     * +wv(5)*vel(ip,jp,kp1)+wv(6)*vel(ip1,jp,kp1)
     * +wv(7)*vel(ip,jp1,kp1)+wv(8)*vel(ip1,jp1,kp1)
      return
c***** end of subroutine vel3 *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine vel3eft(isp,x,y,z,v)
c This version of vel3 has earth-flattening transformation from
c Buland and Chapman (1983)
c  This routine is Cliff Thurber's
c  common block variables:
      common/weight/ wv(8),ip,jp,kp,kpg
      include 'simul2014_common.inc'
c
      re=6371.0
c  use Prothero's intmap here
      call intmap(x,y,z,ip,jp,kp)
c
      ip1=ip+1
      jp1=jp+1
      kp1=kp+1
c
c	write(16,100)x,xl,ip,y,yl,jp,z,zl,kp
c100 format(3(2f7.3,i3))
      xf=(x-xn(ip))/(xn(ip1)-xn(ip))
      yf=(y-yn(jp))/(yn(jp1)-yn(jp))
      zf=(z-zn(kp))/(zn(kp1)-zn(kp))
      xf1=1.0-xf
      yf1=1.0-yf
      zf1=1.0-zf
c
      wv(1)=xf1*yf1*zf1
      wv(2)=xf*yf1*zf1
      wv(3)=xf1*yf*zf1
      wv(4)=xf*yf*zf1
      wv(5)=xf1*yf1*zf
      wv(6)=xf*yf1*zf
      wv(7)=xf1*yf*zf
      wv(8)=xf*yf*zf
c  calculate velocity
c  S-velocity is stored after P-velocity
c  (or V*Q if iuseq=1)
      kpg=kp
      if(isp.eq.1) kp=kp+nz
      kp1=kp+1
      v=wv(1)*vel(ip,jp,kp)+wv(2)*vel(ip1,jp,kp)
     2 +wv(3)*vel(ip,jp1,kp)+wv(4)*vel(ip1,jp1,kp)
     * +wv(5)*vel(ip,jp,kp1)+wv(6)*vel(ip1,jp,kp1)
     * +wv(7)*vel(ip,jp1,kp1)+wv(8)*vel(ip1,jp1,kp1)
      v=(re*v)/(re-z)
      return
c***** end of subroutine vel3eft *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine veladj(nit)
c  routine to set up matrices for velocity inversion
c
c  declaration statements:
      integer nm(5),mbla(5)
      real r2(5),rhm(5),rhsvar(5),sqnorm(5),snrm(5),c1(5),dampc1(5)
      save c1
      character*11 solnam(5)
      parameter(zero=0.0)
c
c  common block variables:
      include 'simul2014_common.inc'
c
      data solnam/'Velocity   ','P-Velocity ','Vp/Vs ratio',
     2 'Sta. Corr. ','Telpath dly'/
      if(iuseq.eq.1) solnam(2)='  Qp * Vp  '
c
      open(96,file ='veladj.print',form='formatted')
c  compute adjustments to velocity model
c
        write(16,385)
 385    format(/,' Compute velocity adjustments to model-Veladj')
c  remove unknowns with fewer than nmin observations
c  mbl was defined in PARSEP as no. rows (ie: no. velocity gridpts.
c    observed), in VELADJ, it's adjusted by the <hitct gridpts.
      kbl=mbl
      if (hitct.eq.0.0) go to 230  !skip next section if hitcutoff not used
c
      WRITE(52,5203) mbla
 5203 format(' veladj start mbla=',5i8)
      do 10 i=1,4
         mbla(i)=0
   10 continue
      mbl=0
      ii=0
      jj=0
c  construct inverse mapping for index
c  index points from nodes(input) to solution arrays (rhs,g), while 
c  jndex points in the other direction.
c       loop over all sampleable nodes.
c  twice as many nodes if both P and S velocity
      do 170 i=1,npari
c  g has no equations for Unsampled nodes
         if (khit(mdexfx(i)).eq.0) go to 170
         jndex(index(i))=i
  170 continue
c       now remove equations from g for poorly sampled nodes
c  jj=counter for old g array, ii=counter for new g array
      do 220 i=1,kbl
         if(khit(mdexfx(jndex(i))).eq.0) goto 180
         if(hit(mdexfx(jndex(i))).ge.hitct) goto 190
  180    jj=jj+i
         go to 220
  190    mbl=mbl+1
c  reorder rhs and index
         rhs(mbl)=rhs(i)
         index(jndex(i))=mbl
         do 210 j=1,i
            jj=jj+1
            if (khit(mdexfx(jndex(j))).eq.0) goto 210
            if(hit(mdexfx(jndex(j))).ge.hitct) goto 200
            go to 210
  200       ii=ii+1
            g(ii)=g(jj)
  210    continue
  220 continue
c  revise jndex array for reduced g array
      do 225 i=1,npari
         if(khit(mdexfx(i)).eq.0) goto 225
         if(hit(mdexfx(i)).lt.hitct) goto 225
         jndex(index(i))=i
  225 continue
c  add damping (different for Vp and Vs) to diagonal elements of g
  230 j=0
      WRITE(52,5202) nparvs, nodes2
 5202 format('nparvs=',i8, ' nodes2=',i8)
      do 240 i=1,mbl
         j=j+i
         kv=(mdexfx(jndex(i))-1)/nodes2+1
      WRITE(52,5201) i,j,jndex(i),mdexfx(jndex(i)),kv
 5201 format('veladj i j jndex(i):',3i8,' mdexfx(jndex(i)),kv:',i8,i4)
c  Correction for vp & sta corr only (w/out vp/vs)
         if((kv.eq.2).and.(iuses.eq.1)) kv=3
         if(kv.gt.3) kv=3
c tele path delay damp in vdamp(4)
         if(mdexfx(jndex(i)).gt.nparvs) kv=4
         g(j)=g(j)+vdamp(kv)
  240 continue
c  if desired store g and rhs for resolution and error calculation
      if (ires.eq.0) go to 272
c  store .rhs
      do 260 n=1,mbl
         rhs1(n)=rhs(n)
  260 continue
      ng=mbl*(mbl+1)/2
      do 270 n=1,ng
         g1(n)=g(n)
  270 continue
c  subtract off vdamp from diagonal of g1
  272 j=0
      do 280 i=1,mbl
         j=j+i
         kv=(mdexfx(jndex(i))-1)/nodes2+1
         if((kv.eq.2).and.(iuses.eq.1)) kv=3
         if(kv.gt.3) kv=3
         if(mdexfx(jndex(i)).gt.nparvs) kv=4
         mbla(kv)=mbla(kv)+1
         if(ires.gt.0) g1(j)=g1(j)-vdamp(kv)
C      PRINT *,'veladj i,jndex(i),mdexfx(jndex(i)),kv',
C     2  i,jndex(i),mdexfx(jndex(i)),kv,' khit, hit ',
C     3  khit(mdexfx(jndex(i))),hit(mdexfx(jndex(i)))
  280 continue
c
  290 ier=0
      write(16,20) kbl,mbl,(solnam(i+1),mbla(i),i=1,4)
      write(36,20) kbl,mbl,(solnam(i+1),mbla(i),i=1,4)
  20  format(' total nodes observed=',i7,', nodes inverted for=',i7,
     2 3(',',1x,a11,':',i7))
cfhdmep 
c  Save number of Vp, Vp/Vs and station correction nodes for later 
c  use with resolution matrix output
      nrowp=mbla(1)
      nrows=mbla(2)
      nrowst=mbla(3)
      nrowte=mbla(4)
c
c  perform lu decomposition
      call ludecp(g,g,mbl,ier)
      if (ier.ne.0) write(16,1001)
 1001 format('  *** error in ludecp ***')
      if(ier.eq.129) write(16,1610)
 1610 format('  ier=129, matrix A is algorithmically not positive ',
     2 'definite')
c  perform lu elimination
      call luelmp(g,rhs,mbl,rhs)
c  inversion complete
c
c  calculate solution vector statistics
c  mean of rhs
      ivadj=0
      if(nparvi.eq.0) goto 36
      do 22 i=1,5
         rhm(i)=0.0
         nm(i)=0
   22 continue
      do 23 n=1,npari
         if(khit(mdexfx(n)).eq.0) goto 23
         if(hit(mdexfx(n)).lt.hitct) goto 23
         k=(mdexfx(n)-1)/nxy2 + 2
         ia=2
         if(k.ge.nz) ia=3
         if(n.gt.nparvi) ia=4
         if(n.gt.nparvsi) ia=5
         rhm(ia)=rhm(ia)+rhs(index(n))
         nm(ia)=nm(ia)+1
C      WRITE(56,*)'  mdexfx(n),index(n),ia,rhs',
C     2  mdexfx(n),index(n),ia,rhs(index(n))
         if(ia.gt.3) goto 23
         rhm(1)=rhm(1)+rhs(index(n))
         nm(1)=nm(1)+1
   23 continue
      do 24 i=1,5
      WRITE(56,*)' i,rhm(sum),nm',i,rhm(i),' ',nm(i)
         if(nm(i).gt.0.0) rhm(i)=rhm(i)/float(nm(i))
         rhsvar(i)=0.0
         sqnorm(i)=0.0
   24 continue
      do 25 n=1,npari
         if(khit(mdexfx(n)).eq.0) goto 25
         if(hit(mdexfx(n)).lt.hitct) goto 25
         k=(mdexfx(n)-1)/nxy2 + 2
         ia=2
         if(k.ge.nz) ia=3
         if(n.gt.nparvi) ia=4
         if(n.gt.nparvsi) ia=5
         sqnorm(ia)=sqnorm(ia)+rhs(index(n))*rhs(index(n))
         r2(ia)=rhs(index(n))-rhm(ia)
         rhsvar(ia)=rhsvar(ia)+r2(ia)*r2(ia)
         if(ia.gt.3) goto 25
         sqnorm(1)=sqnorm(1)+rhs(index(n))*rhs(index(n))
         r2(1)=rhs(index(n))-rhm(1)
         rhsvar(1)=rhsvar(1)+r2(1)*r2(1)
   25 continue
      write(16,1624) nit
      write(36,1624) nit
 1624 format(' iteration',i3,' rhs solution vector:')
      do 27 i=1,5
         if(nm(i).le.0) goto 27 
         snrm(i)=sqrt(sqnorm(i))
         rhsvar(i)=rhsvar(i)/float(nm(i))
         write(16,1625) solnam(i),rhm(i),rhsvar(i),sqnorm(i),snrm(i)
         write(36,1625) solnam(i),rhm(i),rhsvar(i),sqnorm(i),snrm(i)
   27 continue
 1625 format(9x,a11,': mean',f13.9,', variance',
     2 e14.4,', norm squared',e14.4,', norm',e14.4)
c check whether velocity solution norm is below cutoff value
      if(snrm(1).lt.snrmct) then
        nit=99
        write(16,1630) snrm(1),snrmct
 1630   format(/,' **** STOP since snrm',e14.4,' is below',
     2  ' snrmct',e14.4,' ****',/)
        return
      endif
c
c  apply adjustments to velocity model
      do 30 n=1,nparvi
c  gridpoint number in array of all nodes
         nn=mdexfx(n)
         vadj(nn)=zero
         if(khit(nn).eq.0) goto 30
         if(hit(nn).lt.hitct) goto 30
         rh=rhs(index(n))
c  calculate x and z indices of velocity grid
         k=(nn-1)/nxy2+2
         j=2+(nn-1+(2-k)*nxy2)/nx2
         i=1+nn+nx2*(2-j)+nxy2*(2-k)
         if(k.ge.nz) k=k+2      ! if s velocity node
       if(iuseq.eq.0) then
           delm=-vel(i,j,k)*rh/(1.0+rh)
       else
           delm=-vel(i,j,k+nz)*rh/(1.0+rh)
       endif
c**   Modification for S-P times
c    equations solve for delta-Vp/Vs directly
         if(k.ge.nz) delm=rh
         delma=abs(delm)
c  place upper bound on velocity adjustment
         if(k.gt.nz) then
            dvmax=dvsmx
         else
            dvmax=dvpmx
         end if
       if(iuseq.eq.0) then
            if (delma.gt.dvmax) delm=dvmax*delm/delma
            vadj(nn)=delm
       else
          delmaa = delma/vel(i,j,k)
            if (delmaa.gt.dvmax) delm=dvmax*vel(i,j,k)*delm/delma
            vadj(nn)=delm
       end if
          
c  apply adjustment to model
c** S-P CHANGE
         if (k.le.nz) then
          if(iuseq.eq.0) then
            vel(i,j,k)=vel(i,j,k)+delm
            if(vel(i,j,k).lt.cvpmin) vel(i,j,k)=cvpmin
            if(vel(i,j,k).gt.cvpmax) vel(i,j,k)=cvpmax
          else
            vold=vel(i,j,k+nz)
            vel(i,j,k+nz)=vel(i,j,k+nz)+delm
            qold=qval(i,j,k)
            qval(i,j,k)=vel(i,j,k+nz)/vel(i,j,k)
            if(qval(i,j,k).le.0.0) then
                qval(i,j,k) = 1.0
                vel(i,j,k+nz) = vel(i,j,k) * qval(i,j,k)
                vadj(nn) = vel(i,j,k+nz)-vold
            end if
            qadj(nn)=qval(i,j,k)-qold
          endif
       endif
         if(k.ge.nz) then
            vpvs(i,j,k-nz)=vpvs(i,j,k-nz)+delm
         endif
         ivadj=ivadj+1
   30 continue
c
c  start cht 1998
c
c  adjust linked gridpoints
c **** First adjust all nodes linked by constant adjustment
      do 32 n=1,nparv
      if (imerge(n).eq.0) goto 32
c  calculate x and z indices of "master" velocity grid
         nma=jequal(n)
         if(ltype(nma).ne.1) goto 32
c
c      write(16,1616) n
c 1616 format(' applying adjustment to linked node ',i4)
c
c  calculate x-y-z indices of linked velocity grid
         k=(n-1)/nxy2+2
         j=2+(n-1+(2-k)*nxy2)/nx2
         i=1+n+nx2*(2-j)+nxy2*(2-k)
         if(k.ge.nz) k=k+2      ! if s velocity node
c
         km=(nma-1)/nxy2+2
         jm=2+(nma-1+(2-km)*nxy2)/nx2
         im=1+nma+nx2*(2-jm)+nxy2*(2-km)
         if(km.ge.nz) km=km+2      ! if s velocity node
c
c
c  if constant link
c ****      if (ltype(nma).eq.1) then
c      write(16,2626) i,j,k,n,im,jm,km,nma
c 2626 format('  setting value for constant-type linked node: ',
c     &  4i4,2x,4i4)
c
           vadj(n)=vadj(nma)
c           vold=vel(i,j,k)
c **** Let's use the constant perturbation. Maybe better than
c **** constant velocity in some cases eg: vertical linking
c
c **** make constant change to v*Q, then apply to get updated Q
       if(iuseq.eq.0) then
           if (k.lt.nz) then
c ****           vel(i,j,k)=vel(im,jm,km)
c        write(16,2628) vold,vadj(n),vel(i,j,k)
c 2628   format(' Old vel=',f6.2,' vadj=',f6.2,' updated vel=',f6.2)
            vel(i,j,k)=vel(i,j,k)+vadj(n)
            if(vel(i,j,k).lt.cvpmin) vel(i,j,k)=cvpmin
            if(vel(i,j,k).gt.cvpmax) vel(i,j,k)=cvpmax
           else
c ****           vpvs(i,j,k-nz)=vpvs(im,jm,km-nz)
            vpvs(i,j,k-nz)=vpvs(i,j,k-nz)+vadj(n)
           endif
       else
c ****           vel(i,j,k+nz)=vel(im,jm,km+nz)
           vold=vel(i,j,k+nz)
           vel(i,j,k+nz)=vel(i,j,k+nz)+vadj(n)
           qold=qval(i,j,k)
           qval(i,j,k)=vel(i,j,k+nz)/vel(i,j,k)
C           if(qval(i,j,k).le.0.0) then
           if(qval(i,j,k).le.qmin) then
c             qval(i,j,k) = 1.0
             qval(i,j,k) = qmin
             vel(i,j,k+nz) = vel(i,j,k) * qval(i,j,k)
             vadj(n) = vel(i,j,k+nz)-vold
           end if
           qadj(n)=qval(i,j,k)-qold
       endif
   32 continue
      
c
c **** Now that all invert and constant-link have been updated,
c Update the linearly-linked nodes
c  Do iterative update to get correct dv for multiple linked nodes
c   set 10 iterations
      nitup=10
      do 440 itup=1,nitup
      do 33 n=1,nparv
      if (imerge(n).eq.0) goto 33
c  calculate x and z indices of "master" velocity grid
         nma=jequal(n)
         if(ltype(nma).eq.1) goto 33
c
c      write(16,1616) n
c  calculate x-y-z indices of linked velocity grid
         k=(n-1)/nxy2+2
         j=2+(n-1+(2-k)*nxy2)/nx2
         i=1+n+nx2*(2-j)+nxy2*(2-k)
         if(k.ge.nz) k=k+2      ! if s velocity node
c
         km=(nma-1)/nxy2+2
         jm=2+(nma-1+(2-km)*nxy2)/nx2
         im=1+nma+nx2*(2-jm)+nxy2*(2-km)
         if(km.ge.nz) km=km+2      ! if s velocity node
c
c  if linear link
       il=i
       jl=j
       kl=k
c      write(16,2627) i,j,k,n,im,jm,km,nma,il,jl,kl
c 2627 format('  setting value for linear-type linked nodes: ',
c     &  4i4,2x,4i4,2x,3i4)
c ** 17-jun-2010 Do NOT use multiple linked nodes
c Allow multiple linked nodes between inversion nodes
c Note that linear adjustment is linear for the parameter that
c is solved for in the inversion v*Q (not linear for Q)
c
c simul2010 allows diagonal linear linking and simply gives
c  it the average of the individual directions
c check for diagonal.  Do not use multiple for diagonal
      idiag=0
      if((il.ne.im).and.(jl.ne.jm).and.(kl.ne.km)) idiag=1
c
      if (k.lt.nz) then
      nvl=0
      vsum=0.0
      if(iuseq.eq.0) then
        vold=vel(i,j,k)
        kk=k
        kkm=km
        kkl=kl
      else
        kk=k+nz
        kkm=km+nz
        kkl=kl+nz
        vold=vel(i,j,k+nz)
      endif
c x direction
       if (i.ne.im) then
       ia=1
       if(i-im.lt.0) ia= -1
       il=i
  402  il=il+ia
       goto 404
c       if(idiag.eq.1) goto 404
c       if((il.eq.1).or.(il.eq.nx)) goto 404
c       nl=(k-2)*nxy2+(j-2)*nx2+(il-1)
c       if((imerge(nl).eq.1).and.(ltype(jequal(nl)).ne.1)) goto 402
  404  continue
        vnew=vel(im,jm,kkm)+(vel(il,jl,kkl)-vel(im,jm,kkm))*
     &   (xn(i)-xn(im))/(xn(il)-xn(im))
        nvl=nvl+1
        vsum=vsum+vnew
       endif
c
c y direction
       if (j.ne.jm) then
       ja=1
       if(j-jm.lt.0) ja= -1
       jl=j
  406  jl=jl+ja
       goto 408
c       if(idiag.eq.1) goto 408
c       if((jl.eq.1).or.(jl.eq.ny)) goto 408
c       nl=(k-2)*nxy2+(jl-2)*nx2+(i-1)
c       if((imerge(nl).eq.1).and.(ltype(jequal(nl)).ne.1)) goto 406
  408  continue
        vnew=vel(im,jm,kkm)+(vel(il,jl,kkl)-vel(im,jm,kkm))*
     &   (yn(j)-yn(jm))/(yn(jl)-yn(jm))
        nvl=nvl+1
        vsum=vsum+vnew
       endif
c
c z direction
       if (k.ne.km) then
       ka=1
       if(k-km.lt.0) ka= -1
       kl=k
  410  kl=kl+ka
       kkl=kkl+ka
       goto 412
c       if(idiag.eq.1) goto 412
c       if((kl.eq.1).or.(kl.eq.nz)) goto 412
c       nl=(kl-2)*nxy2+(j-2)*nx2+(i-1)
c       if((imerge(nl).eq.1).and.(ltype(jequal(nl)).ne.1)) goto 410
  412  continue
        vnew=vel(im,jm,kkm)+(vel(il,jl,kkl)-vel(im,jm,kkm))*
     &   (zn(k)-zn(km))/(zn(kl)-zn(km))
        nvl=nvl+1
        vsum=vsum+vnew
       endif
c
c updated linearly linked vel is average of directions
       vel(i,j,kk)=vsum/float(nvl)
       if(iuseq.eq.0) then
         if(vel(i,j,kk).lt.cvpmin) vel(i,j,kk)=cvpmin
         if(vel(i,j,kk).gt.cvpmax) vel(i,j,kk)=cvpmax
       endif
       vadj(n)=vel(i,j,kk)-vold
       if(iuseq.ne.0) goto 413
       if(abs(vadj(n)).le.dvpmx) goto 413
        dir=1.0
        if(vadj(n).lt.0.0) dir= -1.0
        vadj(n)=dvpmx*dir
        vel(i,j,kk)=vold+vadj(n)
  413  if(idiag.eq.1) write(96,1161) n,vold,vsum,nvl,vel(i,j,kk),
     &  vadj(n)
 1161  format('for linked node',i5,' vold=',f6.3,', vsum=',
     & f7.2,' nvl=',i1,' vel(i,j,kk)=',f6.3,' vadj=',f6.3)
       if(iuseq.ne.0) then
         qold=qval(i,j,k)
         qval(i,j,k)=vel(i,j,k+nz)/vel(i,j,k)
         if(qval(i,j,k).le.0.0) then
           qval(i,j,k) = 1.0
           vel(i,j,k+nz) = vel(i,j,k) * qval(i,j,k)
         end if
         vadj(n) = vel(i,j,k+nz)-vold
         qadj(n)=qval(i,j,k)-qold
         if(abs(qadj(n)).gt.dvpmx) then
           dir=1.0
           if(qadj(n).lt.0.0) dir= -1.0
           qadj(n)=dvpmx*dir
           qval(i,j,k)=qold+qadj(n)
           vel(i,j,k+nz)=vel(i,j,k)*qval(i,j,k)
         endif
       endif
      endif
c
c Vp/Vs
      if (k.gt.nz) then
      nvl=0
      vsum=0.0
      kk=k-nz
      kkl=kl-nz
      kkm=km-nz
      vold=vpvs(i,j,k-nz)
c x direction
       if (i.ne.im) then
c ***       il=2*i-im
       ia=1
       if(i-im.lt.0) ia= -1
       il=i
  422  il=il+ia
       if((il.eq.1).or.(il.eq.nx)) goto 424
       nl=(k-4)*nxy2+(j-2)*nx2+(il-1)
       if((imerge(nl).eq.1).and.(ltype(jequal(nl)).ne.1)) goto 422
  424  vnew=vpvs(im,jm,kkm)+(vpvs(il,jl,kkl)
     & -vpvs(im,jm,kkm))*(xn(i)-xn(im))/(xn(il)-xn(im))
        nvl=nvl+1
        vsum=vsum+vnew
       endif
c
c y direction
       if (j.ne.jm) then
c ***       jl=2*j-jm
       ja=1
       if(j-jm.lt.0) ja= -1
       jl=j
  426  jl=jl+ja
       if((jl.eq.1).or.(jl.eq.ny)) goto 428
       nl=(k-4)*nxy2+(jl-2)*nx2+(i-1)
       if((imerge(nl).eq.1).and.(ltype(jequal(nl)).ne.1)) goto 426
  428  vnew=vpvs(im,jm,kkm)+(vpvs(il,jl,kkl)
     & -vpvs(im,jm,kkm))*(yn(j)-yn(jm))/(yn(jl)-yn(jm))
        nvl=nvl+1
        vsum=vsum+vnew
c       write(96,*) 'i,j,k   ',i,j,k,' vpvs(i,j,k-nz)',vpvs(i,j,k-nz)
c       write(96,*) 'il,jl,kl',il,jl,kl,'vpvs(il,jl,kkl)',
c     & vpvs(il,jl,kkl)
c      write(96,*) 'im,jm,km',im,jm,km,'vpvs(im,jm,kkm)',
c     & vpvs(im,jm,kkm)
c      write(96,*) 'yn(j),yn(jl),yn(jm)',yn(j),yn(jl),yn(jm),' vnew',
c     & vnew
       endif
c
c z direction
       if (k.ne.km) then
c ***       kl=2*k-km
       ka=1
       if(k-km.lt.0) ka= -1
       kl=k
  430  kl=kl+ka
      kkl=kkl+ka
       if((kl.eq.nz+1).or.(kl.eq.2*nz)) goto 432
       nl=(kl-4)*nxy2+(j-2)*nx2+(i-1)
       if((imerge(nl).eq.1).and.(ltype(jequal(nl)).ne.1)) goto 430
  432  vnew=vpvs(im,jm,kkm)+(vpvs(il,jl,kkl)
     & -vpvs(im,jm,kkm))*(zn(kk)-zn(kkm))/(zn(kkl)-zn(kkm))
        nvl=nvl+1
        vsum=vsum+vnew
c       write(96,*) 'i,j,k   ',i,j,k,' vpvs(i,j,k-nz)',vpvs(i,j,k-nz)
c       write(96,*) 'il,jl,kl',il,jl,kl,'vpvs(il,jl,kkl)',
c     & vpvs(il,jl,kkl)
c      write(96,*) 'im,jm,km',im,jm,km,'vpvs(im,jm,kkm)',
c     & vpvs(im,jm,kkm)
c      write(96,*) 'zn(kk),zn(kkl),zn(kkm)',zn(kk),zn(kkl),zn(kkm),
c     & ' vnew',vnew
       endif
c
c updated linearly linked vel is average of directions
       vpvs(i,j,k-nz)=vsum/float(nvl)
       vadj(n)= vpvs(i,j,k-nz)-vold
      if(abs(vadj(n)).le.dvsmx) goto 433
      dir=1.0
      if(vadj(n).lt.0.0) dir= -1.0
      vadj(n)=dvsmx*dir
      vpvs(i,j,k-nz)=vold+vadj(n)
  433  if(idiag.eq.1) write(96,1162) n,vold,vsum,nvl,vpvs(i,j,k-nz),
     &  vadj(n)
 1162  format('for linked node',i5,' vold=',f6.3,', vsum=',
     & f7.3,' nvl=',i1,' vpvs(i,j,k-nz)=',f6.3,' vadj=',f6.3)
      endif
c
c **** end of linearly-linked node section
   33  continue
  440 continue
c
c  end of linked node section
c
      if(iuses.eq.1) goto 130
c  Compute Vs with new Vp and new Vp/Vs
        do 120 k=1,nz
           ks=k+nz
           do 115 j=1,ny
              do 110 i=1,nx
                 vel(i,j,ks)=vel(i,j,k)/vpvs(i,j,k)
  110         continue
  115      continue
  120   continue
c**
c
c Recalculate damping
c consider c1 to give relative size of 2 terms being minimized
c  min [ sum(data resid squared) + damping*(solution norm squared)]
  130 do 34 i=1,iuses
         if(i.eq.1) ssq=ssqrwp
         if(i.eq.2) ssq=ssqrws
         if(ssq.le.0) goto 34
         if(mbla(i).eq.0) goto 34
         if(nit.eq.1) c1(i)=vdamp(i)*sqnorm(i+1)/ssq
         if(nit.gt.1) dampc1(i)=c1(i)*ssq/sqnorm(i+1)
         write(16,1626) solnam(i+1),vdamp(i),c1(i),dampc1(i)
         write(36,1626) solnam(i+1),vdamp(i),c1(i),dampc1(i)
 1626    format('Damping: ',a11,' damp =',e12.3,', c1 =',e13.4,
     2      ', damp(c1) =',e12.3)
c Replace damping if idmp is turned on
c Do not change damping to be smaller
         if(dampc1(i).lt.vdamp(i)) goto 34
         if((idmp.eq.1).and.(nit.gt.1)) vdamp(i)=dampc1(i)
   34 continue
   36 continue
      isadj=0
c
      if(invdel.eq.0) goto 50
c  Station correction adjustments
      do 40 n=1,nsts
         nvn=n+nparv
         if(khit(nvn).eq.0) goto 40
         if((hit(nvn).lt.hitct).or.(nfixst(n).eq.1)) goto 40
         nv=ndexfx(nvn)
         vadj(nvn)=zero
         nin=index(nv)
         vadj(nvn)=rhs(nin)
         pdl(n)=pdl(n)+vadj(nvn)
         isadj=isadj+1
   40 continue
      if(iuses.eq.1) goto 50
c  S delay adjustments
      do 45 n=1,nsts
         nvn=n+nsts+nparv
         if(khit(nvn).eq.0) goto 45
         if((hit(nvn).lt.hitct).or.(nfixst(n).eq.1)) goto 45
         nv=ndexfx(nvn)
         vadj(nvn)=zero
         nin=index(nv)
         vadj(nvn)=rhs(nin)
         sdl(n)=sdl(n)+vadj(nvn)
         isadj=isadj+1
   45 continue
c
   50 itadj=0
      if(ntel.eq.0) goto 99
c  Teleseismic event path delay term adjustments
      do 52 n=1,ntel
        nvn=n+nparvs
        if(khit(nvn).eq.0) goto 52
        if(hit(nvn).lt.hitct) goto 52
        nv=ndexfx(nvn)
        vadj(nvn)=zero
        nin=index(nv)
        vadj(nvn)=rhs(nin)
        telpad(n)=telpad(n)+vadj(nvn)
        itadj=itadj+1
      WRITE(56,*) 'n,hit,telpad',n,hit(nvn),telpad(n)
   52 continue
c S-P tele paths delay adjustment
      do 55 n=1,ntel
        nvn=n+ntel+nparvs
        if(khit(nvn).eq.0) goto 55
        if(hit(nvn).lt.hitct) goto 55
        nv=ndexfx(nvn)
        vadj(nvn)=zero
        nin=index(nv)
        vadj(nvn)=rhs(nin)
        telpads(n)=telpads(n)+vadj(nvn)
        itadj=itadj+1
      WRITE(56,*) 'n,hit,telpads',n,hit(nvn),telpads(n)
   55 continue
   99 continue
        write(16,1602) ivadj,isadj,itadj
 1602 format(i10,' velocity adjustments,',i10,' station corr. ',
     2 'adj,',i10,' Tele path delay adj.')
c
      close(96)
      return
c***** end of subroutine veladj *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine velbku(istop1)
c  This routine will backup the velocity and station parameters halfway
c  common block variables:
      include 'simul2014_common.inc'
c
c  apply adjustments to velocity model
        ivadj=0
      if(nparvi.eq.0) goto 36
c
c  with Linked nodes, need to do all nodes
c      do 30 n=1,nparvi
      do 30 n=1,nparv
c  gridpoint number in array of all nodes
c         nn=mdexfx(n)
         nn=n
         if(vadj(nn).eq.0.0) goto 30
c  calculate x and z indices of velocity grid
         k=(nn-1)/nxy2+2
         j=2+(nn-1+(2-k)*nxy2)/nx2
         i=1+nn+nx2*(2-j)+nxy2*(2-k)
         if(k.ge.nz) k=k+2   ! if S velocity node
c  compute half adjustment if first backup
         if(istop1.eq.0) vadj(nn)= -0.5*vadj(nn)
c  apply adjustment to model
c** S-P CHANGE
         if (k.le.nz) then
          if(iuseq.eq.0) then
            vel(i,j,k)=vel(i,j,k)+vadj(nn)
          else
            vel(i,j,k+nz)=vel(i,j,k+nz)+vadj(nn)
            qval(i,j,k)=vel(i,j,k+nz)/vel(i,j,k)
            if(qval(i,j,k).le.0.0) then
                qval(i,j,k) = 1.0
                vel(i,j,k+nz) = vel(i,j,k) * qval(i,j,k)
                vadj(nn) = vadj(nn)-vel(i,j,k+nz)
            end if
          endif
       endif
         if(k.ge.nz) then
            vpvs(i,j,k-nz)=vpvs(i,j,k-nz)+vadj(nn)
         endif
         ivadj=ivadj+1
   30 continue
c
      if(iuses.eq.1) goto 36
c  Compute Vs with new Vp and new Vp/Vs
        do 120 k=1,nz
           ks=k+nz
           do 115 j=1,ny
              do 110 i=1,nx
                 vel(i,j,ks)=vel(i,j,k)/vpvs(i,j,k)
  110         continue
  115      continue
  120   continue
c**
c
   36 continue
        isadj=0
      if(invdel.eq.0)goto 50
c  Station correction adjustments
      do 40 n=1,nsts
         nvn=n+nparv
         if(khit(nvn).eq.0) goto 40
         if((hit(nvn).lt.hitct).or.(nfixst(n).eq.1)) goto 40
         if(istop1.eq.0) vadj(nvn)= -0.5*vadj(nvn)
         pdl(n)=pdl(n)+vadj(nvn)
         isadj=isadj+1
   40 continue
      if(iuses.eq.0) goto 50
c  S delay adjustments
      do 45 n=1,nsts
         nvn=n+nsts+nparv
         if(khit(nvn).eq.0) goto 45
         if((hit(nvn).lt.hitct).or.(nfixst(n).eq.1)) goto 45
         if(istop1.eq.0) vadj(nvn)= -0.5*vadj(nvn)
         sdl(n)=sdl(n)+vadj(nvn)
         isadj=isadj+1
   45 continue
c
   50 itadj=0
      if(ntel.eq.0) goto 99
c  Teleseismic event path delay term adjustments
      do 52 n=1,ntel
        nvn=n+nparvs
        if(khit(nvn).eq.0) goto 52
        if(hit(nvn).lt.hitct) goto 52
        if(istop1.eq.0) vadj(nvn)= -0.5*vadj(nvn)
        telpad(n)=telpad(n)+vadj(nvn)
        itadj=itadj+1
   52 continue
c S-P tele paths delay adjustment
      do 55 n=1,nsts
        nvn=n+ntel+nparvs
        if(khit(nvn).eq.0) goto 55
        if(hit(nvn).lt.hitct) goto 55
        if(istop1.eq.0) vadj(nvn)= -0.5*vadj(nvn)
        telpads(n)=telpads(n)+vadj(nvn)
        itadj=itadj+1
   55 continue
   99 continue
        write(16,1602) ivadj,isadj,itadj
 1602 format('Backup',i10,' velocity adjustments,',i10,'station corr. ',
     2 'adj.',i10,' Tele path delay adj.')
c
      return
c***** end of subroutine velbku *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine veld(isp,xx,yy,zz,vx,vy,vz)
c
c*****this routine computes the derivatives of velocity
c     in x, y, and z directions
c
c  common block variables:
      include 'simul2014_common.inc'
c
c  use Prothero's intmap here
      call intmap(xx,yy,zz,ip,jp,kp)
c
      ip1=ip+1
      jp1=jp+1
      kp1=kp+1
c
      xd=xn(ip1)-xn(ip)
      yd=yn(jp1)-yn(jp)
      zd=zn(kp1)-zn(kp)
c
      xf=(xx-xn(ip))/xd
      yf=(yy-yn(jp))/yd
      zf=(zz-zn(kp))/zd
c
      xf1=1.0-xf
      yf1=1.0-yf
      zf1=1.0-zf
c
c  S-velocity is stored in the 2nd half of the velocity array
      if(isp.eq.1) kp=kp+nz
      kp1=kp+1
c
c*****calculate derivatives of velocity
c
      vx=(yf1*zf1*(vel(ip1,jp,kp)-vel(ip,jp,kp))
     *+yf*zf1*(vel(ip1,jp1,kp)-vel(ip,jp1,kp))
     *+yf1*zf*(vel(ip1,jp,kp1)-vel(ip,jp,kp1))
     *+yf*zf*(vel(ip1,jp1,kp1)-vel(ip,jp1,kp1)))/xd
c
      vy=(xf1*zf1*(vel(ip,jp1,kp)-vel(ip,jp,kp))
     *+xf*zf1*(vel(ip1,jp1,kp)-vel(ip1,jp,kp))
     *+xf1*zf*(vel(ip,jp1,kp1)-vel(ip,jp,kp1))
     *+xf*zf*(vel(ip1,jp1,kp1)-vel(ip1,jp,kp1)))/yd
c
      vz=(xf1*yf1*(vel(ip,jp,kp1)-vel(ip,jp,kp))
     *+xf*yf1*(vel(ip1,jp,kp1)-vel(ip1,jp,kp))
     *+xf1*yf*(vel(ip,jp1,kp1)-vel(ip,jp1,kp))
     *+xf*yf*(vel(ip1,jp1,kp1)-vel(ip1,jp1,kp)))/zd
c
      return
c ***** end of subroutine veld *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine wthyp(ne,nobs,ah,rmswt,nwr,w)
c  routine to apply weighting to hypocenter matrix + residuals
c
c  common block variables:
      include 'simul2014_common.inc'
c
c  declaration statements:
      real w(maxobs)
      double precision ah(maxobs,maxobs)
      character*1 phs(2)
      parameter (pi=3.1415926536)
c
      data phs/'P','S'/
c
c  compute weighting
c
      nwr=0
      wnorm=0.0
      do 13 j=1,nobs
c  reading weight
         w(j)=wt(j,ne)
c  change relative weighting of S-P observations
         if(intsp(j,ne).eq.1) w(j)=w(j)*wtsp
c  residual weighting
c  downweighting(linear) 0 to 98% res1 to res2, 98 to 100% res2 to res3
         ares=abs(res(j))
         if(ares.le.res2) then
            wr=1.0-(ares-res1)*dres12
            if (wr.gt.1.0) wr=1.0
         else
            if(res3.gt.res2) then
               wr=0.02-(ares-res2)*dres23
               if (wr.lt.0.0) wr=0.0
            else
               wr=0.0
            endif
         endif
c  distance weighting
         wd=1.0-(dlta(j,ne)-delt1)*ddlt
         if (wd.gt.1.0) wd=1.0
         if (wd.lt.0.0) wd=0.0
c  unnormalized weight
         w(j)=w(j)*wr*wd
         wnorm=wnorm+w(j)
         if (w(j).gt.0.0) nwr=nwr+1
  13  continue
c  check to be sure 4 or more readings with nonzero weight
      if(ne.le.netemp) then
        if(nwr.lt.nparhy) goto 900
      else
        if(nwr.lt.2) goto 900
      endif
c  normalize weights
   12 wfac=nwr/wnorm
c
   20 continue
      wnorm=0.0
      rmswt=0.0
c  normalize weights, apply to hypocenter matrix, and place into
c  a-matrix for inversion
      do 30 i=1,nobs
         w(i)=w(i)*wfac
         wnorm=wnorm+w(i)*w(i)
         rmswt=rmswt+(res(i)*w(i))*(res(i)*w(i))
         do 35 j=1,nparhy
            ah(i,j)=dth(i,j)*w(i)
   35    continue
         ah(i,nparhy+1)=res(i)*w(i)
   30 continue
      rmswt=sqrt(rmswt/wnorm)
cDEP Save weight
c      print *, "WTHYP ne", ne, "w",w
      do 60 j=1,nobs
      wtcomb(j,ne)=w(j)
      if(ne.gt.netemp) wtcomb(j,ne)=wtcomb(j,ne)*wtsht
   60 continue
      return
  900 write(16,1601) ne
 1601 format(' in WTHYP, event',i5,' HAS LESS THAN 4 READINGS WITH ',
     2 '  NON-ZERO WEIGHT, SHOULD BE DELETED')
      if(kttfor.ne.3) then
        write(16,51)
51      format(1x,4('sta ph   wt   res   ttime delta',2x))
        write(16,53) (stn(isto(j,ne)),phs(intsp(j,ne)+1),w(j),res(j),
     2  (secp(j,ne)-seco(ne)),dlta(j,ne),j=1,kobs(ne))
53      format(4(1x,a4,1x,a1,f6.2,f7.3,2f6.2,1x))
      else
        write(16,41)
41      format(1x,4('sta   ph   wt   res   ttime delta',2x))
        write(16,43)(stn6(isto(j,ne)),phs(intsp(j,ne)+1),w(j),res(j),
     2  (secp(j,ne)-seco(ne)),dlta(j,ne),j=1,kobs(ne))
43      format(4(1x,a6,1x,a1,f6.2,f7.3,2f6.2,1x))
      endif
      write(16,1602)
 1602 format(/)
      return
c***** end of subroutine wthyp *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine wthypc(nc,ne,nobs1,nobs,ahc,rmswt,nwr,ssqrwe,wnorme,wc)
c  routine to apply weighting to hypocenter matrix + residuals
c  This version is for an individual event within a cluster
c
c  common block variables:
      include 'simul2014_common.inc'
c
c  declaration statements:
      real wc(mxobsc),w(maxobs)
      double precision ahc(mxobsc,mxobsc)
      character*1 phs(2)
      parameter (pi=3.1415926536)
c
      data phs/'P','S'/
c
c  compute weighting
c
      nwr=0
      wnorm=0.0
      do 13 j=1,nobs
c  reading weight
         w(j)=wtc(j,ne,nc)
c  change relative weighting of S-P observations
         if(intspc(j,ne,nc).eq.1) w(j)=w(j)*wtsp
c  residual weighting
c  downweighting(linear) 0 to 98% res1 to res2, 98 to 100% res2 to res3
         ares=abs(resc(j,ne))
         if(ares.le.res2) then
            wr=1.0-(ares-res1)*dres12
            if (wr.gt.1.0) wr=1.0
         else
            if(res3.gt.res2) then
               wr=0.02-(ares-res2)*dres23
               if (wr.lt.0.0) wr=0.0
            else
               wr=0.0
            endif
         endif
c  distance weighting
         wd=1.0-(dltac(j,ne,nc)-delt1)*ddlt
         if (wd.gt.1.0) wd=1.0
         if (wd.lt.0.0) wd=0.0
c  unnormalized weight
         w(j)=w(j)*wr*wd
         wnorm=wnorm+w(j)
         if (w(j).gt.0.0) nwr=nwr+1
  13  continue
c  check to be sure 4 or more readings with nonzero weight
        if(nwr.lt.nparhy) goto 900
c  normalize weights
   12 wfac=nwr/wnorm
c
   20 continue
      wnorm=0.0
      rmswt=0.0
c  normalize weights, apply to hypocenter matrix, and place into
c  a-matrix for inversion
c  k is counter in all obs for this cluster
      ja1=nparhy*(ne-1)
      jahy=nparhy*ncev(nc)
      do 30 i=1,nobs
        k=nobs1+i
        w(i)=w(i)*wfac
        wc(k)=w(i)
        wnorm=wnorm+w(i)*w(i)
        wres=resc(i,ne)*w(i)
        rmswt=rmswt+(wres*wres)
        do 35 j=1,nparhy
          ahc(k,ja1+j)=dthc(i,j,ne)*w(i)
   35   continue
        ahc(k,jahy+1)=wres
   30 continue
      ssqrwe=rmswt
      wnorme=wnorm
      rmswt=sqrt(rmswt/wnorm)
cDEP Save weight
c      print *, "WTHYP ne", ne, "w",w
      do 60 j=1,nobs
      wtcombc(j,ne,nc)=w(j)
   60 continue
      return
  900 write(16,1601) ne
 1601 format(' ERROR*** -in WTHYPC, event',i5,' HAS LESS THAN 4 ',
     2 '  READINGS WITH NON-ZERO WEIGHT, SHOULD BE DELETED')
      if(kttfor.ne.3) then
        write(16,51)
51      format(1x,4('sta ph   wt   res   ttime delta',2x))
        write(16,53) (stn(isto(j,ne)),phs(intsp(j,ne)+1),w(j),res(j),
     2  (secp(j,ne)-seco(ne)),dlta(j,ne),j=1,kobs(ne))
53      format(4(1x,a4,1x,a1,f6.2,f7.3,2f6.2,1x))
      else
        write(16,41)
41      format(1x,4('sta   ph   wt   res   ttime delta',2x))
        write(16,43)(stn6(isto(j,ne)),phs(intsp(j,ne)+1),w(j),res(j),
     2  (secp(j,ne)-seco(ne)),dlta(j,ne),j=1,kobs(ne))
43      format(4(1x,a6,1x,a1,f6.2,f7.3,2f6.2,1x))
      endif
      write(16,1602)
 1602 format(/)
      return
c***** end of subroutine wthypc *****
      end
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine utm_geo(rlon,rlat,rx,ry,cmerid,UTM_PROJECTION_ZONE,
     2 iway,inorth)
c
c convert geodetic longitude and latitude to UTM, and back
c from geodynamics.org/svn/cig/seismo/3D/SPECFEM3D/trunk/src/shared/utm_geo.f90
c  input values must be double precision
c 25-sept-2014 Donna EP obtained and made small changes.  It is called by 
c  subroutine tmllxy
c changed iway, removed SUPPRESS_UTM_PROJECTION
c added cmerid for user input central meridian, used if UTM_PROJECTION_ZONE=0
c inorth =0 for northern hemisphere, =1 for southern hemisphere
c=====================================================================
c               S p e c f e m 3 D  V e r s i o n  2 . 1
c               ---------------------------------------
c          Main authors: Dimitri Komatitsch and Jeroen Tromp
c    Princeton University, USA and CNRS / INRIA / University of Pau
c (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
c                             July 2012
c This program is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 2 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c=====================================================================
c
c  UTM (Universal Transverse Mercator) projection from the USGS
c=====================================================================
c  subroutine utm_geo(rlon,rlat,rx,ry,UTM_PROJECTION_ZONE,iway,SUPPRESS_UTM_PROJECTION)
c
c DMEP change
c use iway=1 for long/lat to tranverse mercator
c     iway=2 for transverse mercator to lat/lon
c use iway = ILONGLAT2UTM for long/lat to UTM, IUTM2LONGLAT for UTM to lat/long
c a list of UTM zones of the world is available at www.dmap.co.uk/utmworld.htm

      implicit none

c only needed PI (DMEP)
c  include "constants.h"
c
c-----CAMx v2.03
c
c     UTM_GEO performs UTM to geodetic (long/lat) translation, and back.
c
c     This is a Fortran version of the BASIC program "Transverse Mercator
c     Conversion", Copyright 1986, Norman J. Berls (Stefan Musarra, 2/94)
c     Based on algorithm taken from "Map Projections Used by the USGS"
c     by John P. Snyder, Geological Survey Bulletin 1532, USDI.
c
c     Input/Output arguments:
c
c        rlon                  Longitude (deg, negative for West)
c        rlat                  Latitude (deg)
c        rx                    UTM easting (m)
c        ry                    UTM northing (m)
c        UTM_PROJECTION_ZONE  UTM zone
c        iway                  Conversion type
c                              1 (ILONGLAT2UTM) = geodetic to UTM
c                              2 (IUTM2LONGLAT) = UTM to geodetic

      integer UTM_PROJECTION_ZONE,iway,inorth
      double precision rx,ry,rlon,rlat,cmerid
      logical SUPPRESS_UTM_PROJECTION

      double precision PI,degrad,raddeg,semimaj,semimin,scfa
      parameter (PI = 3.141592653589793d0)
      parameter (degrad=PI/180.d0, raddeg=180.d0/PI)
      parameter (semimaj=6378206.4d0, semimin=6356583.8d0)
      parameter (scfa=0.9996d0)

c some extracts about UTM:
c
c There are 60 longitudinal projection zones numbered 1 to 60 starting at 180Â°W.
c Each of these zones is 6 degrees wide, apart from a few exceptions around Norway and Svalbard.
c There are 20 latitudinal zones spanning the latitudes 80Â°S to 84Â°N and denoted
c by the letters C to X, ommitting the letter O.
c Each of these is 8 degrees south-north, apart from zone X which is 12 degrees south-north.
c
c To change the UTM zone and the hemisphere in which the
c calculations are carried out, need to change the fortran code and recompile. The UTM zone is described
c actually by the central meridian of that zone, i.e. the longitude at the midpoint of the zone, 3 degrees
c from either zone boundary.
c To change hemisphere need to change the "north" variable:
c  - north=0 for northern hemisphere and
c  - north=10000000 (10000km) for southern hemisphere. values must be in metres i.e. north=10000000.
c
c Note that the UTM grids are actually Mercators which
c employ the standard UTM scale factor 0.9996 and set the
c Easting Origin to 500,000;
c the Northing origin in the southern
c hemisphere is kept at 0 rather than set to 10,000,000
c
c and this gives a uniform scale across the equator if the
c normal convention of selecting the Base Latitude (origin)
c at the equator (0 deg.) is followed.  Northings are
c positive in the northern hemisphere and negative in the
c southern hemisphere.
      double precision north,east
      parameter (east=500000.d0)

      double precision e2,e4,e6,ep2,xx,yy,dlat,dlon,zone,cm,cmr,delam
      double precision f1,f2,f3,f4,rm,rn,t,c,a,e1,u,rlat1,dlat1,c1,t1,
     2 rn1,r1,d
      double precision rx_save,ry_save,rlon_save,rlat_save

      if(inorth.eq.0) north=0.d0
      if(north.eq.1) north=10000000.d0

c      PRINT *,' utm_geo cmerid ',cmerid
c checks if conversion to utm has to be done
c   if(SUPPRESS_UTM_PROJECTION) then
c   if (iway == ILONGLAT2UTM) then
c          rx = rlon
c          ry = rlat
c        else
c          rlon = rx
c          rlat = ry
c        endif
c        return
c   endif

c save original parameters
      rlon_save = rlon
      rlat_save = rlat
      rx_save = rx
      ry_save = ry

      xx = 0.d0
      yy = 0.d0
      dlat = 0.d0
      dlon = 0.d0

c define parameters of reference ellipsoid
      e2=1.0-(semimin/semimaj)**2.0
      e4=e2*e2
      e6=e2*e4
      ep2=e2/(1.-e2)

c if (iway == IUTM2LONGLAT) then
      if (iway == 2) then
        xx = rx
        yy = ry
      else
        dlon = rlon
        dlat = rlat
      endif
c
c      PRINT *,'UTM_PROJECTION_ZONE ',UTM_PROJECTION_ZONE,' cmerid ',
c     2 cmerid
      if(UTM_PROJECTION_ZONE.eq.0) then
        cm=cmerid
      else
c----- Set Zone parameters
        zone = dble(UTM_PROJECTION_ZONE)
c sets central meridian for this zone
        cm = zone*6.0 - 183.0
      endif
      cmr = cm*degrad
c      PRINT *,'cm',cm,'cmr',cmr
c
c---- Lat/Lon to UTM conversion
c
c  if (iway == ILONGLAT2UTM) then
      if (iway == 1) then

      rlon = degrad*dlon
      rlat = degrad*dlat

      delam = dlon - cm
      if (delam < -180.) delam = delam + 360.
      if (delam > 180.) delam = delam - 360.
      delam = delam*degrad

      f1 = (1. - e2/4. - 3.*e4/64. - 5.*e6/256)*rlat
      f2 = 3.*e2/8. + 3.*e4/32. + 45.*e6/1024.
      f2 = f2*sin(2.*rlat)
      f3 = 15.*e4/256.*45.*e6/1024.
      f3 = f3*sin(4.*rlat)
      f4 = 35.*e6/3072.
      f4 = f4*sin(6.*rlat)
      rm = semimaj*(f1 - f2 + f3 - f4)
      if (dlat == 90. .or. dlat == -90.) then
        xx = 0.
        yy = scfa*rm
      else
        rn = semimaj/sqrt(1. - e2*sin(rlat)**2)
        t = tan(rlat)**2
        c = ep2*cos(rlat)**2
        a = cos(rlat)*delam

        f1 = (1. - t + c)*a**3/6.
        f2 = 5. - 18.*t + t**2 + 72.*c - 58.*ep2
        f2 = f2*a**5/120.
        xx = scfa*rn*(a + f1 + f2)
        f1 = a**2/2.
        f2 = 5. - t + 9.*c + 4.*c**2
        f2 = f2*a**4/24.
        f3 = 61. - 58.*t + t**2 + 600.*c - 330.*ep2
        f3 = f3*a**6/720.
        yy = scfa*(rm + rn*tan(rlat)*(f1 + f2 + f3))
      endif
      xx = xx + east
      yy = yy + north

c
c---- UTM to Lat/Lon conversion
c
      else

      xx = xx - east
      yy = yy - north
      e1 = sqrt(1. - e2)
      e1 = (1. - e1)/(1. + e1)
      rm = yy/scfa
      u = 1. - e2/4. - 3.*e4/64. - 5.*e6/256.
      u = rm/(semimaj*u)

      f1 = 3.*e1/2. - 27.*e1**3./32.
      f1 = f1*sin(2.*u)
      f2 = 21.*e1**2/16. - 55.*e1**4/32.
      f2 = f2*sin(4.*u)
      f3 = 151.*e1**3./96.
      f3 = f3*sin(6.*u)
      rlat1 = u + f1 + f2 + f3
      dlat1 = rlat1*raddeg
      if (dlat1 >= 90. .or. dlat1 <= -90.) then
        dlat1 = dmin1(dlat1,dble(90.) )
        dlat1 = dmax1(dlat1,dble(-90.) )
        dlon = cm
      else
        c1 = ep2*cos(rlat1)**2.
        t1 = tan(rlat1)**2.
        f1 = 1. - e2*sin(rlat1)**2.
        rn1 = semimaj/sqrt(f1)
        r1 = semimaj*(1. - e2)/sqrt(f1**3)
        d = xx/(rn1*scfa)

        f1 = rn1*tan(rlat1)/r1
        f2 = d**2/2.
        f3 = 5.*3.*t1 + 10.*c1 - 4.*c1**2 - 9.*ep2
        f3 = f3*d**2*d**2/24.
        f4 = 61. + 90.*t1 + 298.*c1 + 45.*t1**2. - 252.*ep2 - 3.*c1**2
        f4 = f4*(d**2)**3./720.
        rlat = rlat1 - f1*(f2 - f3 + f4)
        dlat = rlat*raddeg

        f1 = 1. + 2.*t1 + c1
        f1 = f1*d**2*d/6.
        f2 = 5. - 2.*c1 + 28.*t1 - 3.*c1**2 + 8.*ep2 + 24.*t1**2.
        f2 = f2*(d**2)**2*d/120.
        rlon = cmr + (d - f1 + f2)/cos(rlat1)
        dlon = rlon*raddeg
        if (dlon < -180.) dlon = dlon + 360.
        if (dlon > 180.) dlon = dlon - 360.
      endif
      endif

c  if (iway == IUTM2LONGLAT) then
      if (iway == 2) then
        rlon = dlon
        rlat = dlat
        rx = rx_save
        ry = ry_save
      else
        rx = xx
        ry = yy
        rlon = rlon_save
        rlat = rlat_save
      endif

      end subroutine utm_geo
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine llnzmg(latd,xlat,lond,xlon,pnorkm,peaskm)
c
c     programs for conversion to NZ Map Grid coordinates
c     obtained from Chris Pearson, Otago U Surveying Dept.
c     Now save coordinates in KM (instead of M) as pnorkm, peaskm.
c
c      PROGRAM LLNZMG
C     ==============
C     THIS PROGRAM TRANSFORMS GEOGRAPHICAL COORDINATES TO NZMG
C     COORDINATES
       common/nzmg1/ick
      REAL*8 A,N,N0,E,E0,DLAT,DPSI,DLON,LAT,LON,CS,PI
      real xlat,xlon,pnorkm, peaskm
      REAL*8 FI0,LO0,CONV,DRZ,DIZ
      INTEGER*4 LATD,LOND,CD,CM
      COMPLEX*16 Z,ZETA
c
c      print *, 'LLNZMG latd,xlat,lond,xlon',latd,xlat,lond,xlon
      LO0=173.D0
      FI0=41.D0
      A=6378388.D0
      N0=6023150.D0
      E0=2510000.D0
      PI=1.0D0
      PI=DATAN(PI)*4.
C     READ THE GEOGRAPHICAL COORDINATES
C     ---------------------------------
c      WRITE(6,104)
c104   FORMAT('1'//' Enter the LONGITUDE(E) of station in the f',
c     1'orm DDD MM SS.SSSS.'/' Key RETURN.'//)
c      READ(5,*)LOND,LONM,LONS
c      WRITE(6,105)
c105   FORMAT('1'//' Enter the LATITUDE(S) of station in the fo',
c     1'rm DD MM SS.SSSS.'/' Key RETURN.'//)
c      READ(5,*)LATD,LATM,LATS
c      LAT=DBLE(LATD)+DBLE(LATM)/60.+LATS/3600.
c      LON=DBLE(LOND)+DBLE(LONM)/60.+LONS/3600.
c
c  Use absolute values since simul has negative for S
c  and LLNZMG assume all E and S
      LAT=DBLE(abs(LATD))+DBLE(abs(xlat))/60.
      LON=DBLE(abs(LOND))+DBLE(abs(xlon))/60.
c
C     COMPUTE THE DIFFERENCES IN LATITUDE AND LONGITUDE IN DEGREES
C     ------------------------------------------------------------
      DLAT=FI0-LAT
      DLON=LON-LO0
C     COMPUTE THE ISOMETRIC LAT DIFF FROM THE GEODETIC LAT DIFF
C     ---------------------------------------------------------
      CALL CLID(DLAT,DPSI)
C     --------------------------------------------
      DLON=DLON*PI/180.
      ZETA=DCMPLX(DPSI,DLON)
C     COMPUTE THE COMPLEX NUMBER "Z"
C     ------------------------------
      CALL CGCD(ZETA,Z)
C     COMPUTE E AND N ON THE NZMG
C     ---------------------------
      DRZ=DREAL(Z)
      DIZ=DIMAG(Z)
      N=DRZ+N0
      E=DIZ+E0
c
c save in km for simul
      pnorkm=N/1000.
      peaskm=E/1000.
c
C     COMPUTE CONVERGENCE
C     -------------------
      CALL DZDZEE (ZETA,CONV)
      CALL CDDMS (CONV,CD,CM,CS)
C     WRITE HEADINGS
C     --------------
c       if(ick.eq.2) goto 95
c98    WRITE(6,106)
c106   FORMAT(2X,'STATION',4X,'LONGITUDE(E)',4X,'LATITUDE(S) ',6X,
c     1'EASTING',4X,'NORTHING',3X,'CONVERGENCE'/2X,7('='),4X,3('===='),
c     24X,11('='),7X,7('='),4X,8('='),3X,11('='))
cC     WRITE STATION,LONGITUDE,LATITUDE,E,N & CONVERGENCE
c   95 WRITE(6,107)abs(LOND),abs(xlon),abs(LATD),abs(xlat),
c     2 E,N,CD,CM,CS
c107   FORMAT(14x,i3,F8.4,4x,I3,F8.4,2X,2F12.3,2X,2I3,F6.2)
c      ick=2
120   continue
      return
c***** end of subroutine llnzmg *****
      END
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine nzmgll(pnorkm,peaskm,latd,xlat,lond,xlon)
c
c     programs for conversion from NZ Map Grid coordinates 
c     to latitude and longitude
c     obtained from Chris Pearson, Otago U Surveying Dept.
c     Now save coordinates in KM (instead of M) as pnorkm, peaskm.
c
c
c      PROGRAM NZMGLL
C     ==============
C     THIS PROGRAM TRANSFORMS NZMG COORDINATES TO GEOGRAPHICAL
C     COORDINATES
      REAL*8 A,N,N0,E,E0,DLAT,DPSI,DLON,LAT,LON,LATS,LONS,CS,PI
      REAL*8 FI0,LO0,CONV
      real xlat,xlon,pnorkm, peaskm
      INTEGER*4 LATD,LATM,LOND,LONM,CD,CM
      COMPLEX*16 Z,ZETA
      LO0=173.D0
      FI0=41.D0
      A=6378388.D0
      N0=6023150.D0
      E0=2510000.D0
      PI=1.0D0
      PI=DATAN(PI)*4.
C     READ THE NZMG COORDINATES
C     -------------------------
c      WRITE(6,104)
c104   FORMAT('1'//' Enter the EASTING of station and key RETUR',
c     1'N'//)
c      READ(5,*)E
c      WRITE(6,105)
c105   FORMAT('1'//' Enter the NORTHING of station and key RETU',
c     1'RN'//)
c      READ(5,*)N
c
c convert km to metre
      E=peaskm*1000.0
      N=pnorkm*1000.0
c
C     COMPUTE THE COMPLEX NUMBER Z=((N-N0)+i(E-E0))/A
C     -----------------------------------------------
      Z=DCMPLX((N-N0)/A,(E-E0)/A)
C     COMPUTE AN APPROXIMATE VALUE FOR THE COMPLEX NUMBER ZETA
C     --------------------------------------------------------
      CALL AZETA (Z,ZETA)
C     COMPUTE A FINAL VALUE FOR ZETA BY ITERATION
C     -------------------------------------------
      CALL ITZETA (Z,ZETA)
      DPSI=DREAL(ZETA)
      DLON=DIMAG(ZETA)
C     COMPUTE GEODETIC LAT DIFF FROM ISOMETRIC LAT DIFF
C     -------------------------------------------------
      CALL CILD (DPSI,DLAT)
C     COMPUTE LATITUDE FROM LATITUDE DIFFERENCE
C     -----------------------------------------
      LAT=FI0-DLAT+.00005/3600.
      CALL CDDMS (LAT,LATD,LATM,LATS)
C     COMPUTE LONGITUDE FROM LONGITUDE DIFFERENCE
C     -------------------------------------------
      DLON=DLON*180./PI
      LON=LO0+DLON+.00005/3600.
      CALL CDDMS (LON,LOND,LONM,LONS)
C     COMPUTE CONVERGENCE
C     -------------------
      CALL DZDZEE (ZETA,CONV)
      CALL CDDMS (CONV,CD,CM,CS)
c
c      simul2000 use deg and minutes
      xlat=sngl(dfloat(latm)+lats/60.0)
      xlon=sngl(dfloat(lonm)+lons/60.0)
c
C     WRITE HEADINGS
C     --------------
c98    WRITE(6,106)
c106   FORMAT(2X,'STATION',4X,'LONGITUDE(E)',4X,'LATITUDE(S) ',6X,
c     1'EASTING',4X,'NORTHING',3X,'CONVERGENCE'/2X,7('='),4X,3('===='),
c     24X,11('='),7X,7('='),4X,8('='),3X,11('='))
C     WRITE STATION,LONGITUDE,LATITUDE,E,N & CONVERGENCE
C     --------------------------------------------------
c      WRITE(6,107)LOND,xlon,LATD,xlat,E,N,CD,CM,CS
c107   FORMAT(1x,'NZMGLL ',2X,i3,F8.4,4x,i3,F8.4,2X,2F12.3,2X,2I3,F6.2)
 120   continue
      return
c***** end of subroutine nzmgll *****
      END
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c following subroutines are all for NZMG-LL programs
C	-------------------------
      SUBROUTINE AZETA (Z,ZETA)
C	-------------------------
C	THIS SUBROUTINE COMPUTES AN APPROXIMATE VALUE OF "ZETA" FROM A
C	GIVEN VALUE OF "Z" BY MEANS OF A COMPLEX POWER SERIES
       COMPLEX*16 Z,ZETA,B(6)
      B(1)=(1.3231270439D0,0.0D0)
	B(2)=(-.577245789D0,-.007809598D0)
	B(3)=(0.508307513D0,-.112208952D0)
	B(4)=(-.15094762D0,0.18200602D0)
	B(5)=(1.01418179D0,1.64497696D0)
	B(6)=(1.9660549D0,2.5127645D0)
	
	ZETA=(0.0,0.0)
	DO 10 I=1,6
	ZETA=ZETA+B(I)*Z**I
10      CONTINUE
	RETURN
	END
C	-------------------------
	SUBROUTINE ITZETA (Z,ZETA)
C	--------------------------
C	THIS SUBROUTINE COMPUTES A FINAL VALUE FOR THE COMPLEX NUMBER
C	ZETA BY AN ITERATIVE PROCESS USING 3 ITERATIONS
	COMPLEX*16 Z,ZETA,ZU,ZL,B(6)
	B(1)=(0.7557853228D0,0.000000000D0)
	B(2)=(0.249204646D0 ,0.003371507D0)
	B(3)=(-.001541739D0 ,0.041058560D0)
	B(4)=(-.10162907D0  ,0.017276090D0)
	B(5)=(-.26623489D0  ,-.362492180D0)
	B(6)=(-.6870983D0   ,-1.16519670D0)
	DO 3 I=1,3
	ZU=Z
	DO 1 J=1,5
1	ZU=ZU+J*B(J+1)*ZETA**(J+1)
	ZL=B(1)
	DO 2 J=1,5
2	ZL=ZL+(J+1)*B(J+1)*ZETA**J
3	ZETA=ZU/ZL
	RETURN
	END
C	-------------------------
	SUBROUTINE CILD (PSI,FI)
C	------------------------
C	THIS SUBROUTINE CONVERTS AN ISOMETRIC LATITUDE DIFFERENCE TO A
C	GEODETIC LATITUDE DIFFERENCE BY MEANS OF A POWER SERIES.
	REAL*8 PSI,FI,A(9)
	A(1)=1.5627014243D0
	A(2)=0.5185406398D0
	A(3)=-.03333098D0
	A(4)=-.1052906D0
	A(5)=-.0368594D0
	A(6)=0.007317D0
	A(7)=0.01220D0
	A(8)=0.00394D0
	A(9)=-.0013D0
	FI=0.0D0
	DO 10 I=1,9
	FI=FI+A(I)*PSI**I
10      CONTINUE
	FI=FI/.036D0
	RETURN
	END
C       ------------------------
        SUBROUTINE CGCD (ZETA,Z)
C       ------------------------
C       THIS SUBROUTINE CONVERTS THE COMPLEX NUMBER "ZETA" TO "Z"
C               ZETA=(DPSI+iDLON)
C               Z=(N+iE) WHERE N & E ARE NZMG COORDINATES
        COMPLEX*16 ZETA,Z,B(6),A
        A=(6378388.,0.0)
        B(1)=(0.7557853228,0.0)
        B(2)=(0.249204646 ,0.003371507)
        B(3)=(-.001541739 ,0.041058560)
        B(4)=(-.10162907  ,0.017276090)
        B(5)=(-.26623489  ,-.362492180)
        B(6)=(-.6870983   ,-1.16519670)
        Z=(0.0,0.0)
        DO 993 I=1,6
        Z=Z+B(I)*ZETA**I
993     CONTINUE
        Z=Z*A
        RETURN
        END
C     *************************************************************
        SUBROUTINE CLID(FI,PSI)
C       -----------------------
C       THIS SUBROUTINE CONVERTS A GEODETIC LATITUDE DIFFERENCE TO AN
C       ISOMETRIC LATITUDE DIFFERENCE BY MEANS OF A POWER SERIES
        REAL*8 FI,PSI,A(10)
        A(1)=0.6399175073
        A(2)=-.1358797613
        A(3)=0.063294409
        A(4)=-.02526853
        A(5)=0.0117879
        A(6)=-.0055161
        A(7)=0.0026906
        A(8)=-.001333
        A(9)=0.00067
        A(10)=-.00034
        FI=FI*.036
        PSI=0.0
        DO 992 I=1,10
        PSI=PSI+A(I)*FI**I
992     continue
        RETURN
        END
C     **************************************************************
        SUBROUTINE CDDMS (DEGREE,DEG,MIN,SEC)
C       -------------------------------------
C       THIS SUBROUTINE CONVERTS DEGREES TO DEGREES,MINUTES AND SECONDS
C       USING DOUBLE PRECISION (REAL*8)
        REAL*8    DEGREE,SEC
        INTEGER*4 DEG,MIN
        DEG=IDINT(DEGREE)
        RDEG=DEG
        MIN=IDINT((DEGREE-RDEG)*60.)
        RMIN=MIN
        SEC=DEGREE*3600.-RDEG*3600.-RMIN*60.
        IF(DEG.NE.0) GO TO 9
        IF(MIN.NE.0) GO TO 10
        GO TO 11
9       MIN=IABS(MIN)
10      SEC=DABS(SEC)
11      CONTINUE
        RETURN
        END
C     ***************************************************************
        SUBROUTINE DZDZEE (ZETA,CONV)
C       -----------------------------
C       THIS SUBROUTINE COMPUTES THE COMPLEX NUMBER DZ/DZETA FROM A
C       COMPLEX POWER SERIES AND EVALUATES THE CONVERGENCE
        COMPLEX*16 ZETA,DZETA,B(6),A
        REAL*8     CONV,PI
        PI=1.0
        PI=DATAN(PI)*4.
        A=(6378388.,0.0)
        B(1)=(0.7557853228,0.0)
        B(2)=(0.249204646 ,0.003371507)
        B(3)=(-.001541739 ,0.041058560)
        B(4)=(-.10162907  ,0.017276090)
        B(5)=(-.26623489  ,-.362492180)
        B(6)=(-.6870983   ,-1.16519670)
        DZETA=B(1)
        DO 991 I=1,5
        DZETA=DZETA+(I+1)*B(I+1)*ZETA**I
991     continue
        DZETA=DZETA*A
c        CONV=(DATAN(AIMAG(DZETA)/REAL(DZETA)))*180./PI
        CONV=(DATAN(DIMAG(DZETA)/DREAL(DZETA)))*180./PI
        RETURN
        END
C     **
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
c 
       subroutine ll2spc(latd,xlat,lond,xlon,pnorkm,peaskm)
c  This converts from Lat Long to state plane coordinates
c  Use Alaska zone 4 (5004) and don't use whole table from
c  the SPCS83 program
c  This subroutine assembled from various programs downloaded from
c  NOAA site ftp.ngs.noaa.gov/pub/pcsoft/spcs83
c
c  Since this is just to convert back and forth within the SIMUL
c  program, the northings and eastings are not used in an absolute 
c  sense.  So it does not matter whether wgs83 or earlier as long as
c  all the seismic stations are in the same datum.
c  downloaded from noaa website aug2003
c

C     SccsID = "@(#)spcs83.for	1.2	01/28/02"  
*
*     THIS PROGRAM CONVERTS GPS TO PLANE COORDINATES
*     AND VICE VERSA FOR THE NAD83 DATUM.
*     THIS PROGRAM WAS WRITTEN BY E. CARLSON
*     SUBROUTINES TMGRID,OCONST,SKEWD,LCONST, TCONST,
*     TMGEOD,LAMR1,SKEWR,TCONPC,
*     LAMD1 WERE WRITTEN BY T. VINCENTY, NGS, IN JULY 1984
*     AND LAST UPDATED IN FEBUARY 1986.
*
*     VERSION NUMBER  -  1    09-17-87
*
*
**********************************************************************
*                  DISCLAIMER                                         *
*                                                                     *
*   THIS PROGRAM AND SUPPORTING INFORMATION IS FURNISHED BY THE       *
* GOVERNMENT OF THE UNITED STATES OF AMERICA, AND IS ACCEPTED AND     *
* USED BY THE RECIPIENT WITH THE UNDERSTANDING THAT THE UNITED STATES *
* GOVERNMENT MAKES NO WARRANTIES, EXPRESS OR IMPLIED, CONCERNING THE  *
* ACCURACY, COMPLETENESS, RELIABILITY, OR SUITABILITY OF THIS         *
* PROGRAM, OF ITS CONSTITUENT PARTS, OR OF ANY SUPPORTING DATA.       *
*                                                                     *
*   THE GOVERNMENT OF THE UNITED STATES OF AMERICA SHALL BE UNDER NO  *
* LIABILITY WHATSOEVER RESULTING FROM ANY USE OF THIS PROGRAM.  THIS  *
* PROGRAM SHOULD NOT BE RELIED UPON AS THE SOLE BASIS FOR SOLVING A   *
* PROBLEM WHOSE INCORRECT SOLUTION COULD RESULT IN INJURY TO PERSON   *
* OR PROPERTY.                                                        *
*                                                                     *
*   THIS PROGRAM IS PROPERTY OF THE GOVERNMENT OF THE UNITED STATES   *
* OF AMERICA.  THEREFORE, THE RECIPIENT FURTHER AGREES NOT TO ASSERT  *
* PROPRIETARY RIGHTS THEREIN AND NOT TO REPRESENT THIS PROGRAM TO     *
* ANYONE AS BEING OTHER THAN A GOVERNMENT PROGRAM.                    *
*                                                                     *
***********************************************************************


C     SccsID = "@(#)drgppc.for	1.2	01/28/02"  
c      SUBROUTINE DRGPPC(CARDR,ICODE,FILFLAG,EWFLAG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      real xlat,xlon,pnorkm,peaskm
      DIMENSION SPCC(6)
      CHARACTER*5 GDVAL
      DIMENSION GDVAL(1001),GDNUM(1001)
      LOGICAL FILFLAG
      LOGICAL EWFLAG
      CHARACTER*1 AP
      CHARACTER*4 ZN(135),ZONE
      CHARACTER*80 CARDR
      REAL*8 LAM,NORTH,KP,NB,KC,LATC,LONC,LONO,K,KO,NO
      COMMON/TAB/SPCC
      COMMON/CHAR/ZN,AP
      COMMON/CONST/RAD,ER,RF,ESQ,PI
      COMMON/LATLN/LD,LM,SLAT,LOD,LOM,SLON
      COMMON/FILES/I3,I4,I2,ICON
      COMMON/DONUM/ISN
      COMMON/GEODS/GDNUM,GDVAL

      PI=4.D0*DATAN(1.D0)
      RAD=180.D0/PI
      ER=6378137.D0
      RF=298.257222101D0
      F=1.D0/RF
      ESQ=(F+F-F*F)


c      FI=(LD+(LM+SLAT/60.D0)/60.D0)/RAD
c      LAM=(LOD+(LOM+SLON/60.D0)/60.D0)/RAD
      FI=(latd+(xlat)/60.D0)/RAD
      LAM=(lond+(xlon)/60.D0)/RAD

        IF(EWFLAG) THEN
          LAM = (360.0D0/RAD) - LAM
        ENDIF

c Alaska zone 4 
      icode=5004
      iz=6
      AP='T'
****           AK 4
      SPCC(1)=150.D0
      SPCC(2)=500000.D0
      SPCC(3)=54.D0
      SPCC(4)=10000.D0
      SPCC(5)=0.D0
c      ELSEIF(AP(IZ).EQ.'T') THEN
        CM=SPCC(1)/RAD
        FE=SPCC(2)
        OR=SPCC(3)/RAD
        SF=1.D0-1.D0/SPCC(4)
        FN=SPCC(5)

**** FIND ZONE NAME  ********

      ZONE='AK 4'
*
*
*      COMPUTE  ALL CONSTANCES FOR PROJECTION
*
        CALL TCONST (ER,RF,SF,OR,ESQ,EPS,R,A,B,C,U,V,W,SO,
     &                   CM,FE,FN)
*
*       CONVERT LAT AND LONG TO PCS
*
        CALL TMGRID(FI,LAM,NORTH,EAST,CONV,KP,ER,ESQ,EPS,CM,
     &                  FE,FN,SF,SO,R,A,B,C,U,V,W)
*
c convert to km
      pnorkm=NORTH/1000.0
      peaskm=EAST/1000.0
*      print *,'LL2SPC: lat,lon',latd,xlat,lond,xlon,
*     2   ', pnorkm,peaskm',pnorkm,peaskm
*

  10  CONTINUE
*
*
      RETURN
      END
**********************************************************************
C     SccsID = "@(#)tconst.for	1.2	01/28/02"  
*********************************************************************
      SUBROUTINE TCONST (ER,RF,SF,OR,ESQ,EPS,R,A,B,C,U,V,W,SO,
     &                   CM,FE,FN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C***** TRANSVERSE MERCATOR PROJECTION
C      PRECOMPUTATION OF CONSTANTS
C***** Programmed by T. Vincenty, NGS, in July 1984.
C******************** SYMBOLS AND DEFINITIONS  **********************
C   ER is equatorial radius of the ellipsoid (= major semiaxis).
C   RF is reciprocal of flattening of the ellipsoid.
C   SF is scale factor of the central meridian.
C   OR is southernmost parallel of latitude (in radians) for which
C     the northing coordinate is zero at the central meridian.
C   R, A, B, C, U, V, W are ellipsoid constants used for computing
C     meridional distance from latitude and vice versa.
C   SO is meridional distance (multiplied by the scale factor) from
C     the equator to the southernmost parallel of latitude.
C******************************************************************
C
      F=1./RF
      ESQ=(F+F-F**2)
      EPS=ESQ/(1.-ESQ)
      PR=(1.-F)*ER
      EN=(ER-PR)/(ER+PR)
      A=-1.5D0*EN + (9./16.)*EN**3
      B= 0.9375D0*EN**2 - (15./32.)*EN**4
      C=-(35./48.)*EN**3
      U=1.5D0*EN - (27./32.)*EN**3
      V=1.3125D0*EN**2 - (55./32.)*EN**4
      W=(151./96.)*EN**3
      R=ER*(1.-EN)*(1.-EN**2)*(1.+2.25D0*EN**2+(225./64.)*EN**4)
      OMO=OR + A*SIN(2.*OR) + B*SIN(4.*OR) + C*SIN(6.*OR)
      SO=SF*R*OMO
C
      RETURN
      END
**********************************************************************
C     SccsID = "@(#)tmgrid.for	1.2	01/28/02"  
      SUBROUTINE TMGRID(FI,LAM,NORTH,EAST,CONV,KP,ER,ESQ,EPS,CM,
     &                  FE,FN,SF,SO,R,A,B,C,U,V,W)

      IMPLICIT DOUBLE PRECISION(A-H,K-Z)
C
C*****  TRANSVERSE MERCATOR PROJECTION
C       CONVERSION OF GEODETIC COORDINATES TO GRID COORDINATES
C*****  Programmed by T. Vincenty, NGS, in July 1984.
C*****************  SYMBOLS AND DEFINITIONS *************************
C   Latitude positive north, longitude positive west.  All angles are
C     in radian measure.
C   N, E are northing and easting coordinates respectively.
C   LAT, LON are latitude and longitude respectively.
C   CONV is convergence.
C   KP is point scale factor.
C   ER is equatorial radius of the ellipsoid (= major semiaxis).
C   ESQ is the square of first eccentricity of the ellipsoid.
C   EPS is the square of second eccentricity of the ellipsoid.
C   CM is the central meridian of the projection zone.
C   FE is false easting value at the central meridian.
C   FN is "false northing" at the southernmost latitude, usually zero.
C   SF is scale factor at the central meridian.
C   SO is meridional distance (multiplied by the scale factor) from
C     the equator to the southernmost parallel of latitude for the zone.
C   R is the radius of the rectifying sphere (used for computing
C     meridional distance from latitude and vice versa).
C   A, B, C, U, V, W are other precomputed constants for determination
C     of meridional distance from latitude and vice versa.
C
C   The formula used in this subroutine gives geodetic accuracy within
C   zones of 7 degrees in east-west extent.  Within State transverse
C   Mercator projection zones, several minor terms of the equations
C   may be omitted (see a separate NGS publication).  If programmed
C   in full, the subroutine can be used for computations in surveys
C   extending over two zones.
C
C*********************************************************************
      OM=FI + A*SIN(2.*FI) + B*SIN(4.*FI) + C*SIN(6.*FI)
      S=R*OM*SF
      SINFI=SIN(FI)
      COSFI=COS(FI)
      TN=SINFI/COSFI
      TS=TN**2
      ETS=EPS*COSFI**2
      L=(LAM-CM)*COSFI
      LS=L*L
      RN=SF*ER/SQRT(1.-ESQ*SINFI**2)
C
      A2=RN*TN/2.
      A4=(5.-TS+ETS*(9.+4.*ETS))/12.
      A6=(61.+TS*(TS-58.)+ETS*(270.-330.*TS))/360.
      A1=-RN
      A3=(1.-TS+ETS)/6.
      A5=(5.+TS*(TS-18.)+ETS*(14.-58.*TS))/120.
      A7=(61.-479.*TS+179.*TS**2-TS**3)/5040.
      NORTH=S-SO + A2*LS*(1.+LS*(A4+A6*LS)) +FN
      EAST=FE + A1*L*(1.+ LS*(A3+LS*(A5+A7*LS)))
C
C*** CONVERGENCE
      C1=-TN
      C3=(1.+3.*ETS+2.*ETS**2)/3.
      C5=(2.-TS)/15.
      CONV=C1*L*(1.+LS*(C3+C5*LS))
C
C*** POINT SCALE FACTOR
      F2=(1.+ETS)/2.
      F4=(5.-4.*TS+ETS*( 9.-24.*TS))/12.
      KP=SF*(1.+F2*LS*(1.+F4*LS))
C
      RETURN
      END
c
c - - -- - - - - - - - - - - - - - - - - - - - - - - - - -
c
c
      subroutine spc2ll(pnorkm,peaskm,latd,xlat,lond,xlon)
C     SccsID = "@(#)drpcgp.for	1.2	01/28/02"  
c      SUBROUTINE DRPCGP(CARDR,ICODE,FILFLAG,FILPRT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      real xlat,xlon,pnorkm,peaskm
      LOGICAL FILFLAG,FILPRT
      DIMENSION SPCC(6)
      CHARACTER*1 AP
      CHARACTER*4 ZN(135),ZONE
      CHARACTER*80 CARDR
      REAL*8 LAT,LON,NORTH,KP,NB,KC,LATC,LONC,LONO,K,KO,NO
      COMMON/CHAR/ZN,AP
      COMMON/TAB/SPCC
      COMMON/CONST/RAD,ER,RF,ESQ,PI
      COMMON/XY/NORTH,EAST


      PI=4.D0*DATAN(1.D0)
      RAD=180.D0/PI
      ER=6378137.D0
      RF=298.257222101D0
      F=1.D0/RF
      ESQ=(F+F-F*F)

      FILPRT=.TRUE.
      FILFLAG=.TRUE.
      NORTH=pnorkm*1000.0
      EAST=peaskm*1000.0

c Alaska zone 4 
      icode=5004
      iz=6
      AP='T'
      SPCC(1)=150.D0
      SPCC(2)=500000.D0
      SPCC(3)=54.D0
      SPCC(4)=10000.D0
      SPCC(5)=0.D0
c      ELSEIF(AP(IZ).EQ.'T') THEN
        CM=SPCC(1)/RAD
        FE=SPCC(2)
        OR=SPCC(3)/RAD
        SF=1.D0-1.D0/SPCC(4)
        FN=SPCC(5)

**** FIND ZONE NAME  ********
      ZONE='AK 4'

*** PERFORM TRANSVERSE MERCATOR

c      ELSEIF(AP(IZ).EQ.'T') THEN
        CM=SPCC(1)/RAD
        FE=SPCC(2)
        OR=SPCC(3)/RAD
        SF=1.D0-1.D0/SPCC(4)
        FN=SPCC(5)

*
*      COMPUTE  ALL CONSTANCES FOR PROJECTION
*
        CALL TCONPC (SF,OR,EPS,R,SO,V0,V2,V4,V6,ER,ESQ)
*
*       CONVERT PCS TO LAT AND LONG
*
        CALL TMGEOD(NORTH,EAST,LAT,LON,EPS,CM,FE,SF,SO,R,V0,V2,
     &                  V4,V6,FN,ER,ESQ,CONV,KP)
      CALL TODMS(LAT,LD,LM,SLAT)
      CALL TODMS(LON,LOD,LOM,SLON)
      xlat=LM+SLAT/60.0
      xlon=LOM+SLON/60.0
        latd=int(LD)
        lond=int(LOD)
c      print *,'LL2SPC: lat',latd,xlat,', lon',lond,xlon,
c     2   ', pnorkm,peaskm',pnorkm,peaskm
*
  10  CONTINUE
*
*
      RETURN
      END
**************************************************************
C     SccsID = "@(#)tconpc.for	1.2	01/28/02"  
      SUBROUTINE TCONPC(SF,OR,EPS,R,SO,V0,V2,V4,V6,ER,ESQ)

***          TRANSVERSE MERCATOR PROJECTION               ***
*** CONVERSION OF GRID COORDS TO GEODETIC COORDS
*** REVISED SUBROUTINE OF T. VINCENTY  FEB. 25, 1985
************** SYMBOLS AND DEFINITIONS ***********************
*** ER IS THE SEMI-MAJOR AXIS FOR GRS-80
*** SF IS THE SCALE FACTOR AT THE CM
*** SO IS THE MERIDIANAL DISTANCE (TIMES THE SF) FROM THE
***       EQUATOR TO SOUTHERNMOST PARALLEL OF LAT. FOR THE ZONE
*** R IS THE RADIUS OF THE RECTIFYING SPHERE
*** U0,U2,U4,U6,V0,V2,V4,V6 ARE PRECOMPUTED CONSTANTS FOR
***   DETERMINATION OF MERIDIANAL DIST. FROM LATITUDE
*** OR IS THE SOUTHERNMOST PARALLEL OF LATITUDE FOR WHICH THE
***       NORTHING COORD IS ZERO AT THE CM
**************************************************************

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      F=1.D0/298.257222101D0
      EPS=ESQ/(1.D0-ESQ)
      PR=(1.D0-F)*ER
      EN=(ER-PR)/(ER+PR)
      EN2=EN*EN
      EN3=EN*EN*EN
      EN4=EN2*EN2

      C2=-3.D0*EN/2.D0+9.D0*EN3/16.D0
      C4=15.D0*EN2/16.D0-15.D0*EN4/32.D0
      C6=-35.D0*EN3/48.D0
      C8=315.D0*EN4/512.D0
      U0=2.D0*(C2-2.D0*C4+3.D0*C6-4.D0*C8)
      U2=8.D0*(C4-4.D0*C6+10.D0*C8)
      U4=32.D0*(C6-6.D0*C8)
      U6=128.D0*C8

      C2=3.D0*EN/2.D0-27.D0*EN3/32.D0
      C4=21.D0*EN2/16.D0-55.D0*EN4/32.D0
      C6=151.D0*EN3/96.D0
      C8=1097.D0*EN4/512.D0
      V0=2.D0*(C2-2.D0*C4+3.D0*C6-4.D0*C8)
      V2=8.D0*(C4-4.D0*C6+10.D0*C8)
      V4=32.D0*(C6-6.D0*C8)
      V6=128.D0*C8

      R=ER*(1.D0-EN)*(1.D0-EN*EN)*(1.D0+2.25D0*EN*EN+
     &     (225.D0/64.D0)*EN4)
      COSOR=DCOS(OR)
      OMO=OR+DSIN(OR)*COSOR*(U0+U2*COSOR*COSOR+U4*COSOR**4+
     &    U6*COSOR**6)
      SO=SF*R*OMO

      RETURN
      END
**************************************************************
C     SccsID = "@(#)tmgeod.for	1.2	01/28/02"  
      SUBROUTINE TMGEOD(N,E,LAT,LON,EPS,CM,FE,SF,SO,R,V0,V2,
     &                  V4,V6,FN,ER,ESQ,CONV,KP)

***          TRANSVERSE MERCATOR PROJECTION               ***
*** CONVERSION OF GRID COORDS TO GEODETIC COORDS
*** REVISED SUBROUTINE OF T. VINCENTY  FEB. 25, 1985
************** SYMBOLS AND DEFINITIONS ***********************
*** LATITUDE POSITIVE NORTH, LONGITUDE POSITIVE WEST.  ALL
***          ANGLES ARE IN RADIAN MEASURE.
*** LAT,LON ARE LAT. AND LONG. RESPECTIVELY
*** N,E ARE NORTHING AND EASTING COORDINATES RESPECTIVELY
*** K IS POINT SCALE FACTOR
*** ER IS THE SEMI-MAJOR AXIS FOR GRS-80
*** ESQ IS THE SQUARE OF THE 1ST ECCENTRICITY
*** E IS THE 1ST ECCENTRICITY
*** CM IS THE CENTRAL MERIDIAN OF THE PROJECTION ZONE
*** FE IS THE FALSE EASTING VALUE AT THE CM
*** CONV IS CONVERGENCE
*** EPS IS THE SQUARE OF THE 2ND ECCENTRICITY
*** SF IS THE SCALE FACTOR AT THE CM
*** SO IS THE MERIDIANAL DISTANCE (TIMES THE SF) FROM THE
***       EQUATOR TO SOUTHERNMOST PARALLEL OF LAT. FOR THE ZONE
*** R IS THE RADIUS OF THE RECTIFYING SPHERE
*** U0,U2,U4,U6,V0,V2,V4,V6 ARE PRECOMPUTED CONSTANTS FOR
***   DETERMINATION OF MERIDIANAL DIST. FROM LATITUDE
***
*** THE FORMULA USED IN THIS SUBROUTINE GIVES GEODETIC ACCURACY
*** WITHIN ZONES OF 7 DEGREES IN EAST-WEST EXTENT.  WITHIN STATE
*** TRANSVERSE MERCATOR PROJECTION ZONES, SEVERAL MINOR TERMS OF
*** THE EQUATIONS MAY BE OMMITTED (SEE A SEPARATE NGS PUBLICATION).
*** IF PROGRAMMED IN FULL, THE SUBROUTINE CAN BE USED FOR
*** COMPUTATIONS IN SURVEYS EXTENDING OVER TWO ZONES.
***********************************************************************

      IMPLICIT DOUBLE PRECISION(A-H,K-Z)

      OM=(N-FN+SO)/(R*SF)
      COSOM=DCOS(OM)
      FOOT=OM+DSIN(OM)*COSOM*(V0+V2*COSOM*COSOM+V4*COSOM**4+
     &     V6*COSOM**6)
      SINF=DSIN(FOOT)
      COSF=DCOS(FOOT)
      TN=SINF/COSF
      TS=TN*TN
      ETS=EPS*COSF*COSF
      RN=ER*SF/DSQRT(1.D0-ESQ*SINF*SINF)
      Q=(E-FE)/RN
      QS=Q*Q
      B2=-TN*(1.D0+ETS)/2.D0
      B4=-(5.D0+3.D0*TS+ETS*(1.D0-9.D0*TS)-4.D0*ETS*ETS)/12.D0
      B6=(61.D0+45.D0*TS*(2.D0+TS)+ETS*(46.D0-252.D0*TS-
     &    60.D0*TS*TS))/360.D0
      B1=1.D0
      B3=-(1.D0+TS+TS+ETS)/6.D0
      B5=(5.D0+TS*(28.D0+24.D0*TS)+ETS*(6.D0+8.D0*TS))/120.D0
      B7=-(61.D0+662.D0*TS+1320.D0*TS*TS+720.D0*TS**3)/5040.D0
      LAT=FOOT+B2*QS*(1.D0+QS*(B4+B6*QS))
      L=B1*Q*(1.D0+QS*(B3+QS*(B5+B7*QS)))
      LON=-L/COSF+CM
C*********************************************************************
C     COMPUTE CONVERENCE AND SCALE FACTOR
      FI=LAT
      LAM = LON
      SINFI=SIN(FI)
      COSFI=COS(FI)
      L1=(LAM-CM)*COSFI
      LS=L1*L1
C
C*** CONVERGENCE
      C1=-TN
      C3=(1.+3.*ETS+2.*ETS**2)/3.
      C5=(2.-TS)/15.
      CONV=C1*L1*(1.+LS*(C3+C5*LS))
C
C*** POINT SCALE FACTOR
      F2=(1.+ETS)/2.
      F4=(5.-4.*TS+ETS*( 9.-24.*TS))/12.
      KP=SF*(1.+F2*LS*(1.+F4*LS))

      RETURN
      END
C*********************************************************************
C     SccsID = "@(#)todms.for	1.2	01/28/02"  
*********************************************************************
*
*
      SUBROUTINE TODMS(RAD,IDG,MIN,SEC)
C
C     RADIANS TO DEGREES,MINUTES AND SECONDS
C
      REAL*8 RAD,SEC,RHOSEC
      DATA RHOSEC/2.062648062471D05/
      SEC=RAD*RHOSEC
      IDG=SEC/3600.D0
      SEC=SEC-DBLE(IDG*3600)
      MIN=SEC/60.D0
      SEC=SEC-DBLE(MIN*60)
      IF((60.D0-DABS(SEC)).GT.5.D-6) GO TO 100
      SEC=SEC-DSIGN(60.D0,SEC)
      MIN=MIN+ISIGN(1,MIN)
  100 IF(IABS(MIN).LT.60) GO TO 101
      MIN=MIN-ISIGN(60,MIN)
      IDG=IDG+ISIGN(1,IDG)
  101 MIN=IABS(MIN)
      SEC=DABS(SEC)
      IF(RAD.GE.0.D0) GO TO 102
      IF(IDG.EQ.0) MIN=-MIN
      IF(IDG.EQ.0.AND.MIN.EQ.0)SEC=-SEC
  102 RETURN
      END
