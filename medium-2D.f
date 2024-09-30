c--   Modules to replace common blocks

C--   medium parameters
      module MEDPARAM
      implicit none      
      INTEGER NF
      DOUBLE PRECISION:: CENTRMIN,CENTRMAX,BREAL,CENTR,RAU=10.,
     & BMIN, BMAX
      CHARACTER*150 HYDRODIR,NCOLLHISTO    
      save
      end module MEDPARAM

C--   internal medium parameters
      module MEDPARAMINT
      implicit none
      DOUBLE PRECISION:: TAUI,TC,ALPHA,BETA,GAMMA,
     $     D3=0.9d0,ZETA3=1.2d0,D,N0,SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,withflow      
      save
      end module MEDPARAMINT      

C--   temperature grid
      module tempgrid 
      implicit none
      integer ntauvals,nxvals,nyvals,tauAllStat,xAllStat,yAllStat
     $     ,tempsAllStat      
      double precision, allocatable :: taus(:),xs(:),ys(:),temps(:,:,:)
      double precision :: tempmax=0.d0
      save
      end module tempgrid        

C--   velocity grids (ux and uy)      
      module velgrid
      implicit none
       double precision, allocatable :: Uxs(:,:,:), Uys(:,:,:),taus2(:)
     $     ,xs2(:),ys2(:)
      double precision umax,unorm
      integer ntauvals2,nxvals2,nyvals2,tauAllStat2,xAllStat2,yAllStat2
     $     ,uxsAllStat,uysAllStat
      save
      end module velgrid      

C--   max rapidity
      module rapmax
      implicit none
      double precision etamax
      save
      end module rapmax 

C--   max rapidity 2
      module rapmax2
      implicit none
      double precision etamax2
      save
      end module rapmax2        

C--   longitudinal boost of momentum distribution
      module boostmed
      implicit none
      logical boost
      save
      end module boostmed

C--   factor to vary Debye mass
      module MDFAC
      implicit none
      DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      save
      end module MDFac

C--   nuclear thickness function
      module THICKFNC
      implicit none
      DOUBLE PRECISION RMAX,TA(100,2)
      save
      end module THICKFNC

C--   geometrical cross section
      module CROSSEC
      implicit none
      DOUBLE PRECISION IMPMAX,CROSS(200,3)
      save
      end module CROSSEC

C--   identifier of log file
      module logfile
      implicit none
      integer logfid
      save
      end module logfile

C--   variables for integration
      module INTEG
      implicit none
      DOUBLE PRECISION B,R      
      save
      end module INTEG

C--   Ncoll grid
      module ncollgrid 
      implicit none
      integer ncollnxvals,ncollnyvals,ncollNZsize,ncollSum,nzc   
      integer ncollxSTAT,ncollySTAT,ncollncollSTAT,ncollNZSTAT,
     & ncollpDistSTAT
      double precision, allocatable :: ncollxs(:),ncollys(:),ncoll(:,:)
      double precision, allocatable :: ncollNZ(:,:)
      double precision, allocatable :: ncollpDist(:,:)
      save
      end module ncollgrid        

      


      SUBROUTINE MEDINIT(FILE,id,etam,intmass)
      use MEDPARAM
      use MEDPARAMINT
      use tempgrid
      use velgrid
      use rapmax2
      use boostmed
      use MDFAC
      use THICKFNC
      use CROSSEC
      use logfile
      

      IMPLICIT NONE
C--   local variables
      INTEGER I,LUN,POS,IOS,id,intmass
      double precision etam
      CHARACTER*100 LABEL,tempbuf
      CHARACTER*150 vlist,tlist
      Character*150 FILE, buffer
      character(len=:), allocatable :: TlistCommand, VlistCommand
      character firstchar
      logical fileexist

      etamax2 = etam
      logfid = id

      IOS=0
      LUN=77

C--   default settings
      TC=0.17d0
      WOODSSAXON=.TRUE.
      CENTRMIN=0.d0
      CENTRMAX=10.d0
      BMIN=0.d0
      BMAX=4.93d0
      NF=3
      A=208   
      N0=0.17d0
      D=0.54d0
      SIGMANN=6.76
      MDFACTOR=0.45d0
      MDSCALEFAC=0.9d0
      withflow=.true.
      HYDRODIR='/hydro/'
      NCOLLHISTO='/ncollhisto/'

C--   read settings from file
      write(logfid,*)
      inquire(file=FILE,exist=fileexist)
      if(fileexist)then
         write(logfid,*)'Reading medium parameters from ',FILE
         OPEN(unit=LUN,file=FILE,status='old',err=10)
         do 20 i=1,1000
            READ(LUN, '(A)', iostat=ios) BUFFER
	    if (ios.ne.0) goto 30
	    firstchar = buffer(1:1)
	    if (firstchar.eq.'#') goto 20
            POS=SCAN(BUFFER,' ')
            LABEL=BUFFER(1:POS)
            BUFFER=BUFFER(POS+1:)
            IF (LABEL=="TC") THEN
               READ(BUFFER,*,IOSTAT=IOS) TC
            ELSE IF (LABEL=="WOODSSAXON") THEN
               READ(BUFFER,*,IOSTAT=IOS) WOODSSAXON
            ELSE IF (LABEL=="CENTRMIN") THEN
               READ(BUFFER,*,IOSTAT=IOS) CENTRMIN
            ELSE IF (LABEL=="CENTRMAX") THEN
               READ(BUFFER,*,IOSTAT=IOS) CENTRMAX
            ELSE IF (LABEL=="BMIN") THEN
               READ(BUFFER,*,IOSTAT=IOS) BMIN
            ELSE IF (LABEL=="BMAX") THEN
               READ(BUFFER,*,IOSTAT=IOS) BMAX   
            ELSE IF (LABEL=="NF") THEN
               READ(BUFFER,*,IOSTAT=IOS) NF
            ELSE IF (LABEL=="A") THEN
               READ(BUFFER,*,IOSTAT=IOS) A
            ELSE IF (LABEL=="N0") THEN
               READ(BUFFER,*,IOSTAT=IOS) N0
            ELSE IF (LABEL=="D") THEN
               READ(BUFFER,*,IOSTAT=IOS) D
            ELSE IF (LABEL=="SIGMANN") THEN
               READ(BUFFER,*,IOSTAT=IOS) SIGMANN
            ELSE IF (LABEL=="MDFACTOR") THEN
               READ(BUFFER,*,IOSTAT=IOS) MDFACTOR
            ELSE IF (LABEL=="MDSCALEFAC") THEN
               READ(BUFFER,*,IOSTAT=IOS) MDSCALEFAC
            ELSE IF (LABEL=="BOOST") THEN
               READ(BUFFER,*,IOSTAT=IOS) withflow
            elseif(label.eq."HYDRODIR")then
               read(BUFFER,'(A)',iostat=ios) hydrodir
            elseif(label.eq."NCOLLHISTO")then
               read(BUFFER,'(A)',iostat=ios) NCOLLHISTO  
               write(logfid,*)'You provided ncollhisto, please '//
     &           'make sure to also provide bmin, bmax,'//
     &           ' centrmin, centrmax' 
               
	    else
               write(logfid,*)'unknown label ',label
	    endif
 20      continue

 30      close(LUN,status='keep')
         write(logfid,*)'...done'
         goto 40

 10      write(logfid,*)'Could not open medium parameter file, '//
     &        'will run with default settings.'

      else
         write(logfid,*)'No medium parameter file found, '//
     &        'will run with default settings.'
      endif

 40   write(logfid,*)'using parameters:'
      write(logfid,*)'TC         = ',TC
      write(logfid,*)'WOODSSAXON = ',WOODSSAXON
      write(logfid,*)'CENTRMIN   = ',CENTRMIN
      write(logfid,*)'CENTRMAX   = ',CENTRMAX
      write(logfid,*)'BMIN       = ',BMIN
      write(logfid,*)'BMAX       = ',BMAX
      write(logfid,*)'NF         = ',NF
      write(logfid,*)'A          = ',A
      write(logfid,*)'N0         = ',N0
      write(logfid,*)'D          = ',D
      write(logfid,*)'SIGMANN    = ',SIGMANN
      write(logfid,*)'MDFACTOR   = ',MDFACTOR
      write(logfid,*)'MDSCALEFAC = ',MDSCALEFAC
      write(logfid,*)'BOOST      = ',withflow
      write(logfid,*)'HYDRODIR   = ',hydrodir
      write(logfid,*)'NCOLLHISTO = ',ncollhisto
      write(logfid,*)

c--   If Ncoll is given, read it in, otherwise compute TA
c--   and the geometric cross section using a model
      if (NCOLLHISTO.eq.'/ncollhisto/') then
         CALL CALCTA  
         CALL CALCXSECTION   
      else
         CALL READNCOLL(NCOLLHISTO)
      endif


      
      
      
C--   Read in hydro profiles:      
      write(*,*)'The hydro profiles I will read are in ',hydrodir

C--   read in temperature profile
      TlistCommand='ls ' // TRIM(hydrodir) //'/Tc* >'
     &// TRIM(hydrodir) // '/Tlist.dat'

      call execute_command_line(TlistCommand)
      Tlist=TRIM(hydrodir)//'/Tlist.dat'
      call readtemps(Tlist)
      write(*,*)'Done reading temperatures.'
      write(*,*)''

C--   read in fluid rapidity profile
      VlistCommand='ls ' // TRIM(hydrodir) //'/Vc* >'
     &// TRIM(hydrodir) // '/Vlist.dat'

      call execute_command_line(VlistCommand)
      Vlist=TRIM(hydrodir)//'/Vlist.dat'
      call readvelocities(Vlist)
      write(*,*)'Done reading velocities.'
      write(*,*)''

      write(*,*)'I read everything, now test for non-zero grids'
      write(*,*)'uxs (2,1,1): ',uxs(2,1,1)
      write(*,*)'uys (2,1,1): ',uys(2,1,1)
      write(*,*)'xs (10): ',xs(10)
      write(*,*)'temps (1,10,10): ',temps(1,10,10)
      
      End


c -------------------------------------------------------------
c   MEDNEXTEVT function - chooses impact parameter and centrality
c -------------------------------------------------------------
      SUBROUTINE MEDNEXTEVT
      use MEDPARAM
      use MEDPARAMINT
      use CROSSEC
      IMPLICIT NONE
C--   local variables
      integer i,j
      DOUBLE PRECISION PYR,R,b1,b2,gettemp

      if (NCOLLHISTO.eq.'/ncollhisto/') then
c--      No ncoll information, use model         
         r=(pyr(0)*(centrmax-centrmin)+centrmin)/100.
         i=0
         do 130 j=1,200
            if ((r-cross(j,3)/cross(200,3)).ge.0.) then
               i=i+1
            else 
               goto 132
            endif
 130  continue
 132  continue
         b1 = (i-1)*0.1d0
         b2 = i*0.1d0
         breal = (b2*(cross(i-1,3)/cross(200,3)-r)
     &     +b1*(r-cross(i,3)/cross(200,3)))/
     &     (cross(i-1,3)/cross(200,3)-cross(i,3)/cross(200,3))
         centr = r;
         
      else
c--    Ncoll data available, so use bmin, bmax, centrmin and centremax
         centr=(pyr(0)*(centrmax-centrmin)+centrmin)/100.  
         breal=(pyr(0)*(bmax-bmin)+bmin)         
      endif

      
      
    
      END

      double precision function getcentrality()
      use MEDPARAM
      implicit none
      getcentrality=centr
      end

c -------------------------------------------------------------
c   PICKVTXNcoll function.  Choose vertex for initial hard scattering
c -------------------------------------------------------------

      SUBROUTINE PICKVTX(X,Y)
      use MEDPARAM
      use ncollgrid
      IMPLICIT NONE
C--   local variables
      DOUBLE PRECISION X,Y,X1,X2,Y1,Y2,Z,XVAL,YVAL,ZVAL,NTHICK,PYR
      double precision ntest, dx, dy
      integer rndIndex

      if (NCOLLHISTO.eq.'/ncollhisto/') then
c--      No Ncoll data given - use a model      
         X1=BREAL/2.-RAU
         X2=RAU-BREAL/2.
         Y1=-SQRT(4*RAU**2-BREAL**2)/2.
         Y2=SQRT(4*RAU**2-BREAL**2)/2.
 131  XVAL=PYR(0)*(X2-X1)+X1
         YVAL=PYR(0)*(Y2-Y1)+Y1
         IF((NTHICK(XVAL-BREAL/2.,YVAL).EQ.0.d0).OR.
     &     NTHICK(XVAL+BREAL/2.,YVAL).EQ.0.d0) GOTO 131
         ZVAL=PYR(0)*NTHICK(-BREAL/2.,0d0)*NTHICK(BREAL/2.,0d0)
         Z=NTHICK(XVAL-BREAL/2.,YVAL)*NTHICK(XVAL+BREAL/2.,YVAL)
         IF(ZVAL.GT.Z) GOTO 131
         X=XVAL
         Y=YVAL            
      else
c--      Ncoll data given, choose vertex by random sampling of weighted distribution:
c--      Readncoll will have also created a list of all the bins with binary collisions
c--      and duplicates for bins with more than one.  
c--      Select a bin randomly.

         rndIndex = int(pyr(0) * ncollSum)
         
         
c--   	Lastly, choose a random location in the histo square located at nci:
c--	   Remember that readncoll stores the midpoint of the histo edges         
         dx=(ncollxs(2)-ncollxs(1))
         dy=dx
         X = ncollNZ(rndIndex,1)+ (pyr(0)-0.5) * dx
         Y = ncollNZ(rndIndex,2)+ (pyr(0)-0.5) * dy
      endif
      END

      SUBROUTINE SETB(BVAL)
      use MEDPARAM
      IMPLICIT NONE

      DOUBLE PRECISION BVAL
      BREAL=BVAL
      END



      SUBROUTINE GETSCATTERER(X,Y,Z,T,TYPE,PX,PY,PZ,E,MS)
      use MEDPARAM
      use MEDPARAMINT
      use logfile
      
      IMPLICIT NONE
C--   local variables
      DOUBLE PRECISION X,Y,Z,T,MS,PX,PY,PZ,E,MD,TEMP,ux,uy,umag
      INTEGER TYPE
      DOUBLE PRECISION R,PYR,pmax,wt,tau,theta,phi,pi,p,yst,pz2,e2,
     &     epar(2),eperp(2),pf,e3,px2,py2,rp,frap,getfrap,ppar,pperp,
     &     pxb,pyb,pzb,pparb
      DATA PI/3.141592653589793d0/
C--   function calls
      double precision getms,getmd,gettemp,getux,getuy
      R=PYR(0)
      IF(R.LT.(2.*12.*NF*D3/3.)/(2.*12.*NF*D3/3.+3.*16.*ZETA3/2.))THEN
         TYPE=2
         MS=GETMS(X,Y,Z,T)
      ELSE
         TYPE=21
         MS=GETMD(X,Y,Z,T)
      ENDIF
      TEMP=GETTEMP(X,Y,Z,T)
      ux=getux(x,y,z,t)
      uy=getuy(x,y,z,t)
      umag=sqrt(ux**2+uy**2)
      pmax = 10.*temp
      
      if (withflow) then
         yst = 0.5*log((t+z)/(t-z))
         frap = getfrap(ux,uy)
      else 
         yst = 0.d0
         frap = 0.d0
      endif
      if (umag.gt.0.d0) then
        epar(1)=ux/umag           !Vectors for transforming to fluid velocity coordinates
        epar(2)=uy/umag
        eperp(1)=-uy/umag
        eperp(2)=ux/umag
      else
        epar(1)=1.d0           !Vectors for transforming to fluid velocity coordinates
        epar(2)=0.d0
        eperp(1)=0.d0
        eperp(2)=1.d0
      endif
      


      IF(TEMP.LT.1.D-2)THEN
         write(logfid,*)'asking for a scattering centre without medium:'
         write(logfid,*)'at (x,y,z,t)=',X,Y,Z,T
         write(logfid,*)'making one up to continue but '//
     &        'something is wrong!'
         TYPE=21
         PX=0.d0
         PY=0.d0
         PZ=0.d0
         MS=GETMS(0.d0,0.d0,0.d0,0.d0)
         MD=GETMD(0.d0,0.d0,0.d0,0.d0)
         E=SQRT(PX**2+PY**2+PZ**2+MS**2)
         RETURN
      ENDIF

 10   p = pyr(0)**0.3333333*pmax
      E3 = sqrt(p**2+ms**2)
      if (type.eq.2) then
         wt = (exp(ms/temp)-1.)/(exp(E3/temp)-1.)
      else
         wt = (exp(ms/temp)+1.)/(exp(E3/temp)+1.)
      endif
      if (wt.gt.1.) write(logfid,*)'Error in getscatterer: weight = ',wt
      if (wt.lt.0.) write(logfid,*)'Error in getscatterer: weight = ',wt
      if (pyr(0).gt.wt) goto 10
C--   momentum in local fluid rest frame
      phi = pyr(0)*2.*pi
      theta = -acos(2.*pyr(0)-1.)+pi
      px2 = p*sin(theta)*cos(phi)
      py2 = p*sin(theta)*sin(phi)
      pz2 = p*cos(theta)
C--   Transform to fluid velocity coordinates
      ppar=px2*epar(1)+py2*epar(2)
      pperp=px2*eperp(1)+py2*eperp(2)
C--   radial boost
      E2=E3*cosh(frap)+ppar*sinh(frap)
      pparb=E3*sinh(frap) +ppar*cosh(frap)
c--   transform boosted momentum to cartesian coordinates
      pxb=pparb*epar(1)+pperp*eperp(1)
   	pyb=pparb*epar(2)+pperp*eperp(2)
C--   longitudinal boost
      E   = cosh(yst)*E2 + sinh(yst)*pz2
      pzb = sinh(yst)*E2 + cosh(yst)*pz2
      pz=pzb
      px=pxb
      py=pyb
      END


      SUBROUTINE AVSCATCEN(X,Y,Z,T,PX,PY,PZ,E,m)
      use boostmed
      use rapmax2
      IMPLICIT NONE
C--   local variables
      double precision x,y,z,t,px,py,pz,e,getms,m,yst,frap,r,ux,uy
c--   function calls
      double precision getfrap,getux,getuy
      
      if (boost) then
         yst = 0.5*log((t+z)/(t-z))
         if ((z.eq.0.d0).and.(t.eq.0.d0)) yst =0.d0
         if (yst.gt.etamax2) yst=etamax2
         if (yst.lt.-etamax2) yst=-etamax2
         ux=getux(x,y,z,t)
         uy=getuy(x,y,z,t)
         frap = getfrap(ux,uy)
      else
         yst = 0.d0
         frap = 0.d0
      endif
      r = sqrt(x**2+y**2)
      m  = getms(x,y,z,t)
      e  = m*cosh(frap)*cosh(yst)
      px = m*sinh(frap)*x/r
      py = m*sinh(frap)*y/r
      pz = m*sinh(yst)
      end


      SUBROUTINE maxscatcen(PX,PY,PZ,E,m)
      use boostmed
      use rapmax
      use rapmax2
      use velgrid
      IMPLICIT NONE
C--   local variables
      double precision px,py,pz,e,getmsmax,m,yst,frap

      if (boost) then
         yst = etamax2
         frap=atanh(umax)
      else
         yst = 0.d0
         frap = 0.d0
      endif
      m  = getmsmax()
      e  = m*cosh(frap)*cosh(yst)
      px = m*sinh(frap)
      py = 0.d0
      pz = m*sinh(yst)
      end
      


      DOUBLE PRECISION FUNCTION GETMD(X1,Y1,Z1,T1)
      use MDFAC
      IMPLICIT NONE

      DOUBLE PRECISION X1,Y1,Z1,T1,GETTEMP
      GETMD=MDSCALEFAC*3.*GETTEMP(X1,Y1,Z1,T1)
      GETMD=MAX(GETMD,MDFACTOR)
      END



      DOUBLE PRECISION FUNCTION GETMS(X2,Y2,Z2,T2)
      IMPLICIT NONE
      DOUBLE PRECISION X2,Y2,Z2,T2,GETMD
      GETMS=GETMD(X2,Y2,Z2,T2)/SQRT(2.)
      END



      DOUBLE PRECISION FUNCTION GETNEFF(X3,Y3,Z3,T3,PX,PY)
      use MEDPARAM
      use MEDPARAMINT
      use rapmax2
      IMPLICIT NONE

C--   local variables
      DOUBLE PRECISION X3,Y3,Z3,T3,PI,frap,tau,cosheta,ux,uy,eta,
     & PX,PY,PZ,PE,costheta,unorm,pnorm
c--   function calls
      double precision getfrap,getux,getuy,gettemp
      
      DATA PI/3.141592653589793d0/

      if (withflow) then
         tau = sqrt(t3**2-z3**2)
         cosheta = t3/tau
         eta = 0.5*log((t3+z3)/(t3-z3))
         ux=getux(x3,y3,z3,t3)
         uy=getuy(x3,y3,z3,t3)
         frap = getfrap(ux,uy)

C--   Find angle between fluid velocity and parton momentum
         unorm = sqrt(ux**2 + uy**2)
         pnorm = sqrt(px**2 + py**2)
         costheta = (ux*px +uy*py)/(unorm*pnorm)

C--   Compute density and scale by (angled) boost
         if (eta.gt.etamax2) then
            getneff = 0.d0
         else
            GETNEFF=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &        *GETTEMP(X3,Y3,Z3,T3)**3/PI**2
            getneff = getneff/(cosh(frap) - sinh(frap)*costheta)
         endif
      else 
         GETNEFF=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &        *GETTEMP(X3,Y3,Z3,T3)**3/PI**2
      endif

C--   Add longitudinal boost
      getneff=getneff/cosheta
      
      END


          

      DOUBLE PRECISION FUNCTION GETTEMP(X,Y,Z,T)
      use MEDPARAM
      use MEDPARAMINT
      use tempgrid
      use logfile
      use rapmax2
      IMPLICIT NONE
C--   local variables
      integer i,j,k,xline,yline,tauline
      DOUBLE PRECISION X,Y,Z,T,TAU,r,deltatau,deltax,deltay,X1A(2),
     &     X2A(2),X3A(2),YA(2,2,2),YB(2,2),A2,B,c,YC(2),YD,tempi,yst,
     &     MA(2,2),MB(2),MC

      GETTEMP=0.D0
      deltatau=abs(taus(2)-taus(1))
      deltax=abs(xs(2)-xs(1))
      deltay=abs(ys(2)-ys(1))

C--   Fail if outside the lightcone
      IF(ABS(Z).GT.T)RETURN

C--   Fail if outside the rapidity range
      yst = 0.5*log((t+z)/(t-z))
      if (abs(yst).gt.etamax2) return


C--   Compute the proper time and radial position
      TAU=SQRT(abs(T**2-Z**2))


      tauline = min(int((tau-taui)/deltatau)+1,ntauvals-1)
      tauline = max(tauline,1)               


      xline = min(int((x-xs(1))/deltax)+1,nxvals-1)
      xline = max(xline,1)


      yline = min(int((y-ys(1))/deltay)+1,nyvals-1)
      yline = max(yline,1)
      
C--   Need vectors containing tau,x and y values at corners

      do 30 i=1,2
         X1A(i)=taus(tauline-1+i)
         X2A(i)=xs(xline-1+i)
         X3A(i)=ys(yline-1+i)
 30   continue
      
      
C--   If tau is before the initial time, ramp up linearly
      if(tau.le.taui)then

C--   Read in temperatures at x,y corner for lowest tau
         do 310 i=1,2
            do 311 j=1,2
               YB(i,j)=temps(1,xline-1+i,yline-1+j)
 311        continue
 310     continue
         

C--   Find gradients on lines of y, estimate T at x on that line
         do 32 i=1,2
            mB(i)=(YB(2,i)-YB(1,i))/(X2A(2)-X2A(1))
            YC(i)=YB(1,i)+mB(i)*(x-X2A(1))
 32      continue

C--   Find gradient between y vals, estimate T at x,y, and ramp down to tau
         MC=(YC(2)-YC(1))/(X3A(2)-X3A(1))
         YD=YC(1)+MC*(y-X3A(1))
         gettemp=YD*(tau/taui)
      else
C--   For tau > taui, need to interpolate in tau direction as well
         do 340 i=1,2
            do 341 j= 1,2
               do 342 k=1,2
                  YA(i,j,k)=temps(tauline-1+i,xline-1+j,yline-1+k)
 342           continue
 341        continue
 340     continue


C--   Find gradients on lines of varying tau, estimate T at tau vals
         do 350 i=1,2
            do 351 j=1,2
               mA(i,j)=(YA(2,i,j)-YA(1,i,j))/(X1A(2)-X1A(1))
               YB(i,j)= YA(1,i,j)+mA(i,j)*(tau-X1A(1))
 351        continue
 350     continue

C--   Find gradients along lines of varying x, estimate T at x vals

         do 37 i=1,2
            mB(i)=(YB(2,i)-YB(1,i))/(X2A(2)-X2A(1))
            YC(i)=YB(1,i)+mB(i)*(x-X2A(1))
 37      continue

         
C--   Find gradient along line of varying y
         mC=(YC(2)-YC(1))/(X3A(2)-X3A(1))


C--   Estimate temp at (tau,x,y)
         GETTEMP=YC(1)+mC*(y-X3A(1))

      end if
      IF(GETTEMP.LT.TC) GETTEMP=0.d0
      END



      DOUBLE PRECISION FUNCTION GETUX(X,Y,Z,T)
      use MEDPARAM
      use MEDPARAMINT
      use velgrid
      use logfile
      IMPLICIT NONE
C--   local variables
      integer i,j,k,lun,pos,ios,tauline,xline,yline
      double precision tau,r,u,deltatau,deltax,deltay,t,x,y,z,myR,yst,
     &     X1A(2),X2A(2),X3A(2),YA(2,2,2),A2,B,YC(2),C,
     &     mA(2,2),mB(2),mC,YB(2,2),pyr,YD



      
      deltatau=abs(taus2(2)-taus2(1))
      deltax=abs(xs2(2)-xs2(1))
      deltay=abs(ys2(2)-ys2(1))

      IF(ABS(Z).GT.T)RETURN

      if (t.lt.z) write(logfid,*) 'error in GETUX: t < z',t,z
      TAU=SQRT(T**2-Z**2)


      tauline = min(int((tau-taui)/deltatau)+1,ntauvals2-1)
      tauline = max(tauline,1)

      
      xline = min(int((x-xs2(1))/deltax)+1,nxvals2-1)
      xline = max(xline,1)

      
      yline = min(int((y-ys2(1))/deltay)+1,nyvals2-1)
      yline = max(yline,1)




      IF(TAU.LT.TAUI)THEN
         getux = 0.
         return
      ELSE

C--   Need vectors containing tau,x,y and u values at corners
         do 40 i=1,2
            X1A(i)=taus2(tauline-1+i)
            X2A(i)=xs2(xline-1+i)
            X3A(i)=ys2(yline-1+i)

 40      continue
         do 410 i=1,2
            do 411 j= 1,2
               do 412 k=1,2
                  YA(i,j,k)=uxs(tauline-1+i,xline-1+j,yline-1+k)
 412           continue
 411        continue
 410     continue
         

C--   Find gradients on lines parallel with tau and estimate T at tau vals
               do 420 i=1,2
                  do 421 j=1,2
                     mA(i,j)=(YA(2,i,j)-YA(1,i,j))/(X1A(2)-X1A(1))
                     YB(i,j)= YA(1,i,j)+mA(i,j)*(tau-X1A(1))
                     
 421              continue
 420           continue

C--   Find gradients along lines parallel with x and estimate T at x vals
                  do 43 i=1,2
                     mB(i)=(YB(2,i)-YB(1,i))/(X2A(2)-X2A(1))
                     YC(i)=YB(1,i)+mB(i)*(x-X2A(1))
                     
 43               continue
C--   Find gradient along line parallel with y
                  mC=(YC(2)-YC(1))/(X3A(2)-X3A(1))
                  GETUx=YC(1)+mC*(y-X3A(1))

               end if

               END





      
      DOUBLE PRECISION FUNCTION GETUY(X,Y,Z,T)
      use MEDPARAM
      use MEDPARAMINT
      use velgrid
      use logfile
      IMPLICIT NONE
C--   local variables
      integer i,j,k,lun,pos,ios,tauline,xline,yline
      double precision tau,r,u,deltatau,deltax,deltay,t,x,y,z,myR,yst,
     &     X1A(2),X2A(2),X3A(2),YA(2,2,2),A2,B,YC(2),C,
     &     mA(2,2),mB(2),mC,YB(2,2),pyr,YD



      
      deltatau=abs(taus2(2)-taus2(1))
      deltax=abs(xs2(2)-xs2(1))
      deltay=abs(ys2(2)-ys2(1))

      IF(ABS(Z).GT.T)RETURN

      if (t.lt.z) write(logfid,*) 'error in GETUy: t < z',t,z
      TAU=SQRT(T**2-Z**2)


      tauline = min(int((tau-taui)/deltatau)+1,ntauvals2-1)
      tauline = max(tauline,1)

      
      xline = min(int((x-xs2(1))/deltax)+1,nxvals2-1)
      xline = max(xline,1)

      
      yline = min(int((y-ys2(1))/deltay)+1,nyvals2-1)
      yline = max(yline,1)




      IF(TAU.LT.TAUI)THEN
         getuy = 0.
         return
      ELSE

C--   Need vectors containing tau,x,y and u values at corners
         do 40 i=1,2
            X1A(i)=taus2(tauline-1+i)
            X2A(i)=xs2(xline-1+i)
            X3A(i)=ys2(yline-1+i)

 40      continue
         do 410 i=1,2
            do 411 j= 1,2
               do 412 k=1,2
                  YA(i,j,k)=uys(tauline-1+i,xline-1+j,yline-1+k)
 412           continue
 411        continue
 410     continue
         

C--   Find gradients on lines parallel with tau and estimate T at tau vals
               do 420 i=1,2
                  do 421 j=1,2
                     mA(i,j)=(YA(2,i,j)-YA(1,i,j))/(X1A(2)-X1A(1))
                     YB(i,j)= YA(1,i,j)+mA(i,j)*(tau-X1A(1))
                     
 421              continue
 420           continue

C--   Find gradients along lines parallel with x and estimate T at x vals
                  do 43 i=1,2
                     mB(i)=(YB(2,i)-YB(1,i))/(X2A(2)-X2A(1))
                     YC(i)=YB(1,i)+mB(i)*(x-X2A(1))
                     
 43               continue
C--   Find gradient along line parallel with y
                  mC=(YC(2)-YC(1))/(X3A(2)-X3A(1))
                  GETUY=YC(1)+mC*(y-X3A(1))

               end if

               END



      

      DOUBLE PRECISION FUNCTION GETFRAP(UX,UY)
      implicit none
      double precision umag,ux,uy

      umag=sqrt(ux**2+uy**2)
      getfrap=atanh(umag)
      end
      

      DOUBLE PRECISION FUNCTION GETMDMAX()
      use MDFAC
      use tempgrid
      IMPLICIT NONE
      GETMDMAX=MDSCALEFAC*3.*tempmax
      GETMDMAX=MAX(GETMDMAX,MDFACTOR)
      END



      DOUBLE PRECISION FUNCTION GETMDMIN()
      use MEDPARAM
      use MEDPARAMINT
      use MDFAC
      IMPLICIT NONE
      GETMDMIN=MDSCALEFAC*3.*TC
      GETMDMIN=MAX(GETMDMIN,MDFACTOR)
      END



      DOUBLE PRECISION FUNCTION GETMSMAX()
      use tempgrid
      IMPLICIT NONE
      DOUBLE PRECISION SQRT
      GETMSMAX=3.*tempmax/SQRT(2.D0)
      END



      DOUBLE PRECISION FUNCTION GETNATMDMIN()
      use MEDPARAM
      use MEDPARAMINT
      use rapmax2
      use velgrid
      use MDFAC
      IMPLICIT NONE
C--   local variables
      DOUBLE PRECISION T,GETMDMIN,pi
      DATA PI/3.141592653589793d0/
      T=GETMDMIN()/(MDSCALEFAC*3.)
      GETNATMDMIN=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &     *T**3/PI**2
!     getnatmdmin=getnatmdmin*(cosh(rapmax) + sinh(rapmax))
      END



      DOUBLE PRECISION FUNCTION GETLTIMEMAX()
      use MEDPARAM
      use MEDPARAMINT
      use rapmax2
      use velgrid
      use tempgrid
      IMPLICIT NONE
C--   function call
      DOUBLE PRECISION getfrap
c--   local variables
      double precision rapmax
      rapmax=atanh(umax)
      GETLTIMEMAX=TAUI*(tempmax/TC)**3*(cosh(rapmax) + sinh(rapmax))
      END



      DOUBLE PRECISION FUNCTION GETNEFFMAX()
      use MEDPARAM
      use MEDPARAMINT
      use velgrid
      use tempgrid
      IMPLICIT NONE
C--   local variables
      DOUBLE PRECISION PI,neffmaxtemp,rapmax
      DATA PI/3.141592653589793d0/

      neffmaxtemp=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &     *tempmax**3/PI**2
      rapmax=atanh(umax)
      getneffmax=neffmaxtemp/(cosh(rapmax)-sinh(rapmax))
      END
      
      

      DOUBLE PRECISION FUNCTION NPART(XX1,YY1,XX2,YY2)
      use medparamint
      IMPLICIT NONE
C--   local variables
      DOUBLE PRECISION XX1,YY1,XX2,YY2,NTHICK     
      NPART = NTHICK(XX1,YY1)*(1.-EXP(-SIGMANN*NTHICK(XX2,YY2))) +
     &     NTHICK(XX2,YY2)*(1.-EXP(-SIGMANN*NTHICK(XX1,YY1)))
      END
      


      DOUBLE PRECISION FUNCTION NTHICK(X1,Y1)
      use MEDPARAM
      use MEDPARAMINT
      use logfile
      use THICKFNC
      IMPLICIT NONE
C--   local variables
      INTEGER LINE,LMIN,LMAX,I
      DOUBLE PRECISION X1,Y1,XA(4),YA(4),Y,DY,R,C,B,DELTA
      
      R=SQRT(X1**2+Y1**2)
      IF(R.GT.TA(100,1))THEN
	 NTHICK=0.
      ELSE
	 LINE=INT(R*99.d0/TA(100,1)+1)
	 LMIN=MAX(LINE,1)
	 LMIN=MIN(LMIN,99)
	 IF((R.LT.TA(LMIN,1)).OR.(R.GT.TA(LMIN+1,1)))
     &        write(logfid,*)LINE,LMIN,R,TA(LMIN,1),TA(LMIN+1,1)
	 XA(1)=TA(LMIN,1)
	 XA(2)=TA(LMIN+1,1)
	 YA(1)=TA(LMIN,2)
	 YA(2)=TA(LMIN+1,2)
	 C=(YA(2)-YA(1))/(XA(2)-XA(1))
	 B=YA(1)-C*XA(1)
	 NTHICK=C*R+B
      ENDIF
      END



      SUBROUTINE CALCTA()
      use MEDPARAM
      use MEDPARAMINT
      use THICKFNC
      use INTEG
      IMPLICIT NONE
C--   local variables
      INTEGER NSTEPS,I
      DOUBLE PRECISION EPS,HFIRST,Y

      NSTEPS=100
      EPS=1.E-4
      HFIRST=0.1D0

      R=1.12*A**(0.33333)-0.86*A**(-0.33333)
      RMAX=2.*R

      DO 10 I=1,NSTEPS
C--   set transverse position
         B=(I-1)*2.D0*R/NSTEPS
         Y=0.D0
C--   integrate along longitudinal line
         CALL ODEINT(Y,-2*R,2*R,EPS,HFIRST,0.d0,101)
         TA(I,1)=B
         TA(I,2)=Y
 10   CONTINUE
      END



      SUBROUTINE CALCXSECTION()
      use MEDPARAM
      use MEDPARAMINT
      use CROSSEC
      IMPLICIT NONE
C--   local variables
      INTEGER IX,IY,IB
      DOUBLE PRECISION B,P,PROD,X,Y,NTHICK,NPART,pprev

      pprev=0.

      DO 30 IB=1,200
         B=0.1d0*IB
         PROD=1.d0
         DO 10 IX=1,100
            DO 20 IY=1,100
               X=-20.d0+IX*0.4d0
               Y=-20.d0+IY*0.4d0
               PROD=PROD*EXP(-NTHICK(X+B/2.D0,Y)*SIGMANN)
     &              **(0.16d0*NTHICK(X-B/2.D0,Y))
 20         CONTINUE
 10      CONTINUE
         P=(1.D0-PROD)*8.8D0/14.D0*B
         CROSS(IB,1)=B
         CROSS(IB,2)=P
         if (ib.eq.1) then
            cross(ib,3)=0.
         else
            cross(ib,3)=cross(ib-1,3)+(p+pprev)/2.*0.1
         endif
         pprev=p
 30   CONTINUE
      IMPMAX=19.95
      END



      DOUBLE PRECISION FUNCTION MEDDERIV(XVAL,W)
      use MEDPARAMINT
      use INTEG
      IMPLICIT NONE
C--   local variables
      DOUBLE PRECISION XVAL
      INTEGER W

      IF (W.EQ.1) THEN
C--   XVAL corresponds to z-coordinate
         MEDDERIV=N0/(1+EXP((SQRT(B**2+XVAL**2)-R)/D))
      ELSE 
         MEDDERIV=0.D0
      ENDIF
      END

      subroutine readtemps(filename)
      use MEDPARAMINT
      use tempgrid
      use logfile
      implicit none
C--   local variables
      integer i,dunit,tlistun,tauc,xc,yc,iolist,iot,dtauc,dxc
     $,diolist,diot,dataFileLength
      logical listexist,fileexist
      character*80 format,ctau
      character*150 tfileph,dtfileph,filename
      double precision xph,yph,tauph,tph,dxph,dyph,dtph

      
C--   Check that the list of data files exists and open it.
      inquire(file=filename,exist=listexist)
      if(listexist)then
         tlistun=13

         open(unit=tlistun,file=filename,status='old')
         write(*,*)'I have opened ',filename
         iolist=0               !Tracks READ errors and looks for the end of the file    
         tauc=0

C--   Determine how many data files there are
         
         dtauc=0
         do
            read(tlistun,95,iostat=diolist) dtfileph

            if(diolist>0) then
               write(*,*)'Something went wrong reading', dtfileph
               write(*,*)diolist
               exit
            else if (diolist<0) then
               exit
            else
               dtauc=dtauc+1
            end if
            
         end do
         ntauvals=dtauc         !set size of tau array
         write(*,*)'I counted the items in Tlist: ',ntauvals
         rewind(tlistun)
         
         
         
c--   Open the last timestamp file and determine its size
         inquire(file=dtfileph,exist=fileexist)
            if(fileexist)then
               Open(unit=20,file=dtfileph,status='old',err=103)
               diot=0
               dxc=0
               

               do
                  read(20,*,iostat=diot) dxph, dyph, dTph !Check if first line contains number
                  if (diot.eq.0) exit
               end do
               
               dxc=dxc+1
            

               do while(diot.eq.0)                
                  read(20,*,iostat=diot) dxph, dyph, dTph
                  dxc=dxc+1
               end do
               rewind(20)
            else
               write(*,*)"The dummy temperature file doesn't exist"
               write(*,*)'I tried to open this dtfileph: ',dtfileph
            endif
            
c     -Compute the number of x and y positions, assuming that there are the
c     -same number of them in the list. Fortran rounds down to integers
            nxvals=int(sqrt(real(dxc)))+1
            nyvals=nxvals

            allocate (taus(ntauvals), STAT=tauAllStat)
            allocate (xs(nxvals), STAT=xAllStat)
            allocate (ys(nyvals), STAT=yAllStat)
            allocate (temps(ntauvals,nxvals,nyvals), STAT=tempsAllStat)
            
            if (tauAllStat.ne.0) then
               STOP
               write(*,*)'Error allocating time array'
            else if (xAllStat.ne.0) then
               STOP
               write(*,*)'Error allocating x array'
            else if (yAllStat.ne.0) then
               STOP
               write(*,*)'Error allocating y array'
            else if (tempsAllStat.ne.0) then
               STOP
               write(*,*)'Error allocating temps array'
            else
               write(*,*)'I allocated tau,x,y, and temps arrays'
            end if

            
 103     close(20,status='keep') !Close the data file

C--   Determine location of timestamp in filename
         dataFileLength=len(TRIM(dtfileph))
C--   Open data files   
         do while (iolist.eq.0)         
            read(tlistun,95,iostat=iolist) tfileph
            if(iolist.ne.0) go to 102
            tauc=tauc+1          ! Increment tau counter
 95         format(A150)
            ctau=tfileph(dataFileLength-5:dataFileLength)  ! Get value of tau from the (full) filename
            read(ctau,94) tauph  
            taus(tauc)=tauph
 94         format(F5.3)
            inquire(file=tfileph,exist=fileexist)
            if(fileexist)then
               dunit=20+tauc
               Open(unit=dunit,file=tfileph,status='old',err=101)
               

C--   Fill arrays.             
               
               xc=1             !Counters to keep track of x and y
               yc=1

               do               !Check if the first line contains a header
                  read(dunit,*,iostat=iot) xph,yph,Tph
                  if(iot.eq.0) exit
               end do
               xs(1)=xph
               ys(1)=yph
               Temps(tauc,1,1)=Tph
               

               iot=0            !Tracks READ errors and looks for the end of the file
               do while (iot.eq.0) !Iterate over all the rows in the data file
                  read(dunit,*,iostat=iot) xph,yph,Tph !Read into place holders
                  
                  if ((xph.eq.xs(xc)).and.(yph.gt.ys(yc))) then
                     yc=yc+1
                     ys(yc)=yph
                     Temps(tauc,xc,yc)=Tph
                  else if ((xph.gt.xs(xc)).and.(yph.eq.ys(yc))) then
                     xc=xc+1
                     xs(xc)=xph
                     Temps(tauc,xc,yc)=Tph
                  else if ((xph.gt.xs(xc)).and.(yph.lt.ys(yc))) then
                     yc=1
                     xc=xc+1
                     ys(yc)=yph
                     xs(xc)=xph
                     Temps(tauc,xc,yc)=Tph
                  else if ((xph.lt.xs(xc)).and.(yph.gt.ys(yc))) then
                     xc=1
                     yc=yc+1
                     ys(yc)=yph
                     xs(xc)=xph
                     Temps(tauc,xc,yc)=Tph                  
                  end if
                  ntauvals=tauc
                  nxvals=xc
                  nyvals=yc
                  if (Tph.gt.tempmax) then
                     tempmax=Tph
                  end if
C--   End of reading in temperatures from tfileph           
               end do

            else
               write(*,*)'This temperature file  does not exist.'
            endif
 101        close(dunit,status='keep') !Close the data file
         end do
         
         
C--   End of reading all time snapshots
 102     close(tlistun,status='keep') !Close the list of files       
         
      else
         write(*,*) 'There is no list of T snapshots'
      end if
      taui=taus(1)

C--   Give warning if no max temperature.      
      if (tempmax.eq.0.d0) then
         write(*,*)'Warning, max temp = 0.'
      else
         write(*,*)'Maximum temperature: ', tempmax
      end if

      end


      subroutine readvelocities(filename)
      use MEDPARAMINT
      use velgrid
      use logfile
      implicit none
C     --   local variables
      integer i,dunit,vlistun,tauc,xc,yc,iolist,iov,dtauc
     $,diolist,diov,dxc,dataFileLength
      logical listexist,fileexist
      character*80 format,ctau
      character*150 filename,vfileph,dvfileph
      double precision xph,yph,tph,uxph,uyph,dxph,dyph,duxph,duyph


c     -check that the list of data files exists and open it      
      inquire(file=filename,exist=listexist)
      if(listexist)then
         vlistun=14

         open(unit=vlistun,file=filename,status='old')
         write(*,*)'I have opened ',filename

c     -Determine how many data files there are
       

         dtauc=0
         do
            read(vlistun,83,iostat=diolist) dvfileph
            if(diolist>0) then
               write(*,*)'Something went wrong reading', dvfileph
               write(*,*)diolist
               exit
            else if (diolist<0) then
               exit
            else
               dtauc=dtauc+1
            end if
            
         end do
         ntauvals2=dtauc         !set size of tau array
         write(*,*)'I counted the items in Vlist: ',ntauvals2
         rewind(vlistun)

c--   Open the last timestamp file and determine its size
         
         inquire(file=dvfileph,exist=fileexist)
         if(fileexist)then
            Open(unit=30,file=dvfileph,status='old',err=200)
            diov=0
            dxc=0

            do
               read(30,*,iostat=diov) dxph, dyph, duxph,duyph !Check if first line contains number
            if (diov.eq.0) exit
         end do
         
         dxc=dxc+1
            
            do while(diov.eq.0)                
               read(30,*,iostat=diov) dxph, dyph, duxph,duyph
               dxc=dxc+1
            end do
            rewind(30)
         else
            write(*,*)"The dummy velocity file doesn't exist"
         endif
 200     close(dunit,status='keep') !Close the data file



c     -Compute the number of x and y positions, assuming that there are the
c     -same number of them in the list. Fortran rounds _down_ to integers.
            nxvals2=int(sqrt(real(dxc)))+1
            nyvals2=nxvals2

            allocate (taus2(ntauvals2), STAT=tauAllStat2)
            allocate (xs2(nxvals2), STAT=xAllStat2)
            allocate (ys2(nyvals2), STAT=yAllStat2)
            allocate (uxs(ntauvals2,nxvals2,nyvals2), STAT=uxsAllStat)
            allocate (uys(ntauvals2,nxvals2,nyvals2), STAT=uysAllStat)
            
            if (tauAllStat2.ne.0) then
               STOP
               write(*,*)'Error allocating time array'
            else if (xAllStat2.ne.0) then
               STOP
               write(*,*)'Error allocating x array'
            else if (yAllStat2.ne.0) then
               STOP
               write(*,*)'Error allocating y array'
            else if (uxsAllStat.ne.0) then
               STOP
               write(*,*)'Error allocating uxs array'
            else if (uysAllStat.ne.0) then
               STOP
               write(*,*)'Error allocating uyx array'
            else
               write(*,*)'I allocated tau, x, y, uxs, uys arrays'
            end if

 

         iolist=0               !Tracks READ errors and looks for the end of the file    
         tauc=0    
         umax=0
     
C--   Determine location of timestamp value in file name
         dataFileLength=len(TRIM(dvfileph))
C--         Step through each file in the list.  
         do while (iolist.eq.0)         
            read(vlistun,83,iostat=iolist) vfileph
            if(iolist.ne.0) go to 202
            tauc=tauc+1
 83         format(A150)
            ctau=vfileph(dataFileLength-5:dataFileLength)
            read(ctau,82) tph
            taus2(tauc)=tph
 82         format(F5.3)
            inquire(file=vfileph,exist=fileexist)
            if(fileexist)then
               dunit=40+tauc
               Open(unit=dunit,file=vfileph,status='old',err=201)

C--   Fill arrays.           
               xc=1             !Counters to keep track of x and y
               yc=1
               

               
               do               !Check if the first line contains a header
                  read(dunit,*,iostat=iov) xph,yph,uxph,uyph
                  if(iov.eq.0) exit
               end do
              
               xs2(1)=xph
               ys2(1)=yph

               Uxs(tauc,1,1)=uxph
               Uys(tauc,1,1)=uyph 

              
               iov=0            !Tracks READ errors and looks for the end of the file
               do while (iov.eq.0) !Iterate over all the rows in the data file
                  read(dunit,*,iostat=iov) xph,yph,Uxph,Uyph !Place holders           
                  if ((xph.eq.xs2(xc)).and.(yph.gt.ys2(yc))) then
                     yc=yc+1
                     ys2(yc)=yph
                     Uxs(tauc,xc,yc)=uxph
                     Uys(tauc,xc,yc)=uyph
                  else if ((xph.gt.xs2(xc)).and.(yph.eq.ys2(yc))) then
                     xc=xc+1
                     xs2(xc)=xph
                     Uxs(tauc,xc,yc)=uxph
                     Uys(tauc,xc,yc)=uyph
                  else if ((xph.gt.xs2(xc)).and.(yph.lt.ys2(yc))) then
                     yc=1
                     xc=xc+1
                     ys2(yc)=yph
                     xs2(xc)=xph
                     Uxs(tauc,xc,yc)=uxph
                     Uys(tauc,xc,yc)=uyph
                  else if ((xph.lt.xs2(xc)).and.(yph.gt.ys2(yc))) then
                     xc=1
                     yc=yc+1
                     ys2(yc)=yph
                     xs2(xc)=xph
                     Uxs(tauc,xc,yc)=uxph
                     Uys(tauc,xc,yc)=uyph
                  end if

                  unorm=sqrt(uxph**2+uyph**2)
                  if(unorm.gt.umax)umax=unorm
                  ntauvals2=tauc
                  nxvals2=xc
                  nyvals2=yc
C--   End of reading in velocities from vfileph               
               end do

            else
               write(*,*)'This velocity file  does not exist: ',vfileph               
            endif
 201        close(dunit,status='keep') !Close the data file

            
         end do
        

         
C--   End of reading all time snapshots
 202     close(vlistun,status='keep') !Close the lsit of files

         
         
      else
         write(*,*) 'There is no list of V snapshots'
      end if
      end

c -------------------------------------------------------------
c   READNCOLL function - reads Ncollhistogram
c -------------------------------------------------------------      
      SUBROUTINE READNCOLL(filename)
        use MEDPARAM
        use MEDPARAMINT
        use THICKFNC
        use INTEG
		   use ncollgrid
        IMPLICIT NONE
C--   local variables
        INTEGER NSTEPS,I, histoun, iolist, diot, dxc, xc, yc, iot,
     &		nzci, nzcj, ncollint, ncollpi
        DOUBLE PRECISION EPS,HFIRST,Y1,Y2
         logical histoexists
         character*150 filename
         double precision xph, yph, ncollph, dxph, dyph, dncollph

        NSTEPS=100
        EPS=1.E-4
        HFIRST=0.1D0

C--   Check that the NColl histogram file exists and open it.
		inquire(file=filename,exist=histoexists)
		if(histoexists)then
			histoun=30

			open(unit=histoun,file=filename,status='old')
			write(*,*)'I have opened ',filename
			iolist=0               !Tracks READ errors and looks for the end of the file    
			
c--	  See how long this file is
			
			do
				read(histoun,*,iostat=diot) dxph, dyph, dncollph !Check if first line contains number
				if (diot.eq.0) exit
			end do
			
			dxc=dxc+1
		
			do while(diot.eq.0)                
				read(histoun,*,iostat=diot) dxph, dyph, dncollph
				dxc=dxc+1
			end do
			rewind(histoun)

c     -Compute the number of x and y positions, assuming that there are the
c     -same number of them in the list. Fortran rounds down to integers
            ncollnxvals=int(sqrt(real(dxc)))
            ncollnyvals=ncollnxvals

            allocate (ncollxs(ncollnxvals), STAT=ncollxSTAT)
            allocate (ncollys(ncollnyvals), STAT=ncollySTAT)
            allocate (ncoll(ncollnxvals,ncollnyvals), 
     &		 	STAT=ncollncollSTAT)

c	  - Allocate array for non-zero ncoll elements to pick vertex from
	 		ncollNZsize=ncollnxvals*ncollnyvals
			allocate (ncollNZ(ncollNZsize,3), STAT=ncollNZSTAT)	

c	  - Write out errors if anything went wrong
			write(*,*)"ncollhisto.dat has size ",dxc
			write(*,*)"for ncollyvals = ncollxvals: ", ncollnyvals 
			
			if (ncollxSTAT.ne.0) then
               STOP
               write(*,*)'Error allocating x array'
            else if (ncollySTAT.ne.0) then
               STOP
               write(*,*)'Error allocating y array'
            else if (ncollncollSTAT.ne.0) then
               STOP
               write(*,*)'Error allocating ncoll array'
            else
               write(*,*)'I allocated x,y, and ncoll arrays'
            end if

C--   Fill arrays.             
               
			xc=1             !Counters to keep track of x and y
			yc=1
			nzc=1
			ncollSum=0

			do               !Check if the first line contains a header
				read(histoun,*,iostat=diot) xph, yph, ncollph
				if(diot.eq.0) exit
			end do
			ncollxs(1)=xph
			ncollys(1)=yph
			ncoll(1,1)=ncollph


			

			diot=0            !Tracks READ errors and looks for the end of the file
			do while (diot.eq.0) !Iterate over all the rows in the data file
				read(histoun,*,iostat=diot) xph,yph,ncollph !Read into place holders
				
				if ((xph.eq.ncollxs(xc)).and.(yph.gt.ncollys(yc))) then
					yc=yc+1
					ncollys(yc)=yph
					ncoll(xc,yc)=ncollph
				else if ((xph.gt.ncollxs(xc)).and.(yph.eq.ncollys(yc))) then
					xc=xc+1
					ncollxs(xc)=xph
					ncoll(xc,yc)=ncollph
				else if ((xph.gt.ncollxs(xc)).and.(yph.lt.ncollys(yc))) then
					yc=1
					xc=xc+1
					ncollys(yc)=yph
					ncollxs(xc)=xph
					ncoll(xc,yc)=ncollph
				else if ((xph.lt.ncollxs(xc)).and.(yph.gt.ncollys(yc))) then
					xc=1
					yc=yc+1
					ncollys(yc)=yph
					ncollxs(xc)=xph
					ncoll(xc,yc)=ncollph                  
				end if

				if(ncollph.gt.0.) then
					
					ncollNZ(nzc,1) = xph
					ncollNZ(nzc,2) = yph
					ncollNZ(nzc,3) = ncollph
					ncollSum=ncollSum+ncollph
					nzc=nzc+1
				end if
C--   End of reading in Ncoll list 
			end do

C--	  Create ncoll pdist
			allocate (ncollpDist(ncollSum,2), STAT=ncollpDistSTAT)	
			if (ncollpDistSTAT.ne.0) then
               STOP
               write(*,*)'Error allocating ncoll probability array'	
			else
				write(*,*)'I allocated ncoll probability array'
			endif

			nzci = 1
			ncollpi = 1
C--	  Place one entry for each binary collision
			do while (nzci.le.nzc)
				nzcj = 1
				ncollint=int(ncollNZ(nzci,3))
				
				do while (nzcj.le.ncollint)
					ncollpDist(ncollpi,1) = ncollNZ(nzci,1)
					ncollpDist(ncollpi,2) = ncollNZ(nzci,2)
					nzcj = nzcj + 1
					ncollpi = ncollpi +1
				end do
				nzci = nzci +1
			end do         			
			write(*,*)"Number of non-zero entries: " , nzc 
			close(histoun,status='keep')

		else
			write(*,*)"WARNING! I could not open NcollHisto, use a model"
			call CALCTA
		endif

			
10      CONTINUE
      END   

