      program spectextract
      implicit none


      integer i,m,k,ip

      real Distance
      parameter(Distance=10.0e0) !Mpc
      real dist

      integer nn1d
      real x1d(100000),y1d(100000)

      integer nn
      integer nth,nphi,nangle
      parameter(nth=40,nphi=40)
      parameter(nangle=nth*nphi)
      
      real*8 x0(10000),y0(10000)
      real x(10000),y(10000),time,yy(10000,nth)
      real xmin,xmax,ymax,ymin
      real th,phi

      real*8 ydum
      real xl(2),yl(2)
      real*8 epochmin,epochmax

      real*8 ytot(100000),yangle(100000,nth)

      integer kfile
      character*1 kfile1
      character*2 kfile2
      character*3 kfile3
      character*3 fname

c const:
      real Mpc2cm
      parameter(Mpc2cm = 3.08e24)
      real ev2erg
      parameter(ev2erg = 1.60217733e-12) ! 1 [ev] = evtoerg [erg]
      real pi
      parameter(pi=3.141592e0)
      
c begin:
c      write(*,*) "Epoch - min?"
c      read(*,*) epochmin
c      write(*,*) "Epoch - max?"
c      read(*,*) rpochmax

      dist=Distance*Mpc2cm !dist [cm]

c      do i=1,100000
c         ytot(i)=0.0d0
c         do j=1,10
c            yangle(i,j)=0.0d0
c         end do
c      end do


      open(unit=21,status="old",file=
     &     "../input/tmp.d")
      open(unit=22,status="old",file=
     &     "../input/gamtmp.d")

      open(unit=51,status="new",file=
     &     "./spectresult/timesequence.dat")

c      open(unit=31,status="new",file=
c     &     "spect.d")
c      open(unit=32,status="new",file=
c     &     "angle.d")

      read(21,*) nn
      read(21,*) (x0(i),i=1,nn-1)
      read(22,*) 
      read(22,*)

c time loop:
      do 1000 k=1,100000
      read(21,*,end=2000) time,(y0(i),i=1,nn-1)

      
      if(k.le.9) then
         write(kfile1,"(i1)") k
         fname="00"//kfile1
      elseif(k.le.99) then
         write(kfile2,"(i2)") k
         fname="0"//kfile2
      elseif(k.le.999) then
         write(kfile3,"(i3)") k
         fname=kfile3
      else
         write(*,*) "Too large data"
         stop
      end if
c      write(fname,"(i3)") kfile

      write(51,*) k,time,
     &     "./spectresult/spect_"//fname//".d"
      open(unit=31,status="new",file=
     &     "./spectresult/spect_"//fname//".d")
      open(unit=32,status="new",file=
     &     "./spectresult/angle_"//fname//".d")

      ymax=-40.0
      do i=1,nn-1
c input: A, erg /s /A
         y0(i)=y0(i)/4.0/pi/dist/dist !/cm**2.
         y(i)=y0(i)
         x(i)=x0(i)
      end do
c      if(time.ge.epochmin) then
c         write(*,*) "Integrate: ", time
      do i=1,nn-1
         ytot(i)=ytot(i)+y(i)
         write(31,*) x(i),y(i)
      end do

      do i=1,nn-1
         y(i)=0.0d0
      end do
      do m=1,nangle
         read(22,*,end=100) time,th,phi
         read(22,*) (y0(i),i=1,nn-1)
         do i=1,nn-1
c input: A, erg/s/A
            y0(i)=y0(i)/4.0/pi/dist/dist !/cm**2.
            y(i)=y(i)+y0(i)
         end do
c         write(*,*) "day: ",time
c         write(*,*) "th: ",th," phi: ",phi
         if(mod(m,nphi).eq.0) then
            do i=1,nn-1
               ip=m/nphi
               yy(i,ip)=y(i)/real(nphi)
               y(i)=0.0d0
            end do
         end if
c         read(*,*)
      end do
c      if(time.ge.epochmin) then
c      write(*,*) time
      do i=1,nn-1
         write(32,555) x(i),(yy(i,ip),ip=1,nth)
      end do
c      goto 2000
c      end if

      close(31)
      close(32)
 100  continue 
c      read(*,*)
 555  format(2001(1x,1pe12.4))
 1000 continue
 2000 continue

      stop
      end

