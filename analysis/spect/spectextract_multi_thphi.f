      program spectextract
      implicit none


      integer i,m,k,ip,l,np

      real Distance
      parameter(Distance=10.0) !Mpc
      real dist

      integer nn1d
      real x1d(100000),y1d(100000)

      integer nn
      integer nth,nphi,nangle
      parameter(nth=40,nphi=40)
      parameter(nangle=nth*nphi)
      
      integer nthsum,nphisum
      parameter(nthsum=4,nphisum=4)

      real*8 x0(100000),y0(100000)
      real x(100000),y(100000),time,yy0(100000,nangle)
      real yy(100000,nangle/nthsum/nphisum)
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


c      open(unit=21,status="old",file=
c     &     "../../W7smooth2d_inst_all/SNOTRA/output/x0.4/noitex10/tmp.d")
c      open(unit=22,status="old",file=
c     &     "../../W7smooth2d_inst_all/SNOTRA/output/x0.4/noitex10/gamtmp.d")
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
c input: keV, erg /s /keV
         y0(i)=y0(i)/(x0(i)*1.0e3*ev2erg) !photon
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
c input: keV, erg/s/keV
            y0(i)=y0(i)/(x0(i)*1.0e3*ev2erg) !photon
            y0(i)=y0(i)/4.0/pi/dist/dist !/cm**2.
c            y(i)=y(i)+y0(i)
            yy0(i,m)=y0(i)
         end do
         write(*,*) "day: ",time
         write(*,*) "th: ",th," phi: ",phi
c         if(mod(m,nphi).eq.0) then
c            do i=1,nn-1
c               ip=m/nphi
c               yy(i,ip)=y(i)/real(nphi)
c               y(i)=0.0d0
c            end do
c         end if
c         read(*,*)
      end do
c angle binning:
      do ip=1,nangle/nthsum/nphisum
         do i=1,nn-1
            yy(i,ip)=0.0d0
         end do
      end do
c      do ip=1,nangle/nthsum/nphisum
c         do k=1,nthsum
c            do l=1,nphisum
c               do i=1,nn-1
c                  yy(i,ip)=yy(i,ip)+yy0(i,l)
c               end do
c            end do
c         end do
c      end do
      do i=1,nn-1
         do ip=0,nth/nthsum-1
            np=ip*160
            yy(i,10*ip+1)=yy0(i,np+1)+yy0(i,np+2)
     &           +yy0(i,np+3)+yy0(i,np+4)
     &           +yy0(i,np+1+nth*1)+yy0(i,np+2+nth*1)
     &           +yy0(i,np+3+nth*1)+yy0(i,np+4+nth*1)
     &           +yy0(i,np+1+nth*2)+yy0(i,np+2+nth*2)
     &           +yy0(i,np+3+nth*2)+yy0(i,np+4+nth*2)
     &           +yy0(i,np+1+nth*3)+yy0(i,np+2+nth*3)
     &           +yy0(i,np+3+nth*3)+yy0(i,np+4+nth*3)
            yy(i,10*ip+2)=yy0(i,np+5)+yy0(i,np+6)
     &           +yy0(i,np+7)+yy0(i,np+8)
     &           +yy0(i,np+5+nth*1)+yy0(i,np+6+nth*1)
     &           +yy0(i,np+7+nth*1)+yy0(i,np+8+nth*1)
     &           +yy0(i,np+5+nth*2)+yy0(i,np+6+nth*2)
     &           +yy0(i,np+7+nth*2)+yy0(i,np+8+nth*2)
     &           +yy0(i,np+5+nth*3)+yy0(i,np+6+nth*3)
     &           +yy0(i,np+7+nth*3)+yy0(i,np+8+nth*3)
            yy(i,10*ip+3)=yy0(i,np+9)+yy0(i,np+10)
     &           +yy0(i,np+11)+yy0(i,np+12)
     &           +yy0(i,np+9+nth*1)+yy0(i,np+10+nth*1)
     &           +yy0(i,np+11+nth*1)+yy0(i,np+12+nth*1)
     &           +yy0(i,np+9+nth*2)+yy0(i,np+10+nth*2)
     &           +yy0(i,np+11+nth*2)+yy0(i,np+12+nth*2)
     &           +yy0(i,np+9+nth*3)+yy0(i,np+10+nth*3)
     &           +yy0(i,np+11+nth*3)+yy0(i,np+12+nth*3)
            yy(i,10*ip+4)=yy0(i,np+13)+yy0(i,np+14)
     &           +yy0(i,np+15)+yy0(i,np+16)
     &           +yy0(i,np+13+nth*1)+yy0(i,np+14+nth*1)
     &           +yy0(i,np+15+nth*1)+yy0(i,np+16+nth*1)
     &           +yy0(i,np+13+nth*2)+yy0(i,np+14+nth*2)
     &           +yy0(i,np+15+nth*2)+yy0(i,np+16+nth*2)
     &           +yy0(i,np+13+nth*3)+yy0(i,np+14+nth*3)
     &           +yy0(i,np+15+nth*3)+yy0(i,np+16+nth*3)
            yy(i,10*ip+5)=yy0(i,np+17)+yy0(i,np+18)
     &           +yy0(i,np+19)+yy0(i,np+20)
     &           +yy0(i,np+17+nth*1)+yy0(i,np+18+nth*1)
     &           +yy0(i,np+19+nth*1)+yy0(i,np+20+nth*1)
     &           +yy0(i,np+17+nth*2)+yy0(i,np+18+nth*2)
     &           +yy0(i,np+19+nth*2)+yy0(i,np+20+nth*2)
     &           +yy0(i,np+17+nth*3)+yy0(i,np+18+nth*3)
     &           +yy0(i,np+19+nth*3)+yy0(i,np+20+nth*3)
            yy(i,10*ip+6)=yy0(i,np+21)+yy0(i,np+22)
     &           +yy0(i,np+23)+yy0(i,np+24)
     &           +yy0(i,np+21+nth*1)+yy0(i,np+22+nth*1)
     &           +yy0(i,np+23+nth*1)+yy0(i,np+24+nth*1)
     &           +yy0(i,np+21+nth*2)+yy0(i,np+22+nth*2)
     &           +yy0(i,np+23+nth*2)+yy0(i,np+24+nth*2)
     &           +yy0(i,np+21+nth*3)+yy0(i,np+22+nth*3)
     &           +yy0(i,np+23+nth*3)+yy0(i,np+24+nth*3)
            yy(i,10*ip+7)=yy0(i,np+25)+yy0(i,np+26)
     &           +yy0(i,np+27)+yy0(i,np+28)
     &           +yy0(i,np+25+nth*1)+yy0(i,np+26+nth*1)
     &           +yy0(i,np+27+nth*1)+yy0(i,np+28+nth*1)
     &           +yy0(i,np+25+nth*2)+yy0(i,np+26+nth*2)
     &           +yy0(i,np+27+nth*2)+yy0(i,np+28+nth*2)
     &           +yy0(i,np+25+nth*3)+yy0(i,np+26+nth*3)
     &           +yy0(i,np+27+nth*3)+yy0(i,np+28+nth*3)
            yy(i,10*ip+8)=yy0(i,np+29)+yy0(i,np+30)
     &           +yy0(i,np+31)+yy0(i,np+32)
     &           +yy0(i,np+29+nth*1)+yy0(i,np+30+nth*1)
     &           +yy0(i,np+31+nth*1)+yy0(i,np+32+nth*1)
     &           +yy0(i,np+29+nth*2)+yy0(i,np+30+nth*2)
     &           +yy0(i,np+31+nth*2)+yy0(i,np+32+nth*2)
     &           +yy0(i,np+29+nth*3)+yy0(i,np+30+nth*3)
     &           +yy0(i,np+31+nth*3)+yy0(i,np+32+nth*3)
            yy(i,10*ip+9)=yy0(i,np+33)+yy0(i,np+34)
     &           +yy0(i,np+35)+yy0(i,np+36)
     &           +yy0(i,np+33+nth*1)+yy0(i,np+34+nth*1)
     &           +yy0(i,np+35+nth*1)+yy0(i,np+36+nth*1)
     &           +yy0(i,np+33+nth*2)+yy0(i,np+34+nth*2)
     &           +yy0(i,np+35+nth*2)+yy0(i,np+36+nth*2)
     &           +yy0(i,np+33+nth*3)+yy0(i,np+34+nth*3)
     &           +yy0(i,np+35+nth*3)+yy0(i,np+36+nth*3)
            yy(i,10*ip+10)=yy0(i,np+37)+yy0(i,np+38)
     &           +yy0(i,np+39)+yy0(i,np+40)
     &           +yy0(i,np+37+nth*1)+yy0(i,np+38+nth*1)
     &           +yy0(i,np+39+nth*1)+yy0(i,np+40+nth*1)
     &           +yy0(i,np+37+nth*2)+yy0(i,np+38+nth*2)
     &           +yy0(i,np+39+nth*2)+yy0(i,np+40+nth*2)
     &           +yy0(i,np+37+nth*3)+yy0(i,np+38+nth*3)
     &           +yy0(i,np+39+nth*3)+yy0(i,np+40+nth*3)
         end do
      end do
      

c      if(time.ge.epochmin) then
      write(*,*) time

      do i=1,nn-1
         do ip=1,nangle/nthsum/nphisum
            yy(i,ip)=yy(i,ip)/real(nthsum)/real(nphisum)
         end do
      end do

      do i=1,nn-1
         write(32,555) x(i),(yy(i,ip),ip=1,nangle/nthsum/nphisum)
      end do
c      goto 2000
c      end if

      close(31)
      close(32)
 100  continue 
c      read(*,*)
 555  format(201(1x,1pe12.4))
 1000 continue
 2000 continue

      stop
      end

