      program specphot
      implicit none

      integer nn
      parameter(nn=100000)

      real*8 m_zp_u, m_zp_b,m_zp_v
      real*8 m_zp_uvm2, m_zp_uvw1,m_zp_uvw2

      parameter(m_zp_uvw2=17.38)
      parameter(m_zp_uvm2=16.85)
      parameter(m_zp_uvw1=17.44)
      parameter(m_zp_u=18.34)
      parameter(m_zp_b=19.11)
      parameter(m_zp_v=17.89)

      real*8 f_u,f_b, f_v, f_uvm2, f_uvw1, f_uvw2
      real*8 m_u,m_b, m_v, m_uvm2, m_uvw1, m_uvw2

      real*8 lam_obj(nn), f_obj0(nn),f_obj(nn)
      real*8 time

      real*8 lam_u(nn), T_u(nn),minlam_u, maxlam_u
      real*8 lam_b(nn), T_b(nn),minlam_b, maxlam_b
      real*8 lam_v(nn), T_v(nn),minlam_v, maxlam_v
      real*8 lam_uvm2(nn), T_uvm2(nn),minlam_uvm2, maxlam_uvm2
      real*8 lam_uvw1(nn), T_uvw1(nn),minlam_uvw1, maxlam_uvw1
      real*8 lam_uvw2(nn), T_uvw2(nn),minlam_uvw2, maxlam_uvw2

      integer n_u, n_b, n_v, n_uvm2, n_uvw1, n_uvw2
      
      integer nobj

      integer i,j,k,m
      real*8 dum,dw,phn

      real*8 dpc,dmpc,zdst,av,dist
      real*8 dmd,ebmv,redshift
      parameter(dmd=30.0) !distance modulus !10Mpc = 30.
      parameter(ebmv=0.0) !E(B-V)
      parameter(redshift=0.0)

      real*8 lcminpc,pi
      parameter(pi=3.141592d0)

      real*8 hh,cc
      parameter(hh=6.626d-27) ! in cgs, erg/s
      parameter(cc=3.0d10) ! in cgs, cm/s
c E = hc/lam
c photon number = lum/E = lum*lam/h/c
c [erg s-1 cm erg-1 s cm-1 s]=[s-1]

      lcminpc=log10(3.085678d18)
      DPC = 10.D0**(1.D0+0.2D0*DMD)
      DMPC = DPC/1.d+6
      write(96,208)DMD,DMPC
 208  format(/' Dist. mod. =',F8.3,'  Dist (Mpc) = ',F8.3)
      write(96,209)EBMV
 209  format(/' Reddening E(B-V) = ',F8.3)
      ZDST = 1.D0 + 0.2D0*DMD+lcminpc
      DIST = 10.D0**(-2.D0*ZDST)/4.D0/pi

      open(unit=21,status="old",file=
     &     "../data/swift/U_UVOT.txt")
      read(21,*)
      n_u=0
      do i=1,nn
         n_u=n_u+1
         read(21,*,end=10) lam_u(i), T_u(i)
      end do
      write(*,*) "Wrong"
      stop
 10   continue
      n_u=n_u-1
      minlam_u=lam_u(1)
      maxlam_u=lam_u(n_u)
      close(21)

      open(unit=22,status="old",file=
     &     "../data/swift/B_UVOT.txt")
      read(22,*)
      n_b=0
      do i=1,nn
         n_b=n_b+1
         read(22,*,end=20) lam_b(i), T_b(i)
      end do
      write(*,*) "Wrong"
      stop
 20   continue
      n_b=n_b-1
      minlam_b=lam_b(1)
      maxlam_b=lam_b(n_b)
      close(22)

      open(unit=23,status="old",file=
     &     "../data/swift/V_UVOT.txt")
      read(23,*)
      n_v=0
      do i=1,nn
         n_v=n_v+1
         read(23,*,end=30) lam_v(i), T_v(i)
      end do
      write(*,*) "Wrong"
      stop
 30   continue
      n_v=n_v-1
      minlam_v=lam_v(1)
      maxlam_v=lam_v(n_v)
      close(23)

      open(unit=24,status="old",file=
     &     "../data/swift/UVM2.txt")
      read(24,*)
      n_uvm2=0
      do i=1,nn
         n_uvm2=n_uvm2+1
         read(24,*,end=40) lam_uvm2(i), T_uvm2(i)
      end do
      write(*,*) "Wrong"
      stop
 40   continue
      n_uvm2=n_uvm2-1
      minlam_uvm2=lam_uvm2(1)
      maxlam_uvm2=lam_uvm2(n_uvm2)
      close(24)

      open(unit=25,status="old",file=
     &     "../data/swift/UVW1.txt")
      read(25,*)
      n_uvw1=0
      do i=1,nn
         n_uvw1=n_uvw1+1
         read(25,*,end=50) lam_uvw1(i), T_uvw1(i)
      end do
      write(*,*) "Wrong"
      stop
 50   continue
      n_uvw1=n_uvw1-1
      minlam_uvw1=lam_uvw1(1)
      maxlam_uvw1=lam_uvw1(n_uvw1)
      close(25)

      open(unit=26,status="old",file=
     &     "../data/swift/UVW2.txt")
      read(26,*)
      n_uvw2=0
      do i=1,nn
         n_uvw2=n_uvw2+1
         read(26,*,end=60) lam_uvw2(i), T_uvw2(i)
      end do
      write(*,*) "Wrong"
      stop
 60   continue
      n_uvw2=n_uvw2-1
      minlam_uvw2=lam_uvw2(1)
      maxlam_uvw2=lam_uvw2(n_uvw2)
      close(26)

c      open(unit=51,status="old",file=
c     &     "../data/vega/Vega_K94_abs_norm_H85.txt")
c     &     "../data/vega/vega_kurucz.txt")
      open(unit=51,status="old",file=
     &     "../../input/tmp.d")
c      open(unit=52,status="old",file=
c     &     "../input/gamtmp.d")
c      write(*,*) "A"
      open(unit=61,status="new",file="./lcresult/multilc.swift.d")



      read(51,*) nobj
      nobj=nobj-1
      read(51,*) (lam_obj(i),i=1,nobj)
c      read(52,*)
c      read(52,*)
         
      do k=1,10000
         read(51,*,end=5000) time,(f_obj0(i),i=1,nobj)
c         lum=0.0d0
         do m=1,nobj
c input: A, erg /s /A /cm2
            AV=EBMV*(3.10D0+2.002D0
     &           * (1.D4/lam_obj(m)-1.D0/0.55D0))
            f_obj(m)=f_obj0(m)*dist !distance            
            f_obj(m)=f_obj(m)*10.D0**(-0.4D0*AV) !extinction
         end do

         f_uvw2=0.0d0
         f_uvm2=0.0d0
         f_uvw1=0.0d0
         f_u=0.0d0
         f_b=0.0d0
         f_v=0.0d0

c uvw2
         do i=2,nobj-1
            if(lam_obj(i).lt.minlam_uvw2) goto 1100
            if(lam_obj(i).gt.maxlam_uvw2) goto 1100
            dw=0.5*(lam_obj(i+1)-lam_obj(i-1)) !wave bin in A
            dum=lam_obj(i)*1.0d-8 !wavelength in cm
            phn=f_obj(i)*dw*dum/hh/cc !photon s-1 cm-2
            do j=1,n_uvw2
               if(lam_obj(i).le.lam_uvw2(j)) then
                  f_uvw2=f_uvw2+phn*T_uvw2(j)
                  goto 1100
               end if
            end do
 1100       continue
         end do
c uvm2
         do i=2,nobj-1
            if(lam_obj(i).lt.minlam_uvm2) goto 1200
            if(lam_obj(i).gt.maxlam_uvm2) goto 1200
            dw=0.5*(lam_obj(i+1)-lam_obj(i-1)) !wave bin in A
            dum=lam_obj(i)*1.0d-8 !wavelength in cm
            phn=f_obj(i)*dw*dum/hh/cc !photon s-1 cm-2
            do j=1,n_uvm2
               if(lam_obj(i).le.lam_uvm2(j)) then
                  f_uvm2=f_uvm2+phn*T_uvm2(j)
                  goto 1200
               end if
            end do
 1200       continue
         end do
c uvw1
         do i=2,nobj-1
            if(lam_obj(i).lt.minlam_uvw1) goto 1300
            if(lam_obj(i).gt.maxlam_uvw1) goto 1300
            dw=0.5*(lam_obj(i+1)-lam_obj(i-1)) !wave bin in A
            dum=lam_obj(i)*1.0d-8 !wavelength in cm
            phn=f_obj(i)*dw*dum/hh/cc !photon s-1 cm-2
            do j=1,n_uvw1
               if(lam_obj(i).le.lam_uvw1(j)) then
                  f_uvw1=f_uvw1+phn*T_uvw1(j)
                  goto 1300
               end if
            end do
 1300       continue
         end do

c u
         do i=2,nobj-1
            if(lam_obj(i).lt.minlam_u) goto 1400
            if(lam_obj(i).gt.maxlam_u) goto 1400
            dw=0.5*(lam_obj(i+1)-lam_obj(i-1)) !wave bin in A
            dum=lam_obj(i)*1.0d-8 !wavelength in cm
            phn=f_obj(i)*dw*dum/hh/cc !photon s-1 cm-2
            do j=1,n_u
               if(lam_obj(i).le.lam_u(j)) then
                  f_u=f_u+phn*T_u(j)
                  goto 1400
               end if
            end do
 1400       continue
         end do

c b
         do i=2,nobj-1
            if(lam_obj(i).lt.minlam_b) goto 1500
            if(lam_obj(i).gt.maxlam_b) goto 1500
            dw=0.5*(lam_obj(i+1)-lam_obj(i-1)) !wave bin in A
            dum=lam_obj(i)*1.0d-8 !wavelength in cm
            phn=f_obj(i)*dw*dum/hh/cc !photon s-1 cm-2
            do j=1,n_b
               if(lam_obj(i).le.lam_b(j)) then
                  f_b=f_b+phn*T_b(j)
                  goto 1500
               end if
            end do
 1500       continue
         end do

c v
         do i=2,nobj-1
            if(lam_obj(i).lt.minlam_v) goto 1600
            if(lam_obj(i).gt.maxlam_v) goto 1600
            dw=0.5*(lam_obj(i+1)-lam_obj(i-1)) !wave bin in A
            dum=lam_obj(i)*1.0d-8 !wavelength in cm
            phn=f_obj(i)*dw*dum/hh/cc !photon s-1 cm-2
            do j=1,n_v
               if(lam_obj(i).le.lam_v(j)) then
                  f_v=f_v+phn*T_v(j)
                  goto 1600
               end if
            end do
 1600       continue
         end do

         m_uvw2=-2.5*log10(f_uvw2)+m_zp_uvw2
         m_uvm2=-2.5*log10(f_uvm2)+m_zp_uvm2
         m_uvw1=-2.5*log10(f_uvw1)+m_zp_uvw1
         m_u=-2.5*log10(f_u)+m_zp_u
         m_b=-2.5*log10(f_b)+m_zp_b
         m_v=-2.5*log10(f_v)+m_zp_v
         
         write(61,4000) time,
     &        m_uvw2, m_uvm2, m_uvw1, m_u, m_b, m_v
 4000    format(f7.2,6(1x,f6.2))
 
         
      end do
 5000 continue
         
      stop
      end
