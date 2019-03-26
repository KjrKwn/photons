      program specphot
      implicit none

      integer nn
      parameter(nn=100000)

      real*8 fvega_ju, fvega_jb, fvega_jv, fvega_jr, fvega_ji
      real*8 fvega_fu, fvega_fb, fvega_fv, fvega_fr, fvega_fi
      real*8 fvega_jj, fvega_jh, fvega_jk
 
      real*8 lum,mbol

      parameter(fvega_ju=4.2631d-9) !erg/s/cm2/A
      parameter(fvega_jb=6.1855d-9)
      parameter(fvega_jv=3.6021d-9)
      parameter(fvega_jr=2.1422d-9)
      parameter(fvega_ji=1.1117d-9)
      parameter(fvega_jj=3.0949d-10)
      parameter(fvega_jh=1.1534d-10)
      parameter(fvega_jk=4.4936d-11)
      parameter(fvega_fu=3.91476193d-9)
      parameter(fvega_fb=6.18069247d-9)
      parameter(fvega_fv=3.64177726d-9)
      parameter(fvega_fr=2.09980855d-9)
      parameter(fvega_fi=1.12042287d-9)


      real*8 fcal_ju,fcal_jb, fcal_jv, fcal_jr, fcal_ji 
      real*8 mcal_ju,mcal_jb, mcal_jv, mcal_jr, mcal_ji
      real*8 fcal_jj,fcal_jh, fcal_jk
      real*8 mcal_jj,mcal_jh, mcal_jk
      real*8 fcal_fu,fcal_fb, fcal_fv, fcal_fr, fcal_fi 
      real*8 mcal_fu,mcal_fb, mcal_fv, mcal_fr, mcal_fi
      real*8 fcal_fj,fcal_fh, fcal_fk
      real*8 mcal_fj,mcal_fh, mcal_fk
   

      real*8 lam_obj(nn), f_obj0(nn),f_obj(nn)
      real*8 time


      real*8 lam_John_U(nn), T_John_U(nn), minlam_John_U, maxlam_John_U
      real*8 intf_John_U, intT_John_U
      real*8 lam_John_B(nn), T_John_B(nn), minlam_John_B, maxlam_John_B
      real*8 intf_John_B, intT_John_B
      real*8 lam_John_V(nn), T_John_V(nn), minlam_John_V, maxlam_John_V
      real*8 intf_John_V, intT_John_V
      real*8 lam_John_R(nn), T_John_R(nn), minlam_John_R, maxlam_John_R
      real*8 intf_John_R, intT_John_R
      real*8 lam_John_I(nn), T_John_I(nn), minlam_John_I, maxlam_John_I
      real*8 intf_John_I, intT_John_I

      real*8 lam_fcs_U(nn), T_fcs_U(nn), minlam_fcs_U, maxlam_fcs_U
      real*8 intf_fcs_U, intT_fcs_U
      real*8 lam_fcs_B(nn), T_fcs_B(nn), minlam_fcs_B, maxlam_fcs_B
      real*8 intf_fcs_B, intT_fcs_B
      real*8 lam_fcs_V(nn), T_fcs_V(nn), minlam_fcs_V, maxlam_fcs_V
      real*8 intf_fcs_V, intT_fcs_V
      real*8 lam_fcs_R(nn), T_fcs_R(nn), minlam_fcs_R, maxlam_fcs_R
      real*8 intf_fcs_R, intT_fcs_R
      real*8 lam_fcs_I(nn), T_fcs_I(nn), minlam_fcs_I, maxlam_fcs_I
      real*8 intf_fcs_I, intT_fcs_I

      real*8 lam_2m_J(nn), T_2m_J(nn), minlam_2m_J, maxlam_2m_J
      real*8 intf_2m_J, intT_2m_J
      real*8 lam_2m_H(nn), T_2m_H(nn), minlam_2m_H, maxlam_2m_H
      real*8 intf_2m_H, intT_2m_H
      real*8 lam_2m_K(nn), T_2m_K(nn), minlam_2m_K, maxlam_2m_K
      real*8 intf_2m_K, intT_2m_K


      integer nju, njb, njv, njr, nji, n2j, n2h, n2k
      integer nfu, nfb, nfv, nfr, nfi
      
      integer nobj

      integer i,j,k,m
      real*8 dum

      real*8 dpc,dmpc,zdst,av,dist
      real*8 dmd,ebmv,redshift
      parameter(dmd=30.0) !distance modulus !10Mpc = 30.
      parameter(ebmv=0.0) !E(B-V)
      parameter(redshift=0.0)

      real*8 lcminpc,pi
      parameter(pi=3.141592d0)


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
     &     "../data/johnson/U3_res.txt")
      nju=0
      intT_John_U=0.0d0
      do i=1,nn
         nju=nju+1
         read(21,*,end=10) lam_John_U(i), T_John_U(i)
      end do
      write(*,*) "Wrong"
      stop
 10   continue
      nju=nju-1
      do i=2,nju
         intT_John_U=intT_John_U
     &        +T_John_U(i)*0.5*(lam_John_U(i+1)-lam_John_U(i-1))
      end do

      minlam_John_U=lam_John_U(1)
      maxlam_John_U=lam_John_U(nju)


      open(unit=22,status="old",file=
     &     "../data/johnson/B2_res.txt")
      njB=0
      intT_John_B=0.0d0
      do i=1,nn
         njb=njb+1
         read(22,*,end=20) lam_John_B(i), T_John_B(i)
      end do
      write(*,*) "Wrong"
      stop
 20   continue
      njb=njb-1
      do i=2,njb
         intT_John_B=intT_John_B
     &        +T_John_B(i)*0.5*(lam_John_B(i+1)-lam_John_B(i-1))
      end do

      minlam_John_B=lam_John_B(1)
      maxlam_John_B=lam_John_B(njb)

      
      open(unit=23,status="old",file=
     &     "../data/johnson/V_AS_res.txt")
      njv=0
      intT_John_V=0.0d0
      do i=1,nn
         njv=njv+1
         read(23,*,end=30) lam_John_V(i), T_John_V(i)
      end do
      write(*,*) "Wrong"
      stop
 30   continue
      njv=njv-1
      do i=2,njv
         intT_John_V=intT_John_V
     &        +T_John_V(i)*0.5*(lam_John_V(i+1)-lam_John_V(i-1))
      end do

      minlam_John_V=lam_John_V(1)
      maxlam_John_V=lam_John_V(njv)

      
      open(unit=24,status="old",file=
     &     "../data/johnson/R_bessell_res.txt")
      njr=0
      intT_John_R=0.0d0
      do i=1,nn
         njr=njr+1
         read(24,*,end=40) lam_John_R(i), T_John_R(i)
      end do
      write(*,*) "Wrong"
      stop
 40   continue
      njr=njr-1
      do i=2,njr
         intT_John_R=intT_John_R
     &        +T_John_R(i)*0.5*(lam_John_R(i+1)-lam_John_R(i-1))
      end do

      minlam_John_R=lam_John_R(1)
      maxlam_John_R=lam_John_R(njr)

      
      open(unit=25,status="old",file=
     &     "../data/johnson/I_bessell_res.txt")
      nji=0
      intT_John_I=0.0d0
      do i=1,nn
         nji=nji+1
         read(25,*,end=50) lam_John_I(i), T_John_I(i)
      end do
      write(*,*) "Wrong"
      stop
 50   continue
      nji=nji-1
      do i=2,nji
         intT_John_I=intT_John_I
     &        +T_John_I(i)*0.5*(lam_John_I(i+1)-lam_John_I(i-1))
      end do

      minlam_John_I=lam_John_I(1)
      maxlam_John_I=lam_John_I(nji)

      
      open(unit=31,status="old",file=
     &     "../data/focas/focas.txt")
      nfu=0
      intT_fcs_U=0.0d0
      intT_fcs_B=0.0d0
      intT_fcs_V=0.0d0
      intT_fcs_R=0.0d0
      intT_fcs_I=0.0d0
      do i=1,nn
         nfu=nfu+1
         read(31,*,end=100) lam_fcs_U(i), T_fcs_U(i), T_fcs_B(i), 
     &        T_fcs_V(i), T_fcs_R(i), T_fcs_I(i)
      end do
      write(*,*) "Wrong"
      stop
 100  continue
      nfu=nfu-1
      nfb=nfu
      nfv=nfu
      nfr=nfu
      nfi=nfu
      do i=1, nfu+1
         lam_fcs_B(i)=lam_fcs_U(i)
         lam_fcs_V(i)=lam_fcs_U(i)
         lam_fcs_R(i)=lam_fcs_U(i)
         lam_fcs_I(i)=lam_fcs_U(i)
      end do
      do i=2,nfu
         intT_fcs_U=intT_fcs_U
     &        +T_fcs_U(i)*0.5*(lam_fcs_U(i+1)-lam_fcs_U(i-1))
         intT_fcs_B=intT_fcs_B
     &        +T_fcs_B(i)*0.5*(lam_fcs_B(i+1)-lam_fcs_B(i-1))
         intT_fcs_V=intT_fcs_V
     &        +T_fcs_V(i)*0.5*(lam_fcs_V(i+1)-lam_fcs_V(i-1))
         intT_fcs_R=intT_fcs_R
     &        +T_fcs_R(i)*0.5*(lam_fcs_R(i+1)-lam_fcs_R(i-1))
         intT_fcs_I=intT_fcs_I
     &        +T_fcs_I(i)*0.5*(lam_fcs_I(i+1)-lam_fcs_I(i-1))
      end do


      minlam_fcs_U=lam_fcs_U(1)
      maxlam_fcs_U=lam_fcs_U(nfu)
      minlam_fcs_B=lam_fcs_B(1)
      maxlam_fcs_B=lam_fcs_B(nfb)
      minlam_fcs_V=lam_fcs_V(1)
      maxlam_fcs_V=lam_fcs_V(nfv)
      minlam_fcs_R=lam_fcs_R(1)
      maxlam_fcs_R=lam_fcs_R(nfr)
      minlam_fcs_I=lam_fcs_I(1)
      maxlam_fcs_I=lam_fcs_I(nfi)

c 2mass (NIR):

      open(unit=41,status="old",file=
     &     "../data/twomass/jband.2mass.txt")
      n2j=0
      intT_2m_J=0.0d0
      do i=1,nn
         n2j=n2j+1
         read(41,*,end=110) j,lam_2m_J(i), T_2m_J(i)
      end do
      write(*,*) "Wrong"
      stop
 110  continue
      do i=1,n2j
         lam_2m_J(i)=lam_2m_J(i)*1.0d4
      end do
      n2j=n2j-1
      do i=2,n2j
         intT_2m_J=intT_2m_J
     &        +T_2m_J(i)*0.5*(lam_2m_J(i+1)-lam_2m_J(i-1))
      end do

      minlam_2m_J=lam_2m_J(1)
      maxlam_2m_J=lam_2m_J(n2j)
      
      open(unit=42,status="old",file=
     &     "../data/twomass/hband.2mass.txt")
      n2h=0
      intT_2m_H=0.0d0
      do i=1,nn
         n2h=n2h+1
         read(42,*,end=120) j,lam_2m_H(i), T_2m_H(i)
      end do
      write(*,*) "Wrong"
      stop
 120  continue
      do i=1,n2h
         lam_2m_H(i)=lam_2m_H(i)*1.0d4
      end do
      n2h=n2h-1
      do i=2,n2h
         intT_2m_H=intT_2m_H
     &        +T_2m_H(i)*0.5*(lam_2m_H(i+1)-lam_2m_H(i-1))
      end do

      minlam_2m_H=lam_2m_H(1)
      maxlam_2m_H=lam_2m_H(n2h)

      open(unit=43,status="old",file=
     &     "../data/twomass/ksband.2mass.txt")
      n2k=0
      intT_2m_K=0.0d0
      do i=1,nn
         n2k=n2k+1
         read(43,*,end=130) j,lam_2m_K(i), T_2m_K(i)
      end do
      write(*,*) "Wrong"
      stop
 130  continue
      do i=1,n2k
         lam_2m_K(i)=lam_2m_K(i)*1.0d4
      end do
      n2k=n2k-1
      do i=2,n2k
         intT_2m_K=intT_2m_K
     &        +T_2m_K(i)*0.5*(lam_2m_K(i+1)-lam_2m_K(i-1))
      end do

      minlam_2m_K=lam_2m_K(1)
      maxlam_2m_K=lam_2m_K(n2k)


c      open(unit=51,status="old",file=
c     &     "../data/vega/Vega_K94_abs_norm_H85.txt")
c     &     "../data/vega/vega_kurucz.txt")
      open(unit=51,status="old",file=
     &     "../../input/tmp.d")
c      open(unit=52,status="old",file=
c     &     "../input/gamtmp.d")
c      write(*,*) "A"
      open(unit=61,status="new",file="./lcresult/multilc.d")



      read(51,*) nobj
      nobj=nobj-1
      read(51,*) (lam_obj(i),i=1,nobj)
c      read(52,*)
c      read(52,*)
         
      do k=1,10000
         read(51,*,end=5000) time,(f_obj0(i),i=1,nobj)
         lum=0.0d0
         do m=1,nobj
c input: A, erg /s /A
            AV=EBMV*(3.10D0+2.002D0
     &           * (1.D4/lam_obj(m)-1.D0/0.55D0))
            f_obj(m)=f_obj0(m)*dist !distance            
            f_obj(m)=f_obj(m)*10.D0**(-0.4D0*AV) !extinction
            if(m.gt.1) then
               lum=lum+f_obj(m)*(lam_obj(m)-lam_obj(m-1))
            end if
         end do

         intf_John_U=0.0d0
         do i=1,nobj
            if(lam_obj(i).lt.minlam_John_U) goto 1100
            if(lam_obj(i).gt.maxlam_John_U) goto 1100
            do j=1,nju
               if(lam_obj(i).le.lam_John_U(j)) then
                  dum=(T_John_U(j)-T_John_U(j-1))
     &                 /(lam_John_U(j)-lam_John_U(j-1))
     &                 *(lam_obj(i)-lam_John_U(j))+T_John_U(j)
                  intf_John_U=intf_John_U
     &                 +f_obj(i)*dum
     &                 *0.5*(lam_obj(i+1)-lam_obj(i-1))
                  goto 1100
               end if
            end do
 1100       continue
         end do
         
         intf_John_B=0.0d0
         do i=1,nobj
            if(lam_obj(i).lt.minlam_John_B) goto 1200
            if(lam_obj(i).gt.maxlam_John_B) goto 1200
            do j=1,njb
               if(lam_obj(i).le.lam_John_B(j)) then
                  dum=(T_John_B(j)-T_John_B(j-1))
     &                 /(lam_John_B(j)-lam_John_B(j-1))
     &                 *(lam_obj(i)-lam_John_B(j))+T_John_B(j)
                  intf_John_B=intf_John_B
     &                 +f_obj(i)*dum
     &                 *0.5*(lam_obj(i+1)-lam_obj(i-1))
                  goto 1200
               end if
            end do
 1200       continue
         end do
         
         intf_John_V=0.0d0
         do i=1,nobj
            if(lam_obj(i).lt.minlam_John_V) goto 1300
            if(lam_obj(i).gt.maxlam_John_V) goto 1300
            do j=1,njv
               if(lam_obj(i).le.lam_John_V(j)) then
                  dum=(T_John_V(j)-T_John_V(j-1))
     &                 /(lam_John_V(j)-lam_John_V(j-1))
     &                 *(lam_obj(i)-lam_John_V(j))+T_John_V(j)
                  intf_John_V=intf_John_V
     &                 +f_obj(i)*dum
     &                 *0.5*(lam_obj(i+1)-lam_obj(i-1))
                  goto 1300
               end if
            end do
 1300       continue
         end do
         
         intf_John_R=0.0d0
         do i=1,nobj
            if(lam_obj(i).lt.minlam_John_R) goto 1400
            if(lam_obj(i).gt.maxlam_John_R) goto 1400
            do j=1,njr
               if(lam_obj(i).le.lam_John_R(j)) then
                  dum=(T_John_R(j)-T_John_R(j-1))
     &                 /(lam_John_R(j)-lam_John_R(j-1))
     &                 *(lam_obj(i)-lam_John_R(j))+T_John_R(j)
                  intf_John_R=intf_John_R
     &                 +f_obj(i)*dum
     &                 *0.5*(lam_obj(i+1)-lam_obj(i-1))
                  goto 1400
               end if
            end do
 1400       continue
         end do
         
         intf_John_I=0.0d0
         do i=1,nobj
            if(lam_obj(i).lt.minlam_John_I) goto 1500
            if(lam_obj(i).gt.maxlam_John_I) goto 1500
            do j=1,nji
               if(lam_obj(i).le.lam_John_I(j)) then
                  dum=(T_John_I(j)-T_John_I(j-1))
     &                 /(lam_John_I(j)-lam_John_I(j-1))
     &                 *(lam_obj(i)-lam_John_I(j))+T_John_I(j)
                  intf_John_I=intf_John_I
     &                 +f_obj(i)*dum
     &                 *0.5*(lam_obj(i+1)-lam_obj(i-1))
                  goto 1500
               end if
            end do
 1500       continue
         end do
         
         
         intf_fcs_U=0.0d0
         do i=1,nobj
            if(lam_obj(i).lt.minlam_fcs_U) goto 2100
            if(lam_obj(i).gt.maxlam_fcs_U) goto 2100
            do j=1,nfu
               if(lam_obj(i).le.lam_fcs_U(j)) then
                  dum=(T_fcs_U(j)-T_fcs_U(j-1))
     &                 /(lam_fcs_U(j)-lam_fcs_U(j-1))
     &                 *(lam_obj(i)-lam_fcs_U(j))+T_fcs_U(j)
                  intf_fcs_U=intf_fcs_U
     &                 +f_obj(i)*dum
     &                 *0.5*(lam_obj(i+1)-lam_obj(i-1))
                  goto 2100
               end if
            end do
 2100       continue
         end do
         
         intf_fcs_B=0.0d0
         do i=1,nobj
            if(lam_obj(i).lt.minlam_fcs_B) goto 2200
            if(lam_obj(i).gt.maxlam_fcs_B) goto 2200
            do j=1,nfb
               if(lam_obj(i).le.lam_fcs_B(j)) then
                  dum=(T_fcs_B(j)-T_fcs_B(j-1))
     &                 /(lam_fcs_B(j)-lam_fcs_B(j-1))
     &                 *(lam_obj(i)-lam_fcs_B(j))+T_fcs_B(j)
                  intf_fcs_B=intf_fcs_B
     &                 +f_obj(i)*dum
     &                 *0.5*(lam_obj(i+1)-lam_obj(i-1))
                  goto 2200
               end if
            end do
 2200       continue
         end do
         
         intf_fcs_V=0.0d0
         do i=1,nobj
            if(lam_obj(i).lt.minlam_fcs_V) goto 2300
            if(lam_obj(i).gt.maxlam_fcs_V) goto 2300
            do j=1,nfv
               if(lam_obj(i).le.lam_fcs_V(j)) then
                  dum=(T_fcs_V(j)-T_fcs_V(j-1))
     &                 /(lam_fcs_V(j)-lam_fcs_V(j-1))
     &                 *(lam_obj(i)-lam_fcs_V(j))+T_fcs_V(j)
                  intf_fcs_V=intf_fcs_V
     &                 +f_obj(i)*dum
     &                 *0.5*(lam_obj(i+1)-lam_obj(i-1))
                  goto 2300
               end if
            end do
 2300       continue
         end do
         
         intf_fcs_R=0.0d0
         do i=1,nobj
            if(lam_obj(i).lt.minlam_fcs_R) goto 2400
            if(lam_obj(i).gt.maxlam_fcs_R) goto 2400
            do j=1,nfr
               if(lam_obj(i).le.lam_fcs_R(j)) then
                  dum=(T_fcs_R(j)-T_fcs_R(j-1))
     &                 /(lam_fcs_R(j)-lam_fcs_R(j-1))
     &                 *(lam_obj(i)-lam_fcs_R(j))+T_fcs_R(j)
                  intf_fcs_R=intf_fcs_R
     &                 +f_obj(i)*dum
     &                 *0.5*(lam_obj(i+1)-lam_obj(i-1))
                  goto 2400
               end if
            end do
 2400       continue
         end do
         
         intf_fcs_I=0.0d0
         do i=1,nobj
            if(lam_obj(i).lt.minlam_fcs_I) goto 2500
            if(lam_obj(i).gt.maxlam_fcs_I) goto 2500
            do j=1,nfi
               if(lam_obj(i).le.lam_fcs_I(j)) then
                  dum=(T_fcs_I(j)-T_fcs_I(j-1))
     &                 /(lam_fcs_I(j)-lam_fcs_I(j-1))
     &                 *(lam_obj(i)-lam_fcs_I(j))+T_fcs_I(j)
                  intf_fcs_I=intf_fcs_I
     &                 +f_obj(i)*dum
     &                 *0.5*(lam_obj(i+1)-lam_obj(i-1))
                  goto 2500
               end if
            end do
 2500       continue
         end do

         intf_2m_J=0.0d0
         do i=1,nobj
c            write(*,*) lam_obj(i),minlam_2m_J,maxlam_2m_J
            if(lam_obj(i).lt.minlam_2m_J) goto 3100
            if(lam_obj(i).gt.maxlam_2m_J) goto 3100
            do j=1,n2j
c               write(*,*) lam_obj(i),lam_2m_J(j)
               if(lam_obj(i).le.lam_2m_J(j)) then
                  dum=(T_2m_J(j)-T_2m_J(j-1))
     &                 /(lam_2m_J(j)-lam_2m_J(j-1))
     &                 *(lam_obj(i)-lam_2m_J(j))+T_2m_J(j)
                  intf_2m_J=intf_2m_J
     &                 +f_obj(i)*dum
     &                 *0.5*(lam_obj(i+1)-lam_obj(i-1))
                  goto 3100
               end if
            end do
 3100       continue
         end do

         intf_2m_H=0.0d0
         do i=1,nobj
            if(lam_obj(i).lt.minlam_2m_H) goto 3200
            if(lam_obj(i).gt.maxlam_2m_H) goto 3200
            do j=1,n2h
               if(lam_obj(i).le.lam_2m_H(j)) then
                  dum=(T_2m_H(j)-T_2m_H(j-1))
     &                 /(lam_2m_H(j)-lam_2m_H(j-1))
     &                 *(lam_obj(i)-lam_2m_H(j))+T_2m_H(j)
                  intf_2m_H=intf_2m_H
     &                 +f_obj(i)*dum
     &                 *0.5*(lam_obj(i+1)-lam_obj(i-1))
                  goto 3200
               end if
            end do
 3200       continue
         end do

         intf_2m_K=0.0d0
         do i=1,nobj
            if(lam_obj(i).lt.minlam_2m_K) goto 3300
            if(lam_obj(i).gt.maxlam_2m_K) goto 3300
            do j=1,n2k
               if(lam_obj(i).le.lam_2m_K(j)) then
                  dum=(T_2m_K(j)-T_2m_K(j-1))
     &                 /(lam_2m_K(j)-lam_2m_K(j-1))
     &                 *(lam_obj(i)-lam_2m_K(j))+T_2m_K(j)
                  intf_2m_K=intf_2m_K
     &                 +f_obj(i)*dum
     &                 *0.5*(lam_obj(i+1)-lam_obj(i-1))
                  goto 3300
               end if
            end do
 3300       continue
         end do
         
         fcal_ju = intf_John_U/intT_John_U 
         fcal_jb = intf_John_B/intT_John_B 
         fcal_jv = intf_John_V/intT_John_V 
         fcal_jr = intf_John_R/intT_John_R 
         fcal_ji = intf_John_I/intT_John_I 
         fcal_jj = intf_2m_J/intT_2m_J 
         fcal_jh = intf_2m_H/intT_2m_H 
         fcal_jk = intf_2m_K/intT_2m_K 

         fcal_fu = intf_fcs_U/intT_fcs_U 
         fcal_fb = intf_fcs_B/intT_fcs_B 
         fcal_fv = intf_fcs_V/intT_fcs_V 
         fcal_fr = intf_fcs_R/intT_fcs_R 
         fcal_fi = intf_fcs_I/intT_fcs_I 

         
         mcal_ju = -2.5*log10(fcal_ju/fvega_ju)
         mcal_jb = -2.5*log10(fcal_jb/fvega_jb)
         mcal_jv = -2.5*log10(fcal_jv/fvega_jv)
         mcal_jr = -2.5*log10(fcal_jr/fvega_jr)
         mcal_ji = -2.5*log10(fcal_ji/fvega_ji)
         mcal_jj = -2.5*log10(fcal_jj/fvega_jj)
         mcal_jh = -2.5*log10(fcal_jh/fvega_jh)
         mcal_jk = -2.5*log10(fcal_jk/fvega_jk)
         
         mcal_fu = -2.5*log10(fcal_fu/fvega_fu)
         mcal_fb = -2.5*log10(fcal_fb/fvega_fb)
         mcal_fv = -2.5*log10(fcal_fv/fvega_fv)
         mcal_fr = -2.5*log10(fcal_fr/fvega_fr)
         mcal_fi = -2.5*log10(fcal_fi/fvega_fi)

         
         write(*,*) lum
         mbol=-2.5*log10(lum/dist/3.845d33)+4.74
c      write(*,*) "Assuming that Vega is 0 mag in all the bands."

         write(61,4000) time,mbol,  
     &        mcal_ju,mcal_jb,mcal_jv, mcal_jr, mcal_ji,
     &        mcal_jj, mcal_jh, mcal_jk
 4000    format(f7.2,9(1x,f6.2))
 
         
      end do
 5000 continue
         
      stop
      end
