c wavelength bin output
        do i=1,nfreq-1
            freq0=freq(i)
            freq1=freq(i+1)
            freqa=0.5*(freq0+freq1)
            wava=1.0/freqa
            work(i)=wava
         end do
         write(21,*) nfreq
         write(22,*) nfreq
         write(21,200) (work(i),i=nfreq-1,1,-1) !A
         write(22,200) (work(i),i=nfreq-1,1,-1) !A

c direction binning defenision
         thlc(1) = acos(1.0 - 2.0/real(nthlc))
         do i=2,nthlc-1
            thlc(i)=acos(cos(thlc(i-1)) - 2.0/real(nthlc))
         end do
         thlc(nthlc)=pi

         do i=1,nphilc
            philc(i)=2.0*pi*real(i)/real(nphilc)
         end do
         philc(nphilc)=2.0*pi

c spectrum output
        do j=1,tbin
            t1 = t0*(10.0**(dtlog))
            ta=sqrt(t0*t1)
            write(21,200) ta,(fesum(i,j),i=nfreq-1,1,-1) ! tmp.d
            do l=1,nthlc
               do m=1,nphilc
                  write(22,*) ta,thlc(l),philc(m) ! gamtmp.d
                  write(22,200) (feanglesum(i,j,l,m),i=nfreq-1,1,-1) ! gamtmp.d
               end do
            end do
            t0=t1
         end do
 200     format(5000(1x,e11.4))


