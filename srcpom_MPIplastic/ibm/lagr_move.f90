# 1 "lagr_move.F90"

subroutine lagr_move (dt1)
!!  3-D advection-diffusion lagrangian model for plastics
    include 'pom.h'
    include 'lagr.h'

!!  local variables
    double precision :: r1,r2,ra,rx,ry               !random part used for stochastic processes
    integer :: nsp                                   !micro, macro
    integer :: ivariant                                  !micro/macro sub-class
!   1-D Arrays
    real :: dia,dia0,dia02                           !particle diameter
    real :: den,den0                                 !particle density
    real :: rhom                                     !water density
!   wind/wave flags
    logical :: waves,wind
!   velocity
    real :: wmix,w_cur,wt,wb                         !vertical velocity components
    real :: uloc,uwave,uwind,umix                    !horizontal velocity-x components
    real :: vloc,vwave,vwind,vmix                    !horizontal velocity-y components
!   waves
    real :: stx0,sty0,tsk,wk,wl,hsk
!   wind
    real ::windx,windy,AirSurfRatio

    write(*,*)'start lagr'

!   main loop for all SIs
    do 200 k = 1, nsi

        if(variant(k) == 0)goto 200                      !skip for inactive particles
                 
        ivariant = variant(k)
        nsp = species(k)

        if(ref(k) > 0.)dista(k)=0.                   !zero the distance covered for beached particles
        if(ref(k) > 0.) goto 270                     !particles on the beach, goto beaching

!        if(ref(k) .eq. -1.) then
!            xn(k) = xo(k)
!            yn(k) = yo(k)
!            zn(k) = zo(k)
!            goto 240
!        endif

!calculate (i,j) from current position(xo,yo)-------------------------------
        call FCoord (io, jo, alon(1,1)-(stepx/2.), &    
        alat(1,1)-(stepy/2.), stepx, stepy, &
        xo(k), yo(k), pi, dg)

!   Buoyancy velocity, according to critical density.
!   Distinction between small and large particles

!   Get wind, wave etc at the SI location(currently=(io,jo))----------------
!   Waves
        waves=wavesp(ivariant,nsp)                       !wave activation flag (see plasticparams.dat)
        hsk=wavehei(io,jo)                           !significant wave height(m)
        if(hsk > 0.)then
            tsk=waveper(io,jo)                       !wave period(s)
            stx0=stokesx(io,jo)                      !stokes 1-x at surface
            sty0=stokesy(io,jo)                      !stokes 1-y at surface
            wl = grav * tsk  * tsk / (2. * pi)       !wave length
            wk = 2. * pi / wl                        !wave number
        else
            stx0=0.
            sty0=0.
            wk = 0.
        endif
!   Wind
        wind=windp(ivariant,nsp)                         !wind activation flag (see plasticparams.dat)
        windx=uair(io,jo)                            !wind-x at 10m above sea
        windy=vair(io,jo)                            !wind-y at 10m above sea
        AirSurfRatio = A_air(ivariant,nsp)/A_water(ivariant,nsp) ! particle ratio of surface above sea/total
        
!  Get buoyancy related parameters---------------------------------------------
        call vertical_val(floc,io,jo,im,jm,kb,zz,rho,h,zo(k)) ! water density
!        rhom = floc*1025.+1000.
        rhom = 1025.
        dia0 = diam(ivariant,nsp)                        !particle diameter
        dia02= 0.5*dia0                              !particle radius
        den0 = dens(ivariant,nsp)                        !particle initial density

!--get biofilm thickness (hfou, micro) or particle density (macro)
        if(nsp == 1)then          !Micro
            hfou = fou(k)
        else                      !Macro
            hfou = 0.
            if(ivariant == 3)den0 = fou(k)               !get macro density (saved in fou)
            if(ivariant <= 2 .OR. ivariant == 4)then         !sinking of small macro & bags based on Tsat(plasticparams.dat)
                if(icountdays(k) > Tsat(ivariant,nsp))then
                    den0=2300.                       !increase density to sink
                    fou(k)=den0
                endif
            endif
!            den0 = dens(ivariant,nsp)                   !test no sinking
            goto 345                                 !skip biofouling for macro
        endif

                  
! microplastics biofouling/degradation--------------------------------------
      call vertical_val(floc,io,jo,im,jm,kb,zz,bac3d,h,zo(k)) !get bacteria at SIs location
       bac = floc
!      bac = 0.  !no biofouling
      if(bac > bacmax)bac=bacmax

      if(nsp == 1)then
       hth=dia02*(((denf-den0)/(denf-rhom))**0.333-1) !critical thickness (Cubarenko 2016)
       dh0=dh0max50*50.e-6/dia0                       !max biofilm growth dh0=f(1/size) (=dh0max50 for diam=50um)
       def0=0.92*dh0/(hth+3.e-6)                      !varying defouling=f(hth) (Tsiaras 2021)
!       def0=0.05                                     !constant defouling (calculated from Kiorboe 2003)
       dhfou=dh0*bac/bacmean                          !bacteria dependence
       hfou=hfou+(dhfou-def0*hfou)*dt1/86400.         !update biofouling thickness (hfou)
       if(hfou < 0.)hfou=0.
        fou(k)=hfou
       endif
  345 continue

      dia2 = dia02+hfou                              !update radius with biofilm thickness
      dia  = 2.*dia2                                 !particle diameter
         
!--calculate density
        if(nsp == 1)then                             !Micro
            den=(den0*dia02**3.+denf*(dia2**3.-dia02**3.))/(dia2**3.) !Cubarenko 2016 for sphere
            cd=1.                                    !drag coefficient
        else                                         !Macro
            den =den0
            cd=1.
            if(ivariant == 4)then                        !Bags
                cd=0.3
            endif
        endif

!        den=dens(ivariant,nsp)                          !test no biofouling
!        dia = diam(ivariant,nsp)                        !test no biofouling

!---deposition-resuspension
!        if(ref(k).eq.-1)then
        bstress=SQRT(wubot(io,jo)**2+wvbot(io,jo)**2)*1023.
        depstress=0.4
!             erostress=0.025
!             erostress=0.01
!             depstress=100.
!             erostress=0.

        if(bstress < erostress)then
            erosion=0.
        else
            ero = rand(0)
            if(ero < (1.-erostress/bstress))then
                erosion=1.
            else
                erosion=0.
            endif
        endif
!        erosion=0.
!        endif
!        erosion=0.

!  check if at bottom and skip movement when erosion=0
        if(ref(k) == -1. .AND. erosion == 0.)then 
            xn(k) = xo(k)
            yn(k) = yo(k)
            zn(k) = zo(k)
            goto 240
        endif

!  Vertical movement>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!  calculate buoyancy velocity wb = f(dia,den,rhom) based on Elliot (1986)--------------------------------
        dc = 9.52*vk**(2./3.)/(grav**(1./3.)*(1.-den/rhom)** (1./3.)) !critical diameter
        if(dia <= dc) then
            wb=cd*grav*dia**2*(1.-den/rhom)/(18.*vk)                  !vertical velovity (<0 means sinking)
        else
            wb=cd*8./3.*grav*dia*sign(1.,1.-den/rhom)*sqrt(abs(1.-den/rhom))
        end if

!---turbulent vertical velocity wt = f(Kh,waves)------------------------------
        call vertical_val(floc, io, jo, im, jm, kb, zz, &
        kh, h, zo(k))
        if (waves .AND. tsk > 0.) then
            akz = min(-2. * wk * zo(k), 50.)
            dit = 0.028 * (hsk **2) / tsk * exp(-akz)
            difv = dit + floc
        else
            difv = floc
        endif
                    
        wt = sqrt(difv * 2. / dt1) !vertical velocity due to mixing

        ra = rand(0)
        wmix = wt * (2. * ra - 1.) !add stochastic component

!---vertical advection velocity (w_cur)--------------------------------------
        call vertical_val(floc,io,jo,im,jm,kb,zz,w,h,zo(k))
        w_cur = floc

! update position in the vertical
        zn(k) = zo(k) + (wb + wmix + w_cur) * dt1 

!--check if remains at bottom after erosion
        if(ref(k) == -1.)then
            if(zn(k) <= h(io, jo)*zz(kb-1)) then
                zn(k) = h(io, jo)*zz(kb-1)
                xn(k) = xo(k)
                yn(k) = yo(k)
                goto 240
            else
                ref(k) = 0
            endif
        endif

        if(zn(k) > 0) then    !reset when above surface
            zn(k) = zo(k)
            zp1 = zo(k)
            zp2 = 0.
            zp = 0.5 * (zp1 + zp2)
        else
            zp1 = zo(k)
            zp2 = zn(k)
            zp = 0.5 * (zp1 + zp2)
        end if

        if(abs(zn(k)) > h(io, jo)) zn(k) = -h(io, jo)*(-zz(kb-1)) !bottom
        xp = xo(k)
        yp = yo(k)
                     
        fdep(k)=zn(k)

!  Horizontal movement>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!---From xp, yp, zp find local horizontal current velocity (uloc,vloc) at SIs location-------
        call  INTERP_SIGMA(uloc, vloc, io, jo, im, jm, kb, &
        dx(io,jo), dy(io,jo), zz, &
        u, v, h, xp, yp, zp)

        if (waves) then                              !stokes 1
!  stokes 1
            akz=min(-2.*wk*zo(k),50.)
            uwave = stx0*exp(-akz)
            vwave = sty0*exp(-akz)
        else
!---run without waves
            uwave = 0.
            vwave = 0.
        endif


        if (wind .AND. AirSurfRatio > 0.) then       !wind drag
!---leaway 1/windage
            uwind = kwind*sqrt(AirSurfRatio)*windx
            vwind = kwind*sqrt(AirSurfRatio)*windy
        else
            uwind = 0.
            vwind = 0.
        endif

        call vertical_val(floc,io,jo,im,jm,kb,zz,aam,h,zp) !horizontal mixing / use hydrodynamic model coefficient (AAM)
        difh = floc
        uu = sqrt(difh * 3. / dt1)
        rx = rand(0)
        ry = rand(0) 
        umix=uu * (2. * rx - 1.)                    !add stochastic component-x
        vmix=uu * (2. * ry - 1.)                    !add stochastic component-y

!-------update position in the horizontal
        xn(k) = xo(k) + (uloc + uwave + uwind + umix) * dt1 
        yn(k) = yo(k) + (vloc + vwave + vwind + vmix) * dt1 
         
!---calculate distance covered
        dxx=(uloc + uwave + uwind) * dt1 
        dyy=(vloc + vwave + vwind) * dt1 
        if(dxx+dyy < 3*dt1)then
            dista(k)=dista(k)+sqrt(dxx**2.+dyy**2.)*1.e-3 !distance covered (km)
        endif

!---update (i,j) based on new position(xn,yn)
        call FCoord (in, jn, alon(1,1)-(stepx/2.), &
        alat(1,1)-(stepy/2.), stepx, stepy, &
        xn(k), yn(k), pi, dg)

!---check if out of the boundaries
        270 continue
        if(in < 1 .OR. in > im)then
            ref(k) = -2
            variant(k) =  0
            goto 200
        endif
        if(jn < 1 .OR. jn > jm)then
            ref(k) = -2
            variant(k) =  0
            go to 200
        endif
        if(xn(k) < 0. .OR. xn(k) > xm) then
            ref(k) = -2
            variant(k) =  0
            go to 200
        end if
        if(yn(k) < 0. .OR. yn(k) > ym) then
            ref(k) = -2
            variant(k) =  0
            go to 200
        end if

! Beaching>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!---first time on beach
        if(h(in, jn) <= 1. .AND. ref(k) == 0) then
         rrtime=rtime(ivariant,nsp)                       !residence time on beach
            call beaching(rbt(k),ref(k),in,jn,rrtime,lbc,im,jm)
            rbt(k)=rrtime
            ref(k)=77. 
        endif

        if(ref(k) > 0) then
            rbt(k) = rbt(k) - dt1 / 3600. !update time on beach

            if(rbt(k) > 0.)then   !remain on the beach until rbt=0
                call Cart2Geo (xsi(k), ysi(k), alon(1,1)-(stepx/2.), &
                alat(1,1)-(stepy/2.), xo(k), yo(k), pi, dg)
                goto 200
            endif

            icount=0
            50 continue

!--send back to the sea (random on waves direction) when rbt=0
            r1 = rand(0)
            r2 = rand(0)

            if(stx0**2.+sty0**2. > 0.)then
                wdirx=stx0/(sqrt(stx0**2.+sty0**2.))    !waves direction
                wdiry=sty0/(sqrt(stx0**2.+sty0**2.))
            else
                wdirx=0.
                wdiry=0.
            endif
            xn(k) = xo(k) + r1 * rin *wdirx
            yn(k) = yo(k) + r1 * rin *wdiry
!        xn(k) = xo(k) + r1 * rin * cos(2. * pi * r2)
!        yn(k) = yo(k) + r1 * rin * sin(2. * pi * r2)

            call FCoord (in,jn,alon(1,1)-(stepx/2.),alat(1,1)-(stepy/2.), &  !get new (i,j)
                         stepx,stepy,xn(k),yn(k),pi,dg)
            icount=icount+1

            if(h(in, jn) <= 1. .AND. icount < 20) goto 50            !check if at sea or still at coast
!--if send back with waves does not succeeds goback
!        if(h(in, jn) <= 1.)then
!            call goback(rbt(k), ref(k),xn(k),yn(k),xo(k),yo(k),dt1,io,jo,im,jm,kb,dx(io,jo),dy(io,jo),zz,u,v,h)
!        endif

            if(h(in, jn) > 1.)then                                   !back at sea
                ref(k) =0.
                rbt(k) =0.
            else
                call Cart2Geo (xsi(k), ysi(k), alon(1,1)-(stepx/2.), &
                alat(1,1)-(stepy/2.), xo(k), yo(k), pi, dg)
                goto 200
            endif
        endif
        250 continue

!--check deposition
        if(zn(k) <= h(in, jn)*zz(kb-1)) then
            zn(k) = h(in, jn)*zz(kb-1)
! check bottom stress for  deposition
            bstress=SQRT(wubot(in,jn)**2+wvbot(in,jn)**2)*1023.
            if(bstress < depstress)then
                ref(k) = -1
            endif
        endif

        240 continue

!---reinit position sequence
        xo(k) = xn(k)
        yo(k) = yn(k)
        zo(k) = zn(k)

! update lon,lat (xsi,ysi)
        call Cart2Geo (xsi(k), ysi(k), alon(1,1)-(stepx/2.), &
        alat(1,1)-(stepy/2.), xn(k), yn(k), pi, dg)

        fdep(k)=zn(k)
    200 enddo
             
    999 return
end subroutine lagr_move

