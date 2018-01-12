program SiO2_cluster
!*****************************************************************************************!
!Fortran module containing the modules required to do an ab initio calculation
!including the water and the silica framework with a selected silanol.
!
!Author : W.H. Thompson, J.A. Harvey, P.C. Burris, P.N. Wimalasiri
!Copyright Thompson Group March 2017
!
!*****************************************************************************************
  use kinds
  use SiO2_cluster_mod
  implicit none

! scalars
  integer(kind=ip)              :: i,j,k,n,p,q,r,nfrozen_o,n_pore_si,nsil,ngem,iatm,nwat
  integer(kind=ip)              :: io,ih,ihw,iow,joh,out_cnt,in_cnt,iatm_count
  integer(kind=ip)              :: out_cnt_surf,in_cnt_surf,imin_surf,tmp
  integer(kind=ip)              :: out_cnt_hyd,in_cnt_hyd
  integer(kind=ip)              :: nhydroxyl,nframes,ioh,natms,cluster_size 
  integer(kind=ip)              :: ngrab,nstart,nskip,ngrid,fit_type,imin,it
  integer(kind=ip)              :: cluster_cnt,cluster_manual,natms_traj,explicit_cnt
  real(kind=dp)                 :: dt,R_inner,R_outer,rmin,delr
  real(kind=dp)                 :: box,boy,boz,comx,comy,comz,dz
  real(kind=dp)                 :: massO,massH,qO_water,qH_water,qO_frozen
  real(kind=dp)                 :: qSi_frozen,mtot,rSPCE,pi,norm,norm_H1,norm_H2,Eproj_O,Eproj_H
  real(kind=dp)                 :: dist_O,dist_1,dist_2,dist_HO,dist_H1,dist_H2,dist_HSi
  real(kind=dp)                 :: dist_H_O,rOH,dmin,dmin_surf,qO_sil,qH_sil,charge,charge_cnt
  real(kind=dp),parameter       :: angperau=0.52917721092
  character*100                 :: tmp_path,main_path,data_path,scr_path

  real(kind=dp)                 :: hbnd_OO,hbnd_H1O,hbnd_H2O
  integer(kind=ip)              :: choose
  
! arrays
  real(kind=dp),allocatable,dimension(:,:)      :: frozen_o,si,hydroxyl_o,hydroxyl_h
  real(kind=dp),allocatable,dimension(:)        :: xxx_cluster,yyy_cluster,zzz_cluster
  real(kind=dp),allocatable,dimension(:,:)      :: d_conn,surf_Si,surf_O1,surf_O2
  real(kind=dp),allocatable,dimension(:,:)      :: surf_Sis,surf_O1s,surf_O2s  
  real(kind=dp),allocatable,dimension(:,:)      :: hydroxyl_os,hydroxyl_hs
  real(kind=dp),allocatable,dimension(:,:)      :: rO,rH1,rH2,rOs,rH1s,rH2s 
  real(kind=dp),allocatable,dimension(:,:)      :: cluster,cluster_s
  real(kind=dp),dimension(3)                    :: L,eOH,rCM,clus_H,clus_O,clus_Si
  real(kind=dp),dimension(3)                    :: eOH_H1,eOH_H2,rCM_H1,rCM_H2
  real(kind=dp),dimension(3)                    :: clus_froz_O,shift1,E_O,E_H,rtmpO,rtmpH
  real(kind=dp),dimension(3)                    :: shift2,shift3
  real(kind=dp),allocatable,dimension(:)        :: drank,drank_surf,drank_hyd
  integer(kind=ip),dimension(:),allocatable     :: atm_cluster,atm_cluster_plc
  integer(kind=dp),allocatable,dimension(:)     :: irank
  logical,dimension(:),allocatable              :: si_cluster,o_cluster,oh_cluster,inner,inner_hyd
  logical,dimension(:),allocatable              :: outer,inner_surface,outer_surface,outer_hyd
  character*1   :: iat
  character*2   :: tmp_charge   
  character*3   :: ext(900)      
! testing
  character*100   :: test1,test2,test3,test4

! set path                                                                                  
  main_path='/home/p470w698/work/lammps/NWCHEM/SFG/calc_run/run_cluster_charge_corr/sil_11/charge_0_mult_2/'
  data_path='/home/p470w698/work/lammps/NWCHEM/SFG/calc_run/run_cluster_charge_corr/sil_11/charge_0_mult_2/'
  scr_path ='/home/p470w698/scratch/PANFS-SCRATCH-PUBUDU/NWCHEM/temp6/tmp4/'

! set up script file
  open(23,file='sub_script')
  write(23,'(A)') '#!/bin/sh'

! read in simulation parameters

  open(13,file='rand_grab_num.in',status='old')

  read(13,*)
  read(13,*) nwat, nfrozen_o, n_pore_si, nhydroxyl, nframes
  read(13,*)
  read(13,*) natms, nsil, ngem, ioh, cluster_size, cluster_manual
  read(13,*)
  read(13,*) box, boy, boz
  read(13,*)
  read(13,*) ngrab, nstart, nskip, dt, R_inner, R_outer
  read(13,*)
  read(13,*) ngrid, rmin, delr, fit_type
  read(13,*)
  read(13,*) massO, massH
  read(13,*)
  read(13,*) qO_water, qH_water
  read(13,*)
  read(13,*) qO_frozen, qSi_frozen
  read(13,*)
  read(13,*) qO_sil, qH_sil

  close(13)

! prepare for reading files

  do j=1,900
        write(ext(j),'(I0.3)') j
  enddo

! allocate and get the frozen Si and O

  allocate(si(n_pore_si,3))
  allocate(frozen_o(nfrozen_o,3))
  call get_frozen_pore(nfrozen_o,n_pore_si,frozen_o,si)

! write(6,*) 'read in the frozen stuff'
! open trajectory
  open(15,file=trim(data_path)//'traj/traj.xyz',status='old')

! obtain box lengths
  L(1) = box
  L(2) = boy
  L(3) = boz

! calculate quantities
  mtot = massO + massH
  rSPCE = 1.0_dp
  pi = 4_dp*datan(1_dp)

  if (nstart>1) then 
     do j = 1,nstart-1
        do k = 1, 5694
           read(15,*)
        enddo
     enddo
  endif

! read(15,*) test1
! write(6,'(A)') test1
! stop
 
  it = 0
  do i=1,ngrab
! allocate the rest of the arrays
        it = it+1
        allocate(hydroxyl_o(nhydroxyl,3))
        allocate(hydroxyl_h(nhydroxyl,3))
        allocate(hydroxyl_os(nhydroxyl,3))
        allocate(hydroxyl_hs(nhydroxyl,3))
        allocate(rO(nwat,3))
        allocate(rH1(nwat,3))
        allocate(rH2(nwat,3))
        allocate(rOs(nwat,3))
        allocate(rH1s(nwat,3))
        allocate(rH2s(nwat,3))
        allocate(xxx_cluster(cluster_size*3+300))
        allocate(yyy_cluster(cluster_size*3+300))
        allocate(zzz_cluster(cluster_size*3+300))
        allocate(atm_cluster(cluster_size*3+300))
        allocate(atm_cluster_plc(cluster_size*3+300))
        allocate(si_cluster(n_pore_si))
        allocate(o_cluster(nfrozen_o+nhydroxyl))
        allocate(d_conn(n_pore_si,(nfrozen_o+nhydroxyl)))
       
! advance trajectory to nstep start
        io = 0
        ih = 0
        ihw = 0
        iow = 0
  
! read(15,*) test1
! print *, test1

! read in silanol configurations
        read(15,*) !natms_traj
        read(15,*)
        do p = 1, nsil
              read(15,*) iat, (hydroxyl_o(p,k),k=1,3)
              read(15,*) iat, (hydroxyl_h(p,k),k=1,3)
              hydroxyl_h(p,:) = hydroxyl_h(p,:) - L(:)* &
                              dnint((hydroxyl_h(p,:)-hydroxyl_o(p,:))/L(:))
        enddo
  
! read in geminal configurations
        do q = nsil+1, (2*ngem+nsil)
              read(15,*) iat, (hydroxyl_o(q,k),k=1,3)
              read(15,*) iat, (hydroxyl_h(q,k),k=1,3)
              hydroxyl_h(p,:) = hydroxyl_h(p,:) - L(:)* &
                              dnint((hydroxyl_h(p,:)-hydroxyl_o(p,:))/L(:))
        enddo
 
   
! read(15,*) test1,  test2
! print *, test1,  test2
! stop 
 
! read in water configurations
        do r = 1, nwat
              read(15,*) iat, (rO(r,k),k=1,3)
              read(15,*) iat, (rH1(r,k),k=1,3)
              read(15,*) iat, (rH2(r,k),k=1,3)
        
              rH1(r,:) = rH1(r,:) - L(:)*dnint((rH1(r,:)-rO(r,:))/L(:))
              rH2(r,:) = rH2(r,:) - L(:)*dnint((rH2(r,:)-rO(r,:))/L(:))
        enddo

! determine the molecule to grab from a random number

! ngemsi = nsil + 2*ngem
! ngemsi_grab = int(dfloat(ngemsi-1)*ran3(idum))+1
! ig = ngemsi_grab

! calculate connectivity between Si and O in a matrix

        call si_o_connectivity(nfrozen_o,nhydroxyl,n_pore_si,si,                 &
             hydroxyl_o,frozen_o,L,d_conn)

        iatm = 0
        
        do j=1,n_pore_si
                si_cluster(j) = .false.
        enddo

        do j=1,(nfrozen_o+nhydroxyl)
            o_cluster(j) = .false.
        enddo

!       do j=1,nhydroxyl
!           oh_cluster(j) = .false.
!       enddo
        
        call find_central_cluster(n_pore_si,nhydroxyl,nfrozen_o,ioh,             &
             L,si,hydroxyl_o, hydroxyl_h, frozen_o, d_conn, iatm,                &
             xxx_cluster,yyy_cluster, zzz_cluster,cluster_size,                  &
             atm_cluster,atm_cluster_plc,o_cluster,si_cluster)

!       write(6,*) 'made it through find_central'

!       open(unit = 2 , file = 'cluster.xyz')
!       open(unit = 3 , file = 'cluster.in')
!       write(2,"(i10)") iatm
!       write(2,*) "Si Cluster"
        
!       do k=1,iatm
!
!           write(2,"(i10,3f20.10)") atm_cluster(k), xxx_cluster(k), yyy_cluster(k), &
!           zzz_cluster(k)
!           write(3,"(i10,3f20.10,i10)") atm_cluster(i), xxx_cluster(i), yyy_cluster(i), &
!           zzz_cluster(i), atm_cluster_plc(i)
!
!       enddo
!       close(2)
        
!       stop

        call calc_com_cc(cluster_size,xxx_cluster,yyy_cluster,                   &
             zzz_cluster,atm_cluster,box,boy,iatm,comx,comy,comz)

!       write(6,*) 'made it through calc_com'

        call find_cluster(n_pore_si,nhydroxyl,nfrozen_o,ioh,         &
             box,boy,si,hydroxyl_o,hydroxyl_h,frozen_o,d_conn,iatm,  &
             xxx_cluster,yyy_cluster,zzz_cluster,cluster_size,       &
             atm_cluster,atm_cluster_plc,o_cluster,si_cluster,comx,comy,comz)

!       write(6,*) 'made it through find_cluster'

! call add_water(nwat,boz,rO,rH1,rH2,iatm,xxx_cluster,          &
!      yyy_cluster,zzz_cluster,cluster_size,atm_cluster,atm_cluster_plc)

! write(6,*) 'made it through add_water'

! wrap all cluster into the central Si box

!       write(*,*) iatm
!       do i=1,iatm

!               dz = zzz_cluster(i) - zzz_cluster(3)
!               zzz_cluster(i) = zzz_cluster(i) - (boz*anint(dz/boz))

!       enddo

! write cluster

!      open(unit = 2 , file = 'cluster.xyz')
!      open(unit = 3 , file = 'cluster.in')
!      write(2,"(i10)") iatm
!      write(2,*) "Si Cluster"
!      do j=1,iatm
!
!           write(2,"(i10,3f20.10)") atm_cluster(j), xxx_cluster(j), yyy_cluster(j), &
!           zzz_cluster(j)
!           write(3,"(i10,3f20.10,i10)") atm_cluster(i), xxx_cluster(i), yyy_cluster(i), &
!           zzz_cluster(i), atm_cluster_plc(i)
!
!      enddo
!      close(2)
!      close(unit=3)
  
! allocate the arrays to differentiate between inner or outer molecules
        allocate(inner(nwat))
        allocate(outer(nwat))
        allocate(irank(nwat))
        allocate(drank(nwat))

! allocate the arrays to differentiate between inner or outer surface atoms

!       tmp = (iatm-4)/3
!       allocate(surf_Si(tmp,3))
!       allocate(surf_O1(tmp,3))
!       allocate(surf_O2(tmp,3))
!       allocate(surf_Sis(tmp,3))
!       allocate(surf_O1s(tmp,3))
!       allocate(surf_O2s(tmp,3))

        allocate(cluster(iatm,3))
        allocate(cluster_s(2:iatm,3))

! move cluster from scalar to array
 
!       clus_H(1)  = xxx_cluster(1)
!       clus_H(2)  = yyy_cluster(1)
!       clus_H(3)  = zzz_cluster(1)

!       clus_O(1)  = xxx_cluster(2)
!       clus_O(2)  = yyy_cluster(2)
!       clus_O(3)  = zzz_cluster(2)         

!       clus_Si(1) = xxx_cluster(3)
!       clus_Si(2) = yyy_cluster(3)
!       clus_Si(3) = zzz_cluster(3)

!       clus_froz_O(1) = xxx_cluster(4)
!       clus_froz_O(2) = yyy_cluster(4)
!       clus_froz_O(3) = zzz_cluster(4)

! move surface cluster from scalar to array

!       iatm_count = 0
!       do j=1,(iatm-4)/3
!               surf_Si(j,1) = xxx_cluster(j*3+2)
!               surf_Si(j,2) = yyy_cluster(j*3+2)
!               surf_Si(j,3) = zzz_cluster(j*3+2)

!               surf_O1(j,1) = xxx_cluster(j*3+3)
!               surf_O1(j,2) = yyy_cluster(j*3+3)
!               surf_O1(j,3) = zzz_cluster(j*3+3)
                
!               surf_O2(j,1) = xxx_cluster(j*3+4)
!               surf_O2(j,2) = yyy_cluster(j*3+4)
!               surf_O2(j,3) = zzz_cluster(j*3+4)
!               iatm_count = iatm_count+1
!       enddo

        do j=1,iatm
            cluster(j,1) = xxx_cluster(j)
            cluster(j,2) = yyy_cluster(j)
            cluster(j,3) = zzz_cluster(j)
        enddo

!      open(unit = 2 , file = 'cluster.xyz')
!      write(2,"(i10)") iatm
!      write(2,*) "Si Cluster"
!      do j=1,iatm
!
!           write(2,"(i10,3f20.10)") atm_cluster(j),(cluster(j,k),k=1,3)
!
!      enddo
!      close(2)
        
!       write(6,*) ' iatm_count = ', iatm_count
!       write(6,*) ' no of SiO2 clusters = ', (iatm-4)/3
!       stop 

        allocate(inner_surface(10:iatm))
        allocate(outer_surface(10:iatm))
        allocate(inner_hyd(nhydroxyl))
        allocate(outer_hyd(nhydroxyl))
        allocate(drank_surf(10:iatm))
        allocate(drank_hyd(nhydroxyl))
! calculate the OH vector

        norm = 0_dp
        do k=1,3
              eOH(k) = cluster(1,k) - cluster(2,k)
              norm = norm + eOH(k)*eOH(k)
              rCM(k) = (massO*cluster(2,k) + massH*cluster(1,k))/mtot
        enddo

        norm = dsqrt(norm)

        do k = 1, 3
              eOH(k) = eOH(k)/norm
        enddo

! zero the fields

        do k=1,3
              E_O(k) = 0_dp
              E_H(k) = 0_dp
        enddo

        Eproj_O = 0_dp
        Eproj_H = 0_dp

        do j=1,nwat
              inner(j) = .false.
              outer(j) = .false.
              irank(j) = 0
        enddo

        do j=10,iatm
                inner_surface(j) = .false.
                outer_surface(j) = .false.

!               write(6,'(A)') inner_surface(j)
        enddo

!       do j=1,nhydroxyl
!               inner_hyd(j) = .false.
!               outer_hyd(j) = .false.
!       enddo

! center all the positions based on the HO- hydrogen atom

        do j=1,nwat
          do k = 1, 3
              shift1(k)  = L(k)*dnint((cluster(1,k)-rO(j,k))/L(k))
              rOs(j,k)  = rO(j,k)  + shift1(k)
              rH1s(j,k) = rH1(j,k) + shift1(k)
              rH2s(j,k) = rH2(j,k) + shift1(k)
          enddo   
        enddo

        do j=2,iatm
          do k=1,3
              shift2(k)  = L(k)*dnint((cluster(1,k)-cluster(j,k))/L(k))
              cluster_s(j,k) = cluster(j,k) + shift2(k)
          enddo
        enddo    

!       open(unit = 2 , file = 'cluster.xyz')
!       write(2,"(i10)") iatm
!       write(2,*) "Si Cluster"
       
!       write(2,"(A,3f20.10)") 'H',(cluster(1,k),k=1,3)


!       do j=2,iatm
!           if (atm_cluster(j).eq.1) then
!                tmp_charge='Si'

!                else if ((atm_cluster(j).eq.2) .or. (atm_cluster(j).eq.4)) then
!                    tmp_charge='O'

!                    else if (atm_cluster(j).eq.3) then
!                        tmp_charge='H'
!            endif
!            write(2,"(A,3f20.10)") tmp_charge,(cluster_s(j,k),k=1,3)
  
!       enddo
!       close(2)

!       stop
       
!       do j=1,nhydroxyl
!         do k=1,3
!             shift3(k)  = L(k)*dnint((clus_H(k)-hydroxyl_o(j,k))/L(k))
!             hydroxyl_os(j,k) = hydroxyl_o(j,k) + shift3(k)
!             hydroxyl_hs(j,k) = hydroxyl_h(j,k) + shift3(k)
!         enddo
!       enddo    
  
!      stop

! calculate the distance from the OH to other water O's
        out_cnt=0
        do j=1,nwat
              outer(j) = .false.
              dist_O   = 0_dp
              dist_1   = 0_dp
              dist_2   = 0_dp
              dist_HO  = 0_dp
              dist_H1  = 0_dp
              dist_H2  = 0_dp

              do k=1,3
                      dist_O  = dist_O  + (rOs(j,k)  - rCM(k))**2
                      dist_1  = dist_1  + (rH1s(j,k) - rCM(k))**2
                      dist_2  = dist_2  + (rH2s(j,k) - rCM(k))**2

                      dist_HO = dist_HO + (rOs(j,k)  - cluster(1,k))**2
                      dist_H1 = dist_H1 + (rH1s(j,k) - cluster(1,k))**2
                      dist_H2 = dist_H2 + (rH2s(j,k) - cluster(1,k))**2
              enddo

              dist_O  = sqrt(dist_O)
              dist_1  = sqrt(dist_1)
              dist_2  = sqrt(dist_2)

              dist_HO = sqrt(dist_HO)
!             write(6,*) ' dist_HO= ', dist_HO
              dist_H1 = sqrt(dist_H1)
              dist_H2 = sqrt(dist_H2)
              
              drank(j) = dist_HO

              if(dist_HO.le.R_outer) then
                      outer(j) = .true.
                      out_cnt  = out_cnt+1

                      do k=1,3
!                             E_O(k) = E_O(k)  &
!                                     + qO*(clus_O(k) - rOs(j,k)) /dist_OO**3   &
!                                     + qH*(clus_O(k) - rH1s(j,k))/dist_O1**3  &
!                                     + qH*(clus_O(k) - rH2s(j,k))/dist_O2**3
                              E_H(k) = E_H(k)  &
                                      + qO_water*(cluster(1,k) - rOs(j,k)) /dist_HO**3  &
                                      + qH_water*(cluster(1,k) - rH1s(j,k))/dist_H1**3  &
                                      + qH_water*(cluster(1,k) - rH2s(j,k))/dist_H2**3

                      enddo
              endif
        enddo
        write(6,*) ' out_cnt water = ', out_cnt
 
!       stop

        out_cnt_surf = 0
        do j=10,iatm
              outer_surface(j) = .false.
              dist_HSi = 0_dp

              do k=1,3
                dist_HSi = dist_HSi + (cluster_s(j,k) - cluster(1,k))**2
              enddo

              dist_HSi = sqrt(dist_HSi)
              drank_surf(j) = dist_HSi

              if (dist_HSi.le.R_outer) then
                outer_surface(j) = .true.
                out_cnt_surf = out_cnt_surf+1
              endif
        enddo
        write(6,*) ' out_cnt SiO2  = ', out_cnt_surf

!       out_cnt_hyd = 0
!       do j=1,nhydroxyl
!             outer_hyd(j) = .false.
!             dist_H_O = 0_dp

!             do k=1,3
!               dist_H_O = dist_H_O + (hydroxyl_os(j,k) - clus_H(k))**2
!             enddo

!             dist_H_O = sqrt(dist_H_O)
!             drank_hyd(j) = dist_H_O

!             if (dist_H_O.le.R_outer) then
!               outer_hyd(j) = .true.
!               out_cnt_hyd = out_cnt_hyd+1
!             endif
!       enddo
!       write(6,*) ' out_cnt hydroxyl  = ', out_cnt_hyd

!       stop

        do k = 1, 3
              Eproj_O = Eproj_O + E_O(k)*eOH(k)*angperau**2
              Eproj_H = Eproj_H + E_H(k)*eOH(k)*angperau**2
        enddo

! rank the ones that are in the shortest N_inner distances 
! only need to look at those with outer(j) = .true.

        open(20,file='envir_'//ext(i)//'.dat')
        write(20,'(3(A,F8.5))') ' Eproj_O = ',Eproj_O,  &
                                ' Eproj_H = ',Eproj_H
        write(20,'(A,I4)') ' # Outer = ',out_cnt

! first, find the nearest water molecule

        dmin = R_outer
        imin = 0_ip
        imin_surf = 0_ip
        do j=1,nwat
              if(outer(j)) then
                      if(drank(j).lt.dmin) then
                        imin = j
                        dmin = drank(j)
                      endif
              endif
        enddo

        dmin = R_outer
        do j=10,iatm
              if(outer_surface(j)) then
                      if(drank_surf(j).lt.dmin) then
                        imin_surf = j
                        dmin_surf = drank_surf(j)
                      endif
              endif
        enddo

!       write(6,*) ' imin = ',imin
!       write(6,*) ' imin_surf = ',imin_surf
        irank(1) = imin

        norm_H1 = 0_dp
        do k=1,3
           eOH_H1(k) = rH1s(imin,k) - rOs(imin,k)
           norm_H1 = norm_H1 + eOH_H1(k)*eOH_H1(k)
           rCM_H1(k) = (massO*rOs(imin,k) + massH*rH1s(imin,k))/mtot
        enddo

        norm_H1 = dsqrt(norm_H1)

        norm_H2 = 0_dp
        do k=1,3
           eOH_H2(k) = rH1s(imin,k) - rOs(imin,k)
           norm_H2 = norm_H2 + eOH_H2(k)*eOH_H2(k)
           rCM_H2(k) = (massO*rOs(imin,k) + massH*rH1s(imin,k))/mtot
        enddo

        norm_H2 = dsqrt(norm_H2)

! h-bond criteria (check if the water molecule which is very close to silanol is h-bonding)

        hbnd_OO = 0_dp
        hbnd_H1O = 0_dp
        hbnd_H2O = 0_dp        
        do k=1,3
           hbnd_OO  = hbnd_OO  + (rOs(imin,k)-cluster_s(2,k))**2
        enddo

        hbnd_OO = sqrt(hbnd_OO)
        
        if (hbnd_OO.le.3.5) then
           do k=1,3
              hbnd_H1O  = hbnd_H1O  + (rH1s(imin,k)-cluster_s(2,k))**2
              hbnd_H2O  = hbnd_H2O  + (rH2s(imin,k)-cluster_s(2,k))**2              
           enddo

           hbnd_H1O = sqrt(hbnd_H1O)
           hbnd_H2O = sqrt(hbnd_H2O)

           if (hbnd_H1O.le.2.5) then
              choose = 1
           else if (hbnd_H2O.le.2.5) then
              choose = 2
           endif

        else
           write(6,'(A)') 'The closest water molecule does not h-bond to silanol'
           stop
        endif

! calculate SiO-OH vector and SiO-HO vectors 

        e_OO = 0_dp
        e_HO = 0_dp
        dot  = 0_dp
        
        do k=1,3
           e_OO = e_OO + cluster_s(2,k) - rOs(imin,k)

! finally, calculate the angle between SiO-OH and SiO_HO to confirm h-bonding
           if (choose.eq.1) then
              e_OH = e_HO + cluster_s(2,k) - rH1s(imin,k)
              dot = dot + (cluster_s(2,k)*rH1(imin,k))
           else if (choose.eq.2) then
              e_OH = e_HO + cluster_s(2,k) - rH2s(imin,k)
              dot = dot + (cluster_s(2,k)*rH1(imin,k))              
           endif
        enddo 

        if (choose.eq.1) then
           cos_theta = dot/(hbnd_OO*hbnd_H1O)
        else if (choose.eq.2) then
           cos_theta = dot/(hbnd_OO*hbnd_H2O)
        endif

        write(6,'(A,F12.5)') 'cos_theta =',cos_theta
        stop
        
        inner(imin) = .true.
        outer(imin) = .false.
        in_cnt=1_ip
        out_cnt=out_cnt-1_ip        
!       inner_surface(imin_surf) = .true.
!       outer_surface(imin_surf) = .false.

        if (dmin.le.4.0) then
              write(6,*) 'There is at least one water molecule in the inner sphere.'
              else
                      write(6,*) 'There is no water molecule in the inner sphere.'
        endif
                      
        in_cnt_surf=0_ip
!       in_cnt_hyd=0_ip
        dmin = R_inner
        do j=1,nwat
              if (outer(j)) then
                      if (drank(j).lt.dmin) then
                        outer(j) = .false.
                        out_cnt = out_cnt-1_ip
                        inner(j) = .true.
                        in_cnt = in_cnt+1_ip
                      endif
              endif
        enddo

!       do j=10,iatm
!             if (outer_surface(j)) then
!                     if (drank_surf(j).lt.dmin) then
!                       outer_surface(j) = .false.
!                       out_cnt_surf=out_cnt_surf-1
!                       inner_surface(j) = .true.
!                       in_cnt_surf = in_cnt_surf+1
!                     endif
!             endif
!       enddo

!       do j=1,nhydroxyl
!             if (outer_hyd(j)) then
!                     if (drank_hyd(j).lt.dmin) then
!                       outer_hyd(j) = .false.
!                       out_cnt_hyd=out_cnt_hyd-1
!                       inner_hyd(j) = .true.
!                       in_cnt_hyd = in_cnt_hyd+1
!                     endif
!             endif
!       enddo

        write(6,*) ' out_cnt water = ', out_cnt
        write(6,*) ' in_cnt water  = ', in_cnt

        write(6,*) ' out_cnt SiO2 = ', out_cnt_surf
        write(6,*) ' in_cnt SiO2  = ', in_cnt_surf

!       write(6,*) ' out_cnt hydroxyl = ', out_cnt_hyd
!       write(6,*) ' in_cnt hydroxyl  = ', in_cnt_hyd
       
!       stop

! write out configs with wrapped boundary conditions

! check the total charge
          charge_cnt=0_dp
          do j = 10,iatm
                if (atm_cluster(j).eq.1) then
                    charge=qSi_frozen
                    tmp_charge='Si'
                    
                    else if (atm_cluster(j).eq.2) then
                        charge=qO_sil
                        tmp_charge='O'

                        else if (atm_cluster(j).eq.3) then
                            charge=qH_sil
                            tmp_charge='H'

                            else if (atm_cluster(j).eq.4) then
                                charge=qO_frozen
                                tmp_charge='O'
                endif

                if (outer_surface(j)) then
                    charge_cnt=charge_cnt+charge
                endif
          enddo

          write(6,'(A)') ' '
          write(6,'(A,f)') ' total charge= ',charge_cnt

          open(18,file='explicit_'//ext(i)//'.xyz',status='new')
          open(19,file='scan_'//ext(i)//'.xyz',status='new')
          open(21,file='scan_'//ext(i)//'.nw',status='new')
          open(22,file='scan_'//ext(i)//'.sh',status='new')
          
!         write(6,*) 'opened the files'

      explicit_cnt=1_ip
      do n = 1, ngrid
          rOH = 0_dp
          do k=1,3
                  rtmpH(k) = rCM(k) + massO*eOH(k)*(rmin+delr*dfloat(n-1))/mtot
                  rtmpO(k) = rCM(k) - massH*eOH(k)*(rmin+delr*dfloat(n-1))/mtot
                  rOH = rOH + (rtmpH(k) - rtmpO(k))**2
          enddo
          rOH = dsqrt(rOH)

! change the number of atoms to be written in the xyz file every time if you are 
! swiching from cluster choosing to distance criteria


          write(19,'(I4)') 9+out_cnt_surf+in_cnt_surf+3*(in_cnt+out_cnt)
!         write(19,'(I4)') 4+(in_cnt_surf+out_cnt_surf)*3+(in_cnt+out_cnt)*3+(in_cnt_hyd+out_cnt_hyd-1)*2
          write(19,'(A,F8.5,A)') ' rOH = ',rOH,' Angs.'

          write(21,'(A,2(F8.5,A),F7.2,A,F6.3,A)') 'title "SiOH-/H2O E_O = ',  &
          Eproj_O,' au, E_H = ',Eproj_H, ' au, t = ',dt*dfloat(it-1)/1000d0,  &
          ' ps, r_OH = ',rOH,' Angs"'

          write(21,'(A)') 'echo'
          write(21,'(A,F4.2)') 'charge ',charge_cnt
          write(21,'(A)') 'geometry noautoz units angstroms'

          if (explicit_cnt.eq.1) then
             write(18,'(A)') ' '
             write(18,'(A)') ' '
             write(18,'(A,3F12.5)') 'H',(rtmpH(k),k=1,3)
             write(18,'(A,3F12.5)') 'O',(rtmpO(k),k=1,3)
             write(18,'(A,3F12.5)') 'Si',(cluster_s(3,k),k=1,3)
          endif
 
          write(19,'(A,3F12.5)') 'H',(rtmpH(k),k=1,3)
          write(19,'(A,3F12.5)') 'O',(rtmpO(k),k=1,3)
          write(19,'(A,3F12.5)') 'Si',(cluster_s(3,k),k=1,3)
          write(21,'(A,3F12.5)') 'H',(rtmpH(k),k=1,3)
          write(21,'(A,3F12.5)') 'O',(rtmpO(k),k=1,3)
          write(21,'(A,3F12.5)') 'Si',(cluster_s(3,k),k=1,3)

          do j=4,9
              if (atm_cluster(j).eq.1) then
                  tmp_charge='Si'

                  else if ((atm_cluster(j).eq.2) .or. (atm_cluster(j).eq.4)) then
                      tmp_charge='O'

                      else if (atm_cluster(j).eq.3) then
                          tmp_charge='H'
              endif
              
              if (explicit_cnt.eq.1) then
                 write(18,'(A,3F12.5)') trim(tmp_charge),(cluster_s(j,k),k=1,3)
              endif

              write(19,'(A,3F12.5)') trim(tmp_charge),(cluster_s(j,k),k=1,3)
              write(21,'(A,3F12.5)') trim(tmp_charge),(cluster_s(j,k),k=1,3)
          enddo

!         write(19,'(A,3F12.5)') 'Si',(clus_Si(k),k=1,3)
!         write(19,'(A,3F12.5)') 'O',(rtmpO(k),k=1,3)
!         write(19,'(A,3F12.5)') 'H',(rtmpH(k),k=1,3)
!         write(19,'(A,3F12.5)') 'O',(clus_froz_O(k),k=1,3)
!         write(21,'(A,3F12.5)') 'Si',(clus_Si(k),k=1,3)
!         write(21,'(A,3F12.5)') 'O',(rtmpO(k),k=1,3)
!         write(21,'(A,3F12.5)') 'H',(rtmpH(k),k=1,3)
!         write(21,'(A,3F12.5)') 'O',(clus_froz_O(k),k=1,3)
        
          cluster_cnt = 0_ip
!         do j = 1, iatm_count
!               if (inner_surface(j) .and. (cluster_cnt.le.cluster_manual)) then
!                       cluster_cnt = cluster_cnt + 1
!                       write(19,'(A,3F12.5)') 'Si', (surf_Sis(j,k),k=1,3)
!                       write(19,'(A,3F12.5)') 'O', (surf_O1s(j,k),k=1,3)
!                       write(19,'(A,3F12.5)') 'O', (surf_O2s(j,k),k=1,3)
!                       write(21,'(A,3F12.5)') 'Si', (surf_Sis(j,k),k=1,3)
!                       write(21,'(A,3F12.5)') 'O', (surf_O1s(j,k),k=1,3)
!                       write(21,'(A,3F12.5)') 'O', (surf_O2s(j,k),k=1,3)
!               endif
!         enddo

!         do j = 1, nhydroxyl
!               if (inner_hyd(j) .and. (j.ne.ioh)) then
!                       write(19,'(A,3F12.5)') 'O', (hydroxyl_os(j,k),k=1,3)
!                       write(19,'(A,3F12.5)') 'H', (hydroxyl_hs(j,k),k=1,3)
!                       write(21,'(A,3F12.5)') 'O', (hydroxyl_os(j,k),k=1,3)
!                       write(21,'(A,3F12.5)') 'H', (hydroxyl_hs(j,k),k=1,3)
!               endif
!         enddo

          do j = 1, nwat
                  if(inner(j)) then
                     
                     if (explicit_cnt.eq.1) then
                          write(18,'(A,3F12.5)') 'O',(rOs(j,k),k=1,3)
                          write(18,'(A,3F12.5)') 'H',(rH1s(j,k),k=1,3)
                          write(18,'(A,3F12.5)') 'H',(rH2s(j,k),k=1,3)
                     endif

                          write(19,'(A,3F12.5)') 'O',(rOs(j,k),k=1,3)
                          write(19,'(A,3F12.5)') 'H',(rH1s(j,k),k=1,3)
                          write(19,'(A,3F12.5)') 'H',(rH2s(j,k),k=1,3)
                          write(21,'(A,3F12.5)') 'O',(rOs(j,k),k=1,3)
                          write(21,'(A,3F12.5)') 'H',(rH1s(j,k),k=1,3)
                          write(21,'(A,3F12.5)') 'H',(rH2s(j,k),k=1,3)
                  endif
          enddo

!         write(6,*) 'done with inner clusters'
          explicit_cnt = explicit_cnt + 1
          close(18)

          do j = 10,iatm
                if (atm_cluster(j).eq.1) then
                    charge=qSi_frozen
                    tmp_charge='Si'
                    
                    else if (atm_cluster(j).eq.2) then
                        charge=qO_sil
                        tmp_charge='O'

                        else if (atm_cluster(j).eq.3) then
                            charge=qH_sil
                            tmp_charge='H'

                            else if (atm_cluster(j).eq.4) then
                                charge=qO_frozen
                                tmp_charge='O'
                endif

                if (outer_surface(j)) then
                    cluster_cnt = cluster_cnt+1
                    
                    write(21,'(A,3F12.5,A,3F12.5)') 'bq', (cluster_s(j,k),k=1,3), & 
                    ' charge',charge
!                       write(21,'(A,3F12.5,A,3F12.5)') 'bq', (surf_O1s(j,k),k=1,3), &
!                         ' charge',qO_frozen
!                       write(21,'(A,3F12.5,A,3F12.5)') 'bq', (surf_O2s(j,k),k=1,3), &
!                         ' charge',qO_frozen
                    write(19,'(A,3F12.5)') trim(tmp_charge), (cluster_s(j,k),k=1,3)
!                       write(19,'(A,3F12.5)') 'O', (surf_O1s(j,k),k=1,3)
!                       write(19,'(A,3F12.5)') 'O', (surf_O2s(j,k),k=1,3)
                endif
          enddo

          if (cluster_cnt.ne.out_cnt_surf) then 
              write(6,'(A)') ' inconsistency in point charges'
          endif

!         do j = 1, nhydroxyl
!               if (outer_hyd(j)) then
!                       write(21,'(A,3F12.5,A,3F12.5)') 'bq', (hydroxyl_os(j,k),k=1,3), & 
!                         ' charge',qO_sil
!                       write(21,'(A,3F12.5,A,3F12.5)') 'bq', (hydroxyl_hs(j,k),k=1,3), &
!                         ' charge',qH_sil
!                       write(19,'(A,3F12.5)') 'O', (hydroxyl_os(j,k),k=1,3)
!                       write(19,'(A,3F12.5)') 'H', (hydroxyl_hs(j,k),k=1,3)
!               endif
!         enddo

          do j = 1, nwat
                  if(outer(j)) then
                          write(21,'(A,3F12.5,A,F12.5)') 'bq',(rOs(j,k),k=1,3),  &
                          ' charge',qO_water
                          write(21,'(A,3F12.5,A,F12.5)') 'bq',(rH1s(j,k),k=1,3), &
                          ' charge',qH_water
                          write(21,'(A,3F12.5,A,F12.5)') 'bq',(rH2s(j,k),k=1,3), &
                          ' charge',qH_water
                          write(19,'(A,3F12.5)') 'O',(rOs(j,k),k=1,3)
                          write(19,'(A,3F12.5)') 'H',(rH1s(j,k),k=1,3)
                          write(19,'(A,3F12.5)') 'H',(rH2s(j,k),k=1,3)
                  endif
          enddo

          write(21,'(A)') 'end'
                  if(n.eq.1) then
                                
                          write(21,'(A)') 'basis'
                          write(21,'(A)') '* library 3-21G*'
                          write(21,'(A)') 'end'
                          write(21,'(A)') 'dft'
                          write(21,'(A)') ' iterations 5000'
                          write(21,'(A)') ' xc b3lyp'
                          write(21,'(A)') ' decomp'
!                         write(21,'(A)') ' mult 2'
                          write(21,'(A)') 'end'
                          write(21,'(A)') 'property'
                          write(21,'(A)') ' dipole'
!                         write(21,'(A)') ' response 1 7.73178E-2'
                          write(21,'(A)') ' velocity'
                          write(21,'(A)') 'end'
!                         write(21,'(A)') 'esp'
!                         icnt = 0
!                         do j = 1, nmol - 1 
!                                 if(inner(j)) then
!                                 icnt = icnt + 1
!                         write(21,'(A,I3)')' constrain -0.8476',3*icnt + 1
!                         write(21,'(A,I3)')' constrain 0.4238',3*icnt + 2
!                         write(21,'(A,I3)')' constrain 0.4238',3*icnt + 3
!                         endif
!                         enddo
!                         write(21,'(A)') 'end'
                  endif
!                 write(21,'(A)') 'task dft energy'           
                  write(21,'(A)') 'task dft property'
!                 write(21,'(A)') 'task esp'
        enddo
        close(19)
        close(21)

! write out submit script
        write(22,'(A)') '#MSUB -N scan_'//ext(i)
        write(22,'(A)') '#MSUB -q sixhour'
        if(mod(i,25).eq.0) then
              write(22,'(A)') '#MSUB -m ae'
              write(22,'(A)') '#MSUB -M wthompson@ku.edu'
        endif
          write(22,'(A)') '#MSUB -j oe'
          write(22,'(A)') '#MSUB -d ./'
          write(22,'(A)') '#MSUB -l nodes=2:ppn=20:ib,mem=5gb,walltime=6:00:00'
          write(22,*)
          write(22,'(A)') 'module unload pgi/compiler/15'
          write(22,'(A)') 'module unload openmpi/2.0'
          write(22,'(A)') 'module load intel_mpi_intel64/5.1.2.150'
          write(22,'(A)') 'module load nwchem/6.5_intel'
          write(22,'(A)') 'module load scipy'
          write(22,*)
          write(22,'(A)') 'mkdir '//trim(scr_path)//'tmp_mb_'//ext(i)
          write(22,'(A)') 'cp '//trim(main_path)//'scan_'//ext(i)// &
          '.nw '//trim(scr_path)//'tmp_mb_'//ext(i)
                  if(fit_type.eq.1) then
                          write(22,'(A,A)') 'cp '//trim(main_path)//'fit_main',      &
                          ' '//trim(scr_path)//'tmp_mb_'//ext(i)
                          write(22,'(A,A)') 'cp '//trim(main_path)//'fit_guesses.in',&
                          ' '//trim(scr_path)//'tmp_mb_'//ext(i)
                          write(22,'(A,A)') 'cp '//trim(main_path)//'dipole_get.py', &
                          ' '//trim(scr_path)//'tmp_mb_'//ext(i)

                          else if(fit_type.eq.2) then
                                  write(22,'(A,A)') &
                                  'cp '//trim(main_path)//'fit_morse',    &
                                  ' '//trim(scr_path)//'tmp_mb_'//ext(i)
                                   write(22,'(A,A)') &
                                 'cp '//trim(main_path)//'fit_morse.in', &
                                  ' '//trim(scr_path)//'tmp_mb_'//ext(i)
                  endif
          write(22,'(A)') 'cd '//trim(scr_path)//'tmp_mb_'//ext(i)
          write(22,'(A)') "echo Time is `date`"
          write(22,'(A)') "echo Directory is `pwd`"
          write(22,*)


!       write(22,'(A)') "echo This job has allocated `$NPROCS` nodes"
!       write(22,*) 
          write(22,'(A)') 'nwchem scan_'//ext(i)//'.nw > scan_'//ext(i)//'.out'
          write(22,'(A)') "awk '/DFT energy/ {print $5}' scan_"//ext(i) &
          //".out > pot_"//ext(i)//".out"
          write(22,'(A)') 'python dipole_get.py scan_'//ext(i)//'.out'
                  if(fit_type.eq.1) then
                          write(22,'(A)') './fit_main < pot_'//ext(i)   &
                          //'.out > eig_'//ext(i)//'.dat'
                          else if(fit_type.eq.2) then
                                  write(22,'(A)') './fit_morse < pot_'//ext(i)           &
                                  //'.out > eig_'//ext(i)//'.dat'
                  endif
          write(22,*)
          write(22,'(A)') 'mv '//trim(scr_path)//'tmp_mb_'//ext(i)//'/scan_'//ext(i)    &
          //'.out '//trim(main_path)

          write(22,'(A)') 'mv '//trim(scr_path)//'tmp_mb_'//ext(i)//'/pot_'//ext(i)     &
          //'.out '//trim(main_path)

          write(22,'(A)') 'mv '//trim(scr_path)//'tmp_mb_'//ext(i)//'/eig_'//ext(i)     &
          //'.dat '//trim(main_path)

          write(22,'(A,A)') 'mv '//trim(scr_path)//'tmp_mb_'//ext(i)//'/eigs_fit.dat '  &
          //trim(main_path)//'eigs_fit_'//ext(i)//'.dat'

          write(22,'(A,A)') 'mv '//trim(scr_path)//'tmp_mb_'//ext(i)//'/eners_fit.dat ' &
          //trim(main_path)//'eners_fit_'//ext(i)//'.dat'

          write(22,'(A,A)') 'mv '//trim(scr_path)//'tmp_mb_'//ext(i)//'/dipole.dat '    &
          //trim(main_path)//'dipole_'//ext(i)//'.dat'

          write(22,'(A,A)') 'mv '//trim(scr_path)//'tmp_mb_'//ext(i)//'/trans_dip.dat ' &
          //trim(main_path)//'trdip_'//ext(i)//'.dat'


          write(22,'(A)') 'rm scan_'//ext(i)//'.*'
          write(22,'(A)') 'rm dvr_*'
          write(22,'(A)') 'rm dip*'
                  if(fit_type.eq.1) then
                          write(22,'(A)') 'rm fit_main'
                          write(22,'(A)') 'rm fit_guesses.in'
                          else if(fit_type.eq.2) then
                                  write(22,'(A)') 'rm fit_morse'
                                  write(22,'(A)') 'rm fit_morse.in'
                  endif
!         write(22,'(A)') 'rm ./*'
          write(22,'(A)') ' '
          write(22,'(A)') 'cd '//trim(main_path)
          write(22,'(A)') 'rmdir '//trim(scr_path)//'tmp_mb_'//ext(i)
          write(22,*)
          write(22,'(A)') "echo Ending Time is `date`"
          write(22,'(A)') 'exit 0'
          close(22)


! write to script file

        write(23,'(A)') 'msub scan_'//ext(i)//'.sh'

! skip the next nskip - 1 configurations in the trajectory

        do j = 1, nskip - 2
              it = it + 1
!             read(15,*)
!             read(15,*)
                      do k = 1, 5694
                              read(15,*) 
                      enddo
        enddo

!       read(15,*) test1
        write(6,'(A)') ' '
!       stop

! deallocate the changing arrays
        deallocate(hydroxyl_o)
        deallocate(hydroxyl_h)
        deallocate(hydroxyl_os)
        deallocate(hydroxyl_hs)
        deallocate(rO)
        deallocate(rH1)
        deallocate(rH2)
        deallocate(rOs)
        deallocate(rH1s)
        deallocate(rH2s)
        deallocate(xxx_cluster)
        deallocate(yyy_cluster)
        deallocate(zzz_cluster)
        deallocate(atm_cluster)
        deallocate(atm_cluster_plc)
        deallocate(si_cluster)
        deallocate(o_cluster)
        deallocate(d_conn)
        deallocate(inner)
        deallocate(outer)
        deallocate(irank)
        deallocate(drank)
        deallocate(cluster)
        deallocate(cluster_s)
        deallocate(inner_surface)
        deallocate(outer_surface)
        deallocate(inner_hyd)
        deallocate(outer_hyd)
        deallocate(drank_surf)
        deallocate(drank_hyd)
        
! major do loop
  enddo
  close(23)


! deallocate frozen Si and O and close

  deallocate(si)
  deallocate(frozen_o)

end program
