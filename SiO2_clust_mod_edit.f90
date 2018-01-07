module SiO2_cluster_mod
!*****************************************************************************************!
! Fortran module containing the modules required to do an ab initio calculation
! including the water and the silica framework with a selected silanol.
!
! Author : W.H. Thompson, J.A. Harvey, P.N. Wimalasiri
! Copyright Thompson Group March 2017
!
!*****************************************************************************************

  use kinds
  implicit none

  contains

  subroutine get_frozen_pore(nfrozen_o,n_pore_si,frozen_o,si)
  !************************************************************************
  !
  ! Subroutine for geting the positions of the frozen pore. These are 
  ! oxygens with the atom type of '4', and silicas with the atom type of 
  ! '1' in the original lammps data file.
  !
  ! Author - Jacob A. Harvey
  ! Copyright - Thompson Group Jan 2014
  !
  ! Updated - Apr 2016 - finds frozen oxygens AND silicas
  !
  !************************************************************************

	use kinds

	implicit none

        integer(kind=ip)                :: i,j,io,isi,nfrozeno,nsi
        real(kind=dp),dimension(3)      :: si_pos, frozeno_pos
        real(kind=dp), intent(inout)    :: frozen_o(:,:)
        real(kind=dp), intent(inout)    :: si(:,:)
        integer(kind=ip), intent(in)    :: nfrozen_o,n_pore_si


        character*2                     :: temp1,temp2

! write(6,*) 'started with si-o connectivity'
! initialize

        io = 0
        isi = 0
        
	frozen_o = 0.0_dp
        si = 0.0_dp

! open data file

        open(unit=10,file='si.xyz',status='old')
        open(unit=11,file='o.xyz',status='old')

! get different mass types from data file

        read(11,*) nfrozeno
        read(11,*)
        do i=1,nfrozen_o
                read(11,*) temp1, frozeno_pos(1), frozeno_pos(2), frozeno_pos(3)
                        io = io + 1
                        frozen_o(io,:) = frozeno_pos(:)
        end do
        close(11)

                read(10,*) nsi
                read(10,*)
                do j=1,nsi
                        read(10,*) temp2, si_pos(1), si_pos(2), si_pos(3)
                        isi = isi + 1           
                        si(isi,:) = si_pos(:)
                end do
                close(10)
        
!       print *, isi,   n_pore_si
!       print *, io,    nfrozen_o

!       print *, isi,   si(isi,1),      si(isi,2),      si(isi,3)
!       print *, io,    frozen_o(io,1), frozen_o(io,2), frozen_o(io,3)
 
    
        if (isi .ne. n_pore_si) then
            write(*,*) 'Error - number of calculated Si atoms is different than inputed'
            write(*,*) 'calculated = ' , isi , 'inputed = ' , n_pore_si
        endif

       return

  end subroutine get_frozen_pore

  subroutine si_o_connectivity(nfrozen_o,nhydroxyl,n_pore_si,si,  &
             hydroxyl_o,frozen_o,L,d_conn)
  !************************************************************************
  !
  ! Fortran subroutine for calculating the connectivity between Si and O 
  ! atoms with a matrix:
  !
  ! d_conn(i,j) = dist ; if the ith Si and jth O are within 2 A
  !                  0 ; otherwise
  !
  ! The distance 2 A is chosen from the Si-O radial distribution function
  !
  ! Author - Jacob A. Harvey
  ! Copyright - Thompson Group 
  !
  !************************************************************************
  
  use kinds
  
  implicit none
  
  integer(kind=ip),intent(in)           :: nfrozen_o,nhydroxyl,n_pore_si
  integer(kind=ip) :: jatm,i,j
  
  real(kind=dp),intent(in)              :: si(:,:), frozen_o(:,:), L(:)
  real(kind=dp),intent(in)              :: hydroxyl_o(:,:)
  real(kind=dp),dimension(n_pore_si,nhydroxyl+nfrozen_o),intent(out) :: d_conn
  real(kind=dp),dimension(3)            :: d
  real(kind=dp)                         :: dist
  
  ! initialize
   d_conn = 0.0_dp
  
  ! loop over all Si's
  	do i = 1,n_pore_si
  		jatm = 0
  
  ! loop over hydroxyl O's
  		do j = 1,nhydroxyl
  		
  ! calculate Si-O distance
                        d(:) = si(i,:) - hydroxyl_o(j,:)
                  	d(:) = d(:) - (L(:) * anint(d(:)/L(:)))
  			dist = dsqrt(d(1)**2+d(2)**2+d(3)**2)
  
  ! mark if distance is less than 2
  			if (dist .lt. 4.0) then
                 	        d_conn(i,j) = dist
  !                             if (i==104) then
  !                             write(*,'(A,i3,4X,i3,4X,f8.2)') &
  !                             'd_conn= ', i, j, d_conn(i,j)
  !                             endif 
                	endif
  		end do
  
  ! now loop over all frozen O's
  		do j = nhydroxyl + 1 , nhydroxyl + nfrozen_o
  		        jatm = jatm + 1
  
  ! calculate Si-O distance
                        d(:) = si(i,:) - frozen_o(j-nhydroxyl,:)
  		        d(:) = d(:) - (L(:) * anint(d(:)/L(:)))
                        dist = dsqrt(d(1)**2+d(2)**2+d(3)**2)
  
  ! mark if distance is less than 2
                        if (dist .lt. 3.0) then
                                d_conn(i,j) = dist
  !                             if (i==104) then
  !                             write(*,'(A,i3,4X,i3,4X,f8.2)') &
  !                             'd_conn= ', i, j, d_conn(i,j)
  !                             endif
                        endif
                end do
        end do 
  return
  
  end subroutine si_o_connectivity
  
  subroutine find_central_cluster(n_pore_si,nhydroxyl,nfrozen_o,ioh,  &           
             L,si,hydroxyl_o,hydroxyl_h,frozen_o,d_conn,iatm,  &
             xxx_cluster,yyy_cluster,zzz_cluster,cluster_size, &
             atm_cluster,atm_cluster_plc,o_cluster,si_cluster)
  !************************************************************************
  !
  ! Fortran subroutine for locating the central SiO2H for our cluster.
  !
  ! Author - Jacob A. Harvey
  ! Copyright - Thompson Group May 2016
  !
  ! Updated - Nov 2016 - Modified from 2 OH version
  !
  !************************************************************************
  
  use kinds
 
  implicit none
  
  integer(kind=ip),dimension(:),intent(inout)   :: atm_cluster,atm_cluster_plc
  integer(kind=ip),intent(in)   :: n_pore_si,nhydroxyl,nfrozen_o, &  
                                   cluster_size
  integer(kind=ip)              :: i,j,itemp,itemp2,hyd_cnt,fro_cnt,si_cnt,ioh_old
  integer(kind=ip),intent(inout):: iatm,ioh
  
  real(kind=dp),intent(in)      :: hydroxyl_o(:,:),d_conn(:,:),L(:)
  real(kind=dp),intent(in)      :: hydroxyl_h(:,:),frozen_o(:,:),si(:,:)
  real(kind=dp),dimension(:),intent(inout)      :: xxx_cluster, yyy_cluster, &
                                                   zzz_cluster
  real(kind=dp)                 :: dx,dy,dz,dist,boz,min_dist
  real(kind=dp),dimension(nfrozen_o+nhydroxyl)  :: temp_array
  
  logical,dimension(:),intent(inout)            :: si_cluster
  logical,dimension(:),intent(inout)            :: o_cluster
  logical :: found
  
! save silanols to cluster

  fro_cnt=0_ip
  
         iatm = iatm + 1
         xxx_cluster(iatm) = hydroxyl_h(ioh,1)
         yyy_cluster(iatm) = hydroxyl_h(ioh,2)
         zzz_cluster(iatm) = hydroxyl_h(ioh,3)
         atm_cluster(iatm) = 3
         atm_cluster_plc(iatm) = 1

!        write(6,'(A,f20.10)') 'xxx(1)= ',xxx_cluster(iatm)

         iatm = iatm + 1
         xxx_cluster(iatm) = hydroxyl_o(ioh,1)
         yyy_cluster(iatm) = hydroxyl_o(ioh,2)
         zzz_cluster(iatm) = hydroxyl_o(ioh,3)
         atm_cluster(iatm) = 2
         atm_cluster_plc(iatm) = 1
  
!        write(6,'(A,f20.10)') 'xxx(2)= ',xxx_cluster(iatm)

         iatm = iatm + 1
         xxx_cluster(iatm) = si(ioh,1)
         yyy_cluster(iatm) = si(ioh,2)
         zzz_cluster(iatm) = si(ioh,3)
         atm_cluster(iatm) = 1
         atm_cluster_plc(iatm) = 1
  
! write(6,*) ' made it through third set'
  
         si_cluster(ioh) = .true.
         o_cluster(ioh) = .true.
  
  hyd_cnt=1_ip
  si_cnt=1_ip

! find the nearest oxygen
  
  write(6,*) ' find: made it through preamble'
  

  found = .false.
  
  do j=1,20 
      min_dist = 1000.0_dp
      if ((hyd_cnt.eq.2) .and. (fro_cnt.eq.3)) then
          exit
      endif

      do i=1,nhydroxyl
  
!         write(6,*) min_dist 
          if ((d_conn(ioh,i) .gt. 0.0) .and. (.not. o_cluster(i)) .and. &
          (d_conn(ioh,i) .lt. min_dist) .and. (hyd_cnt.lt.2)) then
  
              write(*,'(A,i2,4X,i3,4X,f8.2,4X,f8.2)') "Found new hydroxyl O for Si ", ioh, i, &
              min_dist, d_conn(ioh,i)
              itemp = i
              min_dist = d_conn(ioh,i)
              found = .true.
  	  endif
      enddo

!     write(6,*) ' find: made it through nhydroxyl loop'  
!     if (.not.found) then 
      min_dist = 1000.0_dp
          do i=nhydroxyl+1,nhydroxyl+nfrozen_o
  
              if ((d_conn(ioh,i) .gt. 0.0) .and. (.not. o_cluster(i)) .and. &
              (d_conn(ioh,i) .lt. min_dist) .and. (fro_cnt.lt.3)) then
  
                  write(*,'(A,i2,4X,i3,4X,f8.2,4X,f8.2)') "Found new frozen O for Si ", ioh, i, &
                  min_dist, d_conn(ioh,i)
                  itemp = i
                  min_dist = d_conn(ioh,i)
                  found = .true.
              endif
          enddo
!     endif 
!     write(6,*) ' find: made it through frozen loop'    

      if (.not. found) write(*,*) "Did not find a new oxygen for Si ", ioh

      if ((itemp .gt. nhydroxyl) .and. (fro_cnt.le.3)) then
      iatm = iatm+1
      fro_cnt = fro_cnt+1
      xxx_cluster(iatm) = frozen_o((itemp - nhydroxyl),1)
      xxx_cluster(iatm) = xxx_cluster(iatm) - L(1)* &
      dnint((xxx_cluster(iatm)-xxx_cluster(3))/L(1))
      yyy_cluster(iatm) = frozen_o((itemp - nhydroxyl),2)
      yyy_cluster(iatm) = yyy_cluster(iatm) - L(2)* &
      dnint((yyy_cluster(iatm)-yyy_cluster(3))/L(2))
      zzz_cluster(iatm) = frozen_o((itemp - nhydroxyl),3)
      atm_cluster(iatm) = 4

      o_cluster(itemp) = .true.
      atm_cluster_plc(iatm) = 1
  
      else if ((itemp .le. nhydroxyl) .and. (hyd_cnt.le.2)) then
          iatm = iatm+1
          hyd_cnt = hyd_cnt+1
          xxx_cluster(iatm) = hydroxyl_o(itemp,1)
          xxx_cluster(iatm) = xxx_cluster(iatm) - L(1)* &
          dnint((xxx_cluster(iatm)-xxx_cluster(3))/L(1))
          yyy_cluster(iatm) = hydroxyl_o(itemp,2)
          yyy_cluster(iatm) = yyy_cluster(iatm) - L(2)* &
          dnint((yyy_cluster(iatm)-yyy_cluster(3))/L(2))
          zzz_cluster(iatm) = hydroxyl_o(itemp,3)
          atm_cluster(iatm) = 2

          o_cluster(itemp) = .true.
          atm_cluster_plc(iatm) = 1

          iatm = iatm+1
          xxx_cluster(iatm) = hydroxyl_h(itemp,1)
          xxx_cluster(iatm) = xxx_cluster(iatm) - L(1)* &
          dnint((xxx_cluster(iatm)-xxx_cluster(3))/L(1))
          yyy_cluster(iatm) = hydroxyl_h(itemp,2)
          yyy_cluster(iatm) = yyy_cluster(iatm) - L(2)* &
          dnint((yyy_cluster(iatm)-yyy_cluster(3))/L(2))
          zzz_cluster(iatm) = hydroxyl_h(itemp,3)
          atm_cluster(iatm) = 3
          atm_cluster_plc(iatm) = 1

          if (itemp.le.32) then
              iatm = iatm+1
              xxx_cluster(iatm) = si(itemp,1)
              xxx_cluster(iatm) = xxx_cluster(iatm) - L(1)* &
                   dnint((xxx_cluster(iatm)-xxx_cluster(3))/L(1))
              yyy_cluster(iatm) = si(itemp,2)
              yyy_cluster(iatm) = yyy_cluster(iatm) - L(2)* &
                   dnint((yyy_cluster(iatm)-yyy_cluster(3))/L(2))
              zzz_cluster(iatm) = si(itemp,3)
              atm_cluster(iatm) = 1

              ioh_old = ioh
              ioh = itemp

              si_cluster(itemp) = .true.
              atm_cluster_plc(iatm) = 1

              else if (itemp.gt.32) then

                  itemp2=itemp-32
                  iatm = iatm+1
                  xxx_cluster(iatm) = si((32+(itemp2/2)+mod(itemp2,2)),1)
                  xxx_cluster(iatm) = xxx_cluster(iatm) - L(1)* &
                       dnint((xxx_cluster(iatm)-xxx_cluster(3))/L(1))
                  yyy_cluster(iatm) = si((32+(itemp2/2)+mod(itemp2,2)),2)
                  yyy_cluster(iatm) = yyy_cluster(iatm) - L(2)* &
                       dnint((yyy_cluster(iatm)-yyy_cluster(3))/L(2))
                  zzz_cluster(iatm) = si((32+(itemp2/2)+mod(itemp2,2)),3)
                  atm_cluster(iatm) = 1

                  ioh_old = ioh
                  ioh = 32+(itemp2/2)+mod(itemp2,2)

                  si_cluster(32+(itemp2/2)+mod(itemp2,2)) = .true.
                  atm_cluster_plc(iatm) = 1
              endif
                    
      endif

  write(6,'(A,i)') 'ioh= ',ioh
  enddo

! write(6,'(A,I)') ' hydroxyl= ', hyd_cnt
! write(6,'(A,I)') ' frozen= ', fro_cnt

  ioh = ioh_old
! write(*,'(i3,4x,f8.2,4X,f8.2,4X,f8.2)') iatm, xxx_cluster(iatm), &
! yyy_cluster(iatm), zzz_cluster(iatm)
  return  
  end subroutine find_central_cluster

  subroutine calc_com_cc(cluster_size,xxx_cluster,yyy_cluster, &
  	     zzz_cluster,atm_cluster,box,boy,iatm,comx,comy,comz)
  !************************************************************************
  !
  ! Fortran subroutine for calculating the COM of the central cluster.
  !
  ! Author - Jacob A. Harvey
  ! Copyright - Thompson Group May 2016
  !
  !************************************************************************
  
  use kinds
  implicit none
  
  integer(kind=ip)              :: cluster_size,iatm,i
  integer(kind=ip),intent(inout),dimension(:) :: atm_cluster
  real(kind=dp)                 :: dy,dx,mass,mass_tot
  real(kind=dp),intent(out)     :: comx,comy,comz
  real(kind=dp),intent(in)      :: box,boy    
  real(kind=dp),intent(inout),dimension(:)      :: xxx_cluster,yyy_cluster,zzz_cluster
  
  comx = 0.0_dp
  comy = 0.0_dp
  comz = 0.0_dp
  mass_tot = 0.0_dp
  
  ! loop over all atoms in cluster
  
  do i = 1 , iatm
  
  ! set mass
  
        if (atm_cluster(i) .eq. 1) then
  		mass = 28.0855_dp
  	else if ((atm_cluster(i) .eq. 2) .or. (atm_cluster(i) .eq. 4)) then
  		mass = 15.999_dp
  	else if (atm_cluster(i) .eq. 3) then
  		mass = 1.008_dp
  	endif
  
  ! wrap all atoms to Si box
  
  	dx = xxx_cluster(i) - xxx_cluster(3)
  	comx = comx + (mass * (xxx_cluster(i) - (box * anint(dx/box))))
  	
  	dy = yyy_cluster(i) - yyy_cluster(3)
  	comy = comy + (mass * (yyy_cluster(i) - (boy * anint(dy/boy))))
  	comz = comz + (mass * zzz_cluster(i))
  	
  	mass_tot = mass_tot + mass
  enddo
  
  comx = comx / mass_tot
  comy = comy / mass_tot
  comz = comz / mass_tot

  return
  end subroutine calc_com_cc
  
  subroutine find_cluster(n_pore_si,nhydroxyl,nfrozen_o,ioh,         &
             box,boy,si,hydroxyl_o,hydroxyl_h,frozen_o,d_conn,iatm,  &
             xxx_cluster,yyy_cluster,zzz_cluster,cluster_size,       &
             atm_cluster,atm_cluster_plc,o_cluster,si_cluster,comx,comy,comz) 
  !************************************************************************
  !
  ! Fortran subroutine for finding the atoms for a cluster centered around
  ! a central Si2O5H2 cluster.
  !
  ! Author - Jacob A. Harvey
  ! Copyright - Thompson Group Apr 2016
  !
  !************************************************************************  
  use kinds
  implicit none
                     
  integer(kind=ip),intent(in)   :: n_pore_si,nhydroxyl,nfrozen_o,ioh, &
                                   cluster_size
  integer(kind=ip)              :: i,itemp,temp_si,j,temp_o     
  real(kind=dp),intent(in)      :: hydroxyl_o(:,:), d_conn(:,:)
  real(kind=dp),intent(in)      :: hydroxyl_h(:,:),frozen_o(:,:),si(:,:)
  real(kind=dp),intent(inout),dimension(:)      :: xxx_cluster,yyy_cluster,zzz_cluster
  integer(kind=ip),intent(inout),dimension(:)   :: atm_cluster,atm_cluster_plc
  integer(kind=ip),intent(inout)                :: iatm
  logical,intent(inout),dimension(:)            :: si_cluster,o_cluster !oh_cluster
  real(kind=dp),intent(in)      :: box,boy,comx,comy,comz 
  real(kind=dp)                 :: dx,dy,dz,dist,min_dist
  logical                       :: found

! loop over cluster size

  do i=1,cluster_size-1
! first find the next Si

        min_dist = 1000.0_dp

        do j = 1 , n_pore_si

! only look at Si's that aren't already in the cluster
                if (.not. si_cluster(j)) then
! calculate distance
                dx = si(j,1) - comx
                dy = si(j,2) - comy
                dz = si(j,3) - comz

                dx = dx - (box * anint(dx/box))
                dy = dy - (boy * anint(dy/boy))

                dist = dsqrt(dx**2+dy**2+dz**2)
                        if (dist .lt. min_dist) then
                                min_dist = dist
                                temp_si = j
                        endif
                endif
        enddo

! save Si

        iatm = iatm + 1
        xxx_cluster(iatm) = si(temp_si,1)
        xxx_cluster(iatm) = xxx_cluster(iatm) - box*  &
        dnint((xxx_cluster(iatm) - comx)/box)
        yyy_cluster(iatm) = si(temp_si,2)
        yyy_cluster(iatm) = yyy_cluster(iatm) - boy*  &
        dnint((yyy_cluster(iatm) - comy)/boy)
        zzz_cluster(iatm) = si(temp_si,3)
        atm_cluster(iatm) = 1
        atm_cluster_plc(iatm) = 1
        si_cluster(temp_si) = .true.

! now find the 2 closest oxygens to that Si
        min_dist = 1000.0_dp
        found = .false.

        do j = 1 , nhydroxyl
                if ((d_conn(temp_si,j) .gt. 0.0) .and. (.not. o_cluster(j)) .and. &
                   (d_conn(temp_si,j) .lt. min_dist)) then
                        temp_o = j
                        min_dist = d_conn(temp_si,j)
                        found = .true.
                endif
        enddo
        
        do j = nhydroxyl + 1 , nhydroxyl + nfrozen_o
                if ((d_conn(temp_si,j) .gt. 0.0) .and. (.not. o_cluster(j)) .and. &
                (d_conn(temp_si,j) .lt. min_dist)) then
                        temp_o = j
                        min_dist = d_conn(temp_si,j)
                        found = .true.
                endif
        enddo

        if (.not. found) then 
                write(*,*) "Did not find a closest oxygen for Si: " , temp_si
                write(*,'(i3,4x,f8.2,4X,f8.2,4X,f8.2)') temp_si, si(temp_si,1), &
                si(temp_si,2), si(temp_si,3)
!               iatm = iatm - 1        
!               cycle
        endif        

        if (temp_o .gt. nhydroxyl) then
                
                iatm = iatm + 1
                xxx_cluster(iatm) = frozen_o((temp_o - nhydroxyl),1)
                xxx_cluster(iatm) = xxx_cluster(iatm) - box* &
                dnint((xxx_cluster(iatm)-xxx_cluster(iatm-1))/box)
                yyy_cluster(iatm) = frozen_o((temp_o - nhydroxyl),2)
                yyy_cluster(iatm) = yyy_cluster(iatm) - boy* &
                dnint((yyy_cluster(iatm)-yyy_cluster(iatm-1))/boy)
                zzz_cluster(iatm) = frozen_o((temp_o - nhydroxyl),3)
                atm_cluster(iatm) = 4
                atm_cluster_plc(iatm) = 1
        
                else if (temp_o .le. nhydroxyl) then
                        
                        iatm = iatm + 1
                        xxx_cluster(iatm) = hydroxyl_o((temp_o),1)
                        xxx_cluster(iatm) = xxx_cluster(iatm) - box* &
                        dnint((xxx_cluster(iatm)-xxx_cluster(iatm-1))/box)
                        yyy_cluster(iatm) = hydroxyl_o((temp_o),2)
                        yyy_cluster(iatm) = yyy_cluster(iatm-1) - boy* &
                        dnint((yyy_cluster(iatm)-yyy_cluster(temp_si))/boy)
                        zzz_cluster(iatm) = hydroxyl_o((temp_o),3)
                        atm_cluster(iatm) = 2
                        atm_cluster_plc(iatm) = 1

!                       oh_cluster(temp_o) = .true.

                        iatm = iatm + 1
                        xxx_cluster(iatm) = hydroxyl_h((temp_o),1)
                        xxx_cluster(iatm) = xxx_cluster(iatm) - box* &
                        dnint((xxx_cluster(iatm)-xxx_cluster(iatm-1))/box)
                        yyy_cluster(iatm) = hydroxyl_h((temp_o),2)
                        yyy_cluster(iatm) = yyy_cluster(iatm-1) - boy* &
                        dnint((yyy_cluster(iatm)-yyy_cluster(temp_si))/boy)
                        zzz_cluster(iatm) = hydroxyl_h((temp_o),3)
                        atm_cluster(iatm) = 3
                        atm_cluster_plc(iatm) = 1

        endif

        o_cluster(temp_o) = .true.
        min_dist = 1000.0_dp
        found = .false.

        do j = 1 , nhydroxyl
                if ((d_conn(temp_si,j) .gt. 0.0) .and. (.not. o_cluster(j)) .and. &
                (d_conn(temp_si,j) .lt. min_dist)) then

                temp_o = j
                min_dist = d_conn(temp_si,j)
                found = .true.
                endif
        enddo

        do j = nhydroxyl + 1 , nhydroxyl + nfrozen_o
                if ((d_conn(temp_si,j) .gt. 0.0) .and. (.not. o_cluster(j)) .and. &
                (d_conn(temp_si,j) .lt. min_dist)) then
        
                temp_o = j
                min_dist = d_conn(temp_si,j)
                found = .true.

                endif
        enddo

        if (.not. found) then
                write(*,*) "Did not find a closest oxygen for Si: " , temp_si
                write(*,'(i3,4x,f8.2,4X,f8.2,4X,f8.2)') temp_si, si(temp_si,1), &
                si(temp_si,2), si(temp_si,3)
!               iatm = iatm - 2
!               cycle
        end if

        if (temp_o .gt. nhydroxyl) then
                
                iatm = iatm + 1
                xxx_cluster(iatm) = frozen_o((temp_o - nhydroxyl),1)
                xxx_cluster(iatm) = xxx_cluster(iatm) - box* &
                dnint((xxx_cluster(iatm)-xxx_cluster(iatm-2))/box)
                yyy_cluster(iatm) = frozen_o((temp_o - nhydroxyl),2)
                yyy_cluster(iatm) = yyy_cluster(iatm) - boy* &
                dnint((yyy_cluster(iatm)-yyy_cluster(iatm-2))/boy)
                zzz_cluster(iatm) = frozen_o((temp_o - nhydroxyl),3)
                atm_cluster(iatm) = 4
                atm_cluster_plc(iatm) = 1
        
                else if (temp_o .le. nhydroxyl) then
                
                        iatm = iatm + 1
                        xxx_cluster(iatm) = hydroxyl_o((temp_o),1)
                        xxx_cluster(iatm) = xxx_cluster(iatm) - box* &
                        dnint((xxx_cluster(iatm)-xxx_cluster(iatm-2))/box)
                        yyy_cluster(iatm) = hydroxyl_o((temp_o),2)
                        yyy_cluster(iatm) = yyy_cluster(iatm) - boy* &
                        dnint((yyy_cluster(iatm)-yyy_cluster(iatm-2))/boy)
                        zzz_cluster(iatm) = hydroxyl_o((temp_o),3)
                        atm_cluster(iatm) = 2
                        atm_cluster_plc(iatm) = 1

!                       oh_cluster(temp_o) = .true.

                        iatm = iatm + 1
                        xxx_cluster(iatm) = hydroxyl_h((temp_o),1)
                        xxx_cluster(iatm) = xxx_cluster(iatm) - box* &
                        dnint((xxx_cluster(iatm)-xxx_cluster(iatm-2))/box)
                        yyy_cluster(iatm) = hydroxyl_h((temp_o),2)
                        yyy_cluster(iatm) = yyy_cluster(iatm) - boy* &
                        dnint((yyy_cluster(iatm)-yyy_cluster(iatm-2))/boy)
                        zzz_cluster(iatm) = hydroxyl_h((temp_o),3)
                        atm_cluster(iatm) = 3
                        atm_cluster_plc(iatm) = 1

        endif
        o_cluster(temp_o) = .true.
  
  enddo
! write(*,'(i3,4x,f8.2,4X,f8.2,4X,f8.2)') iatm, xxx_cluster(iatm), &
! yyy_cluster(iatm), zzz_cluster(iatm)

! do j=1,(nhydroxyl+nfrozen_o)
!       if (o_cluster(j)) then
!       write (6,*) j
!       endif
! enddo
! stop
  return
  end subroutine find_cluster

  subroutine add_water(nwat,boz,rO,rH1,rH2,iatm,xxx_cluster,          &
             yyy_cluster,zzz_cluster,cluster_size,atm_cluster,atm_cluster_plc)
  !************************************************************************
  !
  ! Fortran subroutine for adding waters to the cluster. Waters are chosen
  ! such that the water O is within 7.831 A of the silanol H.
  !
  ! Author - Jacob A. Harvey
  ! Copyright - Thompson Group Nov 2016
  !
  !************************************************************************

  use kinds
  implicit none

  real(kind=dp),intent(inout),dimension(:)      :: xxx_cluster,yyy_cluster,zzz_cluster
  integer(kind=ip),intent(in)                   :: cluster_size,nwat
  integer(kind=ip),intent(inout)                :: iatm
  real(kind=dp),intent(in)                      :: rO(:,:),rH1(:,:),rH2(:,:),boz
  integer(kind=ip),intent(inout),dimension(:)   :: atm_cluster,atm_cluster_plc
  integer(kind=ip)      :: i

  real(kind=dp)         :: dx,dy,dz,dist

! loop over waters

  do i=1,nwat  
  
! determine distance to silanol H
        
        dx = rO(i,1) - xxx_cluster(1)
        dy = rO(i,2) - yyy_cluster(1)
        dz = rO(i,3) - zzz_cluster(1)

        dz = dz - (boz*anint(dz/boz))
        
        dist = dsqrt(dx**2+dy**2+dz**2)
  
! save water to cluster if within 7.831 A of silanol H

        if ((dist .lt. 7.831_dp) .and. (dist .gt. 3.831_dp)) then
                iatm = iatm + 1
                xxx_cluster(iatm) = rO(i,1)
                yyy_cluster(iatm) = rO(i,2)
                zzz_cluster(iatm) = rO(i,3)
                atm_cluster(iatm) = 5
                atm_cluster_plc(iatm) = 2

                iatm = iatm + 1
                xxx_cluster(iatm) = rH1(i,1)
                yyy_cluster(iatm) = rH1(i,2)
                zzz_cluster(iatm) = rH1(i,3)
                atm_cluster(iatm) = 6
                atm_cluster_plc(iatm) = 2
                
                iatm = iatm + 1
                xxx_cluster(iatm) = rH2(i,1)
                yyy_cluster(iatm) = rH2(i,2)
                zzz_cluster(iatm) = rH2(i,3)
                atm_cluster(iatm) = 7
                atm_cluster_plc(iatm) = 2
                
                else if (dist .lt. 3.831_dp) then
                        iatm = iatm + 1
                        xxx_cluster(iatm) = rO(i,1)
                        yyy_cluster(iatm) = rO(i,2)
                        zzz_cluster(iatm) = rO(i,3)
                        atm_cluster(iatm) = 5
                        atm_cluster_plc(iatm) = 1

                        iatm = iatm + 1
                        xxx_cluster(iatm) = rH1(i,1)
                        yyy_cluster(iatm) = rH1(i,2)
                        zzz_cluster(iatm) = rH1(i,3)
                        atm_cluster(iatm) = 6
                        atm_cluster_plc(iatm) = 1

                        iatm = iatm + 1
                        xxx_cluster(iatm) = rH2(i,1)
                        yyy_cluster(iatm) = rH2(i,2)
                        zzz_cluster(iatm) = rH2(i,3)
                        atm_cluster(iatm) = 7
                        atm_cluster_plc(iatm) = 1
                endif 

  enddo
  
  write(*,'(i3,4x,f8.2,4X,f8.2,4X,f8.2)') iatm, xxx_cluster(iatm), &
  yyy_cluster(iatm), zzz_cluster(iatm)

  write(6,*)    iatm

  return
  end subroutine add_water
 
end module SiO2_cluster_mod
