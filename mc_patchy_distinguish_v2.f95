!	Monte-Carlo PROGRAM TO SIMULATE patchy spheres with crowding agents
!	modified date - 18/12/2022
!	need to check the update position of patch positions when it goes beyond the boundary
!


	program mc

	implicit none

	real*8,allocatable :: uh(:,:,:),ah(:,:,:),ap(:,:),a_p(:,:),ac(:,:),a_c(:,:),x(:),patch_rsqrt(:,:)
	real*8,allocatable :: uhn(:,:,:),ahn(:,:,:),a_h(:,:,:),rot(:,:),apn(:,:)
	real*8		   :: boxl,yt,OEp,OEc,OEpc,rsqrt,rhex,rtwelve, R, T, RT, p_r
	real*8		   :: Ei, Ef, dE, d_m, theta11, theta12, theta21, theta22, two_pie, sin_theta, sin_psi, sin_phi, rand, psi, phi, onemcospsi
	real*8 		   :: ny,nxz,nxy,nxsq,nx,nzsq, nz, nyz, nysq, dpsi_max, detrot, cos_theta, cos_psi, cos_phi,rsqrt_old
	real*8 		   :: eps11, eps12, eps22, sgmahex11, sgmahex12, sgmahex22, sgma
	real*8		   :: accp_patch,rej_patch,accp_crowd,rej_crowd,Etotal_final
	integer*8	   :: j, iatom, icrowd, natom, ncrowd, m, jatom, jcrowd, imc, nmc, b, c, ipatch, npatch, jpatch

!	OPEN FILES
	open(unit=1,file='input_mc.dat')
	open(unit=2,file='ini_simu.xyz')
	open(unit=10,file='output.dat')
	open(unit=11,file='accep_prob.dat')
	open(unit=13,file='traj.xyz')

!	READ THE INPUT DATA
!
	read(1,*)
    	read(1,*) boxl, nmc, natom, ncrowd, d_m, dpsi_max, R, T
	read(1,*)
	read(1,*) eps11, eps12, eps22
	read(1,*)
	read(1,*) sgma, sgmahex11, sgmahex12, sgmahex22

!	ALLOCATE DATA
	allocate ( ah(1:natom,1:2,1:3),ap(1:natom,1:3),a_p(1:natom,1:3),ac(1:ncrowd,1:3),a_c(1:ncrowd,1:3),x(1:3),patch_rsqrt(1:natom,1:2))
	allocate ( uh(1:natom,1:2,1:3),uhn(1:natom,1:2,1:3),ahn(1:natom,1:2,1:3),a_h(1:natom,1:2,1:3),rot(1:3,1:3),apn(1:natom,1:3) )
	RT = R * T; npatch = natom
    
! 	WRITE IMPORTANT PARAMTERS OF THESE MC SIMULATION

	write(10,*) 'PARAMTERS OF MC SIMULATION'
	write(10,*) 'boxlength	No of MC Steps	No. of Particles	No. of Crowding agents	d_m	R	T'
	write(10,*) boxl, nmc, natom, ncrowd, d_m, R, T

! 	GENERATE THE INITIAL CONFIGURATION FILE

!	GENERATE THE ATOM POSITIONS RANDOMLY

20	do iatom = 1, natom
		call random_number(yt); ap(iatom,1) = yt*boxl
		call random_number(yt);	ap(iatom,2) = yt*boxl
		call random_number(yt);	ap(iatom,3) = yt*boxl
	enddo !i

!	GENERATE THE PATCH POSITIONS 

	do iatom = 1,natom
		ah(iatom,1,1) = ap(iatom,1) + sgma; ah(iatom,1,2) = ap(iatom,2) ; ah(iatom,1,3) = ap(iatom,3)
		ah(iatom,2,1) = ap(iatom,1) + sgma*cos(1.5708) ; 
		ah(iatom,2,2) = ap(iatom,2) + sgma*sin(1.5708) ; ah(iatom,2,3) = ap(iatom,3)
	enddo !iatom

	do iatom = 1,natom
		if(ah(iatom,1,1).ge.boxl) goto 20
		if(ah(iatom,1,1).lt.0) goto 20
		if(ah(iatom,1,2).ge.boxl) goto 20
		if(ah(iatom,1,2).lt.0) goto 20
		if(ah(iatom,1,3).ge.boxl) goto 20
		if(ah(iatom,1,3).lt.0) goto 20

		if(ah(iatom,2,1).ge.boxl) goto 20
		if(ah(iatom,2,1).lt.0) goto 20 	
		if(ah(iatom,2,2).ge.boxl) goto 20
		if(ah(iatom,2,2).lt.0) goto 20 	
		if(ah(iatom,2,3).ge.boxl) goto 20
		if(ah(iatom,2,3).lt.0) goto 20 	
	enddo !iatom		

!	GENERATE THE CROWD POSITIONS RANDOMLY

	do icrowd = 1, ncrowd
		call random_number(yt); ac(icrowd,1) = yt*boxl
		call random_number(yt);	ac(icrowd,2) = yt*boxl
		call random_number(yt);	ac(icrowd,3) = yt*boxl
	enddo !i

!	CALCULATE ANGLE BETWEEN PATCHES OF DIFFERENT PARTICLES

	do iatom = 1,natom
		patch_rsqrt(iatom,:) = 0
		uh(iatom,:,:) = 0

			uh(iatom,1,:) = ah(iatom,1,:)-ap(iatom,:)
			uh(iatom,1,:) = uh(iatom,1,:)-boxl*anint(uh(iatom,1,:)/boxl)
			patch_rsqrt(iatom,1) = sqrt(uh(iatom,1,1)*uh(iatom,1,1)+uh(iatom,1,2)*uh(iatom,1,2)+uh(iatom,1,3)*uh(iatom,1,3)) 
			write(*,*) iatom,uh(iatom,1,:),patch_rsqrt(iatom,1)
		
			uh(iatom,1,:) = uh(iatom,1,:)/patch_rsqrt(iatom,1)
			patch_rsqrt(iatom,1) = sqrt(uh(iatom,1,1)*uh(iatom,1,1)+uh(iatom,1,2)*uh(iatom,1,2)+uh(iatom,1,3)*uh(iatom,1,3)) 
			write(*,*) iatom,uh(iatom,1,:),patch_rsqrt(iatom,1)
			
			uh(iatom,2,:) = ah(iatom,2,:)-ap(iatom,:)
			uh(iatom,2,:) = uh(iatom,2,:)-boxl*anint(uh(iatom,2,:)/boxl)
			patch_rsqrt(iatom,2) = sqrt(uh(iatom,2,1)*uh(iatom,2,1)+uh(iatom,2,2)*uh(iatom,2,2)+uh(iatom,3,3)*uh(iatom,3,3))
			write(*,*) iatom,uh(iatom,2,:),patch_rsqrt(iatom,2)

			uh(iatom,2,:) = uh(iatom,2,:)/patch_rsqrt(iatom,2)
			patch_rsqrt(iatom,2) = sqrt(uh(iatom,2,1)*uh(iatom,2,1)+uh(iatom,2,2)*uh(iatom,2,2)+uh(iatom,2,3)*uh(iatom,2,3)) 
			write(*,*) iatom,uh(iatom,2,:),patch_rsqrt(iatom,2)
			

			write(*,*) patch_rsqrt(iatom,:)

			do jatom = 1,natom
			if (iatom.ne.jatom) then
				x(:)=ap(jatom,:)-ap(iatom,:) !x() is the distance between ith and jth particle
				rsqrt_old=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
				write(*,*) x(:), rsqrt_old
				!x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention
				rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
				write(*,*) x(:), rsqrt

				x(:) = x(:)/rsqrt
				rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
				write(*,*) x(:), rsqrt

				uh(jatom,1,:) = ah(jatom,1,:)-ap(jatom,:)
				uh(jatom,1,:) = uh(jatom,1,:)-boxl*anint(uh(jatom,1,:)/boxl)
				patch_rsqrt(jatom,1) = sqrt(uh(jatom,1,1)*uh(jatom,1,1)+uh(jatom,1,2)*uh(jatom,1,2)+uh(jatom,1,3)*uh(jatom,1,3)) 
				write(*,*) jatom,uh(jatom,1,:),patch_rsqrt(jatom,1)
			
				uh(jatom,2,:) = ah(jatom,2,:)-ap(jatom,:)
				uh(jatom,2,:) = uh(jatom,2,:)-boxl*anint(uh(jatom,2,:)/boxl)
				patch_rsqrt(jatom,2) = sqrt(uh(jatom,2,1)*uh(jatom,2,1)+uh(jatom,2,2)*uh(jatom,2,2)+uh(jatom,3,3)*uh(jatom,3,3))
				write(*,*) jatom,uh(jatom,2,:),patch_rsqrt(jatom,2)

				uh(jatom,1,:) = uh(jatom,1,:)/patch_rsqrt(jatom,1)
				patch_rsqrt(jatom,1) = sqrt(uh(jatom,1,1)*uh(jatom,1,1)+uh(jatom,1,2)*uh(jatom,1,2)+uh(jatom,1,3)*uh(jatom,1,3)) 
				write(*,*) jatom,uh(jatom,1,:),patch_rsqrt(jatom,1)

				uh(jatom,2,:) = uh(jatom,2,:)/patch_rsqrt(jatom,2)
				patch_rsqrt(jatom,2) = sqrt(uh(jatom,2,1)*uh(jatom,2,1)+uh(jatom,2,2)*uh(jatom,2,2)+uh(jatom,2,3)*uh(jatom,2,3)) 
				write(*,*) jatom,uh(jatom,2,:),patch_rsqrt(jatom,2)

				theta11 = dacos ( (x(1)*uh(iatom,1,1) + x(2)*uh(iatom,1,2) + x(3)*uh(iatom,1,3)) / (rsqrt * patch_rsqrt(iatom,1)))
				theta12 = dacos ( (x(1)*uh(iatom,2,1) + x(2)*uh(iatom,2,2) + x(3)*uh(iatom,2,3)) / (rsqrt * patch_rsqrt(iatom,2)))
				theta21 = dacos ( (-x(1)*uh(jatom,1,1) - x(2)*uh(jatom,1,2) - x(3)*uh(jatom,1,3)) / (rsqrt * patch_rsqrt(jatom,1)))
				theta22 = dacos ( (-x(1)*uh(jatom,2,1) - x(2)*uh(jatom,2,2) - x(3)*uh(jatom,2,3)) / (rsqrt * patch_rsqrt(jatom,2)))

				write(*,*) theta11,((x(1)*uh(iatom,1,1) + x(2)*uh(iatom,1,2) + x(3)*uh(iatom,1,3)) / (rsqrt * patch_rsqrt(iatom,1)))
				write(*,*) theta12,((x(1)*uh(iatom,2,1) + x(2)*uh(iatom,2,2) + x(3)*uh(iatom,2,3)) / (rsqrt * patch_rsqrt(iatom,2)))
				write(*,*) theta21,((-x(1)*uh(jatom,1,1) - x(2)*uh(jatom,1,2) - x(3)*uh(jatom,1,3)) / (rsqrt * patch_rsqrt(jatom,1)))
				write(*,*) theta22,((-x(1)*uh(jatom,2,1) - x(2)*uh(jatom,2,2) - x(3)*uh(jatom,2,3)) / (rsqrt * patch_rsqrt(jatom,2)))

				if(rsqrt_old.gt.boxl*0.5d0) then
					write(*,*) iatom, jatom, theta11*57.2958, theta12*57.2958, theta21*57.2958, theta22*57.2958
					write(*,*) iatom, jatom, rsqrt_old,rsqrt, patch_rsqrt(iatom,:),patch_rsqrt(jatom,:)
						theta11 = 3.14159-theta11; theta12 = 3.14159-theta12; theta21 = 3.14159-theta21; theta22 = 3.14159-theta22
				endif

				write(*,*) iatom, jatom, theta11*57.2958, theta12*57.2958, theta21*57.2958, theta22*57.2958

			endif
			enddo !jatom
		enddo !ipatch

!	WRITE THE INITIAL CONFIGURATION FILE

	write(2,*) natom*3+ncrowd
	write(2,*)

	do iatom = 1,natom
		write(2,*) 'C',ap(iatom,1),ap(iatom,2),ap(iatom,3)
		write(2,*) 'P',ah(iatom,1,1),ah(iatom,1,2),ah(iatom,1,3)
		write(2,*) 'P',ah(iatom,2,1),ah(iatom,2,2),ah(iatom,2,3)
	enddo !i
	do icrowd = 1,ncrowd
		write(2,*) 'Ar',ac(icrowd,1),ac(icrowd,2),ac(icrowd,3)
	enddo !i



							! ---------------- MC simulation ------------------!


	accp_crowd = 0; rej_crowd = 0
	accp_patch = 0; rej_patch = 0

	do imc = 1, nmc   ! nmc is the no of simulation
		write(*,*) imc
		!write(2,*) '		6'
		!write(2,*)

		do iatom = 1,natom
			patch_rsqrt(ipatch,:) = 0
			Ei = 0; Ef = 0


!	RANDOM MOTION

			two_pie = 2.0d0 * 3.141592654d0

		!	rotational
	
			rand=0.0d0
			call random_number(rand)
			phi 	   = rand * two_pie
			cos_phi = dcos( phi )
			sin_phi = dsin( phi )

			rand=0.0d0
			call random_number(rand)
			cos_theta  = 2.0d0 * rand - 1.0d0
			sin_theta  = sqrt(1.0d0 - cos_theta * cos_theta)

			rand=0.0d0
			call random_number(rand)
			psi 	   = (2.0d0 * rand - 1.0d0) * dpsi_max
			cos_psi    = dcos( psi )
			sin_psi    = dsin( psi ) 
			onemcospsi = 1 - cos_psi

			nx = sin_theta * cos_phi
			ny = sin_theta * sin_phi
			nz = cos_theta

			nxy  = nx * ny; nyz  = ny * nz; nxz  = nx * nz   
			nxsq = nx * nx; nysq = ny * ny; nzsq = nz * nz

!	Rotational Matrix
! 	rot(i,j) = rotational matrix element at ith row & jth column
 
			rot(1,1) = cos_psi + nxsq * onemcospsi
			rot(1,2) = nxy * onemcospsi + nz * sin_psi
			rot(1,3) = nxz * onemcospsi - ny * sin_psi

			rot(2,1) = nxy * onemcospsi - nz * sin_psi
			rot(2,2) = cos_psi + nysq*onemcospsi
			rot(2,3) = nyz * onemcospsi + nx * sin_psi

			rot(3,1) = nxz * onemcospsi + ny * sin_psi
			rot(3,2) = nyz * onemcospsi - nx * sin_psi
			rot(3,3) = cos_psi + nzsq*onemcospsi

			
!	calculation of determinant

			detrot=rot(1,1)*(rot(2,2)*rot(3,3)-rot(2,3)*rot(3,2))
			detrot=detrot+rot(1,2)*(rot(3,1)*rot(2,3)-rot(2,1)*rot(3,3))
			detrot=detrot+rot(1,3)*(rot(2,1)*rot(3,2)-rot(3,1)*rot(2,2))

			!write(*,*) imc,detrot

			if(detrot.lt.0.9999.or.detrot.gt.1.0001) then
				write(*,*) 'blah!!','determinant is not 1',detrot
				stop
			endif

!	 MOVE iTH PARTICLE RANDOMLY   !

30			do j=1,3
				call random_number(yt)
				apn(iatom,j)=ap(iatom,j)!+(2*yt-1)*d_m  ! move ith particle randomly
				if(apn(iatom,j).gt.boxl) apn(iatom,j)=apn(iatom,j)-boxl*1.0d0
				if(apn(iatom,j).lt.0.0d0) apn(iatom,j)=apn(iatom,j)+boxl*1.0d0
			enddo !j

			a_h(iatom,1,:) = ah(iatom,1,:) +  (apn(iatom,:) - ap(iatom,:))
			a_h(iatom,2,:) = ah(iatom,2,:) +  (apn(iatom,:) - ap(iatom,:))
			
!	 MULTIPLY WITH ROTATIONAL MATRIX	!

			ahn(iatom,1,1)=rot(1,1)*(a_h(iatom,1,1)-apn(iatom,1))+rot(1,2)*(a_h(iatom,1,2)-apn(iatom,2))&
			+rot(1,3)*(a_h(iatom,1,3)-apn(iatom,3))
			ahn(iatom,1,2)=rot(2,1)*(a_h(iatom,1,1)-apn(iatom,1))+rot(2,2)*(a_h(iatom,1,2)-apn(iatom,2))&
			+rot(2,3)*(a_h(iatom,1,3)-apn(iatom,3))
			ahn(iatom,1,3)=rot(3,1)*(a_h(iatom,1,1)-apn(iatom,1))+rot(3,2)*(a_h(iatom,1,2)-apn(iatom,2))&
			+rot(3,3)*(a_h(iatom,1,3)-apn(iatom,3))

			ahn(iatom,2,1)=rot(1,1)*(a_h(iatom,2,1)-apn(iatom,1))+rot(1,2)*(a_h(iatom,2,2)-apn(iatom,2))&
			+rot(1,3)*(a_h(iatom,2,3)-apn(iatom,3))
			ahn(iatom,2,2)=rot(2,1)*(a_h(iatom,2,1)-apn(iatom,1))+rot(2,2)*(a_h(iatom,2,2)-apn(iatom,2))&
			+rot(2,3)*(a_h(iatom,2,3)-apn(iatom,3))
			ahn(iatom,2,3)=rot(3,1)*(a_h(iatom,2,1)-apn(iatom,1))+rot(3,2)*(a_h(iatom,2,2)-apn(iatom,2))&
			+rot(3,3)*(a_h(iatom,2,3)-apn(iatom,3))
 				
			ahn(iatom,1,:) = ahn(iatom,1,:) + apn(iatom,:); ahn(iatom,2,:) = ahn(iatom,2,:) + apn(iatom,:)

			!if(ahn(iatom,1,1).ge.boxl) goto 30
			!if(ahn(iatom,1,1).lt.0) goto 30
			!if(ahn(iatom,1,2).ge.boxl) goto 30
			!if(ahn(iatom,1,2).lt.0) goto 30
			!if(ahn(iatom,1,3).ge.boxl) goto 30
			!if(ahn(iatom,1,3).lt.0) goto 30

			!if(ahn(iatom,2,1).ge.boxl) goto 30
			!if(ahn(iatom,2,1).lt.0) goto 30 	
			!if(ahn(iatom,2,2).ge.boxl) goto 30
			!if(ahn(iatom,2,2).lt.0) goto 30 	
			!if(ahn(iatom,2,3).ge.boxl) goto 30
			!if(ahn(iatom,2,3).lt.0) goto 30 	

			
			write(2,*) 'C',apn(iatom,1),apn(iatom,2),apn(iatom,3)
			write(2,*) 'P',ahn(iatom,1,1),ahn(iatom,1,2),ahn(iatom,1,3)
			write(2,*) 'P',ahn(iatom,2,1),ahn(iatom,2,2),ahn(iatom,2,3)


!	CALCULATE ANGLE BETWEEN PATCHES OF DIFFERENT PARTICLES

	
			patch_rsqrt(iatom,:) = 0
			uh(iatom,:,:) = 0

			uh(iatom,1,:) = ahn(iatom,1,:)-apn(iatom,:)
			uh(iatom,1,:) = uh(iatom,1,:)-boxl*anint(uh(iatom,1,:)/boxl)
			patch_rsqrt(iatom,1) = sqrt(uh(iatom,1,1)*uh(iatom,1,1)+uh(iatom,1,2)*uh(iatom,1,2)+uh(iatom,1,3)*uh(iatom,1,3)) 
			
			uh(iatom,2,:) = ahn(iatom,2,:)-apn(iatom,:)
			uh(iatom,2,:) = uh(iatom,2,:)-boxl*anint(uh(iatom,2,:)/boxl)
			patch_rsqrt(iatom,2) = sqrt(uh(iatom,2,1)*uh(iatom,2,1)+uh(iatom,2,2)*uh(iatom,2,2)+uh(iatom,3,3)*uh(iatom,3,3))

			uh(iatom,1,:) = uh(iatom,1,:)/patch_rsqrt(iatom,1)
			patch_rsqrt(iatom,1) = sqrt(uh(iatom,1,1)*uh(iatom,1,1)+uh(iatom,1,2)*uh(iatom,1,2)+uh(iatom,1,3)*uh(iatom,1,3)) 

			uh(iatom,2,:) = uh(iatom,2,:)/patch_rsqrt(iatom,2)
			patch_rsqrt(iatom,2) = sqrt(uh(iatom,2,1)*uh(iatom,2,1)+uh(iatom,2,2)*uh(iatom,2,2)+uh(iatom,2,3)*uh(iatom,2,3)) 
			

			do jatom = 1,natom
				if (iatom.ne.jatom) then
					x(:)=ap(jatom,:)-apn(iatom,:) !x() is the distance between ith and jth particle
					rsqrt_old = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
					!x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention
					rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
					x(:)=x(:)/rsqrt
					rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))

					uh(jatom,1,:) = ah(jatom,1,:)-ap(jatom,:)
					uh(jatom,1,:) = uh(jatom,1,:)-boxl*anint(uh(jatom,1,:)/boxl)
			
					uh(jatom,2,:) = ah(jatom,2,:)-ap(jatom,:)
					uh(jatom,2,:) = uh(jatom,2,:)-boxl*anint(uh(jatom,2,:)/boxl)

					patch_rsqrt(jatom,1) = sqrt(uh(jatom,1,1)*uh(jatom,1,1)+uh(jatom,1,2)*uh(jatom,1,2)+uh(jatom,1,3)*uh(jatom,1,3)) 
					patch_rsqrt(jatom,2) = sqrt(uh(jatom,2,1)*uh(jatom,2,1)+uh(jatom,2,2)*uh(jatom,2,2)+uh(jatom,3,3)*uh(jatom,3,3))

			uh(jatom,1,:) = uh(jatom,1,:)/patch_rsqrt(jatom,1)
			patch_rsqrt(jatom,1) = sqrt(uh(jatom,1,1)*uh(jatom,1,1)+uh(jatom,1,2)*uh(jatom,1,2)+uh(jatom,1,3)*uh(jatom,1,3)) 

			uh(jatom,2,:) = uh(jatom,2,:)/patch_rsqrt(jatom,2)
			patch_rsqrt(jatom,2) = sqrt(uh(jatom,2,1)*uh(jatom,2,1)+uh(jatom,2,2)*uh(jatom,2,2)+uh(jatom,2,3)*uh(jatom,2,3)) 
			

					theta11 = dacos ( (x(1)*uh(iatom,1,1) + x(2)*uh(iatom,1,2) + x(3)*uh(iatom,1,3)) / (rsqrt * patch_rsqrt(iatom,1)))
					theta12 = dacos ( (x(1)*uh(iatom,2,1) + x(2)*uh(iatom,2,2) + x(3)*uh(iatom,2,3)) / (rsqrt * patch_rsqrt(iatom,2)))
					theta21 = dacos ( (-x(1)*uh(jatom,1,1) - x(2)*uh(jatom,1,2) - x(3)*uh(jatom,1,3)) / (rsqrt * patch_rsqrt(jatom,1)))
					theta22 = dacos ( (-x(1)*uh(jatom,2,1) - x(2)*uh(jatom,2,2) - x(3)*uh(jatom,2,3)) / (rsqrt * patch_rsqrt(jatom,2)))

					if(rsqrt_old.gt.boxl*0.5d0) then
						theta11 = 3.14159-theta11; theta12 = 3.14159-theta12; theta21 = 3.14159-theta21; theta22 = 3.14159-theta22
					endif

					write(*,*) iatom, jatom, theta11*57.2958, theta12*57.2958, theta21*57.2958, theta22*57.2958
	
				endif
			enddo !jatom
		

		enddo !iatom
		

	!---- write no of nth cycle,Energy at nth cycle,Energy of each particle after nth cycle  in 'energy.dat' file to plot---!

		write(11,*) imc,accp_patch/(accp_patch+rej_patch),accp_crowd/(accp_crowd+rej_crowd), d_m, dpsi_max


	enddo !imc

	end program mc
