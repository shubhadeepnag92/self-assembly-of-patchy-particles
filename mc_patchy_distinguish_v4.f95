!	Monte-Carlo PROGRAM TO SIMULATE patchy spheres with crowding agents
!	modified date - 18/12/2022
!	need to check the update position of patch positions when it goes beyond the boundary
!


	program mc

	implicit none

	real,allocatable :: uh(:,:,:),ah(:,:,:),ap(:,:),a_p(:,:),ac(:,:),a_c(:,:),x(:),patch_rsqrt(:,:)
	real,allocatable :: uhn(:,:,:),ahn(:,:,:),a_h(:,:,:),rot(:,:),apn(:,:)
	real		   :: boxl,yt,OEp,OEc,OEpc,rsqrt,rhex,rtwelve, R, T, RT, p_r, theta0
	real		   :: Ei, Ef, dE, d_m, theta11, theta12, theta21, theta22, two_pie, sin_theta, sin_psi, sin_phi, rand, psi, phi, onemcospsi
	real		   :: ny,nxz,nxy,nxsq,nx,nzsq, nz, nyz, nysq, dpsi_max, detrot, cos_theta, cos_psi, cos_phi,rsqrt_old
	real		   :: eps11, eps12, eps22, sgmahex11, sgmahex12, sgmahex22, sgma, rsqr, sgmatwelve
	real		   :: accp_patch,rej_patch,accp_crowd,rej_crowd,E_final, E_ini
	integer	   :: j, iatom, icrowd, natom, ncrowd, m, jatom, jcrowd, imc, nmc, b, c, ipatch, npatch, jpatch

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
	read(1,*) sgma, sgmatwelve

!	ALLOCATE DATA
	allocate ( ah(1:natom,1:2,1:3),ap(1:natom,1:3),a_p(1:natom,1:3),ac(1:ncrowd,1:3),a_c(1:ncrowd,1:3),x(1:3),patch_rsqrt(1:natom,1:2))
	allocate ( uh(1:natom,1:2,1:3),uhn(1:natom,1:2,1:3),ahn(1:natom,1:2,1:3),a_h(1:natom,1:2,1:3),rot(1:3,1:3),apn(1:natom,1:3) )
	RT = R * T; npatch = natom; theta0 = 0.26
    
! 	WRITE IMPORTANT PARAMTERS OF THESE MC SIMULATION

	write(10,*) 'PARAMTERS OF MC SIMULATION'
	write(10,*) 'boxlength	No of MC Steps	No. of Particles	No. of Crowding agents	d_m	R	T'
	write(10,*) boxl, nmc, natom, ncrowd, d_m, R, T

! 	GENERATE THE INITIAL CONFIGURATION FILE

!	GENERATE THE ATOM POSITIONS RANDOMLY

	do iatom = 1, natom
		call random_number(yt); ap(iatom,1) = yt*boxl
		call random_number(yt);	ap(iatom,2) = yt*boxl
		call random_number(yt);	ap(iatom,3) = yt*boxl
	enddo !i

!	GENERATE THE PATCH POSITIONS 

	do iatom = 1,natom
		ah(iatom,1,1) = ap(iatom,1) + sgma*0.5; ah(iatom,1,2) = ap(iatom,2) ; ah(iatom,1,3) = ap(iatom,3)
		ah(iatom,2,1) = ap(iatom,1) + sgma*0.5*cos(2.0944) ; 
		ah(iatom,2,2) = ap(iatom,2) + sgma*0.5*sin(2.0944) ; ah(iatom,2,3) = ap(iatom,3)
	enddo !iatom	

!	GENERATE THE CROWD POSITIONS RANDOMLY

	do icrowd = 1, ncrowd
		call random_number(yt); ac(icrowd,1) = yt*boxl
		call random_number(yt);	ac(icrowd,2) = yt*boxl
		call random_number(yt);	ac(icrowd,3) = yt*boxl
	enddo !i

!	CALCULATE ANGLE BETWEEN PATCHES OF DIFFERENT PARTICLES
	E_ini = 0
	do iatom = 1,natom
		patch_rsqrt(iatom,:) = 0
		uh(iatom,:,:) = 0

		uh(iatom,1,:) = ah(iatom,1,:)-ap(iatom,:)
		uh(iatom,1,:) = uh(iatom,1,:)-boxl*anint(uh(iatom,1,:)/boxl)
		patch_rsqrt(iatom,1) = sqrt(uh(iatom,1,1)*uh(iatom,1,1)+uh(iatom,1,2)*uh(iatom,1,2)+uh(iatom,1,3)*uh(iatom,1,3)) 
		uh(iatom,1,:) = uh(iatom,1,:)/patch_rsqrt(iatom,1)
			
		uh(iatom,2,:) = ah(iatom,2,:)-ap(iatom,:)
		uh(iatom,2,:) = uh(iatom,2,:)-boxl*anint(uh(iatom,2,:)/boxl)
		patch_rsqrt(iatom,2) = sqrt(uh(iatom,2,1)*uh(iatom,2,1)+uh(iatom,2,2)*uh(iatom,2,2)+uh(iatom,3,3)*uh(iatom,3,3))
		uh(iatom,2,:) = uh(iatom,2,:)/patch_rsqrt(iatom,2)
			
		do jatom = 1,natom
			if (iatom.ne.jatom) then
			x(:)=ap(jatom,:)-ap(iatom,:) !x() is the distance between ith and jth particle
			x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention
			rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
			x(:) = x(:)/rsqrt

			uh(jatom,1,:) = ah(jatom,1,:)-ap(jatom,:)
			uh(jatom,1,:) = uh(jatom,1,:)-boxl*anint(uh(jatom,1,:)/boxl)
			patch_rsqrt(jatom,1) = sqrt(uh(jatom,1,1)*uh(jatom,1,1)+uh(jatom,1,2)*uh(jatom,1,2)+uh(jatom,1,3)*uh(jatom,1,3)) 
			uh(jatom,1,:) = uh(jatom,1,:)/patch_rsqrt(jatom,1)
			
			uh(jatom,2,:) = ah(jatom,2,:)-ap(jatom,:)
			uh(jatom,2,:) = uh(jatom,2,:)-boxl*anint(uh(jatom,2,:)/boxl)
			patch_rsqrt(jatom,2) = sqrt(uh(jatom,2,1)*uh(jatom,2,1)+uh(jatom,2,2)*uh(jatom,2,2)+uh(jatom,3,3)*uh(jatom,3,3))
			uh(jatom,2,:) = uh(jatom,2,:)/patch_rsqrt(jatom,2)

			theta11 = acos ( (x(1)*uh(iatom,1,1) + x(2)*uh(iatom,1,2) + x(3)*uh(iatom,1,3)) )
			theta12 = acos ( (x(1)*uh(iatom,2,1) + x(2)*uh(iatom,2,2) + x(3)*uh(iatom,2,3)) )
			theta21 = acos ( (-x(1)*uh(jatom,1,1) - x(2)*uh(jatom,1,2) - x(3)*uh(jatom,1,3)) )
			theta22 = acos ( (-x(1)*uh(jatom,2,1) - x(2)*uh(jatom,2,2) - x(3)*uh(jatom,2,3)) )

			!write(*,*) iatom, jatom, theta11*57.2958, theta12*57.2958, theta21*57.2958, theta22*57.2958

			if(rsqrt.le.10) then
				rsqr = rsqrt*rsqrt
				rhex = rsqr*rsqr*rsqr
				rtwelve = rhex*rhex
				E_ini =E_ini + eps11*(sgmatwelve/rtwelve)
			endif
			
			if(rsqrt.ge.sgma.and.rsqrt.le.sgma+1.0) then
				write(*,*) 'Eureka'
				if(theta11.le.0.3.and.theta21.le.0.3) E_ini = E_ini - eps11
				if(theta11.le.0.3.and.theta22.le.0.3) E_ini = E_ini - eps11
				if(theta12.le.0.3.and.theta21.le.0.3) E_ini = E_ini - eps11
				if(theta12.le.0.3.and.theta22.le.0.3) E_ini = E_ini - eps11
			endif

			endif

			write(30,*) iatom,jatom,E_ini,rtwelve,rsqrt
		enddo !jatom
	enddo !iatom

	write(30,*) 'Initial Energy',E_ini

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
			if(mod(imc,100).eq.0)  write(2,*) natom*3+ncrowd
			if(mod(imc,100).eq.0)  write(2,*)

		do iatom = 1,natom
			
			Ei = 0; Ef = 0

	!	CALCULATE ANGLE BETWEEN PATCHES OF DIFFERENT PARTICLES
	
			patch_rsqrt(iatom,:) = 0
			uh(iatom,:,:) = 0

			uh(iatom,1,:) = ah(iatom,1,:)-ap(iatom,:)
			uh(iatom,1,:) = uh(iatom,1,:)-boxl*anint(uh(iatom,1,:)/boxl)
			patch_rsqrt(iatom,1) = sqrt(uh(iatom,1,1)*uh(iatom,1,1)+uh(iatom,1,2)*uh(iatom,1,2)+uh(iatom,1,3)*uh(iatom,1,3)) 
			uh(iatom,1,:) = uh(iatom,1,:)/patch_rsqrt(iatom,1)
			
			uh(iatom,2,:) = ah(iatom,2,:)-ap(iatom,:)
			uh(iatom,2,:) = uh(iatom,2,:)-boxl*anint(uh(iatom,2,:)/boxl)
			patch_rsqrt(iatom,2) = sqrt(uh(iatom,2,1)*uh(iatom,2,1)+uh(iatom,2,2)*uh(iatom,2,2)+uh(iatom,3,3)*uh(iatom,3,3))
			uh(iatom,2,:) = uh(iatom,2,:)/patch_rsqrt(iatom,2)

			do jatom = 1,natom
				if (iatom.ne.jatom) then

					uh(jatom,1,:) = ah(jatom,1,:)-ap(jatom,:)
					uh(jatom,1,:) = uh(jatom,1,:)-boxl*anint(uh(jatom,1,:)/boxl)
					patch_rsqrt(jatom,1) = sqrt(uh(jatom,1,1)*uh(jatom,1,1)+uh(jatom,1,2)*uh(jatom,1,2)+uh(jatom,1,3)*uh(jatom,1,3)) 
					uh(jatom,1,:) = uh(jatom,1,:)/patch_rsqrt(jatom,1)
			
					uh(jatom,2,:) = ah(jatom,2,:)-ap(jatom,:)
					uh(jatom,2,:) = uh(jatom,2,:)-boxl*anint(uh(jatom,2,:)/boxl)
					patch_rsqrt(jatom,2) = sqrt(uh(jatom,2,1)*uh(jatom,2,1)+uh(jatom,2,2)*uh(jatom,2,2)+uh(jatom,3,3)*uh(jatom,3,3))
					uh(jatom,2,:) = uh(jatom,2,:)/patch_rsqrt(jatom,2)

					x(1)=ap(jatom,1)-ap(iatom,1); x(1)=x(1)-boxl*anint(x(1)/boxl)
					x(2)=ap(jatom,2)-ap(iatom,2); x(2)=x(2)-boxl*anint(x(2)/boxl)
					x(3)=ap(jatom,3)-ap(iatom,3); x(3)=x(3)-boxl*anint(x(3)/boxl)
					rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
					x(1)=x(1)/rsqrt; x(2)=x(2)/rsqrt; x(3)=x(3)/rsqrt
			
					theta11 =  acos ( (x(1)*uh(iatom,1,1) + x(2)*uh(iatom,1,2) + x(3)*uh(iatom,1,3)) )
					theta12 =  acos ( (x(1)*uh(iatom,2,1) + x(2)*uh(iatom,2,2) + x(3)*uh(iatom,2,3)) )

					x(1)=ap(iatom,1)-ap(jatom,1); x(1)=x(1)-boxl*anint(x(1)/boxl)
					x(2)=ap(iatom,2)-ap(jatom,2); x(2)=x(2)-boxl*anint(x(2)/boxl)
					x(3)=ap(iatom,3)-ap(jatom,3); x(3)=x(3)-boxl*anint(x(3)/boxl)
					rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
					x(1)=x(1)/rsqrt; x(2)=x(2)/rsqrt; x(3)=x(3)/rsqrt

					theta21 =  acos ( (x(1)*uh(jatom,1,1) + x(2)*uh(jatom,1,2) + x(3)*uh(jatom,1,3)) )
					theta22 =  acos ( (x(1)*uh(jatom,2,1) + x(2)*uh(jatom,2,2) + x(3)*uh(jatom,2,3)) )

					!write(32,*) 'before',iatom, jatom, theta11*57.2958, theta12*57.2958, theta21*57.2958, theta22*57.2958

					rsqr = rsqrt*rsqrt
					rhex = rsqr*rsqr*rsqr
					rtwelve = rhex*rhex

					if(rsqrt.gt.10) Ei = Ei 
					if(rsqrt.gt.sgma+0.895.and.rsqrt.le.10) Ei = Ei + eps11*(sgmatwelve/rtwelve)
					if(rsqrt.gt.sgma.and.rsqrt.le.sgma+0.895) Ei = Ei
					if(rsqrt.le.sgma) Ei = Ei + eps11*(sgmatwelve/rtwelve)

					!if(rsqrt.le.10) then
					!	rsqr = rsqrt*rsqrt
					!	rhex = rsqr*rsqr*rsqr
					!	rtwelve = rhex*rhex
					!	Ei =Ei + eps11*(sgmatwelve/rtwelve)
					!endif
			
					if(rsqrt.ge.sgma.and.rsqrt.le.sgma+0.895) then
				write(32,*) 'ini',imc,iatom,jatom,rsqrt,theta11*57.29,theta12*57.29,theta21*57.29,theta22*57.29, theta0,Ei
				if(theta11.le.theta0.and.theta21.le.theta0) Ei = Ei - eps11*5!0.5*((theta11-theta0)**2+(theta21-theta0)**2)-5
				if(theta11.le.theta0.and.theta22.le.theta0) Ei = Ei - eps11*5!0.5*((theta11-theta0)**2+(theta22-theta0)**2)-5
				if(theta12.le.theta0.and.theta21.le.theta0) Ei = Ei - eps11*5!0.5*((theta12-theta0)**2+(theta21-theta0)**2)-5
				if(theta12.le.theta0.and.theta22.le.theta0) Ei = Ei - eps11*5!0.5*((theta12-theta0)**2+(theta22-theta0)**2)-5
				write(32,*) imc,iatom,jatom,rsqrt,'Eureka',theta11*57.29,theta12*57.29,theta21*57.29,theta22*57.29,theta0,Ei
					endif

					!write(*,*) rsqrt,Ef
	
				endif
			enddo !jatom

			
!	RANDOM MOTION

30			two_pie = 2.0d0 * 3.141592654d0

		!	rotational
	
			rand=0.0d0
			call random_number(rand)
			phi 	   = rand * two_pie
			cos_phi = cos( phi )
			sin_phi = sin( phi )

			rand=0.0d0
			call random_number(rand)
			cos_theta  = 2.0d0 * rand - 1.0d0
			sin_theta  = sqrt(1.0d0 - cos_theta * cos_theta)

			rand=0.0d0
			call random_number(rand)
			psi 	   = (2.0d0 * rand - 1.0d0) * dpsi_max
			cos_psi    = cos( psi )
			sin_psi    = sin( psi ) 
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

			!write(*,*) detrot

			if(detrot.lt.0.9999.or.detrot.gt.1.0001) then
				write(*,*) 'blah!!','determinant is not 1',detrot
				stop
			endif

			!	 MOVE iTH PARTICLE RANDOMLY   !

			do j=1,3
				call random_number(yt)
				apn(iatom,j)=ap(iatom,j)+(2*yt-1)*d_m  ! move ith particle randomly
				do while(apn(iatom,j).gt.boxl)
				call random_number(yt)
				apn(iatom,j)=ap(iatom,j)+(2*yt-1)*d_m
				end do!while
				do while(apn(iatom,j).lt.0.0)
				call random_number(yt)
				apn(iatom,j)=ap(iatom,j)+(2*yt-1)*d_m
				end do!while
				!if(apn(iatom,j).gt.boxl) apn(iatom,j)=apn(iatom,j)-boxl*1.0d0
				!if(apn(iatom,j).lt.0.0d0) apn(iatom,j)=apn(iatom,j)+boxl*1.0d0
			enddo !j

			a_h(iatom,1,1) = ah(iatom,1,1) + (apn(iatom,1)-ap(iatom,1))
			a_h(iatom,1,2) = ah(iatom,1,2) + (apn(iatom,2)-ap(iatom,2))
			a_h(iatom,1,3) = ah(iatom,1,3) + (apn(iatom,3)-ap(iatom,3))
			!write(34,*) imc,iatom,sqrt((a_h(iatom,1,1)-apn(iatom,1))*(a_h(iatom,1,1)-apn(iatom,1))&
			!+(a_h(iatom,1,2)-apn(iatom,2))*(a_h(iatom,1,2)-apn(iatom,2))+(a_h(iatom,1,3)-apn(iatom,3))*(a_h(iatom,1,3)-apn(iatom,3)))
		
			ahn(iatom,1,1)=rot(1,1)*(a_h(iatom,1,1)-apn(iatom,1))+rot(1,2)*(a_h(iatom,1,2)-apn(iatom,2))&
			+rot(1,3)*(a_h(iatom,1,3)-apn(iatom,3))
			ahn(iatom,1,2)=rot(2,1)*(a_h(iatom,1,1)-apn(iatom,1))+rot(2,2)*(a_h(iatom,1,2)-apn(iatom,2))&
			+rot(2,3)*(a_h(iatom,1,3)-apn(iatom,3))
			ahn(iatom,1,3)=rot(3,1)*(a_h(iatom,1,1)-apn(iatom,1))+rot(3,2)*(a_h(iatom,1,2)-apn(iatom,2))&
			+rot(3,3)*(a_h(iatom,1,3)-apn(iatom,3))
			ahn(iatom,1,1)=ahn(iatom,1,1)+apn(iatom,1)
			ahn(iatom,1,2)=ahn(iatom,1,2)+apn(iatom,2)
			ahn(iatom,1,3)=ahn(iatom,1,3)+apn(iatom,3)
			!write(34,*) imc,iatom,sqrt((ahn(iatom,1,1)-apn(iatom,1))*(ahn(iatom,1,1)-apn(iatom,1))&
			!+(ahn(iatom,1,2)-apn(iatom,2))*(ahn(iatom,1,2)-apn(iatom,2))+(ahn(iatom,1,3)-apn(iatom,3))*(ahn(iatom,1,3)-apn(iatom,3)))

			a_h(iatom,2,1) = ah(iatom,2,1) + (apn(iatom,1)-ap(iatom,1))
			a_h(iatom,2,2) = ah(iatom,2,2) + (apn(iatom,2)-ap(iatom,2))
			a_h(iatom,2,3) = ah(iatom,2,3) + (apn(iatom,3)-ap(iatom,3))
			!write(34,*) imc,iatom,sqrt((a_h(iatom,2,1)-apn(iatom,1))*(a_h(iatom,2,1)-apn(iatom,1))&
			!+(a_h(iatom,2,2)-apn(iatom,2))*(a_h(iatom,2,2)-apn(iatom,2))+(a_h(iatom,2,3)-apn(iatom,3))*(a_h(iatom,2,3)-apn(iatom,3)))
		
			ahn(iatom,2,1)=rot(1,1)*(a_h(iatom,2,1)-apn(iatom,1))+rot(1,2)*(a_h(iatom,2,2)-apn(iatom,2))&
			+rot(1,3)*(a_h(iatom,2,3)-apn(iatom,3))
			ahn(iatom,2,2)=rot(2,1)*(a_h(iatom,2,1)-apn(iatom,1))+rot(2,2)*(a_h(iatom,2,2)-apn(iatom,2))&
			+rot(2,3)*(a_h(iatom,2,3)-apn(iatom,3))
			ahn(iatom,2,3)=rot(3,1)*(a_h(iatom,2,1)-apn(iatom,1))+rot(3,2)*(a_h(iatom,2,2)-apn(iatom,2))&
			+rot(3,3)*(a_h(iatom,2,3)-apn(iatom,3))
			ahn(iatom,2,1)=ahn(iatom,2,1)+apn(iatom,1)
			ahn(iatom,2,2)=ahn(iatom,2,2)+apn(iatom,2)
			ahn(iatom,2,3)=ahn(iatom,2,3)+apn(iatom,3)
			!write(34,*) imc,iatom,sqrt((ahn(iatom,2,1)-apn(iatom,1))*(ahn(iatom,2,1)-apn(iatom,1))&
			!+(ahn(iatom,2,2)-apn(iatom,2))*(ahn(iatom,2,2)-apn(iatom,2))+(ahn(iatom,2,3)-apn(iatom,3))*(ahn(iatom,2,3)-apn(iatom,3)))

			
			if(mod(imc,100).eq.0) then
				write(2,*) 'C',apn(iatom,1),apn(iatom,2),apn(iatom,3)
				write(2,*) 'P',ahn(iatom,1,1),ahn(iatom,1,2),ahn(iatom,1,3)
				write(2,*) 'P',ahn(iatom,2,1),ahn(iatom,2,2),ahn(iatom,2,3)
			endif

!	CALCULATE ANGLE BETWEEN PATCHES OF DIFFERENT PARTICLES
			patch_rsqrt(iatom,:) = 0
			uh(iatom,:,:) = 0

			uh(iatom,1,:) = ahn(iatom,1,:)-apn(iatom,:)
			uh(iatom,1,:) = uh(iatom,1,:)-boxl*anint(uh(iatom,1,:)/boxl)
			patch_rsqrt(iatom,1) = sqrt(uh(iatom,1,1)*uh(iatom,1,1)+uh(iatom,1,2)*uh(iatom,1,2)+uh(iatom,1,3)*uh(iatom,1,3)) 
			uh(iatom,1,:) = uh(iatom,1,:)/patch_rsqrt(iatom,1)
			
			uh(iatom,2,:) = ahn(iatom,2,:)-apn(iatom,:)
			uh(iatom,2,:) = uh(iatom,2,:)-boxl*anint(uh(iatom,2,:)/boxl)
			patch_rsqrt(iatom,2) = sqrt(uh(iatom,2,1)*uh(iatom,2,1)+uh(iatom,2,2)*uh(iatom,2,2)+uh(iatom,2,3)*uh(iatom,2,3))
			uh(iatom,2,:) = uh(iatom,2,:)/patch_rsqrt(iatom,2)

			do jatom = 1,natom
				if (iatom.ne.jatom) then

					uh(jatom,1,1) = ah(jatom,1,1)-ap(jatom,1); uh(jatom,1,1) = uh(jatom,1,1)-boxl*anint(uh(jatom,1,1)/boxl)
					uh(jatom,1,2) = ah(jatom,1,2)-ap(jatom,2); uh(jatom,1,2) = uh(jatom,1,2)-boxl*anint(uh(jatom,1,2)/boxl)
					uh(jatom,1,3) = ah(jatom,1,3)-ap(jatom,3); uh(jatom,1,3) = uh(jatom,1,3)-boxl*anint(uh(jatom,1,3)/boxl)
					patch_rsqrt(jatom,1) = sqrt(uh(jatom,1,1)*uh(jatom,1,1)+uh(jatom,1,2)*uh(jatom,1,2)+uh(jatom,1,3)*uh(jatom,1,3)) 
					uh(jatom,1,1) = uh(jatom,1,1)/patch_rsqrt(jatom,1)
					uh(jatom,1,2) = uh(jatom,1,2)/patch_rsqrt(jatom,1)
					uh(jatom,1,3) = uh(jatom,1,3)/patch_rsqrt(jatom,1)
					!write(33,*) imc,iatom,jatom, uh(jatom,1,:), uh(jatom,1,:)* uh(jatom,1,:)

					uh(jatom,2,1) = ah(jatom,2,1)-ap(jatom,1); uh(jatom,2,1) = uh(jatom,2,1)-boxl*anint(uh(jatom,2,1)/boxl)
					uh(jatom,2,2) = ah(jatom,2,2)-ap(jatom,2); uh(jatom,2,2) = uh(jatom,2,2)-boxl*anint(uh(jatom,2,2)/boxl)
					uh(jatom,2,3) = ah(jatom,2,3)-ap(jatom,3); uh(jatom,2,3) = uh(jatom,2,3)-boxl*anint(uh(jatom,2,3)/boxl)
					patch_rsqrt(jatom,2) = sqrt(uh(jatom,2,1)*uh(jatom,2,1)+uh(jatom,2,2)*uh(jatom,2,2)+uh(jatom,2,3)*uh(jatom,2,3))
					uh(jatom,2,1) = uh(jatom,2,1)/patch_rsqrt(jatom,2)
					uh(jatom,2,2) = uh(jatom,2,2)/patch_rsqrt(jatom,2)
					uh(jatom,2,3) = uh(jatom,2,3)/patch_rsqrt(jatom,2)
					!write(33,*) imc,iatom,jatom,uh(jatom,2,:), uh(jatom,2,:)* uh(jatom,2,:)

					x(1)=ap(jatom,1)-apn(iatom,1); x(1)=x(1)-boxl*anint(x(1)/boxl)
					x(2)=ap(jatom,2)-apn(iatom,2); x(2)=x(2)-boxl*anint(x(2)/boxl)
					x(3)=ap(jatom,3)-apn(iatom,3); x(3)=x(3)-boxl*anint(x(3)/boxl)
					rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
					x(1)=x(1)/rsqrt; x(2)=x(2)/rsqrt; x(3)=x(3)/rsqrt

					theta11 = (x(1)*uh(iatom,1,1) + x(2)*uh(iatom,1,2) + x(3)*uh(iatom,1,3))
					theta12 = (x(1)*uh(iatom,2,1) + x(2)*uh(iatom,2,2) + x(3)*uh(iatom,2,3))

					x(1)=ap(jatom,1)-apn(iatom,1); x(1)=x(1)-boxl*anint(x(1)/boxl)
					x(2)=ap(jatom,2)-apn(iatom,2); x(2)=x(2)-boxl*anint(x(2)/boxl)
					x(3)=ap(jatom,3)-apn(iatom,3); x(3)=x(3)-boxl*anint(x(3)/boxl)
					rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
					x(1)=x(1)/rsqrt; x(2)=x(2)/rsqrt; x(3)=x(3)/rsqrt

					theta21 = (x(1)*uh(jatom,1,1) + x(2)*uh(jatom,1,2) + x(3)*uh(jatom,1,3))
					theta22 = (x(1)*uh(jatom,2,1) + x(2)*uh(jatom,2,2) + x(3)*uh(jatom,2,3))
	
					if(abs(theta11).gt.1.0.or.abs(theta12).gt.1.0.or.abs(theta21).gt.1.0.or.abs(theta22).gt.1.0) then
						write(*,*) imc,theta11,theta12,theta21,theta22
						goto 30
					endif

					theta11 =  acos (theta11)! (x(1)*uh(iatom,1,1) + x(2)*uh(iatom,1,2) + x(3)*uh(iatom,1,3)) )
					theta12 =  acos (theta12)! (x(1)*uh(iatom,2,1) + x(2)*uh(iatom,2,2) + x(3)*uh(iatom,2,3)) )
					theta21 =  acos (theta21)! (-x(1)*uh(jatom,1,1) - x(2)*uh(jatom,1,2) - x(3)*uh(jatom,1,3)) )
					theta22 =  acos (theta22)! (-x(1)*uh(jatom,2,1) - x(2)*uh(jatom,2,2) - x(3)*uh(jatom,2,3)) )

					!write(32,*) imc,'after', iatom, jatom
	!write(32,*) theta11*57.2958,(x(1)*uh(iatom,1,1) + x(2)*uh(iatom,1,2) + x(3)*uh(iatom,1,3)),x(:),uh(iatom,1,:)
	!write(32,*) theta12*57.2958,(x(1)*uh(iatom,2,1) + x(2)*uh(iatom,2,2) + x(3)*uh(iatom,2,3)),x(:),uh(iatom,2,:)
	!write(32,*) theta21*57.2958,(-x(1)*uh(jatom,1,1) - x(2)*uh(jatom,1,2) - x(3)*uh(jatom,1,3)),x(:),uh(jatom,1,:)
	!write(32,*) theta22*57.2958,(-x(1)*uh(jatom,2,1) - x(2)*uh(jatom,2,2) - x(3)*uh(jatom,2,3)),x(:),uh(jatom,2,:)

					rsqr = rsqrt*rsqrt
					rhex = rsqr*rsqr*rsqr
					rtwelve = rhex*rhex
					
					if(rsqrt.gt.10) Ef = Ef 
					if(rsqrt.gt.sgma+0.895.and.rsqrt.le.10) Ef = Ef + eps11*(sgmatwelve/rtwelve)
					if(rsqrt.gt.sgma.and.rsqrt.le.sgma+0.895) Ef = Ef 
					if(rsqrt.le.sgma) Ef = Ef + eps11*(sgmatwelve/rtwelve)

					!if(rsqrt.le.10) then
					!	rsqr = rsqrt*rsqrt
					!	rhex = rsqr*rsqr*rsqr
					!	rtwelve = rhex*rhex
					!	Ef =Ef + eps11*(sgmatwelve/rtwelve)
					!endif

					if(rsqrt.gt.sgma.and.rsqrt.le.sgma+0.895) then
					write(32,*) 'final',imc,iatom,jatom,rsqrt,theta11*57.29,theta12*57.29,theta21*57.29,theta22*57.29, theta0,Ef
						if(theta11.le.theta0.and.theta21.le.theta0) then
							Ef = Ef - eps11*5!0.5*((theta11-theta0)**2+(theta21-theta0)**2)-5.0!; write(35,*) imc,'Ef',Ef
						endif
						if(theta11.le.theta0.and.theta22.le.theta0) then
							Ef = Ef - eps11*5!0.5*((theta11-theta0)**2+(theta22-theta0)**2)-5.0!; write(35,*) imc,'Ef',Ef
						endif
						if(theta12.le.theta0.and.theta21.le.theta0) then
							 Ef = Ef - eps11*5!0.5*((theta12-theta0)**2+(theta21-theta0)**2)-5.0!; write(35,*) imc,'Ef',Ef
						endif
						if(theta12.le.theta0.and.theta22.le.theta0) then
							 Ef = Ef - eps11*5!0.5*((theta12-theta0)**2+(theta22-theta0)**2)-5.0!; write(35,*) imc,'Ef',Ef
						endif

					write(32,*) imc,iatom,jatom,rsqrt,'Eureka',theta11*57.29,theta12*57.29,theta21*57.29,theta22*57.29,theta0,Ef

					endif

					!write(*,*) rsqrt,Ef
	
				endif
			enddo !jatom

			dE=Ef-Ei   ! dE is the difference between new and initial energy
			if (dE.gt.0) then
				call random_number(yt)
				p_r=exp(-dE/(RT)) !Kb is the Boltzmann factor and T is the temperature
			endif
			if (dE.lt.0.or.yt.lt.p_r) then
				ap(iatom,1) = apn(iatom,1); ap(iatom,2) = apn(iatom,2); ap(iatom,3) = apn(iatom,3)
				ah(iatom,1,1) = ahn(iatom,1,1); ah(iatom,1,2) = ahn(iatom,1,2); ah(iatom,1,3) = ahn(iatom,1,3)
				ah(iatom,2,1) = ahn(iatom,2,1); ah(iatom,2,2) = ahn(iatom,2,2); ah(iatom,2,3) = ahn(iatom,2,3)
				accp_patch=accp_patch+1
			else
				rej_patch=rej_patch+1
			endif
		write(31,*) imc,iatom,Ei,Ef,dE,yt,p_r,accp_patch

		enddo !iatom
		

	!---- write no of nth cycle,Energy at nth cycle,Energy of each particle after nth cycle  in 'energy.dat' file to plot---!

		if(mod(imc,100).eq.0) write(11,*) imc,accp_patch/(accp_patch+rej_patch),accp_crowd/(accp_crowd+rej_crowd), d_m, dpsi_max


! --------- WRITE THE TRAJECTORY AT INTERVAL OF 1000 STEPS ------------!
	if(mod(imc,1).eq.0) then
		write(13,*) natom*3+ncrowd
		write(13,*)
		do iatom=1,natom ! natom is the no of atom in simulation box
			write(13,*) 'Xe',ap(iatom,1),ap(iatom,2),ap(iatom,3)
			write(13,*) 'P',ah(iatom,1,1),ah(iatom,1,2),ah(iatom,1,3)
			write(13,*) 'P',ah(iatom,2,1),ah(iatom,2,2),ah(iatom,2,3)
		enddo !i
		do icrowd=1,ncrowd ! natom is the no of atom in simulation box
			write(13,*) 'Ar',ac(icrowd,1),ac(icrowd,2),ac(icrowd,3)
		enddo !i
	endif

	E_final = 0
	do iatom = 1,natom
		patch_rsqrt(iatom,:) = 0
		uh(iatom,:,:) = 0

			uh(iatom,1,:) = ah(iatom,1,:)-ap(iatom,:)
			uh(iatom,1,:) = uh(iatom,1,:)-boxl*anint(uh(iatom,1,:)/boxl)
			patch_rsqrt(iatom,1) = sqrt(uh(iatom,1,1)*uh(iatom,1,1)+uh(iatom,1,2)*uh(iatom,1,2)+uh(iatom,1,3)*uh(iatom,1,3)) 
			uh(iatom,1,:) = uh(iatom,1,:)/patch_rsqrt(iatom,1)
			
			uh(iatom,2,:) = ah(iatom,2,:)-ap(iatom,:)
			uh(iatom,2,:) = uh(iatom,2,:)-boxl*anint(uh(iatom,2,:)/boxl)
			patch_rsqrt(iatom,2) = sqrt(uh(iatom,2,1)*uh(iatom,2,1)+uh(iatom,2,2)*uh(iatom,2,2)+uh(iatom,3,3)*uh(iatom,3,3))
			uh(iatom,2,:) = uh(iatom,2,:)/patch_rsqrt(iatom,2)
			
			do jatom = 1,natom
			if (iatom.ne.jatom) then
				x(:)=ap(jatom,:)-ap(iatom,:) !x() is the distance between ith and jth particle
				x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention
				rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
				x(:) = x(:)/rsqrt

				uh(jatom,1,:) = ah(jatom,1,:)-ap(jatom,:)
				uh(jatom,1,:) = uh(jatom,1,:)-boxl*anint(uh(jatom,1,:)/boxl)
				patch_rsqrt(jatom,1) = sqrt(uh(jatom,1,1)*uh(jatom,1,1)+uh(jatom,1,2)*uh(jatom,1,2)+uh(jatom,1,3)*uh(jatom,1,3)) 
				uh(jatom,1,:) = uh(jatom,1,:)/patch_rsqrt(jatom,1)
			
				uh(jatom,2,:) = ah(jatom,2,:)-ap(jatom,:)
				uh(jatom,2,:) = uh(jatom,2,:)-boxl*anint(uh(jatom,2,:)/boxl)
				patch_rsqrt(jatom,2) = sqrt(uh(jatom,2,1)*uh(jatom,2,1)+uh(jatom,2,2)*uh(jatom,2,2)+uh(jatom,3,3)*uh(jatom,3,3))
				uh(jatom,2,:) = uh(jatom,2,:)/patch_rsqrt(jatom,2)


				theta11 = ( (x(1)*uh(iatom,1,1) + x(2)*uh(iatom,1,2) + x(3)*uh(iatom,1,3)) )
				theta12 = ( (x(1)*uh(iatom,2,1) + x(2)*uh(iatom,2,2) + x(3)*uh(iatom,2,3)) )

				x(:)=ap(jatom,:)-ap(iatom,:) !x() is the distance between ith and jth particle
				x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention
				rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
				x(:) = x(:)/rsqrt

				theta21 = ( (x(1)*uh(jatom,1,1) + x(2)*uh(jatom,1,2) + x(3)*uh(jatom,1,3)) )
				theta22 = ( (x(1)*uh(jatom,2,1) + x(2)*uh(jatom,2,2) + x(3)*uh(jatom,2,3)) )

				!write(*,*) iatom, jatom, theta11*57.2958, theta12*57.2958, theta21*57.2958, theta22*57.2958

				rsqr = rsqrt*rsqrt
				rhex = rsqr*rsqr*rsqr
				rtwelve = rhex*rhex

				if(rsqrt.gt.10) E_final = E_final 
				if(rsqrt.gt.sgma+0.895.and.rsqrt.le.10) E_final = E_final + eps11*(sgmatwelve/rtwelve)
				if(rsqrt.gt.sgma.and.rsqrt.le.sgma+0.895) E_final = E_final
				if(rsqrt.le.sgma) E_final = E_final + eps11*(sgmatwelve/rtwelve)

				!if(rsqrt.le.10) then
				!	rsqr = rsqrt*rsqrt
				!	rhex = rsqr*rsqr*rsqr
				!	rtwelve = rhex*rhex
				!	E_final =E_final + eps11*(sgmatwelve/rtwelve)
				!endif
			
				if(rsqrt.ge.sgma.and.rsqrt.le.sgma+0.895) then
			write(33,*) iatom,jatom,rsqrt,theta11*57.29,theta12*57.29,theta21*57.29,theta22*57.29, theta0*57.29,E_final
			if(theta11.le.theta0.and.theta21.le.theta0) E_final = E_final - eps11*5!0.5*((theta11-theta0)**2+(theta21-theta0)**2)-5
			if(theta11.le.theta0.and.theta22.le.theta0) E_final = E_final - eps11*5!0.5*((theta11-theta0)**2+(theta22-theta0)**2)-5
			if(theta12.le.theta0.and.theta21.le.theta0) E_final = E_final - eps11*5!0.5*((theta12-theta0)**2+(theta21-theta0)**2)-5
			if(theta12.le.theta0.and.theta22.le.theta0) E_final = E_final - eps11*5!0.5*((theta12-theta0)**2+(theta22-theta0)**2)-5
				endif

			endif
			enddo !jatom
		enddo !iatom

		write(30,*) imc, E_final

	enddo !imc

	end program mc