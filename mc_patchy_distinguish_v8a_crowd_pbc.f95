!	Monte-Carlo PROGRAM TO SIMULATE patchy spheres with crowding agents
!	modified date - 18/12/2022
!	need to check the update position of patch positions when it goes beyond the boundary
!


	program mc

	implicit none

	real,allocatable :: uh(:,:,:),ah(:,:,:),ap(:,:),a_p(:,:),ac(:,:),a_c(:,:),x(:),patch_rsqrt(:,:), acn(:,:)
	real,allocatable :: uhn(:,:,:),ahn(:,:,:),a_h(:,:,:),rot(:,:),apn(:,:), bonds(:), before_bond(:,:), after_bond(:,:)

	real		   :: boxl,yt,OEp,OEc,OEpc,rsqrt,rhex,rtwelve, R, T, RT, p_r, theta0, maxtheta, rad_2_deg, const
	real		   :: Ei, Ef, dE, d_m, theta11, theta12, theta21, theta22, two_pie, sin_theta, sin_psi, sin_phi, rand, psi, phi, onemcospsi
	real		   :: ny,nxz,nxy,nxsq,nx,nzsq, nz, nyz, nysq, dpsi_max, detrot, cos_theta, cos_psi, cos_phi,rsqrt_old,com1,com2,com3
	real		   :: eps11, eps12, eps22, sgmahex11, sgmahex12, sgmahex22, sgma, rsqr, sgmatwelve,sgmahex, E_final_patch
	real		   :: accp_patch,rej_patch,accp_crowd,rej_crowd,E_final, E_ini, patch_width, eps_patch, eps_drive, full_patch
	real              :: zero_prob,one_prob,two_prob,three_prob,four_prob,five_prob,six_prob,first_occurance, Elj_final
	real		   :: sgma12, sgma22, sgmatwelve12, sgmatwelve22,order_parameter
	integer	   :: j, iatom, icrowd, natom, ncrowd, m, jatom, jcrowd, imc, nmc, b, c
	integer		:: ipatch, npatch, jpatch, stable_at_target, total_bonds, nbonds
	integer,allocatable :: count_bonds(:), bond1(:),bond2(:)
	character(2) :: s

!	OPEN FILES
	open(unit=1,file='input_mc.dat')
	open(unit=2,file='ini_simu.xyz')
	open(unit=3,file='initial_simu.xyz')
	open(unit=10,file='output.dat')
	open(unit=11,file='accep_prob.dat')
	open(unit=12,file='structure_dist.dat')
	open(unit=13,file='traj.xyz')

!	READ THE INPUT DATA
!
	read(1,*)
    	read(1,*) boxl, nmc, natom, ncrowd, d_m, dpsi_max, R, T
	read(1,*)
	read(1,*) theta0, maxtheta, eps_patch, eps_drive, patch_width
	read(1,*)
	read(1,*) eps11, eps12, eps22
	read(1,*)
	read(1,*) sgma, sgmatwelve, sgmahex
	read(1,*)
	read(1,*) sgma12, sgmatwelve12, sgmahex12
	read(1,*) 
	read(1,*) sgma22, sgmatwelve22, sgmahex22
 
!	ALLOCATE DATA
	allocate ( ah(1:natom,1:2,1:3),ap(1:natom,1:3),a_p(1:natom,1:3),ac(1:ncrowd,1:3),a_c(1:ncrowd,1:3),x(1:3),patch_rsqrt(1:natom,1:2))
	allocate ( uh(1:natom,1:2,1:3),uhn(1:natom,1:2,1:3),ahn(1:natom,1:2,1:3),a_h(1:natom,1:2,1:3),rot(1:3,1:3),apn(1:natom,1:3) )
	allocate ( acn(1:ncrowd,1:3) )
	RT = R * T; npatch = natom; rad_2_deg = 57.29; full_patch = 0
	zero_prob = 0;one_prob = 0;two_prob = 0;three_prob = 0;four_prob = 0;five_prob = 0;six_prob = 0;first_occurance = 0
	stable_at_target = 0; total_bonds = natom*2
	allocate ( bond1(1:natom),bond2(1:natom),bonds(1:natom),count_bonds(0:total_bonds) )
	allocate ( before_bond(1:natom,1:natom), after_bond(1:natom,1:natom) )
	before_bond = 0; after_bond = 0
	do nbonds = 0,total_bonds
		count_bonds(nbonds) = 0
	enddo !nbonds
	    
! 	WRITE IMPORTANT PARAMTERS OF THESE MC SIMULATION

	write(10,*) 'PARAMTERS OF MC SIMULATION'
	write(10,*) 'boxlength	No of MC Steps	No. of Particles	No. of Crowding agents	d_m	R	T	Theta0		maxtheta'
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
		ah(iatom,2,1) = ap(iatom,1) + sgma*0.5*cos(maxtheta) ; 
		ah(iatom,2,2) = ap(iatom,2) + sgma*0.5*sin(maxtheta) ; ah(iatom,2,3) = ap(iatom,3)
	enddo !iatom	

	read(2,*)
	read(2,*)

	do iatom = 1,0!natom
		read(2,*) s,ap(iatom,1),ap(iatom,2),ap(iatom,3)
	enddo !iatom
		
	do iatom = 1,0!natom
		read(2,*) s,ah(iatom,1,1),ah(iatom,1,2),ah(iatom,1,3)
		read(2,*) s,ah(iatom,2,1),ah(iatom,2,2),ah(iatom,2,3)
	enddo !iatom

!	GENERATE THE CROWD POSITIONS RANDOMLY

	do icrowd = 1, ncrowd
		call random_number(yt); ac(icrowd,1) = yt*boxl
		call random_number(yt);	ac(icrowd,2) = yt*boxl
		call random_number(yt);	ac(icrowd,3) = yt*boxl
	enddo !i

!	CALCULATE THE INITINAL ENERGY

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
		patch_rsqrt(iatom,2) = sqrt(uh(iatom,2,1)*uh(iatom,2,1)+uh(iatom,2,2)*uh(iatom,2,2)+uh(iatom,2,3)*uh(iatom,2,3))
		uh(iatom,2,:) = uh(iatom,2,:)/patch_rsqrt(iatom,2)
			
		do jatom = 1,natom
			
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
			patch_rsqrt(jatom,2) = sqrt(uh(jatom,2,1)*uh(jatom,2,1)+uh(jatom,2,2)*uh(jatom,2,2)+uh(jatom,2,3)*uh(jatom,2,3))
			uh(jatom,2,:) = uh(jatom,2,:)/patch_rsqrt(jatom,2)

			theta11 = (x(1)*uh(iatom,1,1) + x(2)*uh(iatom,1,2) + x(3)*uh(iatom,1,3))
			theta12 = (x(1)*uh(iatom,2,1) + x(2)*uh(iatom,2,2) + x(3)*uh(iatom,2,3))

			if(theta11.gt.1.0.and.theta11.lt.1.01) theta11 = theta11 - (theta11-1)
			if(theta12.gt.1.0.and.theta12.lt.1.01) theta12 = theta12 - (theta12-1)

			theta11 = acos ( theta11 )
			theta12 = acos ( theta12 )

			x(:)=ap(iatom,:)-ap(jatom,:) !x() is the distance between ith and jth particle
			x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention
			rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
			x(:) = x(:)/rsqrt

			theta21 = (x(1)*uh(jatom,1,1) + x(2)*uh(jatom,1,2) + x(3)*uh(jatom,1,3))
			theta22 = (x(1)*uh(jatom,2,1) + x(2)*uh(jatom,2,2) + x(3)*uh(jatom,2,3))

			if(theta21.gt.1.0.and.theta21.lt.1.01) theta21 = theta21 - (theta21-1)
			if(theta22.gt.1.0.and.theta22.lt.1.01) theta22 = theta22 - (theta22-1)

			theta21 = acos ( theta21 )
			theta22 = acos ( theta22 )

			if(rsqrt.ge.sgma.and.rsqrt.le.sgma+patch_width) then
				!write(*,*) 'Eureka'
				
				if(cos(theta11).ge.cos(theta0).and.cos(theta21).ge.cos(theta0)) then
                                E_ini = E_ini - eps_patch
				endif
				if(cos(theta11).ge.cos(theta0).and.cos(theta22).ge.cos(theta0)) then
                                E_ini = E_ini - eps_patch
				endif
				if(cos(theta12).ge.cos(theta0).and.cos(theta21).ge.cos(theta0)) then
                                E_ini = E_ini - eps_patch
				endif
				if(cos(theta12).ge.cos(theta0).and.cos(theta22).ge.cos(theta0)) then
                                E_ini = E_ini - eps_patch
				endif
			endif

			!write(*,*) iatom,jatom,E_ini,rsqrt

			

		enddo !jatom
	enddo !iatom
	
        
	write(32,*) 'Initial Energy',E_ini
        write(32,*) 'MC Cycle','        Final Energy'

!	WRITE THE INITIAL CONFIGURATION FILE

	write(3,*) natom*3+ncrowd
	write(3,*)

	do iatom = 1,natom
		write(3,*) 'Xe',ap(iatom,1),ap(iatom,2),ap(iatom,3)
		write(3,*) 'P',ah(iatom,1,1),ah(iatom,1,2),ah(iatom,1,3)
		write(3,*) 'P',ah(iatom,2,1),ah(iatom,2,2),ah(iatom,2,3)
	enddo !i
	do icrowd = 1,ncrowd
		write(2,*) 'Ar',ac(icrowd,1),ac(icrowd,2),ac(icrowd,3)
	enddo !i

								! ---------------- MC simulation ------------------!


	accp_patch = 0; rej_patch = 0


!	STARTING MC CYCLE

	do imc = 1, nmc   ! nmc is the no of simulation


!	STARTING iatom LOOP

		do icrowd = 1,ncrowd

			Ei = 0; Ef = 0
			do jcrowd = 1,ncrowd

				if(icrowd.ne.jcrowd) then

					x(1) = 0; x(2) = 0; x(3) = 0
					x(1) = ac(icrowd,1) - ac(jcrowd,1)
					x(2) = ac(icrowd,2) - ac(jcrowd,2)
					x(3) = ac(icrowd,3) - ac(jcrowd,3) 
				        x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention	
					rsqrt = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
					x(1) = x(1)/rsqrt; x(2) = x(2)/rsqrt; x(3) = x(3)/rsqrt

					if(rsqrt.lt.5) then
						rsqr = rsqrt*rsqrt
						rhex = rsqr*rsqr*rsqr
						rtwelve = rhex*rhex
						Ei =Ei + eps11*((sgmatwelve22/rtwelve) )
					endif

				endif

			enddo !jcrowd

			do iatom = 1,natom

				x(1) = 0; x(2) = 0; x(3) = 0
				x(1) = ac(icrowd,1) - ap(iatom,1)
				x(2) = ac(icrowd,2) - ap(iatom,2)
				x(3) = ac(icrowd,3) - ap(iatom,3)
                                x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention
				rsqrt = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
				x(1) = x(1)/rsqrt; x(2) = x(2)/rsqrt; x(3) = x(3)/rsqrt

				if(rsqrt.lt.5) then
					rsqr = rsqrt*rsqrt
					rhex = rsqr*rsqr*rsqr
					rtwelve = rhex*rhex
					Ei =Ei + eps11*((sgmatwelve12/rtwelve))
				endif

			enddo !iatom

			do j=1,3
				call random_number(yt)
				acn(icrowd,j)=ac(icrowd,j)+(2*yt-1)*d_m  ! move ith particle randomly

!	REFLECTION WALL			
				!do while(acn(icrowd,j).gt.boxl)
				!call random_number(yt)
				!acn(icrowd,j)=ac(icrowd,j)+(2*yt-1)*d_m
				!end do!while
				!do while(acn(icrowd,j).lt.0.0)
				!call random_number(yt)
				!acn(icrowd,j)=ac(icrowd,j)+(2*yt-1)*d_m
				!end do!while


                        if(acn(icrowd,j).gt.boxl) acn(icrowd,j)=acn(icrowd,j)-boxl*1.0d0
                        if(acn(icrowd,j).lt.0.0d0) acn(icrowd,j)=acn(icrowd,j)+boxl*1.0d0

			enddo !j


			do jcrowd = 1,ncrowd

				if(icrowd.ne.jcrowd) then

					x(1) = 0; x(2) = 0; x(3) = 0
					x(1) = acn(icrowd,1) - ac(jcrowd,1)
					x(2) = acn(icrowd,2) - ac(jcrowd,2)
					x(3) = acn(icrowd,3) - ac(jcrowd,3) 
					x(:)=x(:)-boxl*anint(x(:)/boxl)  !Minimum image convention
					rsqrt = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
					x(1) = x(1)/rsqrt; x(2) = x(2)/rsqrt; x(3) = x(3)/rsqrt

					if(rsqrt.lt.5) then
						rsqr = rsqrt*rsqrt
						rhex = rsqr*rsqr*rsqr
						rtwelve = rhex*rhex
						Ef =Ef + eps11*((sgmatwelve22/rtwelve))
					endif

				endif

			enddo !jcrowd

			do iatom = 1,natom

				x(1) = 0; x(2) = 0; x(3) = 0
				x(1) = acn(icrowd,1) - ap(iatom,1)
				x(2) = acn(icrowd,2) - ap(iatom,2)
				x(3) = acn(icrowd,3) - ap(iatom,3)
                                x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention
				rsqrt = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
				x(1) = x(1)/rsqrt; x(2) = x(2)/rsqrt; x(3) = x(3)/rsqrt

				if(rsqrt.lt.5) then
					rsqr = rsqrt*rsqrt
					rhex = rsqr*rsqr*rsqr
					rtwelve = rhex*rhex
					Ef =Ef + eps11*((sgmatwelve12/rtwelve))
				endif

			enddo !iatom

			dE = 0; yt = 0; p_r = 0
			dE=Ef-Ei   ! dE is the difference between new and initial energy

			if (dE.gt.0) then
				call random_number(yt)
				p_r=exp(-dE/(RT)) !Kb is the Boltzmann factor and T is the temperature
			endif

			if (dE.lt.0.or.yt.lt.p_r) then
				ac(icrowd,1) = acn(icrowd,1); ac(icrowd,2) = acn(icrowd,2); ac(icrowd,3) = acn(icrowd,3)
			endif


		enddo !icrowd


		do iatom = 1,natom

		Ei = 0; Ef = 0

		do icrowd = 1,ncrowd

			x(1) = 0; x(2) = 0; x(3) = 0
			x(1) = ac(icrowd,1) - ap(iatom,1)
			x(2) = ac(icrowd,2) - ap(iatom,2)
			x(3) = ac(icrowd,3) - ap(iatom,3)
                        x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention
			rsqrt = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
			x(1) = x(1)/rsqrt; x(2) = x(2)/rsqrt; x(3) = x(3)/rsqrt

			if(rsqrt.lt.5) then
				rsqr = rsqrt*rsqrt
				rhex = rsqr*rsqr*rsqr
				rtwelve = rhex*rhex
				Ei =Ei + eps11*((sgmatwelve12/rtwelve))
			endif

		enddo !icrowd 

!	CALCULATE INITIAL ENERGY FOR iatom
			
		patch_rsqrt(iatom,:) = 0
		uh(iatom,:,:) = 0

		uh(iatom,1,:) = ah(iatom,1,:)-ap(iatom,:)
		uh(iatom,1,:) = uh(iatom,1,:)-boxl*anint(uh(iatom,1,:)/boxl)
		patch_rsqrt(iatom,1) = sqrt(uh(iatom,1,1)*uh(iatom,1,1)+uh(iatom,1,2)*uh(iatom,1,2)+uh(iatom,1,3)*uh(iatom,1,3)) 
		uh(iatom,1,:) = uh(iatom,1,:)/patch_rsqrt(iatom,1)
			
		uh(iatom,2,:) = ah(iatom,2,:)-ap(iatom,:)
		uh(iatom,2,:) = uh(iatom,2,:)-boxl*anint(uh(iatom,2,:)/boxl)
		patch_rsqrt(iatom,2) = sqrt(uh(iatom,2,1)*uh(iatom,2,1)+uh(iatom,2,2)*uh(iatom,2,2)+uh(iatom,2,3)*uh(iatom,2,3))
		uh(iatom,2,:) = uh(iatom,2,:)/patch_rsqrt(iatom,2)
		
                before_bond(iatom,:) = 0; after_bond(iatom,:) = 0

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
			patch_rsqrt(jatom,2) = sqrt(uh(jatom,2,1)*uh(jatom,2,1)+uh(jatom,2,2)*uh(jatom,2,2)+uh(jatom,2,3)*uh(jatom,2,3))
			uh(jatom,2,:) = uh(jatom,2,:)/patch_rsqrt(jatom,2)

			theta11 = (x(1)*uh(iatom,1,1) + x(2)*uh(iatom,1,2) + x(3)*uh(iatom,1,3))
			theta12 = (x(1)*uh(iatom,2,1) + x(2)*uh(iatom,2,2) + x(3)*uh(iatom,2,3))

			if(theta11.gt.1.0.and.theta11.lt.1.01) theta11 = 1.00
			if(theta12.gt.1.0.and.theta12.lt.1.01) theta12 = 1.00

			theta11 = acos ( theta11 )
			theta12 = acos ( theta12 )

			x(:)=ap(iatom,:)-ap(jatom,:) !x() is the distance between ith and jth particle
			x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention
			rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
			x(:) = x(:)/rsqrt

			theta21 = (x(1)*uh(jatom,1,1) + x(2)*uh(jatom,1,2) + x(3)*uh(jatom,1,3))
			theta22 = (x(1)*uh(jatom,2,1) + x(2)*uh(jatom,2,2) + x(3)*uh(jatom,2,3))

			if(theta21.gt.1.0.and.theta21.lt.1.01) theta21 = 1.00
			if(theta22.gt.1.0.and.theta22.lt.1.01) theta22 = 1.00

			theta21 = acos ( theta21 )
			theta22 = acos ( theta22 )

			!write(*,*) iatom,jatom,'theta11',theta11*rad_2_deg,'theta12',theta12*rad_2_deg,'theta21', &
			!		theta21*rad_2_deg,'theta22',theta22*rad_2_deg

			if(rsqrt.lt.5) then
				rsqr = rsqrt*rsqrt
				rhex = rsqr*rsqr*rsqr
				rtwelve = rhex*rhex
				Ei =Ei + eps11*((sgmatwelve/rtwelve))	!-(sgmahex/rhex))
			endif
			
			if(rsqrt.ge.sgma.and.rsqrt.le.sgma+patch_width) then
				!write(*,*) 'Initial Energy', Ei
				!write(*,*) 'initial',iatom,jatom,'theta11',theta11*rad_2_deg,'theta12',theta12*rad_2_deg,'theta21', &
				!	theta21*rad_2_deg,'theta22',theta22*rad_2_deg,'theta0',theta0*rad_2_deg
				!write(*,*) 'initial',iatom,jatom,'theta11',cos(theta11),'theta12',cos(theta12),'theta21', &
				!	cos(theta21),'theta22',cos(theta22),'theta0',cos(theta0)

				if(cos(theta11).ge.cos(theta0).and.cos(theta21).ge.cos(theta0)) then
						Ei = Ei - eps_patch; before_bond(iatom,jatom) = 1
				endif

				if(cos(theta11).ge.cos(theta0).and.cos(theta22).ge.cos(theta0)) then
						Ei = Ei - eps_patch; before_bond(iatom,jatom) = 1
				endif

				if(cos(theta12).ge.cos(theta0).and.cos(theta21).ge.cos(theta0)) then
						Ei = Ei - eps_patch; before_bond(iatom,jatom) = 1
				endif

				if(cos(theta12).ge.cos(theta0).and.cos(theta22).ge.cos(theta0)) then
						Ei = Ei - eps_patch; before_bond(iatom,jatom) = 1
				endif

				!write(*,*) 'Initial Energy', Ei
			endif

			!write(*,*) 'MC Step',imc,iatom,jatom,Ei,rsqrt

			endif

		enddo !jatom

		!write(40,*) iatom,bias(iatom,:)



!	RANDOM DISPLACEMENT

!	ROTATION

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

!	REFLECTION WALL			
				!do while(apn(iatom,j).gt.boxl)
				!call random_number(yt)
				!apn(iatom,j)=ap(iatom,j)+(2*yt-1)*d_m
				!end do!while
				!do while(apn(iatom,j).lt.0.0)
				!call random_number(yt)
				!apn(iatom,j)=ap(iatom,j)+(2*yt-1)*d_m
				!end do!while

!	PERIODIC BOUNDARY
				if(apn(iatom,j).gt.boxl) apn(iatom,j)=apn(iatom,j)-boxl*1.0d0
				if(apn(iatom,j).lt.0.0d0) apn(iatom,j)=apn(iatom,j)+boxl*1.0d0
			enddo !j

!	CHANGE THE PATCHES BY TRANSLATION FROM THE OLD PATCH POSITION	#PATCH INDEX 1
			a_h(iatom,1,1) = ah(iatom,1,1) + (apn(iatom,1)-ap(iatom,1))
			a_h(iatom,1,2) = ah(iatom,1,2) + (apn(iatom,2)-ap(iatom,2))
			a_h(iatom,1,3) = ah(iatom,1,3) + (apn(iatom,3)-ap(iatom,3))
			
!	ROTATING THE PATCHES		#PATCH INDEX 1
	ahn(iatom,1,1)=rot(1,1)*(a_h(iatom,1,1)-apn(iatom,1))+rot(1,2)*(a_h(iatom,1,2)-apn(iatom,2))+rot(1,3)*(a_h(iatom,1,3)-apn(iatom,3))
	ahn(iatom,1,2)=rot(2,1)*(a_h(iatom,1,1)-apn(iatom,1))+rot(2,2)*(a_h(iatom,1,2)-apn(iatom,2))+rot(2,3)*(a_h(iatom,1,3)-apn(iatom,3))
	ahn(iatom,1,3)=rot(3,1)*(a_h(iatom,1,1)-apn(iatom,1))+rot(3,2)*(a_h(iatom,1,2)-apn(iatom,2))+rot(3,3)*(a_h(iatom,1,3)-apn(iatom,3))

!	FINAL POSITION OF NEW PATCH		#PATCH INDEX 1
			ahn(iatom,1,1)=ahn(iatom,1,1)+apn(iatom,1)
			ahn(iatom,1,2)=ahn(iatom,1,2)+apn(iatom,2)
			ahn(iatom,1,3)=ahn(iatom,1,3)+apn(iatom,3)


!	CHANGE THE PATCHES BY TRANSLATION FROM THE OLD PATCH POSITION	#PATCH INDEX 2
			a_h(iatom,2,1) = ah(iatom,2,1) + (apn(iatom,1)-ap(iatom,1))
			a_h(iatom,2,2) = ah(iatom,2,2) + (apn(iatom,2)-ap(iatom,2))
			a_h(iatom,2,3) = ah(iatom,2,3) + (apn(iatom,3)-ap(iatom,3))
			
!	ROTATING THE PATCHES		#PATCH INDEX 2		
	ahn(iatom,2,1)=rot(1,1)*(a_h(iatom,2,1)-apn(iatom,1))+rot(1,2)*(a_h(iatom,2,2)-apn(iatom,2))+rot(1,3)*(a_h(iatom,2,3)-apn(iatom,3))
	ahn(iatom,2,2)=rot(2,1)*(a_h(iatom,2,1)-apn(iatom,1))+rot(2,2)*(a_h(iatom,2,2)-apn(iatom,2))+rot(2,3)*(a_h(iatom,2,3)-apn(iatom,3))
	ahn(iatom,2,3)=rot(3,1)*(a_h(iatom,2,1)-apn(iatom,1))+rot(3,2)*(a_h(iatom,2,2)-apn(iatom,2))+rot(3,3)*(a_h(iatom,2,3)-apn(iatom,3))

!	FINAL POSITION OF NEW PATCH		#PATCH INDEX 2
			ahn(iatom,2,1)=ahn(iatom,2,1)+apn(iatom,1)
			ahn(iatom,2,2)=ahn(iatom,2,2)+apn(iatom,2)
			ahn(iatom,2,3)=ahn(iatom,2,3)+apn(iatom,3)
			

			
!	CALCULATE FINAL ENERGY FOR iatom
			
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
			x(:)=ap(jatom,:)-apn(iatom,:) !x() is the distance between ith and jth particle
			x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention
			rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
			x(:) = x(:)/rsqrt

			uh(jatom,1,:) = ah(jatom,1,:)-ap(jatom,:)
			uh(jatom,1,:) = uh(jatom,1,:)-boxl*anint(uh(jatom,1,:)/boxl)
			patch_rsqrt(jatom,1) = sqrt(uh(jatom,1,1)*uh(jatom,1,1)+uh(jatom,1,2)*uh(jatom,1,2)+uh(jatom,1,3)*uh(jatom,1,3)) 
			uh(jatom,1,:) = uh(jatom,1,:)/patch_rsqrt(jatom,1)
			
			
			uh(jatom,2,:) = ah(jatom,2,:)-ap(jatom,:)
			uh(jatom,2,:) = uh(jatom,2,:)-boxl*anint(uh(jatom,2,:)/boxl)
			patch_rsqrt(jatom,2) = sqrt(uh(jatom,2,1)*uh(jatom,2,1)+uh(jatom,2,2)*uh(jatom,2,2)+uh(jatom,2,3)*uh(jatom,2,3))
			uh(jatom,2,:) = uh(jatom,2,:)/patch_rsqrt(jatom,2)
			if(uh(jatom,2,1).gt.1.0.or.uh(jatom,2,2).gt.1.0.or.uh(jatom,2,3).gt.1.0) stop

			theta11 = (x(1)*uh(iatom,1,1) + x(2)*uh(iatom,1,2) + x(3)*uh(iatom,1,3))
			theta12 = (x(1)*uh(iatom,2,1) + x(2)*uh(iatom,2,2) + x(3)*uh(iatom,2,3))

			if(theta11.gt.1.0.and.theta11.lt.1.01) theta11 = 1.00
			if(theta12.gt.1.0.and.theta12.lt.1.01) theta12 = 1.00

			theta11 = acos ( theta11 )
			theta12 = acos ( theta12 )

			x(:)=apn(iatom,:)-ap(jatom,:) !x() is the distance between ith and jth particle
			x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention
			rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
			x(:) = x(:)/rsqrt

			theta21 = (x(1)*uh(jatom,1,1) + x(2)*uh(jatom,1,2) + x(3)*uh(jatom,1,3))
			theta22 = (x(1)*uh(jatom,2,1) + x(2)*uh(jatom,2,2) + x(3)*uh(jatom,2,3))

			if(theta21.gt.1.0.and.theta21.lt.1.01) theta21 = 1.00
			if(theta22.gt.1.0.and.theta22.lt.1.01) theta22 = 1.00

			theta21 = acos ( theta21 )
			theta22 = acos ( theta22 )

			!write(*,*) iatom,jatom,'theta11',theta11*rad_2_deg,'theta12',theta12*rad_2_deg,'theta21', &
			!		theta21*rad_2_deg,'theta22',theta22*rad_2_deg


			if(rsqrt.lt.5) then
				rsqr = rsqrt*rsqrt
				rhex = rsqr*rsqr*rsqr
				rtwelve = rhex*rhex
				Ef =Ef + eps11*((sgmatwelve/rtwelve))	!-(sgmahex/rhex))
			endif

				!write(*,*) 
			if(rsqrt.ge.sgma.and.rsqrt.le.sgma+patch_width) then

				!write(*,*) 'Final Energy', Ef
				!write(*,*) 'final',iatom,jatom,'theta11',theta11*rad_2_deg,'theta12',theta12*rad_2_deg,'theta21', &
				!	theta21*rad_2_deg,'theta22',theta22*rad_2_deg,'theta0',theta0*rad_2_deg
				!write(*,*) 'final',iatom,jatom,'theta11',cos(theta11),'theta12',cos(theta12),'theta21', &
				!	cos(theta21),'theta22',cos(theta22),'theta0',cos(theta0)

				if(cos(theta11).ge.cos(theta0).and.cos(theta21).ge.cos(theta0)) then
						Ef = Ef - eps_patch; after_bond(iatom,jatom) = 1
				endif

				if(cos(theta11).ge.cos(theta0).and.cos(theta22).ge.cos(theta0)) then
						Ef = Ef - eps_patch; after_bond(iatom,jatom) = 1
				endif

				if(cos(theta12).ge.cos(theta0).and.cos(theta21).ge.cos(theta0)) then
						Ef = Ef - eps_patch; after_bond(iatom,jatom) = 1
				endif

				if(cos(theta12).ge.cos(theta0).and.cos(theta22).ge.cos(theta0)) then
						Ef = Ef - eps_patch; after_bond(iatom,jatom) = 1
				endif

				!write(*,*) 'Final Energy', Ef
			endif

			!write(*,*) 'MC Step',imc,iatom,jatom,Ef,rsqrt

			endif

		if (before_bond(iatom,jatom).eq.0.and.after_bond(iatom,jatom).eq.0) then
			Ef = Ef
			!write(39,*) imc,iatom,jatom,before_bond(iatom,jatom),after_bond(iatom,jatom),'no change',Ef
		endif
		if (before_bond(iatom,jatom).eq.0.and.after_bond(iatom,jatom).eq.1) then
			!write(40,*) Ef
			Ef = Ef - eps_drive
			!write(40,*) imc,iatom,jatom,before_bond(iatom,jatom),after_bond(iatom,jatom),'drive', Ef
		endif
		if (before_bond(iatom,jatom).eq.1.and.after_bond(iatom,jatom).eq.0) then
			!write(41,*) Ef
			Ef = Ef + eps_drive
			!write(41,*) imc,iatom,jatom,before_bond(iatom,jatom),after_bond(iatom,jatom),'drive', Ef
		endif

		enddo !jatom

		do icrowd = 1,ncrowd

			x(1) = 0; x(2) = 0; x(3) = 0
			x(1) = ac(icrowd,1) - apn(iatom,1)
			x(2) = ac(icrowd,2) - apn(iatom,2)
			x(3) = ac(icrowd,3) - apn(iatom,3)
                        x(:) = x(:)-boxl*anint(x(:)/boxl) !Minimum image  convention
			rsqrt = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
			x(1) = x(1)/rsqrt; x(2) = x(2)/rsqrt; x(3) = x(3)/rsqrt

			if(rsqrt.lt.5) then
				rsqr = rsqrt*rsqrt
				rhex = rsqr*rsqr*rsqr
				rtwelve = rhex*rhex
				Ef =Ef + eps11*((sgmatwelve12/rtwelve))
			endif

		enddo !icrowd 


!	SAMPLING


		dE = 0; yt = 0; p_r = 0
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


		!write(31,*) imc,iatom,Ei,Ef,dE,yt,p_r,RT,accp_patch


		enddo !iatom


!---- write no of nth cycle,Energy at nth cycle,Energy of each particle after nth cycle  in 'energy.dat' file to plot---!

		if(mod(imc,1000).eq.0) then
                         !if(accp_patch/(accp_patch+rej_patch).gt.0.65) d_m=d_m*1.05
                         !if(accp_patch/(accp_patch+rej_patch).lt.0.65) d_m=d_m*0.95
                         
			 do jatom = 1,natom
			 	com1 = com1 + ap(jatom,1)
				com2 = com2 + ap(jatom,2)
				com3 = com3 + ap(jatom,3)
			enddo !jatom
			com1 = com1/natom; com2 = com2/natom; com3 = com3/natom
                         write(11,15) imc,accp_patch/(accp_patch+rej_patch), &
                                        com1,com2,com3,sqrt(com1*com1+com2*com2+com3*com3)
15              format(2x,i9,2x,f7.2,2x,f7.2,2x,f7.2,2x,f7.2,2x,f7.2)

                        accp_patch = 0; rej_patch = 0; com1 = 0; com2 = 0; com3 = 0
                endif


! --------- WRITE THE TRAJECTORY AT INTERVAL OF 1000 STEPS ------------!

		if(mod(imc,5000).eq.0) then
			write(13,*) natom*3+ncrowd
			write(13,*)
			do iatom=1,natom ! natom is the no of atom in simulation box
                                write(13,*) 'A',ap(iatom,1),ap(iatom,2),ap(iatom,3)
				write(13,*) 'P',ah(iatom,1,1),ah(iatom,1,2),ah(iatom,1,3)
				write(13,*) 'P',ah(iatom,2,1),ah(iatom,2,2),ah(iatom,2,3)
			enddo !i
			
			do icrowd=1,ncrowd ! natom is the no of atom in simulation box
				write(13,*) 'M',ac(icrowd,1),ac(icrowd,2),ac(icrowd,3)
			enddo !i
		endif
	
14      format(2x,A,2x,f7.2,2x,f7.2,2x,f7.2)	

!	CALCULATE THE FINAL ENERGY
	if(mod(imc,20).eq.0) then
        order_parameter = 0; nbonds = 0
	E_final = 0; Elj_final = 0
	do iatom = 1,natom
		patch_rsqrt(iatom,:) = 0
		uh(iatom,:,:) = 0

		uh(iatom,1,:) = ah(iatom,1,:)-ap(iatom,:)
		uh(iatom,1,:) = uh(iatom,1,:)-boxl*anint(uh(iatom,1,:)/boxl)
		patch_rsqrt(iatom,1) = sqrt(uh(iatom,1,1)*uh(iatom,1,1)+uh(iatom,1,2)*uh(iatom,1,2)+uh(iatom,1,3)*uh(iatom,1,3)) 
		uh(iatom,1,:) = uh(iatom,1,:)/patch_rsqrt(iatom,1)
			
		uh(iatom,2,:) = ah(iatom,2,:)-ap(iatom,:)
		uh(iatom,2,:) = uh(iatom,2,:)-boxl*anint(uh(iatom,2,:)/boxl)
		patch_rsqrt(iatom,2) = sqrt(uh(iatom,2,1)*uh(iatom,2,1)+uh(iatom,2,2)*uh(iatom,2,2)+uh(iatom,2,3)*uh(iatom,2,3))
		uh(iatom,2,:) = uh(iatom,2,:)/patch_rsqrt(iatom,2)

		E_final_patch = 0; bond1(iatom) = 0; bond2(iatom) = 0
			
		do jatom = 1,natom
			if(iatom.ne.jatom) then
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
			patch_rsqrt(jatom,2) = sqrt(uh(jatom,2,1)*uh(jatom,2,1)+uh(jatom,2,2)*uh(jatom,2,2)+uh(jatom,2,3)*uh(jatom,2,3))
			uh(jatom,2,:) = uh(jatom,2,:)/patch_rsqrt(jatom,2)

			theta11 = (x(1)*uh(iatom,1,1) + x(2)*uh(iatom,1,2) + x(3)*uh(iatom,1,3))
			theta12 = (x(1)*uh(iatom,2,1) + x(2)*uh(iatom,2,2) + x(3)*uh(iatom,2,3))

			if(theta11.gt.1.0.and.theta11.lt.1.01) theta11 = 1.00
			if(theta12.gt.1.0.and.theta12.lt.1.01) theta12 = 1.00

			theta11 = acos ( theta11 )
			theta12 = acos ( theta12 )

			x(:)=ap(iatom,:)-ap(jatom,:) !x() is the distance between ith and jth particle
			x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention
			rsqrt=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
			x(:) = x(:)/rsqrt

			theta21 = (x(1)*uh(jatom,1,1) + x(2)*uh(jatom,1,2) + x(3)*uh(jatom,1,3))
			theta22 = (x(1)*uh(jatom,2,1) + x(2)*uh(jatom,2,2) + x(3)*uh(jatom,2,3))

			if(theta21.gt.1.0.and.theta21.lt.1.01) theta21 = 1.00
			if(theta22.gt.1.0.and.theta22.lt.1.01) theta22 = 1.00

			theta21 = acos ( theta21 )
			theta22 = acos ( theta22 )
			
			!write(33,*) imc,iatom,jatom,rsqrt,cos(theta11),cos(theta12),cos(theta21),cos(theta22),cos(theta0)
			!write(33,*) imc,iatom,jatom,rsqrt,theta11,theta12,theta21,theta22,theta0

			if(rsqrt.lt.5) then
				rsqr = rsqrt*rsqrt
				rhex = rsqr*rsqr*rsqr
				rtwelve = rhex*rhex
				Elj_final = Elj_final + eps11*((sgmatwelve/rtwelve))	!-(sgmahex/rhex))
			endif

			if(rsqrt.ge.sgma.and.rsqrt.le.sgma+patch_width) then
				!write(*,*) 'Eureka'
				if(cos(theta11).ge.cos(theta0).and.cos(theta21).ge.cos(theta0)) then
					!write(33,*) imc,iatom,jatom,E_final,'11','21',rsqrt
					E_final_patch = E_final_patch - eps_patch
                                        bond1(iatom) = bond1(iatom) + 1
					!write(33,*) imc,iatom,jatom,E_final,'11','21',rsqrt
				endif
				if(cos(theta11).ge.cos(theta0).and.cos(theta22).ge.cos(theta0)) then
					!write(33,*) imc,iatom,jatom,E_final,'11','22',rsqrt
					E_final_patch = E_final_patch - eps_patch
                                        bond1(iatom) = bond1(iatom) + 1
					!write(33,*) imc,iatom,jatom,E_final,'11','22',rsqrt
				endif
				if(cos(theta12).ge.cos(theta0).and.cos(theta21).ge.cos(theta0)) then
					!write(33,*) imc,iatom,jatom,E_final,'12','21',rsqrt
					E_final_patch = E_final_patch - eps_patch
                                        bond2(iatom) = bond2(iatom) + 1
					!write(33,*) imc,iatom,jatom,E_final,'12','21',rsqrt
				endif
				if(cos(theta12).ge.cos(theta0).and.cos(theta22).ge.cos(theta0)) then
					!write(33,*) imc,iatom,jatom,E_final,'12','22',rsqrt
					E_final_patch = E_final_patch - eps_patch
                                        bond2(iatom) = bond2(iatom) + 1
					!write(33,*) imc,iatom,jatom,E_final,'12','22',rsqrt
				endif
			endif

	!		write(*,*) iatom,jatom,E_final,rsqrt

			endif			

		enddo !jatom

        	if(bond1(iatom).le.1.and.bond2(iatom).le.1) then
        		bonds(iatom) = int(abs(E_final_patch/eps_patch)*1.0)
			E_final = E_final + E_final_patch
        	else
        		bonds(iatom) = 0
        		E_final = E_final
        	endif

	enddo !iatom
        
	E_final = E_final/(natom*eps_patch*2)
	
        do iatom = 1,natom
		nbonds = nbonds + bonds(iatom)
        enddo !iatom
        
        count_bonds(nbonds) = count_bonds(nbonds) + 1
                
        order_parameter = real(nbonds*1.0/total_bonds)

         !write(32,13) imc,E_final, Elj_final, Elj_final +E_final*(natom*eps_patch*2), nbonds, order_parameter
!13      format(2x,i8,2x,f8.2,2x,f8.2,2x,f8.2,i2,2x,f8.2)

        if(abs(E_final).ge.0.0.and.abs(E_final).le.0.1) zero_prob  = zero_prob + 1
        if(abs(E_final).ge.0.1.and.abs(E_final).le.0.2) one_prob = one_prob + 1
        if(abs(E_final).ge.0.3.and.abs(E_final).le.0.4) two_prob = two_prob + 1
        if(abs(E_final).ge.0.4.and.abs(E_final).le.0.6) three_prob = three_prob + 1
        if(abs(E_final).ge.0.6.and.abs(E_final).le.0.7) four_prob  = four_prob + 1
        if(abs(E_final).ge.0.8.and.abs(E_final).le.0.9) five_prob = five_prob + 1
        if(abs(E_final).ge.0.98.and.abs(E_final).le.1.1) six_prob = six_prob + 1
        if(int(six_prob).eq.1.and.first_occurance.eq.0.and.order_parameter.eq.1) first_occurance = imc
	if(int(six_prob)*20.eq.imc) stable_at_target = imc
	!if(stable_at_target.ne.imc) goto 12
	

	endif

	!if(stable_at_target.ne.imc) goto 12
        
	enddo !imc

12 	write(12,*) 'zero_prob','       one_prob','     two_prob','     three_prob', &
            '   four_prob','    five_prob','    six_prob','     first_occurance'
        write(12,*) zero_prob*2/(nmc*0.1),one_prob*2/(nmc*0.1), &
                two_prob*2/(nmc*0.1),three_prob*2/(nmc*0.1),four_prob*2/(nmc*0.1), &
                five_prob*2/(nmc*0.1),six_prob*2/(nmc*0.1),first_occurance,stable_at_target

        write(12,*) count_bonds
        
	end program mc
