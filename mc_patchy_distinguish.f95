!	Monte-Carlo PROGRAM TO SIMULATE patchy spheres with crowding agents

	program mc

	implicit none

	real*8,allocatable :: uh(:,:,:),ah(:,:,:),ap(:,:),a_p(:,:),ac(:,:),a_c(:,:),x(:),patch_rqsrt(:,:)
	real*8,allocatable :: uhn(:,:,:),ahn(:,:,:),a_h(:,:,:),rot(:,:),apn(:,:)
	real*8		   :: boxl,yt,OEp,OEc,OEpc,rsqr,rhex,rtwelve, R, T, RT, p_r
	real*8		   :: Ei, Ef, dE, d_m, theta11, theta12, theta21, theta22, two_pie, sin_theta, sin_psi, sin_phi, rand, psi, phi, onemcospsi
	real*8 		   :: ny,nxz,nxy,nxsq,nx,nzsq, nz, nyz, nysq, dpsi_max, detrot, cos_theta, cos_psi, cos_phi
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
	allocate ( ah(1:natom,1:2,1:3),ap(1:natom,1:3),a_p(1:natom,1:3),ac(1:ncrowd,1:3),a_c(1:ncrowd,1:3),x(1:3),patch_rqsrt(1:natom,1:2))
	allocate ( uh(1:natom,1:2,1:3),uhn(1:natom,1:2,1:3),ahn(1:natom,1:2,1:3),a_h(1:natom,1:2,1:3),rot(1:3,1:3),apn(1:natom,1:3) )
	RT = R * T; npatch = natom
    
! 	WRITE IMPORTANT PARAMTERS OF THESE MC SIMULATION

	write(10,*) 'PARAMTERS OF MC SIMULATION'
	write(10,*) 'boxlength	No. of Particles	No. of Crowding agents'
	write(10,*) boxl, nmc, natom, ncrowd, d_m, R, T, RT

! 	 READ THE INITIAL CONFIGURATION FILE

	do iatom = 1, natom
		call random_number(yt); ap(iatom,1) = yt*boxl
		call random_number(yt);	ap(iatom,2) = yt*boxl
		call random_number(yt);	ap(iatom,3) = yt*boxl
	enddo !i

	do ipatch = 1,npatch
			ah(ipatch,1,1) = ap(ipatch,1) + sgma; ah(ipatch,1,2) = ap(ipatch,2) ; ah(ipatch,1,3) = ap(ipatch,3)
			ah(ipatch,2,1) = ap(ipatch,1) + sgma*cos(2.0944) ; 
			ah(ipatch,2,2) = ap(ipatch,2) + sgma*sin(2.0944) ; ah(ipatch,2,3) = ap(ipatch,3)
	enddo !ipatch

!	CALCULATE ANGLE BETWEEN PATCHES OF DIFFERENT PARTICLES

	do ipatch = 1,npatch
		patch_rqsrt(ipatch,1:2) = 0

			uh(ipatch,1,1) = (ah(ipatch,1,1)-ap(ipatch,1))
			uh(ipatch,1,2) = (ah(ipatch,1,2)-ap(ipatch,2))
			uh(ipatch,1,3) = (ah(ipatch,1,3)-ap(ipatch,3))
			
			uh(ipatch,2,1) = (ah(ipatch,2,1)-ap(ipatch,1))
			uh(ipatch,2,2) = (ah(ipatch,2,2)-ap(ipatch,2))
			uh(ipatch,2,3) = (ah(ipatch,2,3)-ap(ipatch,3)) 

			patch_rqsrt(ipatch,1) = sqrt(uh(ipatch,1,1)*uh(ipatch,1,1)+uh(ipatch,1,2)*uh(ipatch,1,2)+uh(ipatch,1,3)*uh(ipatch,1,3)) 
			patch_rqsrt(ipatch,2) = sqrt(uh(ipatch,2,1)*uh(ipatch,2,1)+uh(ipatch,2,2)*uh(ipatch,2,2)+uh(ipatch,3,3)*uh(ipatch,3,3))

			do jpatch = 1, npatch
			if (ipatch.ne.jpatch) then
					x(:)=ap(ipatch,:)-ap(jpatch,:) !x() is the distance between ith and jth particle
					x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention

				rsqr=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))

				uh(jpatch,1,1) = (ah(jpatch,1,1)-ap(jpatch,1))
				uh(jpatch,1,2) = (ah(jpatch,1,2)-ap(jpatch,2))
				uh(jpatch,1,3) = (ah(jpatch,1,3)-ap(jpatch,3))
			
				uh(jpatch,2,1) = (ah(jpatch,2,1)-ap(jpatch,1))
				uh(jpatch,2,2) = (ah(jpatch,2,2)-ap(jpatch,2))
				uh(jpatch,2,3) = (ah(jpatch,2,3)-ap(jpatch,3)) 

				patch_rqsrt(jpatch,1) = sqrt(uh(jpatch,1,1)*uh(jpatch,1,1)+uh(jpatch,1,2)*uh(jpatch,1,2)+uh(jpatch,1,3)*uh(jpatch,1,3)) 
				patch_rqsrt(jpatch,2) = sqrt(uh(jpatch,2,1)*uh(jpatch,2,1)+uh(jpatch,2,2)*uh(jpatch,2,2)+uh(jpatch,3,3)*uh(jpatch,3,3)) 


				theta11 = dacos ( (x(1)*uh(ipatch,1,1) + x(2)*uh(ipatch,1,2) + x(3)*uh(ipatch,1,3)) / (rsqr * patch_rqsrt(ipatch,1)))
				theta12 = dacos ( (x(1)*uh(ipatch,2,1) + x(2)*uh(ipatch,2,2) + x(3)*uh(ipatch,2,3)) / (rsqr * patch_rqsrt(ipatch,2)))
				theta21 = dacos ( (-x(1)*uh(jpatch,1,1) - x(2)*uh(jpatch,1,2) - x(3)*uh(jpatch,1,3)) / (rsqr * patch_rqsrt(jpatch,1)))
				theta22 = dacos ( (-x(1)*uh(jpatch,2,1) - x(2)*uh(jpatch,2,2) - x(3)*uh(jpatch,2,3)) / (rsqr * patch_rqsrt(jpatch,2)))

				write(*,*) ipatch, jpatch, theta11*57.2958, theta12*57.2958, theta21*57.2958, theta22*57.2958

			endif
			enddo !jatom
		enddo !ipatch

	do icrowd = 1, ncrowd
		call random_number(yt); ac(icrowd,1) = yt*boxl
		call random_number(yt);	ac(icrowd,2) = yt*boxl
		call random_number(yt);	ac(icrowd,3) = yt*boxl
	enddo !i

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

!	COMPUTE THE INITIAL ENERGY OF THE SYSTEM

	OEp = 0.0; OEc = 0.0; OEpc = 0.0 !Set initial energy at zero
	do iatom=1,natom
		do jatom=iatom+1,natom                     !natom is the no of atoms in the box
				do m=1,3
					x(m)=ap(iatom,m)-ap(jatom,m) !x() is the distance between ith and jth particle
					x(m)=x(m)-boxl*anint(x(m)/boxl) !Minimum image convention
				enddo !m
			rsqr=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
			if(rsqr.lt.144) then
				rhex=rsqr*rsqr*rsqr
				rtwelve=rhex*rhex
				OEp = OEp + 4*eps11*((sgmahex11/rhex))	! eo is the interaction energy between ith and jth particle	
				!write(20,*) iatom, jatom, ap(iatom,:), ap(jatom,:)
				!write(20,*) x(:)
				!write(20,*) rsqr		
			endif
		enddo !j
	enddo !i	!OEp is the total initial energy

	do iatom=1,natom
		do jcrowd=1,ncrowd                     !natom is the no of atoms in the box
			do m=1,3
				x(m)=ap(iatom,m)-ac(jcrowd,m) !x() is the distance between ith and jth particle
				x(m)=x(m)-boxl*anint(x(m)/boxl) !Minimum image convention
			enddo !m
			rsqr=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
			if(rsqr.lt.144) then
				rhex=rsqr*rsqr*rsqr
				rtwelve=rhex*rhex
				OEc=OEc+4*eps12*((sgmahex12/rhex))	! eo is the interaction energy between ith and jth particle
			endif
		enddo !j
	enddo !i	!OEc is the total initial energy

	do icrowd=1,ncrowd
		do jcrowd=icrowd+1,ncrowd                     !natom is the no of atoms in the box
			do m=1,3
				x(m)=ac(icrowd,m)-ac(jcrowd,m) !x() is the distance between ith and jth particle
				x(m)=x(m)-boxl*anint(x(m)/boxl) !Minimum image convention
			enddo !m
			rsqr=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
			if(rsqr.lt.144) then
				rhex=rsqr*rsqr*rsqr
				rtwelve=rhex*rhex
				OEpc=OEpc+4*eps22*((sgmahex22/rhex))	! eo is the interaction energy between ith and jth particle
			endif
		enddo !j
	enddo !i	!OEpc is the total initial energy



! 	WRITE INITIAL ENERGY IN OUTPUT

	write(10,*) 'Initial energy of the system in Kj/mole'
	write(10,*) OEp, OEc, OEpc
	write(10,*) 'Initial energy of each particle in Kj/mole'
	write(10,*) OEp/natom, OEp/ncrowd

	write(10,*) '----------------------------------------'


									! ---------------- MC simulation ------------------!


	accp_crowd = 0; rej_crowd = 0
	accp_patch = 0; rej_patch = 0

	do imc = 1, nmc   ! nmc is the no of simulation
		write(*,*) imc
		write(2,*) '	6'
		write(2,*) 

		do ipatch = 1,npatch
			patch_rqsrt(ipatch,1:2) = 0
			Ei = 0; Ef = 0

			!write(2,*) ipatch
			!write(2,*) 'Xe',ap(ipatch,1),ap(ipatch,2),ap(ipatch,3)
			!write(2,*) 'P',ah(ipatch,1,1),ah(ipatch,1,2),ah(ipatch,1,3)
			!write(2,*) 'P',ah(ipatch,2,1),ah(ipatch,2,2),ah(ipatch,2,3)

!	BEFORE RANDOM DISPLACEMENT

			uh(ipatch,1,1) = (ah(ipatch,1,1)-ap(ipatch,1))
			uh(ipatch,1,2) = (ah(ipatch,1,2)-ap(ipatch,2))
			uh(ipatch,1,3) = (ah(ipatch,1,3)-ap(ipatch,3))
			
			uh(ipatch,2,1) = (ah(ipatch,2,1)-ap(ipatch,1))
			uh(ipatch,2,2) = (ah(ipatch,2,2)-ap(ipatch,2))
			uh(ipatch,2,3) = (ah(ipatch,2,3)-ap(ipatch,3))

			patch_rqsrt(ipatch,1) = sqrt(uh(ipatch,1,1)*uh(ipatch,1,1)+uh(ipatch,1,2)*uh(ipatch,1,2)+uh(ipatch,1,3)*uh(ipatch,1,3)) 
			patch_rqsrt(ipatch,2) = sqrt(uh(ipatch,2,1)*uh(ipatch,2,1)+uh(ipatch,2,2)*uh(ipatch,2,2)+uh(ipatch,3,3)*uh(ipatch,3,3)) 

			do jpatch = 1, npatch
				if (ipatch.ne.jpatch) then
				
				x(:)=ap(ipatch,:)-ap(jpatch,:) !x() is the distance between ith and jth particle
				x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention

				rsqr=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))

				uh(jpatch,1,1) = (ah(jpatch,1,1)-ap(jpatch,1))
				uh(jpatch,1,2) = (ah(jpatch,1,2)-ap(jpatch,2))
				uh(jpatch,1,3) = (ah(jpatch,1,3)-ap(jpatch,3))
			
				uh(jpatch,2,1) = (ah(jpatch,2,1)-ap(jpatch,1))
				uh(jpatch,2,2) = (ah(jpatch,2,2)-ap(jpatch,2))
				uh(jpatch,2,3) = (ah(jpatch,2,3)-ap(jpatch,3)) 

				patch_rqsrt(jpatch,1) = sqrt(uh(jpatch,1,1)*uh(jpatch,1,1)+uh(jpatch,1,2)*uh(jpatch,1,2)+uh(jpatch,1,3)*uh(jpatch,1,3)) 
				patch_rqsrt(jpatch,2) = sqrt(uh(jpatch,2,1)*uh(jpatch,2,1)+uh(jpatch,2,2)*uh(jpatch,2,2)+uh(jpatch,3,3)*uh(jpatch,3,3)) 


				theta11 = dacos ( (x(1)*uh(ipatch,1,1) + x(2)*uh(ipatch,1,2) + x(3)*uh(ipatch,1,3)) / (rsqr * patch_rqsrt(ipatch,1)))
				theta12 = dacos ( (x(1)*uh(ipatch,2,1) + x(2)*uh(ipatch,2,2) + x(3)*uh(ipatch,2,3)) / (rsqr * patch_rqsrt(ipatch,2)))
				theta21 = dacos ( (-x(1)*uh(jpatch,1,1) - x(2)*uh(jpatch,1,2) - x(3)*uh(jpatch,1,3)) / (rsqr * patch_rqsrt(jpatch,1)))
				theta22 = dacos ( (-x(1)*uh(jpatch,2,1) - x(2)*uh(jpatch,2,2) - x(3)*uh(jpatch,2,3)) / (rsqr * patch_rqsrt(jpatch,2)))

				!write(40,*) imc
				!write(40,*) 'Before'
				!write(40,*) x(:), uh(ipatch,1,:), uh(jpatch,2,:)

                                if(rsqr.gt.11.6.and.rsqr.lt.15) then
				if(theta11.le.0.24.and.theta21.le.0.24) Ei = Ei - eps11
				if(theta11.le.0.24.and.theta22.le.0.24) Ei = Ei - eps11
				if(theta12.le.0.24.and.theta21.le.0.24) Ei = Ei - eps11
				if(theta12.le.0.24.and.theta22.le.0.24) Ei = Ei - eps11
                                endif

                                if(rsqr.lt.11.6) Ei = Ei + eps11*(sgma**12/rsqr**6)

				write(20,*) imc, ipatch, jpatch, theta11*57.2958, theta12*57.2958,theta21*57.2958, theta22*57.2958,Ei
				!(x(1)*uh(ipatch,1,1) + x(2)*uh(ipatch,1,2) + x(3)*uh(ipatch,1,3)) / (rsqr * patch_rqsrt(ipatch,1)) &
				!, (x(1)*uh(ipatch,2,1) + x(2)*uh(ipatch,2,2) + x(3)*uh(ipatch,2,3)) / (rsqr * patch_rqsrt(ipatch,2)) &
				!, (-x(1)*uh(jpatch,1,1) - x(2)*uh(jpatch,1,2) - x(3)*uh(jpatch,1,3)) / (rsqr * patch_rqsrt(jpatch,1)) &
				!, (-x(1)*uh(jpatch,2,1) - x(2)*uh(jpatch,2,2) - x(3)*uh(jpatch,2,3)) / (rsqr * patch_rqsrt(jpatch,2))

				endif
			enddo !jpatch


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

			!write(*,*) detrot

			if(detrot.lt.0.9999.or.detrot.gt.1.0001) then
				write(*,*) 'blah!!','determinant is not 1',detrot
				stop
			endif

			!	 MOVE iTH PARTICLE RANDOMLY   !

			do j=1,3
				call random_number(yt)
				apn(ipatch,j)=ap(ipatch,j)+(2*yt-1)*d_m  ! move ith particle randomly
				if(apn(ipatch,j).gt.boxl) apn(ipatch,j)=apn(ipatch,j)-boxl*1.0d0
				if(apn(ipatch,j).lt.0.0d0) apn(ipatch,j)=apn(ipatch,j)+boxl*1.0d0
			enddo !j

			a_h(ipatch,1,:) = ah(ipatch,1,:) +  (apn(ipatch,:) - ap(ipatch,:))
			a_h(ipatch,2,:) = ah(ipatch,2,:) +  (apn(ipatch,:) - ap(ipatch,:))
			
			ahn(ipatch,1,1)=rot(1,1)*(a_h(ipatch,1,1)-apn(ipatch,1))+rot(1,2)*(a_h(ipatch,1,2)-apn(ipatch,2))&
			+rot(1,3)*(a_h(ipatch,1,3)-apn(ipatch,3))
			ahn(ipatch,1,2)=rot(2,1)*(a_h(ipatch,1,1)-apn(ipatch,1))+rot(2,2)*(a_h(ipatch,1,2)-apn(ipatch,2))&
			+rot(2,3)*(a_h(ipatch,1,3)-apn(ipatch,3))
			ahn(ipatch,1,3)=rot(3,1)*(a_h(ipatch,1,1)-apn(ipatch,1))+rot(3,2)*(a_h(ipatch,1,2)-apn(ipatch,2))&
			+rot(3,3)*(a_h(ipatch,1,3)-apn(ipatch,3))

			ahn(ipatch,2,1)=rot(1,1)*(a_h(ipatch,2,1)-apn(ipatch,1))+rot(1,2)*(a_h(ipatch,2,2)-apn(ipatch,2))&
			+rot(1,3)*(a_h(ipatch,2,3)-apn(ipatch,3))
			ahn(ipatch,2,2)=rot(2,1)*(a_h(ipatch,2,1)-apn(ipatch,1))+rot(2,2)*(a_h(ipatch,2,2)-apn(ipatch,2))&
			+rot(2,3)*(a_h(ipatch,2,3)-apn(ipatch,3))
			ahn(ipatch,2,3)=rot(3,1)*(a_h(ipatch,2,1)-apn(ipatch,1))+rot(3,2)*(a_h(ipatch,2,2)-apn(ipatch,2))&
			+rot(3,3)*(a_h(ipatch,2,3)-apn(ipatch,3))
 				
			ahn(ipatch,1,:) = ahn(ipatch,1,:) + apn(ipatch,:); ahn(ipatch,2,:) = ahn(ipatch,2,:) + apn(ipatch,:)

			uhn(ipatch,1,1) = (ahn(ipatch,1,1)-apn(ipatch,1))
			uhn(ipatch,1,2) = (ahn(ipatch,1,2)-apn(ipatch,2))
			uhn(ipatch,1,3) = (ahn(ipatch,1,3)-apn(ipatch,3))
			
			uhn(ipatch,2,1) = (ahn(ipatch,2,1)-apn(ipatch,1))
			uhn(ipatch,2,2) = (ahn(ipatch,2,2)-apn(ipatch,2))
			uhn(ipatch,2,3) = (ahn(ipatch,2,3)-apn(ipatch,3)) 

			patch_rqsrt(ipatch,1) = sqrt (uh(ipatch,1,1)*uh(ipatch,1,1)+uh(ipatch,1,2)*uh(ipatch,1,2)+uh(ipatch,1,3)*uh(ipatch,1,3))
			patch_rqsrt(ipatch,2) = sqrt (uh(ipatch,2,1)*uh(ipatch,2,1)+uh(ipatch,2,2)*uh(ipatch,2,2)+uh(ipatch,2,3)*uh(ipatch,2,3))

			do jpatch = 1, npatch
				if (ipatch.ne.jpatch) then
				
				x(:)=apn(ipatch,:)-ap(jpatch,:) !x() is the distance between ith and jth particle
				x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention
				
				rsqr=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))

				uh(jpatch,1,1) = (ah(jpatch,1,1)-ap(jpatch,1))
				uh(jpatch,1,2) = (ah(jpatch,1,2)-ap(jpatch,2))
				uh(jpatch,1,3) = (ah(jpatch,1,3)-ap(jpatch,3))
			
				uh(jpatch,2,1) = (ah(jpatch,2,1)-ap(jpatch,1))
				uh(jpatch,2,2) = (ah(jpatch,2,2)-ap(jpatch,2))
				uh(jpatch,2,3) = (ah(jpatch,2,3)-ap(jpatch,3)) 

				patch_rqsrt(jpatch,1) = sqrt(uh(jpatch,1,1)*uh(jpatch,1,1)+uh(jpatch,1,2)*uh(jpatch,1,2)+uh(jpatch,1,3)*uh(jpatch,1,3)) 
				patch_rqsrt(jpatch,2) = sqrt(uh(jpatch,2,1)*uh(jpatch,2,1)+uh(jpatch,2,2)*uh(jpatch,2,2)+uh(jpatch,3,3)*uh(jpatch,3,3)) 


				theta11 = dacos ((x(1)*uh(ipatch,1,1) + x(2)*uh(ipatch,1,2) + x(3)*uh(ipatch,1,3)) / (rsqr * patch_rqsrt(ipatch,1)))
				theta12 = dacos ((x(1)*uh(ipatch,2,1) + x(2)*uh(ipatch,2,2) + x(3)*uh(ipatch,2,3)) / (rsqr * patch_rqsrt(ipatch,2)))
				theta21 = dacos ((-x(1)*uh(jpatch,1,1) - x(2)*uh(jpatch,1,2) - x(3)*uh(jpatch,1,3)) / (rsqr * patch_rqsrt(jpatch,1)))
				theta22 = dacos ((-x(1)*uh(jpatch,2,1) - x(2)*uh(jpatch,2,2) - x(3)*uh(jpatch,2,3)) / (rsqr * patch_rqsrt(jpatch,2)))
				
				!write(40,*) 'After'
				!write(40,*) x(:), uh(ipatch,1,:), uh(jpatch,2,:)
        
                                if(rsqr.gt.11.6.and.rsqr.lt.15) then                               
				if(theta11.le.0.24.and.theta21.le.0.24) Ef = Ef - eps11
				if(theta11.le.0.24.and.theta22.le.0.24) Ef = Ef - eps11
				if(theta12.le.0.24.and.theta21.le.0.24) Ef = Ef - eps11
				if(theta12.le.0.24.and.theta22.le.0.24) Ef = Ef - eps11
                                endif
                                if(rsqr.lt.11.6) Ef = Ef + eps11*(sgma**12/rsqr**6)

				write(20,*) imc, ipatch, jpatch, theta11*57.2958, theta12*57.2958, theta21*57.2958, theta22*57.2958,Ef,rsqr
				!(x(1)*uh(ipatch,1,1) + x(2)*uh(ipatch,1,2) + x(3)*uh(ipatch,1,3)) / (rsqr * patch_rqsrt(ipatch,1)) &
				!, (x(1)*uh(ipatch,2,1) + x(2)*uh(ipatch,2,2) + x(3)*uh(ipatch,2,3)) / (rsqr * patch_rqsrt(ipatch,2)) &
				!, (-x(1)*uh(jpatch,1,1) - x(2)*uh(jpatch,1,2) - x(3)*uh(jpatch,1,3)) / (rsqr * patch_rqsrt(jpatch,1)) &
				!, (-x(1)*uh(jpatch,2,1) - x(2)*uh(jpatch,2,2) - x(3)*uh(jpatch,2,3)) / (rsqr * patch_rqsrt(jpatch,2))
				endif
			enddo !jatom
	
			!write(2,*) 'after displacement'
			write(2,*) 'C',apn(ipatch,1),apn(ipatch,2),apn(ipatch,3)
			write(2,*) 'P',ahn(ipatch,1,1),ahn(ipatch,1,2),ahn(ipatch,1,3)
			write(2,*) 'P',ahn(ipatch,2,1),ahn(ipatch,2,2),ahn(ipatch,2,3)
	

			write(20,*) ipatch,Ei,Ef,eps11*(sgma**12/rsqr**6)

			dE=Ef-Ei   ! dE is the difference between new and initial energy
			if (dE.gt.0) then
				call random_number(yt)
				p_r=exp(-dE/(RT)) !Kb is the Boltzmann factor and T is the temperature
			endif
			if (dE.lt.0.or.yt.lt.p_r) then
				ap(ipatch,:) = apn(ipatch,:)
				ah(ipatch,1,:) = ahn(ipatch,1,:)
				ah(ipatch,2,:) = ahn(ipatch,2,:)
				accp_patch=accp_patch+1
			else
				rej_patch=rej_patch+1
			endif
		
			!write(2,*) 'after sampling', dE, Ef, Ei, yt, p_r, accp_patch, rej_patch
			!write(2,*) 'Xe',ap(ipatch,1),ap(ipatch,2),ap(ipatch,3)
			!write(2,*) 'P',ah(ipatch,1,1),ah(ipatch,1,2),ah(ipatch,1,3)
			!write(2,*) 'P',ah(ipatch,2,1),ah(ipatch,2,2),ah(ipatch,2,3)

		enddo !ipatch

		do icrowd = 1,ncrowd

!	BEFORE RANDOM DISPLACEMENT
			do jcrowd = 1, ncrowd
			if (icrowd.ne.jcrowd) then
				do m=1,3
					x(m)=ac(icrowd,m)-ac(jcrowd,m) !x() is the distance between ith and jth particle
					x(m)=x(m)-boxl*anint(x(m)/boxl) !Minimum image convention
				enddo !m
				rsqr=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
				if(rsqr.lt.144) then
					rhex=rsqr*rsqr*rsqr
					rtwelve=rhex*rhex
					Ei = Ei + 4*eps22*((sgmahex22/rhex))	! eo is the interaction energy between ith and jth particle
				endif
			endif
			enddo !jatom

			do jatom = 1,natom                     !natom is the no of atoms in the box
				do m=1,3
					x(m)=ac(icrowd,m)-ap(jatom,m) !x() is the distance between ith and jth particle
					x(m)=x(m)-boxl*anint(x(m)/boxl) !Minimum image convention
				enddo !m
				rsqr=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
				if(rsqr.lt.144) then
					rhex=rsqr*rsqr*rsqr
					rtwelve=rhex*rhex
					Ei = Ei + 4*eps12*((sgmahex12/rhex))	! eo is the interaction energy between ith and jth particle
				endif
			enddo !jcrowd

!	 MOVE iTH PARTICLE RANDOMLY   !

			do j=1,3
				call random_number(yt)
				a_c(icrowd,j)=ac(icrowd,j)+(2*yt-1)*d_m  ! move ith particle randomly
				if(a_c(icrowd,j).gt.boxl) a_c(icrowd,j)=a_c(icrowd,j)-boxl*1.0d0
				if(a_c(icrowd,j).lt.0.0d0) a_c(icrowd,j)=a_c(icrowd,j)+boxl*1.0d0
			enddo !j

!	AFTER RANDOM DISPLACEMENT
			do jcrowd = 1, ncrowd
			if (icrowd.ne.jcrowd) then
				do m=1,3
					x(m)=a_c(icrowd,m)-ac(jcrowd,m) !x() is the distance between ith and jth particle
					x(m)=x(m)-boxl*anint(x(m)/boxl) !Minimum image convention
				enddo !m
				rsqr=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
				if(rsqr.lt.144) then
					rhex=rsqr*rsqr*rsqr
					rtwelve=rhex*rhex
					Ef = Ef + 4*eps22*((sgmahex22/rhex))	! eo is the interaction energy between ith and jth particle
				endif
			endif
			enddo !jatom

			do jatom=1,natom                     !natom is the no of atoms in the box
				do m=1,3
					x(m)=a_c(icrowd,m)-ap(jatom,m) !x() is the distance between ith and jth particle
					x(m)=x(m)-boxl*anint(x(m)/boxl) !Minimum image convention
				enddo !m
				rsqr=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
				if(rsqr.lt.144) then
					rhex=rsqr*rsqr*rsqr
					rtwelve=rhex*rhex
					Ef = Ef + 4*eps12*((sgmahex12/rhex))	! eo is the interaction energy between ith and jth particle
				endif
			enddo !jcrowd

!	 SAMPLING  !
			dE=Ef-Ei   ! dE is the difference between new and initial energy
			if (dE.gt.0) then
				call random_number(yt)
				p_r=exp(-dE/(RT)) !Kb is the Boltzmann factor and T is the temperature
			endif
			if (dE.lt.0.or.yt.lt.p_r) then
				do j=1,3
					ac(icrowd,j)=a_c(icrowd,j)	!accept the new configuration
				enddo !j
				accp_crowd=accp_crowd+1
			else
				rej_crowd=rej_crowd+1
			endif

		enddo !icrowd

	!---- write no of nth cycle,Energy at nth cycle,Energy of each particle after nth cycle  in 'energy.dat' file to plot---!

	write(11,*) imc,accp_patch/(accp_patch+rej_patch),accp_crowd/(accp_crowd+rej_crowd), d_m, dpsi_max

!        if((accp_patch/(accp_patch+rej_patch).le.0.6)) then
!                  d_m = d_m*1.05; dpsi_max = dpsi_max1.05
!        else
!                  d_m = d_m*0.95; dpsi_max = dpsi_max * 0.95
!        endif

! --------- WRITE THE TRAJECTORY AT INTERVAL OF 1000 STEPS ------------!
	if(mod(imc,10).eq.0) then
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

	Etotal_final = 0
	do ipatch = 1,npatch
		patch_rqsrt(ipatch,1:2) = 0

!	BEFORE RANDOM DISPLACEMENT

			uh(ipatch,1,1) = (ah(ipatch,1,1)-ap(ipatch,1))
			uh(ipatch,1,2) = (ah(ipatch,1,2)-ap(ipatch,2))
			uh(ipatch,1,3) = (ah(ipatch,1,3)-ap(ipatch,3))
			
			uh(ipatch,2,1) = (ah(ipatch,2,1)-ap(ipatch,1))
			uh(ipatch,2,2) = (ah(ipatch,2,2)-ap(ipatch,2))
			uh(ipatch,2,3) = (ah(ipatch,2,3)-ap(ipatch,3))

			patch_rqsrt(ipatch,1) = sqrt(uh(ipatch,1,1)*uh(ipatch,1,1)+uh(ipatch,1,2)*uh(ipatch,1,2)+uh(ipatch,1,3)*uh(ipatch,1,3)) 
			patch_rqsrt(ipatch,2) = sqrt(uh(ipatch,2,1)*uh(ipatch,2,1)+uh(ipatch,2,2)*uh(ipatch,2,2)+uh(ipatch,3,3)*uh(ipatch,3,3)) 

			do jpatch = 1, npatch
			if (ipatch.ne.jpatch) then
				
				x(:)=ap(ipatch,:)-ap(jpatch,:) !x() is the distance between ith and jth particle
				x(:)=x(:)-boxl*anint(x(:)/boxl) !Minimum image convention

				rsqr=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))

				uh(jpatch,1,1) = (ah(jpatch,1,1)-ap(jpatch,1))
				uh(jpatch,1,2) = (ah(jpatch,1,2)-ap(jpatch,2))
				uh(jpatch,1,3) = (ah(jpatch,1,3)-ap(jpatch,3))
			
				uh(jpatch,2,1) = (ah(jpatch,2,1)-ap(jpatch,1))
				uh(jpatch,2,2) = (ah(jpatch,2,2)-ap(jpatch,2))
				uh(jpatch,2,3) = (ah(jpatch,2,3)-ap(jpatch,3)) 

				patch_rqsrt(jpatch,1) = sqrt(uh(jpatch,1,1)*uh(jpatch,1,1)+uh(jpatch,1,2)*uh(jpatch,1,2)+uh(jpatch,1,3)*uh(jpatch,1,3)) 
				patch_rqsrt(jpatch,2) = sqrt(uh(jpatch,2,1)*uh(jpatch,2,1)+uh(jpatch,2,2)*uh(jpatch,2,2)+uh(jpatch,3,3)*uh(jpatch,3,3)) 


				theta11 = dacos ( (x(1)*uh(ipatch,1,1) + x(2)*uh(ipatch,1,2) + x(3)*uh(ipatch,1,3)) / (rsqr * patch_rqsrt(ipatch,1)))
				theta12 = dacos ( (x(1)*uh(ipatch,2,1) + x(2)*uh(ipatch,2,2) + x(3)*uh(ipatch,2,3)) / (rsqr * patch_rqsrt(ipatch,2)))
				theta21 = dacos ( (-x(1)*uh(jpatch,1,1) - x(2)*uh(jpatch,1,2) - x(3)*uh(jpatch,1,3)) / (rsqr * patch_rqsrt(jpatch,1)))
				theta22 = dacos ( (-x(1)*uh(jpatch,2,1) - x(2)*uh(jpatch,2,2) - x(3)*uh(jpatch,2,3)) / (rsqr * patch_rqsrt(jpatch,2)))

				if(theta11.le.1.0472.and.theta21.le.1.0472) Etotal_final = Etotal_final - eps11
				if(theta11.le.1.0472.and.theta22.le.1.0472) Etotal_final = Etotal_final - eps11
				if(theta12.le.1.0472.and.theta21.le.1.0472) Etotal_final = Etotal_final - eps11
				if(theta12.le.1.0472.and.theta22.le.1.0472) Etotal_final = Etotal_final - eps11

			endif
			enddo !jpatch

		enddo !ipatch

		write(30,*) imc, Etotal_final

	enddo !imc

	end program mc
