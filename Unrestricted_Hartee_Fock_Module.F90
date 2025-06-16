module plain_hartree_fock

!---------------------------------------------------------------------------------
! This module contains an in-core Hartree-Fock implementation using the talsh 
! library to provide the one and two electron integrals. 
!
! This file was part of a branch of the DIRAC relativistic computational software 
! program to practise and understand unrestricted Hartree-Fock in combination 
! with Fortran software development wihtin a team environment using GitLab.
! As result of the relativistic nature of the computational program, tensors contained
! complex parts and therefore requiring a slightly different approach. Compared to
! normal Hartree-Fock.
!---------------------------------------------------------------------------------
! Author: Melle de Groot, Master Student Chemistry, track Molecular Sciences.
! Course: Scientific Software Development with Fortran / VU University.
!---------------------------------------------------------------------------------
! 
! Starting date of the project: March 1st 2019
! End date of the project: 	April 1st 2019
!
!---------------------------------------------------------------------------------

       use exacorr_datatypes

        implicit none
        save 
        private
        public hartree_fock_driver, diagonalize_complex_matrix

        contains
        !----------------------------------------------------------------------------------------------------------------------------
        subroutine hartree_fock_driver(oneint,twoint)
	! This Hartree Fock driver constructs the HF wavefunction based upon user input
	! by calculation of the Fock and Density matrices until a self-consitent field is
	! achieved. Upon request during input of a DIIS accelerated algorithm, provided 
	! with a maximum size of the DIIS spacem during the SCF iterations, an 
	! extrapolated list of weighted Fock matrices is contstructed to speed up convergence.
	! The one- and two body tensors are already provided by the program. The user has to 
	! specify the amount of electrons, their cartesian coordinates, the maximum number of 
	! SCF iterators and the type of HF calculation: restricted or unrestricted. 
	! All results will be logged to an output file which provides information with 
	! respect to the type of calculation and whether the thresholds are met within set 
	! boundaries. 
			  
	! The program starts by allocating space for serveral tensors and will additionally
	! call the user input reader, followed by a check of input and logging of results.
	! Then when everything is in order, the SCF loop is started and altered by the fact
	! whether the DIIS algorithm has been requested in the input file. 
	! At the end all memory is deallocated and the used is provided with a formatted 
	! output textfile. The example on GitHub presents an unrestricted calculation of H2
	! with the DIIS algorithm being used. The one and two body parts have been provided 
	! by another part of the DIRAC program. 
        
         implicit none
					
	 ! One- and Two-electron Part
         type(one_el_t)                                              :: oneint 				       		! one body part (h)
         complex(8)                                                  :: twoint(:,:,:,:)					! two body part (G)
   
         ! Initialization 
         complex(8), dimension(oneint%n_spinor, oneint%n_spinor)     :: P_total, F_matrix, C                              ! Density, Fock and expansion coeficients
         complex(8), dimension(oneint%n_spinor * oneint%n_spinor)    :: DIIS_e_vector			        	  ! DIIS error vector
         complex(8), allocatable, dimension(:,:,:)                   :: DIIS_F_list                     		  ! DIIS list of extrapolated and weighted Fock matrices
         complex(8), allocatable, dimension(:,:)                     :: DIIS_e_list, DIIS_B                     	  ! DIIS error vector list, DIIS B Matrix
         real(8), allocatable, dimension(:)                          :: DIIS_weights, eigenvalues			  ! weights for DIIS algorithm
         real(8)                                                     :: norm, two_el_energy, one_el_energy                ! vector norm, one- and two electron energies
         integer                                                     :: scf_stepnumber = 0                                ! stepnumber during SCF iteration loop
         logical                                                     :: converged = .false. 				  ! SCF convergence achieved
         
         ! Input reader variables
         integer                                                     :: MAX_DIIS , MAX_ITER, N_ELECTRONS		  ! Maximum size of the DIIS space, number of iterations and electrons
         real(8)                                                     :: SCF_TRESHOLD					  ! Self consistent field threshold in Hartree
         character(8)                                                :: HF_TYPE, DIIS 					  ! Type of Hartree-Fock calculation (RHF or UHF), DIIS method (ON or OFF)
      
      
         ! Reads input
         call input_reader (HF_TYPE, N_ELECTRONS, MAX_ITER, SCF_TRESHOLD, DIIS, MAX_DIIS) 

         ! Logs output
         call logger (scf_stepnumber, (one_el_energy + two_el_energy + oneint%e_core) , norm, converged, HF_TYPE &
                     ,N_ELECTRONS, MAX_ITER,SCF_TRESHOLD, DIIS, MAX_DIIS, one_el_energy, two_el_energy, &
                     oneint%e_core, eigenvalues, oneint%n_spinor)
        
         ! Allocate memory for DIIS alorithm when requested.
         if ( trim ( DIIS ) .eq. 'DIIS-ON' ) then
            allocate(DIIS_F_list(scf_stepnumber, oneint%n_spinor, oneint%n_spinor))
            allocate(DIIS_e_list(scf_stepnumber, oneint%n_spinor * oneint%n_spinor)) 
         endif 
         
         ! SCF calculation
         do scf_stepnumber = 1, MAX_ITER
   
            ! Initializes the coefficient matrix 
            if ( scf_stepnumber .eq. 1 ) then 
               call diagonalize_complex_matrix (oneint%h_core, C, eigenvalues)        
            
            ! Diagonalizes the Fock matrix
            elseif ( scf_stepnumber .gt. 1 ) then
              call diagonalize_complex_matrix (F_matrix, C, eigenvalues)
            endif
               
            ! When unrestriced Hartree-Fock is requested 
            if ( trim ( HF_TYPE ) .eq. 'UHF' ) then 
            
               ! Construct unrestriced density and Fock matrices
               call construct_density_uhf (N_ELECTRONS, oneint%n_spinor, C, P_total)  
               call construct_fock_uhf (N_ELECTRONS, oneint%n_spinor, P_total, oneint%h_core, twoint, F_matrix &
                    , one_el_energy, two_el_energy)
            
            ! When restricted Hartree-Fock is requested                 
            elseif ( trim ( HF_TYPE ) .eq. 'RHF' ) then
            
            	 ! Construct restriced density and Fock matrices
               call construct_density_rhf (N_ELECTRONS, oneint%n_spinor, C, P_total)
               call construct_fock_rhf (N_ELECTRONS, oneint%n_spinor, P_total, oneint%h_core, twoint, F_matrix &
                    , one_el_energy, two_el_energy) 
            endif

            ! when DIIS is requested
            if  ( trim ( DIIS ) .eq. 'DIIS-ON' ) then
               call DIIS_compute_fock (oneint%n_spinor, scf_stepnumber ,P_total, MAX_DIIS &
                    , SCF_TRESHOLD, DIIS_e_list , DIIS_F_list, F_matrix, converged, norm) 
  
 		! and convergence criteria are met
                if ( converged .eqv. .true. ) then
                  
                  ! log the output 
                  call logger (scf_stepnumber, (one_el_energy + two_el_energy + oneint%e_core) , norm, converged, HF_TYPE &
                              , N_ELECTRONS, MAX_ITER, SCF_TRESHOLD, DIIS, MAX_DIIS, one_el_energy, two_el_energy, &
                              oneint%e_core, eigenvalues, oneint%n_spinor)
                  exit

               endif 

            ! when DIIS acceleration is off
            elseif ( trim ( DIIS ) .eq. 'DIIS-OFF' ) then             
               ! Check if converged
               call check_convergence (oneint%n_spinor, SCF_TRESHOLD, F_matrix, P_total &
                    , norm, DIIS_e_vector, converged)
               ! If convergence is achieved, log output
               if ( converged .eqv. .true. ) then
                  call logger (scf_stepnumber, (one_el_energy + two_el_energy + oneint%e_core) , norm, converged ,HF_TYPE &
                              ,N_ELECTRONS ,MAX_ITER ,SCF_TRESHOLD ,DIIS ,MAX_DIIS , one_el_energy, two_el_energy, &
                              oneint%e_core, eigenvalues, oneint%n_spinor)     
                  exit

               endif 

            endif
            
            ! Log output at every SCF step number until maximum number of iterations 
            call logger (scf_stepnumber, (one_el_energy + two_el_energy + oneint%e_core) , norm, converged ,HF_TYPE &
            ,N_ELECTRONS ,MAX_ITER ,SCF_TRESHOLD ,DIIS ,MAX_DIIS , one_el_energy, two_el_energy, &
            oneint%e_core, eigenvalues, oneint%n_spinor)
          
          enddo

         ! Free memory when DIIS was requested
         if ( (allocated(DIIS_F_LIST)) .AND. (allocated(DIIS_e_list)) ) then
            
            deallocate(DIIS_F_list) 
            deallocate(DIIS_e_list)

         endif

        end subroutine hartree_fock_driver

 	!-----------------------------------------------------------------------------------------------------------------------------------------------------
 				
        subroutine logger (scf_stepnumber,hf_energy, vector_norm, converged, HF_TYPE, N_ELECTRONS, &
                         MAX_ITER, SCF_TRESHOLD, DIIS, MAX_DIIS, one_el_energy, two_el_energy,V_NN, &
                           eigenvalues, N_spinor)
        
        ! This subroutine writes all the in- and output of the program to the DIRAC output file. The input parameters
        ! are: the SCF step number, the Hartree-Fock energy, the Frobenius norm of [F,P], whether the density matrix 
        ! is converged or not, the requested type of HF calculation, the number of electrons, the maximum number of
        ! SCF iterations, the SCF convergence treshold, whether the DIIS method is used, the maximum size of the DIIS 
        ! space, the one electron energy, the two electron energy, the nuclear repulsion term, the orbital energies 
        ! and the number of molecular spinors.  

        ! Input 
        integer, intent(in)           		        :: scf_stepnumber, N_ELECTRONS, MAX_ITER, MAX_DIIS, N_spinor      
        real(8), intent(in)             		:: hf_energy, vector_norm, SCF_TRESHOLD, one_el_energy, two_el_energy, V_NN ! V_NN = Nuclear Replusion
        logical, intent(in)             		:: converged
        character(8), intent(in)        		:: HF_TYPE, DIIS
        real(8), allocatable, dimension(:) 	        :: eigenvalues
        
        ! Local variables 
        integer                         		:: i  
        real(8)                         		:: x 
        
        allocate ( eigenvalues (N_spinor) )     


	! Print output when looping)
          if ( scf_stepnumber .eq. 0 ) then
             print '(/,22x,a)' ,'======================================'
             print '(22x,a)'   ,'         HARTREE-FOCK PROGRAM         '          
             print '(22x,a,/)' ,'======================================'
             print '(14x,a)'    ,'Input Reader'
             print '(14x,a)', repeat ('-',43)
             print '(14x,a,14x,a,2x,a)','Calculation type:','|', HF_TYPE
             print '(14x,a,11x,a,x,i2)','Number of electrons:','|',N_ELECTRONS
             print '(14x,a,12x,a,2x,i2)','Maximum SCF cycles:','|',MAX_ITER 
             print '(14x,a,6x,a,x,es10.3)','SCF convergence treshold:','|',SCF_TRESHOLD
             print '(14x,a,7x,a,2x,a)','Convergence accelerator:', '|',DIIS
             print '(14x,a,8x,a,2x,i2)','Size of the DIIS space:','|',MAX_DIIS
             print '(14x,a,3/)', repeat ('-',43)
          endif 
          
          if ( scf_stepnumber .gt. 0 ) then
             ! Head of scf table 
             print *, repeat ('-', 76)
             print '(1x,a,2x,a,2x,a,6x,a,6x,a,4x,a,4x,a,4x,a,4x,a,3x,a,3x,a)' &
                ,'|','Iter','|','Energy','|','[F,P]','|','|Error|','|','Conv. Acc.','|'        
             print *, repeat ('-', 76)

             print '(5x,i2,5x,f14.10,5x,es10.3,5x,es10.3,10x,a)' &
                   ,scf_stepnumber,hf_energy,vector_norm,vector_norm,'DIIS'
             print *, repeat ('-', 76)
          
          endif

          if ( converged .eqv. .true. ) then
                  
             print '(/,22x,a)', '*******************************'
             print '(22x,a)'  , '            SUCCES             '
             print '(23x,a,x,i2,x,a)','SCF CONVERGED AFTER',scf_stepnumber,'CYCLES' 
             print '(22x,a,/)', '*******************************'
             
             ! SCF ENERGY
             print '(11x,a,20x,a)'    ,'Total SCF Energy','[a.u.]'
             print '(11x,a)', repeat ('-',49)
             print '(11x,a,14x,a,2x,f19.16)','Total Energy:','|', hf_energy
             print '(11x,a)', repeat ('-',49)
             print '(/,11x,a,19x,a)'    ,'Energy Components','[a.u.]'
             print '(11x,a)', repeat ('-',49)
             print '(11x,a,7x,a,2x,f19.16)','One electron energy:','|', one_el_energy
             print '(11x,a,7x,a,2x,f19.16)','Two electron energy:','|', two_el_energy
             print '(11x,a,9x,a,2x,f19.16)','Electronic energy:','|', (one_el_energy + two_el_energy)
             print '(11x,a,9x,a,2x,f19.16)','Nuclear repulsion:','|', V_NN
             print '(11x,a)', repeat ('-',49)

             ! Print orbital energies
             print '(28x,a)', repeat ('-',16)
             print '(28x,a)' , 'ORBITAL ENERGIES'           
             print '(28x,a,/)', repeat ('-',16) 
             print '(24x,a,3x,a,10x,a)','NO','OCC','[a.u.]' 

             do i = 1, N_spinor
                 
		if ( ( i .le. N_ELECTRONS ) .AND. (trim ( HF_TYPE ) .eq. 'UHF' ) )  then
                
			x = 1.D0

                 elseif ( ( i .gt. N_ELECTRONS ) .AND. (trim ( HF_TYPE ) .eq. 'UHF' ) ) then

			x = 0.D0
                 
		elseif ( ( i .le. N_ELECTRONS / 2 ) .AND. (trim ( HF_TYPE ) .eq. 'RHF' ) ) then
                
			x = 2.D0
                 
		elseif ( ( i .gt. N_ELECTRONS / 2 ) .AND. (trim ( HF_TYPE ) .eq. 'RHF' ) ) then
                
			x = 0.D0
                endif      
               
		print '(21x,i4,3x,f6.4,6x,f9.6)',i,x, eigenvalues(i)
             
	     enddo

             print '(/,22x,a,/)'   ,'*** Terminating Normally ***'

          elseif ( ( converged .eqv. .false. ) .AND. ( scf_stepnumber .eq. MAX_ITER ) ) then
             print '(/,22x,a)', '*******************************'
             print '(22x,a)'  , '            WARNING            '
             print '(23x,a)'  , '      SCF DID NOT CONVERGE     '
             print '(22x,a,/)', '*******************************'
             print '(22x,a,/)'   ,'*** Terminating with Errors ***'
          endif

        deallocate(eigenvalues)
       
        end subroutine logger 
        
        !----------------------------------------------------------------------------------------------------------------------------------------------------
        subroutine input_reader (HF_TYPE, N_ELECTRONS, MAX_ITER, SCF_TRESHOLD, DIIS, MAX_DIIS)
                       
        ! This subroutine reads input from a textfile supplied by the user of the program. Input variables are:
        ! 1. the type of Hartree Fock calculation 
        ! 2. number of electrons 
        ! 3. maximum number of SCF cycles 
        ! 4. the final SCF convergence criteria 
        ! 5. whether the DIIS algorithm is used and if so:
        ! 6. the maximum number of DIIS vectors to be used. 
       
        character(8) :: HF_TYPE, DIIS
        integer      :: N_ELECTRONS, MAX_ITER, MAX_DIIS
        real (8)     :: SCF_TRESHOLD

        open (9 , file = "input.txt", status = "unknown", form = "formatted") !, IOSTAT = ios, ERR = err  )
        
        read (9 ,* )  HF_TYPE        ! Calculation type: RHF/UHF
        read (9 ,* )  N_ELECTRONS    ! Number of electrons
        read (9 ,* )  MAX_ITER	     ! Maximum number of iterations
        read (9 ,* )  SCF_TRESHOLD   ! SCF Treshold value in Hartree 
        read (9 ,* )  DIIS, MAX_DIIS ! DIIS accelation and size of the DIIS space
             
        close (9) 
        
        end subroutine input_reader 
        
        !----------------------------------------------------------------------------------------------------------------------------------------------------
        subroutine construct_fock_rhf(N_elec, N_spinor, P_total, h_matrix, G_tensor, F_matrix, one_e_energy, two_e_energy)

        ! This subroutine constructs the Fock matrix when a restricted calculation is requested, taking the number
        ! of electrons, the number of molecular spinors, the total density matrix (P), the one-body part and two electron
        ! integrals as input. As output it provides the Fock matrix (F), the one- and two electron energy.
        
        ! Input
        integer, intent(in) 								:: N_spinor, N_elec   ! Number of spinors
        complex(8), dimension (N_spinor, N_spinor), intent(in) 				:: P_total, h_matrix  ! Density matrix and one-body part
        complex(8), dimension (N_spinor, N_spinor, N_spinor, N_spinor), intent(in) 	:: G_tensor   	      ! two-body, fourth order tensor (G)

        ! Intermediate variables
        integer :: i, j, k, l
        complex(8), dimension (N_spinor, N_spinor) 					:: G_matrix  	       ! Contracted G-tensor with Density Matrix (P)
        complex(8), PARAMETER 							        :: ZERO = (0.D0, 0.D0)   ! Parameter for ZERO

        ! Output
        complex(8), dimension (N_spinor, N_spinor), intent(out) 			:: F_matrix    	       ! Fock matrix
        real(8), intent(out) 							        :: one_e_energy   ! One electron energy
        real(8), intent(out)						                :: two_e_energy   ! Two electron energy

        ! nullify G_matrix
        G_matrix = ZERO

        ! Contraction of the G-tensor with the total density matrix (P) in column-major order.
        do l = 1, N_elec / 2
        
           do k = 1, N_elec / 2
           
              do j = 1, N_spinor
              
                 do i = 1, N_spinor
                 
                    G_matrix(i,j) = G_matrix(i,j) + P_total(k,l) * ( G_tensor(i,j,k,l) - 0.5D0 * G_tensor(i,l,k,j) ) 
                 
                 enddo
              
              enddo
           
           enddo
        
        enddo

        ! Formation of the Fock matrix by matrix addititon of contracted 2nd order G-tensor and one body part
        F_matrix = h_matrix + G_matrix

        ! nullify tne electron energy
        one_e_energy = ZERO

        ! Calculate one electron energy in column-major order.
        do j = 1, N_elec / 2 
  
           do i = 1, N_elec / 2
              
              one_e_energy = one_e_energy + P_total(i,j) * h_matrix(i,j)
           
           enddo
        
        enddo

        ! nullify two electron energy
        two_e_energy = ZERO

        ! Calculate two electron energy in column major order.
        do j = 1, N_elec / 2
        
           do i = 1, N_elec / 2
           
              two_e_energy = two_e_energy + 0.5D0 * P_total(i,j) * G_matrix(i,j)
              
           enddo
           
        enddo
        
        end subroutine construct_fock_rhf

        !----------------------------------------------------------------------------------------------------------------------------------------------------
        subroutine construct_fock_uhf (N_elec, N_spinor, P_total, h_matrix, G_tensor, F_matrix , one_e_energy, two_e_energy)
       
        ! This subroutine constructs the Fock matrix when an unrestricted calculation is requested, taking the number
        ! of electrons, the number of molecular spinors, the total density matrix, the one-body part and two electron
        ! integrals as input. As output it provides the Fock matrix, the one electron energy and two electron energy. 
       
        ! Input
        integer, intent(in) 								:: N_spinor, N_elec   ! Number of spinors and electrons
        complex(8), dimension (N_spinor, N_spinor), intent(in) 				:: P_total, h_matrix  ! Density matrix and one-body
        complex(8), dimension (N_spinor ,N_spinor, N_spinor, N_spinor), intent(in) 	:: G_tensor   	      ! two-body, fourth order tensor (G)

        ! Intermediate variables
        integer :: i, j, k, l 
        complex(8), dimension (N_spinor, N_spinor) 			                :: G_matrix  	            ! Contracted G-tensor with density matrix (P)
        complex(8), PARAMETER 							        :: ZERO = (0.D0, 0.D0)      ! Parameter ZERO
      
        ! Output  
        complex(8), dimension (N_spinor, N_spinor), intent(out) 	                :: F_matrix     ! Fock matrix (F)
        real(8), intent(out) 								:: one_e_energy ! One electron energy
        real(8), intent(out) 						                :: two_e_energy ! Two electron energy  

        ! nullify G_matrix
        G_matrix = ZERO
        
        ! Contraction of the 4rh order G-tensor with the total density matrix (P) in column-major order.
        do l = 1, N_elec
         
           do k = 1, N_elec
            
              do j = 1, N_spinor
               
                 do i = 1, N_spinor
              
                    G_matrix(i,j) = G_matrix(i,j) + P_total(k,l) * (G_tensor(i,j,k,l) - G_tensor(i,l,k,j))
              
                 enddo 
              
              enddo 
           
           enddo 
        
        enddo  
        
        ! Formation of the Fock matrix by matrix addititon of contracted 2nd order G-tensor and one body part
        F_matrix = h_matrix + G_matrix     
        
        ! nullify one electron energy
        one_e_energy = ZERO

	! Calculate one electron energy in column major order.
        do j = 1, N_elec
           
           do i = 1, N_elec
              
              one_e_energy = one_e_energy + P_total(i,j) * h_matrix(i,j)
           
           enddo
        
        enddo

        ! nullify two electron energy 
        two_e_energy = ZERO

	! Calculate two electron energy in column major order.
        do j = 1, N_elec
        
           do i = 1, N_elec
              
              two_e_energy = two_e_energy + 0.5D0 * P_total(i,j) * G_matrix(i,j)
              
           enddo
           
        enddo

        end subroutine construct_fock_uhf

        !----------------------------------------------------------------------------------------------------------------------------------------------------
        subroutine density_initial (N_elec, N_spinor, P_total)
        
        ! This subroutine constructs the initial RHF density matrix (P), taking the total
        ! number of electrons and spinors as input. 
        
        ! Input
        integer, intent(in) 								:: N_elec     ! Number of electrons
        integer, intent(in) 								:: N_spinor   ! Number of spinors
        
        ! Intermediate variables
        integer :: i, j                             
        complex(8), PARAMETER 						                :: ZERO = (0.D0, 0.D0)   ! Parameter ZERO
        complex(8), PARAMETER 							        :: TWO = (2.D0, 0.D0)    ! Parameter TWO
        
        ! Output
        complex(8), dimension(N_spinor, N_spinor), intent(out)                          :: P_total   ! Density matrix (P)

	! nullify density matrix
        P_total = ZERO 
        
        ! Fill initial density matrix diagonally with number of electrons (RHF)
        do i = 1, N_elec / 2
        
           P_total (i,i) = TWO 
        
        enddo
        
        end subroutine density_initial
        
        !----------------------------------------------------------------------------------------------------------------------------------------------------
        subroutine construct_density_rhf (N_elec, N_spinor, C, P_total)

        ! This subroutine constructs the RHF density matrix and takes the
        ! expansion coefficients of the diagonalized Fock matrix as input.

        ! Input
        integer, intent(in) 									:: N_elec, N_spinor	! Number of electrons, spinors
        complex(8), dimension(N_spinor, N_spinor), intent(in) 	                                :: C           		! Coeffient matrix from diagonal Fock matrix

        ! Intermediate variables
        integer :: i, j, k	
        complex(8), PARAMETER 					            			:: ZERO = (0.D0, 0.D0)   ! Parameter ZERO
	
        ! Output
        complex(8), dimension(N_spinor, N_spinor), intent(out)                                  :: P_total   ! Density matrix (P)

        ! nullify density matrix 
        P_total = ZERO               

        ! Construction of the RHF total density matrix in column major order with double precision.
        do k = 1, N_elec / 2 
           
           do j = 1, N_spinor
              
              do i = 1, N_spinor
                 
                 P_total(i,j) = P_total(i,j) + 2.D0 * C(i,k) * DCONJG(C(j,k))  
        
              enddo
        
           enddo
        
        enddo

        end subroutine construct_density_rhf
        
        !----------------------------------------------------------------------------------------------------------------------------------------------------
        subroutine construct_density_uhf (N_elec, N_spinor, C, P_total)
        ! This subroutine constructs the UHF density matrix and takes the
        ! expansion coefficients of the diagonalized Fock matrix as input. 
        
        ! Input
        integer, intent(in) 						:: N_elec, N_spinor     ! Number of electrons, spinors 
        complex(8), dimension(N_spinor, N_spinor), intent(in)   	:: C            	! Coeffient matrix from diagonal Fock matrix  
        
        ! Intermediate variables 
        integer :: i, j, k          
        complex(8), PARAMETER 						:: ZERO = (0.D0, 0.D0)   	! Parameter ZERO  
        
        ! Output 
        complex(8), dimension(N_spinor, N_spinor), intent(out)  	:: P_total ! Density matrix (P) 
       
        ! nullify density matrix 
        P_total = ZERO 
        
        ! Construction of the UHF density matrix in column major order with double precision. 
        do k = 1, N_elec 
        
           do j = 1, N_spinor 
        
              do i = 1, N_spinor               
                 
                 P_total(i,j) = P_total(i,j) + C(i,k) * DCONJG(C(j,k))
              
              enddo 
           
           enddo 
        
        enddo 
        
        end subroutine construct_density_uhf
       
        !---------------------------------------------------------------------------------------------------------------------------------------------------- 
        subroutine check_convergence (N_spinor, TRESHOLD, F, P, vector_norm, error_vector, converged)
        
        ! This subroutine checks whether the new density matrix (P) is converged by taking the vector norm of
        ! commutator of the Fock and previous density matrix [F,P]. It requires as input: a treshold value,
        ! the number of spinors, the Fock matrix (F), and the density matrix (P). 
        ! In addition, the commutator matrix is reshaped into an 1D Error array for the DIIS algorithm.  

        ! Input 
        integer, intent(in)                                  	:: N_spinor ! Number of spinors
        real(8), intent(in)                                  	:: TRESHOLD ! Convergence treshold in Hartree 
        complex(8), dimension(N_spinor, N_spinor), intent(in) 	:: F        ! Fock matrix 
        complex(8), dimension(N_spinor, N_spinor), intent(in) 	:: P        ! Total density matrix (P, n-1)
                 
        ! Intermediate variables 
        integer                                               	:: i, j                                 
        complex(8), dimension(N_spinor, N_spinor)              	:: commutator ! Fock and Density matrices [F(n),P(n-1)](ij)
        
        ! Output 
        complex(8), dimension(N_spinor * N_spinor), intent(out) :: error_vector      ! The error vector
        real(8), intent(out)                                  	:: vector_norm       ! |Error| or Frobenius norm
        logical, intent(out)                                  	:: converged         ! SCF converged or not

	! nullify convergence and Frobenius norm
        converged = .false.
        vector_norm = 0.D0
        
        ! Calculate the commutatior of the Fock and Density matrices [F,P]  
        commutator = MATMUL(F,P) - MATMUL(P,F)    
       
        ! Reshape commutator into a 1D array for DIIS
        error_vector = reshape(commutator, (/N_spinor * N_spinor/)) 

        ! Calculate the Frobenius norm by squaring the commutator [F(n), P(n-1)](ij) in column major order with double precision
        do j = 1, N_spinor 
           
           do i = 1, N_spinor 
              
              vector_norm = vector_norm + DSQRT((real(commutator(i,j), 8))**2 + (DIMAG(commutator(i,j)))**2) 
           
           enddo 
        
        enddo    

        ! Check for SCF Convergence
        if ( vector_norm < TRESHOLD ) then

               converged = .true.

        endif 
        
        end subroutine check_convergence 
        
        !----------------------------------------------------------------------------------------------------------------------------------------------------
        ! Start of DIIS algorithm implementation.
        ! Sources:
        !   DIIS procedure implemented by R. Meli: https://github.com/RMeli/Hartree-Fock 
        !   Helgaker, P. Jørgensen and J. Olsen - Molecular Electronic-Structure Theory, 2000
        !----------------------------------------------------------------------------------------------------------------------------------------------------
      
        subroutine DIIS_add_vector(N_spinor,scf_step,e_vector, MAX_DIIS_SPACE, e_list)
        
        ! This subroutine adds the constucted error vector (e_vector) to an error vector list (e_list). 
        
        ! Input 
        integer, intent(in) 				        	  :: N_spinor   		! Number of spinors
        integer, intent(in) 						  :: scf_step   		! Current SCF step 
        complex(8), dimension(N_spinor*N_spinor), intent(in) 	          :: e_vector 			! 1D error array 
        integer, intent(in) 				                  :: MAX_DIIS_SPACE             ! Maximum size of the DIIS space

        ! Intermediate
        integer                                                           :: i
        complex(8), allocatable, dimension(:,:) 	       		  :: copy_e_list 	   ! Temp. copy of complex error list 
        
        ! I/O
        complex(8), allocatable, dimension(:,:), intent(inout) 	          :: e_list        	   ! Complex error list 
        
        ! If size DIIS space not met, allocate memory
        if ( scf_step .le. MAX_DIIS_SPACE ) then

              ! Allocate temporary list
              allocate(copy_e_list(scf_step, N_spinor * N_spinor))
           
              ! Copy the old error list into temporary
              do i = 1, scf_step - 1
              
                 copy_e_list (i, :) = e_list (i, :)  
                 
              enddo 
           
              ! Add the new vector to the copy of the error list.
              copy_e_list(scf_step, :) = e_vector 
           
              ! Update error list and free memory 
              deallocate(e_list)
              call move_alloc(copy_e_list,e_list)
				
        ! If size DIIS space is met:        
        ! Discard the oldest vector and add new one to the error list  
        elseif ( scf_step .gt. MAX_DIIS_SPACE ) then 

           allocate(copy_e_list (MAX_DIIS_SPACE, N_spinor * N_spinor ) )
              
              ! Discard the oldest vector while copying
              copy_e_list ( 1 : ( MAX_DIIS_SPACE - 1 ), : ) = e_list ( 2 : MAX_DIIS_SPACE , : )
              
              ! Add the newest vector 
              copy_e_list ( MAX_DIIS_SPACE, : ) = e_vector

              ! Update the error list and free memory
              deallocate (e_list)
              call move_alloc(copy_e_list,e_list)         
                
        endif

        end subroutine DIIS_add_vector
        
        !----------------------------------------------------------------------------------------------------------------------------------------------------
        subroutine DIIS_add_Fock (N_spinor, scf_step, F, MAX_DIIS_SPACE, Fock_list)
        
        ! This subroutine adds the current Fock matrix to a list of Fock matrices.
        
        ! input 
        integer, intent(in) 							:: N_spinor, scf_step ! Number spinors and scf stepnumber
        complex(8), dimension(N_spinor,N_spinor), intent(in) 			:: F                  ! Current Fock matrix 
        integer, intent(in) :: MAX_DIIS_SPACE 

        ! Intermediate
        integer 			 					:: i
        complex(8), allocatable, dimension(:,:,:) 		        	:: copy_Fock_list      ! copy of list of Fock matrices 
        
        ! Output 
        complex(8), allocatable, dimension(:,:,:), intent(inout) 	        :: Fock_list 	       ! List of extrapolated weighted Fock matrices 
        
        ! if size DIIS space not met
        if ( scf_step .le. MAX_DIIS_SPACE ) then    

              ! Allocate temporary list 
              allocate(copy_Fock_list(scf_step,N_spinor,N_spinor))  
           
              ! copy the old Fock list into temporary
              do i = 1, scf_step - 1

                 copy_Fock_list(i,:,:) = Fock_list (i,:,:)  

              enddo 
           
              ! Copy the new Fock matrix into the temp. Fock list 
              copy_Fock_list(scf_step,:,:) = F  
           
              ! Update the Fock list and free memory 
              deallocate(Fock_list)
              call move_alloc(copy_Fock_list,Fock_list)
        
        ! Size DIIS space met:
        elseif ( scf_step .gt. MAX_DIIS_SPACE ) then    

           allocate (copy_Fock_list(MAX_DIIS_SPACE,N_spinor,N_spinor))
           
           ! Replace the oldest Fock matrix
           copy_Fock_list ( 1 : ( MAX_DIIS_SPACE - 1 ), : , : ) = Fock_list ( 2 : ( MAX_DIIS_SPACE ) , : , : )

           copy_Fock_list (MAX_DIIS_SPACE ,: ,: ) = F

           deallocate(Fock_list)
           call move_alloc(copy_Fock_list,Fock_list)

        endif

        end subroutine DIIS_add_Fock
        
        !----------------------------------------------------------------------------------------------------------------------------------------------------
        subroutine DIIS_construct_B(N_spinor,diis_space,e_list,B)
        
        ! This subroutine constructs the B matrix of the DIIS algorithm.

        ! input
        integer, intent(in)                                             	:: N_spinor , diis_space
        complex(8), dimension(diis_space, N_spinor * N_spinor), intent(in) 	:: e_list   

        ! intermediate 
        integer 								:: i, j

        ! output
        complex(8), dimension(diis_space+1, diis_space+1), intent(out)    	:: B ! matrix B

        ! Calculate B matrix in column major bij taking the dot products of the error lists in column major order.
        do j = 1, diis_space
          
           do i = 1, diis_space 
          
              B (i ,j ) = dot_product( e_list (i , : ) , e_list( j , : ) )          ! check if this is correct  
          
           enddo 
       
        enddo
        
        ! Fill rest of matrix B 

        B (: , diis_space + 1 ) 					 = -1.D0
        B  (diis_space + 1 ,: )   					 = -1.D0
        B (diis_space + 1 ,diis_space + 1 )                              = 0.0D0

        end subroutine DIIS_construct_B 

        !----------------------------------------------------------------------------------------------------------------------------------------------------
        subroutine DIIS_compute_weights (diis_space,matrix_B,weights,info) 
        
        ! input
        integer, intent(in)                                          		:: diis_space  ! number of DIIS vectors
        complex(8), dimension(diis_space + 1 ,diis_space + 1 ), intent(in) 	:: matrix_B  ! Matrix B        
        
        ! intermediate
        real(8), allocatable, dimension(:,:)                         		:: reordered_RHS ! B * SOL = RHS
        real(8), allocatable, dimension(:,:,:)                      	        :: reordered_B
        real(8), parameter                                          	        :: ZERO = (0.D0), MINUS_ONE = (-1.D0)
        integer, dimension(diis_space+1)                             		:: IPIV ! pivot indices

        ! output
        real(8), dimension(diis_space), intent(out)                 		:: weights ! weights 
        integer, intent(out)                                         	        :: info ! DGESV
        
        allocate (reordered_RHS(diis_space+1,2))
        allocate (reordered_B(diis_space+1,diis_space+1,2))

        ! Fill right hand side of linear system with real values
        reordered_RHS(: ,: )              = ZERO 
        reordered_RHS(diis_space + 1 ,1 ) = MINUS_ONE  ! real part of right-hand-side = -1

        ! Reorder with DP 
        reordered_B(: ,: ,1 )             = REAL(matrix_B,8)
        reordered_B(: ,: ,2 )             = DIMAG(matrix_B)  ! complex part not used

        ! Solve linear system: B * solution = RHS (using only the real parts and RHS is replaced by sol).  
        call DGESV (diis_space + 1 ,1 ,reordered_B(: ,: ,1 ) ,diis_space + 1 ,IPIV ,reordered_RHS(: ,1 ),diis_space + 1 ,info )
                
        ! Discard the lagrange multiplier 
        weights = reordered_RHS(1 : diis_space, 1)
        
        ! Free memory
        deallocate(reordered_RHS)
        deallocate(reordered_B)

        end subroutine DIIS_compute_weights

        !----------------------------------------------------------------------------------------------------------------------------------------------------
        subroutine DIIS_compute_fock (N_spinor,DIIS_space,P,MAX_DIIS_SPACE,TRESHOLD,e_list,Fock_list,F, converged,vector_norm)
        ! This subroutine construct an extrapolated Fock matrix using the DIIS algorithm.
 
        ! input 
        integer, intent(in)                                		:: N_spinor		! Size of MOs and DIIS spaces
        complex(8), dimension(N_spinor, N_spinor), intent (in) 		:: P    		! Density matrix 
        integer, intent (in)                                  		:: MAX_DIIS_SPACE       ! maximum number of DIIS vectors
        real(8), intent (in)                                  		:: TRESHOLD             ! treshold value for convergence in Hartree  

        ! intermediate 
        complex(8), dimension(N_spinor*N_spinor)              		:: error_vector 
        real(8), allocatable, dimension(:)                    		:: weights 
        complex(8), allocatable, dimension(:,:)               		:: matrix_B 
        integer                                               		:: i, j, info = -1
        real(8)                                               		:: vector_norm          ! Frobenius norm
        logical                                               		:: converged
        
        ! i/o 
        complex(8), allocatable, dimension(:,:), intent(inout)    	:: e_list               ! error list  
        complex(8), allocatable, dimension(:,:,:), intent(inout)  	:: Fock_list            ! list of extrapolated and weighted Fock matrices 
        complex(8), dimension(N_spinor, N_spinor), intent(inout)  	:: F                    ! Fock matrix
        integer, intent(inout) :: DIIS_space
        
        ! Allocate arrays    
        if ( DIIS_space .gt. MAX_DIIS_SPACE ) then

           allocate(matrix_B(MAX_DIIS_SPACE+1,MAX_DIIS_SPACE+1))
           allocate(weights(MAX_DIIS_SPACE))      
        
        elseif ( DIIS_space .le. MAX_DIIS_SPACE ) then 

           allocate(matrix_B(DIIS_space+1,DIIS_space+1))
           allocate(weights(DIIS_space))
        
        endif 

        ! Check convergence
        call check_convergence (N_spinor, TRESHOLD, F, P, vector_norm, error_vector, converged)            
        
        ! Add error vector to error list
        call DIIS_add_vector(N_spinor,DIIS_space,error_vector, MAX_DIIS_SPACE, e_list) 
        
        ! Add Fock matrix to Fock matrix list
        call DIIS_add_Fock(N_spinor,DIIS_space,F,MAX_DIIS_SPACE ,Fock_list)
        
        ! Take care of size of DIIS space
        if ( DIIS_space .gt. MAX_DIIS_SPACE ) then

           DIIS_space = MAX_DIIS_SPACE
           
        endif 

        ! Calculate B matrix 
        call DIIS_construct_B(N_spinor, DIIS_space, e_list, matrix_B)
        ! Compute the weights
        call DIIS_compute_weights(DIIS_space, matrix_B, weights, info) 
        
        ! nullify Fock Matrix          
        F = (0.D0, 0.D0)
        
        ! Generate weighted list of Fock matrices according to max DIIS space set
        do i = 1, DIIS_space
         
           F = F + DCMPLX ( weights(i) ) * Fock_list (i, : , : ) ! temp. convert real weights to complex(8)  
        
        enddo 
       
        ! Free memory
        deallocate(matrix_B)
        deallocate(weights)  
 
        end subroutine DIIS_compute_fock

        !----------------------------------------------------------------------------------------------------------------------------------------------------  
        
        subroutine diagonalize_complex_matrix (the_matrix, C, eigenvalues)

!       Illustrates the diagonalization of a complex Hermitian matrix
!       Note that this routines requires the matrix in a different
!       format, with the real/imaginary dimension as the last instead
!       of the first.

        complex(8), intent(in) 			:: the_matrix(:,:)
        real(8), allocatable 			:: reordered_matrix(:,:,:)
        real(8), allocatable 			:: eigenvectors(:,:,:)
        real(8), allocatable, intent(out) 	:: eigenvalues(:)
        complex(8), intent(out) 		:: C(:,:)
       
        integer 				:: ierr, i, n
        
        n = size(the_matrix ,1 )

        if (size(the_matrix,2) /= n) call quit ('Matrix has to be square')

        ! allocate the required memory
        allocate( reordered_matrix(n,n,2) )
        allocate( eigenvectors(n,n,2) )
        allocate( eigenvalues(n) )

        reordered_matrix(:,:,1) = REAL(the_matrix,8)
        reordered_matrix(:,:,2) = DIMAG(the_matrix)
        call qdiag90(2,n,reordered_matrix,n,n,eigenvalues,1,eigenvectors,n,n,ierr)

        if (ierr.ne.0) then
            call quit('qdiag90 failed ')
        endif  
        
        ! Reorder with DP
        C = CMPLX(eigenvectors(:,:,1), eigenvectors(:,:,2),8)

        ! free the memory
        deallocate( reordered_matrix )
        deallocate( eigenvectors )
        deallocate( eigenvalues )

        end subroutine diagonalize_complex_matrix
    
end module plain_hartree_fock
