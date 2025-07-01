!-----------------------------------------------------------------------
! Subroutine: sample_latent_statuses
! Author: Caleb Frankenberger
!-----------------------------------------------------------------------

subroutine sample_latent_statuses(p, Yt_mat, Z_mat, N, SeSp, Ycols, Zrows, Zcols, U, Ztil) 
    implicit none
    
    !-- Arguments --!
    integer, intent(in) :: N        ! Number of individuals
    integer, intent(in) :: Ycols    ! Number of columns in Yt_mat
    integer, intent(in) :: Zrows    ! Number of rows in Z_mat (number of pools)
    integer, intent(in) :: Zcols    ! Number of columns in Z_mat
    integer, intent(out) :: Ztil(Zrows, Zcols) ! Holds the true pool status
    real*8, intent(in) :: p(N)                 ! Prior probability of disease for each individual
    integer, intent(inout) :: Yt_mat(N, Ycols) ! Individual status and pool membership
    integer, intent(inout) :: Z_mat(Zrows, Zcols) ! Pool test results and composition
    real*8, intent(in) :: SeSp(Zrows, 2)       ! Sensitivity and Specificity for each pool [Se, Sp]
    real*8, intent(in) :: U(N)                 ! Vector of uniform random numbers for sampling

    !-- Local Vars --!
    integer :: i, j, k, s        ! Counters and indices
    integer :: p_id                 ! Current pool ID
    integer :: p_sz                 ! Current pool size
    integer :: p_res                ! Current pool's observed test result
    integer :: num_pos              ! Number of other positive individuals in the current pool
    integer :: n_pools              ! Number of pools the current individual 'i' belongs to
    real*8 :: prob0                 ! Conditional probability the individual is negative
    real*8 :: prob1                 ! Unnormalized conditional probability the individual is positive
    real*8 :: Se, Sp                ! Sensitivity and Specificity for the current pool
    real*8 :: RSe, RSp              ! Se/Sp based probabilities for the observed test result
    
    ! Loop each individual to sample their latent status
    do i=1, N
        ! Exclude the current individual
        Yt_mat(i, 1) = 0
        
        ! Initialize likelihoods
        prob0 = 1.0D0
        prob1 = 1.0D0
        
        ! Number of pools this individual belongs to
        n_pools = Yt_mat(i, 2)
        
        ! Iterate over each pool the current individual is a member of
        do j=1, n_pools
            p_id  = Yt_mat(i, 2 + j)
            p_res = Z_mat(p_id, 1)
            p_sz  = Z_mat(p_id, 2)
            Se    = SeSp(p_id, 1)
            Sp    = SeSp(p_id, 2)
            
            ! Determine effective RSe and RSp based on observed result of current pool
            ! - If the pool tested positive:
            if (p_res == 1) then
                RSe = Se
                RSp = 1.0D0 - Sp 
            ! - If the pool tested negative:
            else
                RSe = 1.0D0 - Se  
                RSp = Sp        
            end if
            
            ! Count the other positives in the pool
            num_pos = 0
            do k=1, p_sz
                s = Z_mat(p_id, 2 + k)
                num_pos = num_pos + Yt_mat(s, 1)
            end do
            
            ! Update prob0 based off whether the pool contained any other positives
            ! - At least one other member tested positive
            if(num_pos > 0 ) then
                prob0 = prob0 * RSe
            ! - All other members are negative
            else
                prob0 = prob0 * RSp
            end if
            
            ! Update prob1
            prob1 = prob1*RSe
        end do
        
        ! Multiply by prior probs
        prob0 = (1.0 - p(i)) * prob0
        prob1 = p(i) * prob1
        
        ! Normalize
        prob0 = prob0/(prob0 + prob1)
        
        ! Sample new status for current individual
        if(U(i) > prob0) then
            Yt_mat(i, 1) = 1
        else
            Yt_mat(i, 1) = 0
        end if
        
    end do
    
   ! Construct Ztil
    do i = 1, Zrows
        p_sz = Z_mat(i, 2)
        
        ! Store the test result 
        Ztil(i, 1) = 0
        do j = 1, p_sz
            s = Z_mat(i, 2 + j)
            if (Yt_mat(s, 1) == 1) then
                Ztil(i, 1) = 1
                exit
            end if
        end do
        
        ! Store the pool size
        Ztil(i, 2) = p_sz

        ! Default the rest of the row to 0
        do j = 3, Zcols
            Ztil(i, j) = 0
        end do

        ! Now store the status of each pool member
        do j = 1, p_sz
            s = Z_mat(i, 2 + j)            ! Individual ID
            Ztil(i, 2 + j) = Yt_mat(s, 1)  ! Store their current sampled status
        end do
    end do
    
end subroutine