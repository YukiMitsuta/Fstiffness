! Copyright 2023/01/31 Yuki Mitsuta (mitsutay[at]omu.ac.jp)
!
! Distributed under terms of the MIT license.
  subroutine getFstiffness(nimage,ndim,BandWidth,kappa,image, Eholo,Fstiffness)
    ! subroutine to calculate the stress of stiffness 
    ! nimage:     the number of images
    ! ndim  :     the dimensiuon of each image
    ! BandWidth : the band width of stiffness term. 
    !             The larger this parameter is, the higher the stress of stiffness.
    ! kappa :     the spring constant of NEB
    ! image :     the array of the coordinates of images
    !             first variable is the number of image (max is nimage)
    !             second variable is the number of dimension (max is ndim)
    ! Eholo :     tangent vector of each images
    ! Fstiffness: the return array of the stress of stiffness
    integer :: j,k,nimage,ndim
    double precision :: Bandwidth,kappa
    double precision, dimension(nimage,ndim) :: image,Eholo,Eholo_vert
    double precision, dimension(nimage,ndim) :: image_inner,image_outer
    double precision, dimension(nimage,ndim) :: v_inner, v_outer, Fstiffness
    double precision :: tdeltaInner,tdeltaOuter,tdeltadelta
    call calcEholo_vert(nimage,ndim,image,Eholo,Eholo_vert)
    do k = 1, nimage
      do j = 1,ndim
        image_inner(k,j) = image(k,j) - Eholo_vert(k,j)*BandWidth*0.5
        image_outer(k,j) = image(k,j) + Eholo_vert(k,j)*BandWidth*0.5
        Fstiffness(k,j) = 0.0
        Eholo_vert(k,j) = 0.0
      end do
    end do
    do k = 2, nimage-1
      do j = 1,ndim
        v_inner(k,j) = image_inner(k,j) - image_inner(k+1,j)
        v_outer(k,j) = image_outer(k,j) - image_outer(k+1,j)
      end do
      call dot(v_inner(k,:),v_inner(k,:),ndim,tdeltainner)
      tdeltainner = sqrt(tdeltainner)
      call dot(v_outer(k,:),v_outer(k,:),ndim,tdeltaouter)
      tdeltaouter = sqrt(tdeltaouter)
      tdeltadelta = tdeltainner - tdeltaouter
      do j = 1,ndim
        Fstiffness(k,j) = Fstiffness(k,j) - Eholo_vert(k,j)*tdeltadelta*kappa*0.5
      end do
      do j = 1,ndim
        v_inner(k,j) = image_inner(k,j) - image_inner(k-1,j)
        v_outer(k,j) = image_outer(k,j) - image_outer(k-1,j)
      end do
      call dot(v_inner(k,:),v_inner(k,:),ndim,tdeltainner)
      tdeltainner = sqrt(tdeltainner)
      call dot(v_outer(k,:),v_outer(k,:),ndim,tdeltaouter)
      tdeltaouter = sqrt(tdeltaouter)
      tdeltadelta = tdeltainner - tdeltaouter
      do j = 1,ndim
        Fstiffness(k,j) = Fstiffness(k,j) - Eholo_vert(k,j)*tdeltadelta*kappa*0.5
      end do
    end do
  end subroutine getFstiffness
  subroutine calcEholo_vert(nimage,ndim,image, Eholo,Eholo_vert)  !
    ! subroutine to calculate the verical vector of the tangent vector
    ! nimage:     the number of images
    ! ndim  :     the dimensiuon of each image
    ! image :     the array of the coordinates of images
    !             first variable is the number of image (max is nimage)
    !             second variable is the number of dimension (max is ndim)
    ! Eholo :     the tangent vector of each images
    ! Eholo_vert: teh return array of the vertical vecors of the tangent vectors
    integer :: i,j,k,nimage,ndim,whileN
    double precision, dimension(nimage,ndim) :: Eholo,Eholo_vert,image
    double precision, dimension(ndim) :: tau,tau_nei,v1,v2
    logical :: findvertlist(nimage)
    double precision :: a,x,xx
    double precision :: v1taudot, v2taudot, calcEholovertTh
    calcEholovertTh = 0.001
    do i = 1,nimage
      findvertlist(i) = .false.
    end do
    do whileN = 1, 1000
      if (all(findvertlist)) cycle
      do k = 1, nimage
        if (.not.findvertlist(k)) then
          tau = Eholo(k,:)
          if (k == 1) then
            if (findvertlist(2)) then
              tau_nei = Eholo_vert(2,:)
              call dot(tau_nei,tau,ndim,x)
              call dot(tau,tau,ndim,xx)
              a = - x/xx
              do j = 1,ndim
                Eholo_vert(k,j) = a*tau(j)+tau_nei(j)
              end do
              findvertlist(k) = .true.
            end if
          else if (k == nimage) then
            if (findvertlist(nimage-1)) then
              tau_nei = Eholo_vert(nimage-1,:)
              call dot(tau_nei,tau,ndim,x)
              call dot(tau,tau,ndim,xx)
              a = - x/xx
                do j = 1,ndim
                  Eholo_vert(k,j) = a*tau(j)+tau_nei(j)
                end do
              findvertlist(k) = .true.
            end if
          else
            x = 0.0
            do j = 1,ndim
              v1(j) = image(k-1,j)-image(k,j)
              x = x + v1(j)*v1(j)
            end do
            x = sqrt(x)
            do j = 1,ndim
              v1(j) = v1(j)/x
            end do
            call dot(v1,tau,ndim,v1taudot)
            x = 0.0
            do j = 1,ndim
              v2(j) = image(k+1,j)-image(k,j)
              x = x + v2(j)*v2(j)
            end do
            x = sqrt(x)
            do j = 1,ndim
              v2(j) = v2(j)/x
            end do
            call dot(v2,tau,ndim,v2taudot)
            !print *, "v1taudot = ", v1taudot
            !print *, "v2taudot = ", v2taudot
            if (1.0-calcEholovertTh<v1taudot .and. 1.0-calcEholovertTh<v2taudot) then
              if (findvertlist(k-1)) then
                tau_nei = Eholo_vert(k-1,:)
                call dot(tau_nei,tau,ndim,x)
                call dot(tau,tau,ndim,xx)
                a = - x/xx
                do j = 1,ndim
                  Eholo_vert(k,j) = a*tau(j)+tau_nei(j)
                end do
                findvertlist(k) = .true.
              else if (findvertlist(k+1)) then
                tau_nei = Eholo_vert(k+1,:)
                call dot(tau_nei,tau,ndim,x)
                call dot(tau,tau,ndim,xx)
                a = - x/xx
                do j = 1,ndim
                  Eholo_vert(k,j) = a*tau(j)+tau_nei(j)
                end do
                findvertlist(k) = .true.
              end if
            else if (calcEholovertTh <= v1taudot) then
              call dot(v2,tau,ndim,x)
              call dot(v1,tau,ndim,xx)
              a = - x/xx
              do j = 1,ndim
                Eholo_vert(k,j) = a*v1(j) + v2(j)
              end do
              call dot(Eholo_vert(k,:),Eholo_vert(k,:),ndim,x)
              x = sqrt(x)
              do j = 1,ndim
                Eholo_vert(k,j) =  Eholo_vert(k,j)/x
              end do
              findvertlist(k) = .true.
            else if (calcEholovertTh <= v2taudot) then
              call dot(v1,tau,ndim,x)
              call dot(v2,tau,ndim,xx)
              a = - x/xx
              do j = 1,ndim
                Eholo_vert(k,j) = a*v2(j) + v1(j)
              end do
              call dot(Eholo_vert(k,:),Eholo_vert(k,:),ndim,x)
              x = sqrt(x)
              do j = 1,ndim
                Eholo_vert(k,j) =  Eholo_vert(k,j)/x
              end do
              findvertlist(k) = .true.
            else
              do j = 1,ndim
                Eholo_vert(k,j) = v1(j)
              end do
              findvertlist(k) = .true.
            end if
          end if
        end if
      end do
    end do
    do k = 2, nimage
      call dot(Eholo_vert(k-1,:),Eholo_vert(k,:),ndim,x)
      if (x < 0.0) then
        do j = 1,ndim
          Eholo_vert(k,j) = -Eholo_vert(k,j)
        end do
      end if
    end do
  end subroutine calcEholo_vert
  subroutine dot(a,b,ndim,x)  !
    ! dot of two vector
    ! a,b : the vectors to evaluate dot
    ! ndim : the dimension of the vectors
    ! x : the return of the dot calculation
    double precision :: x
    double precision, dimension(ndim) :: a,b
    integer :: i,ndim
    x = 0.0
    do i = 1,ndim
      x = x + a(i)*b(i)
    end do
  end subroutine dot
