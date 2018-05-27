module cvgk
    use iso_c_binding
    use f95_precision
    implicit none
    contains
    
    !convergence number with normalized square summation
    real(c_double) function cvgnss(x0, dx, n) bind(c, &
      name = 'cvg_dnssv')
        integer(c_int), intent(in):: n[value]
        real(c_double), intent(in):: x0(n), dx(n)
        real(c_double):: ssx0, ssdx
        
    ! Body of dcvgnss
        ssx0 = 0.0;
        ssdx = 0.0;
        ssx0 = SUM(x0 * x0);
        ssdx = SUM(dX * dX);
        cvgnss = ssdx / (1 + ssx0);
        
        return
    end function cvgnss
    
    subroutine strain3(F, E, ty) bind(c, name='dstrain3')
        integer(c_int),intent(in):: ty[value]
        real(c_double),intent(in):: F(3,3)
        real(c_double),intent(out):: E(6)
        real(c_double):: I(3,3) = 0.0_dp, Em(3,3) = 0.0_dp
        
        I(1,1) = 1.0_dp; I(2,2) = 1.0_dp; I(3,3) = 1.0_dp;
        select case (ty)
            case (0)    !lagrange
        	    Em = (Matmul(F,transpose(F))-I)/2.0_dp;
            case (1)    !eulerian
                Em = (I-Matmul(F,transpose(F)))*2.0_dp;
            case (2)    !Infinitesimal or Engineering
                Em = (F + transpose(F))/2.0_dp - I;
            case (3)    !Engineering
                Em = (F + transpose(F))/2.0_dp - I;
        end select
        
        E(1) = Em(1,1); E(2) = Em(2,2); E(3) = Em(3,3);
        E(4) = Em(1,2); E(5) = Em(2,3); E(6) = Em(3,1);
        if (ty == 3) E(4:6) = 2.0_dp*E(4:6);
        
    end subroutine strain3
    subroutine strain2(F, E, ty) bind(c, name='dstrain2')
        integer(c_int),intent(in):: ty[value]
        real(c_double),intent(in):: F(2,2)
        real(c_double),intent(out):: E(3)
        real(c_double):: I(2,2) = 0.0_dp, Em(2,2) = 0.0_dp
        
        I(1,1) = 1.0_dp; I(2,2) = 1.0_dp;
        select case (ty)
            case (0)    !lagrange
        	    Em = (Matmul(F,transpose(F))-I)/2.0_dp;
            case (1)    !eulerian
                Em = (I-Matmul(F,transpose(F)))*2.0_dp;
            case (2)    !Infinitesimal or Engineering
                Em = (F + transpose(F))/2.0_dp - I;
            case (3)    !Engineering
                Em = (F + transpose(F))/2.0_dp - I;
        end select
        
        E(1) = Em(1,1); E(2) = Em(2,2); E(3) = Em(1,2);
        if (ty == 3) E(3) = 2.0_dp*E(3);
        
    end subroutine strain2
    
end module cvgk

module DataEx
    use iso_c_binding
    use f95_precision
    implicit none
    contains
    
    ! get a sub array from input A according to given center(x,y)
    subroutine SubArrd(A, Sub, m, n, x, y, rx, ry, dim) bind(c, name='dsubm')
        integer(c_int),intent(in):: m[value], n[value], x[value], y[value]
        integer(c_int),intent(in):: ry[value], rx[value], dim[value]
        real(c_double),intent(in):: A(dim,n,m)
        real(c_double),intent(out):: Sub(dim,-rx:rx,-ry:ry)
        integer:: br, er, bc, ec
        
        br = y-ry; er = y+ry;
        bc = x-rx; ec = x+rx;
        Sub = A(:,bc:ec,br:er);
        
    end subroutine SubArrd
    subroutine SubArrf(A, Sub, m, n, x, y, rx, ry, dim) bind(c, name='ssubm')
        integer(c_int),intent(in):: m[value], n[value], x[value], y[value]
        integer(c_int),intent(in):: ry[value], rx[value], dim[value]
        real(c_float),intent(in):: A(dim,n,m)
        real(c_float),intent(out):: Sub(dim,-rx:rx,-ry:ry)
        integer:: br, er, bc, ec
        
        br = y-ry; er = y+ry;
        bc = x-rx; ec = x+rx;
        Sub = A(:,bc:ec,br:er);
        
    end subroutine SubArrf
    
    ! get a sub array from A according to range: r0->r1, c0->c1
    subroutine SubArrd_R(A, Sub, m, n, r0, r1, c0, c1) bind(c, name = 'dsubm_r')
        integer(c_int),intent(in):: m[value], n[value]
        integer(c_int),intent(in):: r0[value], r1[value]
        integer(c_int),intent(in):: c0[value], c1[value]
        real(c_double),intent(in):: A(n,m)
        real(c_double),intent(out):: Sub(c1-c0+1, r1-r0+1)
        
        Sub = A(c0:c1, r0:r1);
        
    end subroutine SubArrd_R
    subroutine SubArrf_R(A, Sub, m, n, r0, r1, c0, c1) bind(c, name = 'ssubm_r')
        integer(c_int),intent(in):: m[value], n[value]
        integer(c_int),intent(in):: r0[value], r1[value]
        integer(c_int),intent(in):: c0[value], c1[value]
        real(c_float), intent(in):: A(n,m)
        real(c_float),intent(out):: Sub(c1-c0+1, r1-r0+1)
        
        Sub = A(c0:c1, r0:r1);
        
    end subroutine SubArrf_R
    
    ! copy a block from A to B
    subroutine copydp(A, m1, n1, B, m2, n2, r0, r1, c0, c1) bind(c, name = '_dcopy')
        integer(c_int),intent(in):: m1[value], n1[value]
        integer(c_int),intent(in):: m2[value], n2[value]
        integer(c_int),intent(in):: r0[value], r1[value]
        integer(c_int),intent(in):: c0[value], c1[value]
        real(c_double),intent(in):: A(n1,m1)
        real(c_double),intent(out):: B(n2,m2)
        
        if (r1 > m1 .OR. c1 > n1) then
            B(c0:c1,r0:r1) = A;
        elseif (r1 > m2 .OR. c1 > n2) then
            return;
        else 
            B(c0:c1,r0:r1) = A(c0:c1,r0:r1);
        endif
        
    end subroutine copydp
    ! copy a block in A to a block in B
    subroutine cpydp(A, m1, n1, B, m2, n2, r1, r2) bind(c, name = '_dcpy')
        integer(c_int),intent(in):: m1[value], n1[value]
        integer(c_int),intent(in):: m2[value], n2[value]
        integer(c_int),intent(in), optional:: r1(2,2), r2(2,2)
        real(c_double),intent(in):: A(n1,m1)
        real(c_double),intent(out):: B(n2,m2)
        integer(c_int):: ar1, ar2, ac1, ac2
        integer(c_int):: br1, br2, bc1, bc2
        
        !copy a block from A to a block in B
        if(present(r1) .AND. present(r2)) then
            Ar1 = r1(1,1); Ar2 = r1(2,1);
            Ac1 = r1(1,2); Ac2 = r1(2,2);
            Br1 = r2(1,1); Br2 = r2(2,1);
            Bc1 = r2(1,2); Bc2 = r2(2,2);
            B(Bc1:Bc2, Br1:Br2) = A(Ac1:Ac2, Ar1:Ar2);
            return;
        end if
        !copy a block from A to full B
        if(present(r1)) then
            Ar1 = r1(1,1); Ar2 = r1(2,1);
            Ac1 = r1(1,2); Ac2 = r1(2,2);
            B = A(Ac1:Ac2, Ar1:Ar2);
            return;
        endif
        !copy full A to a block of B
        if(present(r2)) then
            Br1 = r2(1,1); Br2 = r2(2,1);
            Bc1 = r2(1,2); Bc2 = r2(2,2);
            B(Bc1:Bc2, Br1:Br2) = A;
        else
            B = A;
        endif
        
    end subroutine cpydp
    subroutine copysp(A, m1, n1, B, m2, n2, r0, r1, c0, c1) bind(c, name = '_scopy')
        integer(c_int),intent(in):: m1[value], n1[value]
        integer(c_int),intent(in):: m2[value], n2[value]
        integer(c_int),intent(in):: r0[value], r1[value]
        integer(c_int),intent(in):: c0[value], c1[value]
        real(c_float),intent(in):: A(n1,m1)
        real(c_float),intent(out):: B(n2,m2)
        
        if (r1 > m1 .OR. c1 > n1) then
            B(c0:c1,r0:r1) = A;
        elseif (r1 > m2 .OR. c1 > n2) then
            return;
        else 
            B(c0:c1,r0:r1) = A(c0:c1,r0:r1);
        endif
        
    end subroutine copysp
    
    ! deconvolution by using element-wise division
    !-- A is a complex array {(r1,i1), (r2,i2), ...}
    !-- K is the kernel of FFT
    subroutine deconv(A, K, n) bind(c, name = 'deconv')
        integer(c_int),intent(in):: n[value]
        real(c_double),intent(in):: K(2,n)
        real(c_double),intent(inout)::A(2,n)
        real(c_double):: ca(2), ck(2), norm
        integer:: i
        
        !$OMP DO
        DO i = 1, n
            ca = A(:,i); ck = K(:,i);
            norm = SUM(ck*ck);
            A(1,i) = SUM(ca*ck)/norm;
            A(2,i) = (ca(2)*ck(1)-ca(1)*ck(2))/norm;
        END DO
        !$OMP END DO
        
    end subroutine deconv
    
    subroutine dpad(A, B, m, n, b1, b2, ty) bind(c,name='dpad')
        integer(c_int),intent(in):: m[value], n[value]
        integer(c_int),intent(in):: b1[value], b2[value]
        integer(c_int),intent(in), optional:: ty[value]
        real(c_double),intent(in):: A(n,m)
        real(c_double),intent(out):: B(1-b2:n+b2,1-b1:m+b1)
        integer:: i
        
        B(1:n,1:m) = A;
        if (present(ty)) then
            select case (ty)
                case (0)    !padding with zero value
            	    goto 100;
                case (1)    !padding with edge value
                    B(1-b2:0,1-b1:0) = A(1,1);
                    B(1-b2:0,m+1:m+b1) = A(1,m);
                    B(n+1:n+b2,1-b1:0) = A(n,1);
                    B(n+1:n+b2,m+1:m+b1) = A(n,m);
                    do i = 1, b1
                        B(1:n,i-b1) = A(:,1);   !left pad
                        B(1:n,m+i) = A(:,m);    !right pad
                    end do
                    do i = 1, b2
                        B(i-b2,1:m) = A(1,:);   !top pad
                        B(n+i,1:m) = A(m,:); !bottom pad
                    end do
                case default
                    goto 100;
            end select
        else
100         B(1-b2:0,:) = 0.0_dp;
            B(n+1:n+b2,:) = 0.0_dp;
            B(1:n,1-b1:0) = 0.0_dp;
            B(1:n,m+1:m+b1) = 0.0_dp;
        end if
        
    end subroutine dpad
    
    subroutine dpad_u8(Ai, B, m, n, b1, b2, ty) bind(c,name='dpad_u8')
        integer(c_int),intent(in):: m[value], n[value]
        integer(c_int),intent(in):: b1[value], b2[value]
        integer(c_int),intent(in), optional:: ty[value]
        integer(kind = 1),intent(in):: Ai(n,m)
        real(c_double),intent(out):: B(1-b2:n+b2,1-b1:m+b1)
        !integer(kind = 2) :: A(n,m)
        integer:: i, j
        
        !where (Ai < 0)
        !    B(1:n,1:m) = Ai + 256;
        !elsewhere
        !    B(1:n,1:m) = Ai;
        !end where
        !$OMP DO
        do i = 1, n
            do j = 1, m
                if (Ai(i,j) < 0) then
                    B(i,j) = Ai(i,j) + 256;
                else
                    B(i,j) = Ai(i,j);
                endif
            enddo
        enddo
        !$OMP END DO
        !B(1:n,1:m) = A;
        if (present(ty)) then
            select case (ty)
                case (0)    !padding with zero value
            	    goto 100;
                case (1)    !padding with edge value
                    B(1-b2:0,1-b1:0) = B(1,1);
                    B(1-b2:0,m+1:m+b1) = B(1,m);
                    B(n+1:n+b2,1-b1:0) = B(n,1);
                    B(n+1:n+b2,m+1:m+b1) = B(n,m);
                    do i = 1, b1
                        B(1:n,i-b1) = B(1:n,1);   !left pad
                        B(1:n,m+i) = B(1:n,m);    !right pad
                    end do
                    do i = 1, b2
                        B(i-b2,1:m) = B(1,1:m);   !top pad
                        B(n+i,1:m) = B(n,1:m);    !bottom pad
                    end do
                case default
                    goto 100;
            end select
        else
100         B(1-b2:0,:) = 0.0_dp;
            B(n+1:n+b2,:) = 0.0_dp;
            B(1:n,1-b1:0) = 0.0_dp;
            B(1:n,m+1:m+b1) = 0.0_dp;
        end if
        
    end subroutine dpad_u8
    
    subroutine dfft_col(A, m, n, isign) bind(c,name='_dfft_row')
    USE utl, ONLY : swaprow
        integer, parameter :: DPC = KIND((1.0D0,1.0D0))
        integer(c_int),intent(in):: m[value], n[value]
        integer(c_int),intent(in):: isign[value]
        complex(c_double_complex), intent(inout):: A(n,m)
        integer:: i, j, istep, im, mmax, n2
        real(c_double) :: theta, pi
        complex(c_double_complex):: temp(m), w, wp, ws
        
        pi = 3.141592653589793238462643383279502884197_dp;
        n2 = n/2; j = n2;
        !$OMP DO
        DO i = 1, n-2
            if (j > i) then
                temp = A(j+1,:);
                A(j+1,:) = A(i+1,:);
                A(i+1,:) = temp;
            end if
            im = n2;
            DO
                if (im < 2 .OR. j < im) exit;
                j = j - im;
                im = im / 2;
            ENDDO
            j = j + m;
        ENDDO
        !$OMP END DO
        mmax = 1;
        DO  !outer loop executed log2(N) times
            if (n <= mmax) exit;
            istep = 2*mmax;
            theta = pi/(isign*mmax);
            wp = cmplx(-2.0_dp*sin(0.5_dp*theta)**2, sin(theta), kind = dpc);
            w = cmplx(1.0_dp, 0.0_dp, kind = dpc);
            !$OMP DO
                DO im = 1, mmax
                    ws  = w;
                    DO i = im, n, istep
                        j = i + mmax;
                        temp = ws*A(j,:);
                        A(j,:) = A(i,:) - temp;
                        A(i,:) = A(i,:) + temp;
                    ENDDO
                    w = w*wp + w;
                ENDDO
            !$OMP END DO
            mmax = istep;
        ENDDO
    end subroutine dfft_col
    
    !--- complex operations
    subroutine dcmplx_opter(A, optr, B, n) bind(c, name = 'dcmplx')
        integer(c_int),intent(in):: n[value]
        character(c_char), intent(in):: optr[value]
        complex(c_double_complex), intent(inout):: A(n)
        complex(c_double_complex), intent(inout):: B(n)
        
        if (optr == '+') A = A + B;
        if (optr == '-') A = A - B;
        if (optr == '*') A = A * B;
        if (optr == '/') A = A / B;
        
    end subroutine dcmplx_opter
    subroutine scmplx_opter(A, optr, B, n) bind(c, name = 'scmplx')
        integer(c_int),intent(in):: n[value]
        character(c_char), intent(in):: optr[value]
        complex(c_float_complex), intent(inout):: A(n)
        complex(c_float_complex), intent(inout):: B(n)
        
        if (optr == '+') A = A + B;
        if (optr == '-') A = A - B;
        if (optr == '*') A = A * B;
        if (optr == '/') A = A / B;
        
    end subroutine scmplx_opter
    
    function dssdmd(m, n, I, avg, MASK) bind(c, name = '_dssdmd')
        integer(c_int),intent(in):: m[value], n[value]
        real(c_double),intent(in):: avg[value]
        real(c_double),intent(in):: I(n,m)
        logical(c_bool),intent(in), optional:: MASK(n,m)
        real(c_double) :: dssdmd
        
        if (present(MASK)) then
            dssdmd = SUM((I - avg)*(I - avg), MASK);
        else
            dssdmd = SUM((I - avg)*(I - avg));
        endif
        return
    end function dssdmd
    
end module DataEx