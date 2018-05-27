!------------------------------------------------
!--Description: a module of functions or routines
!  to execut basic operation of vector and matrix
!--Author: Dgelom Su
!------------------------------------------------   
module blak
    include "mkl_direct_call.fi"
    use f95_precision
    use iso_c_binding
    use lapack95
    use blas95
    implicit none
    contains

    !-Purp: det(a(n, n))
    !-Type: double precision
    !-Note: a is transposed if it's pased from C/C++ routine,
    !       but this doesn't change the returned determinant
    !       since det(a) = det(a^T).
    real(c_double) function ddetm(a, n) bind(c, name = '_ddetm')
      USE les, ONLY : dlu
        integer(c_int), intent(in):: n[value]
        real(c_double), intent(in):: a(n, n)
        real(c_double):: b(n, n), sign
        integer:: lda, ipiv(n), info, k
        !DIR$ ATTRIBUTES ALIGN : 32 :: b
        
        ! Body of ddet:
        !--1. exect PLU factorization by using MKL
        lda = MAX(1, n); b = a;
        info = dlu(b, ipiv, lda);
        if (info == 0) then
            ddetm = 0; return;
        endif
        
        ddetm = 1.0; sign = info;
        do k = 1, n, 1
            ddetm = ddetm * b(k, k);
        enddo
        ddetm = sign * ddetm;
        return
    end function ddetm
    real(c_float) function sdetm(a, n) bind(c, name = '_sdetm')
      USE les, ONLY : slu
        integer(c_int), intent(in):: n[value]
        real(c_float), intent(in):: a(n, n)
        real(c_float):: b(n, n), sign
        integer:: lda, ipiv(n), info, k
        !DIR$ ATTRIBUTES ALIGN : 32 :: b
        
        !--1. exect PLU factorization by using MKL
        lda = MAX(1, n); b = a;
        info = slu(b, ipiv, lda);
        if (info == 0) then
            sdetm = 0; return;
        endif
        
        sdetm = 1.0; sign = info;
        do k = 1, n, 1
            sdetm = sdetm * b(k, k);
        enddo
        sdetm = sign * sdetm;
        return
    end function sdetm

    !-Purp: dot product of vector x and y
    !-Type: double precision
    real(c_double) function ddotv(x, y, n) bind(c, name = '_ddotv')
        integer(c_int), intent(in):: n[value]
        real(c_double), intent(in):: x(n), y(n)
        
        ! Body of ddotv
        ddotv = dot_product(x, y);
        
        return
      end function ddotv
      real(c_float) function sdotv(x, y, n) bind(c, name = '_sdotv')
        integer(c_int), intent(in):: n[value]
        real(c_float), intent(in):: x(n), y(n)
        
        ! Body of ddotv
        sdotv = dot_product(x, y);
        
        return
    end function sdotv
    
    !-Purp: assign a scalar to array: a(m,n) = s
    !-Type: single precision
    subroutine sufa(a, s, m, n) bind(c, name = '_sufa') 
        integer(c_int), intent(in) :: m[value], n[value]
        real(c_float), intent(in) :: s[value]
        real(c_float), intent(inout) :: a(n, m)
        
        !-- body:
        a = s;
        
    end subroutine sufa
    !-Type: double precision
    subroutine dufa(a, s, m, n) bind(c, name = '_dufa')
        integer(c_int), intent(in) :: m[value], n[value]
        real(c_double), intent(in) :: s[value]
        real(c_double), intent(inout) :: a(n, m)
        
        !-- body:
        a = s;
        
    end subroutine dufa
    
    !-Purp: scalar times array: a ¡û s ¡Á a(m, n) 
    !-Type: single precision
    subroutine ssxa(s, a, m, n) bind(c, name = '_ssxa')
        integer(c_int), intent(in) :: m[value], n[value]
        real(c_float), intent(in) :: s[value]
        real(c_float), intent(inout) :: a(n, m)

        !-- body:
        a = s * a;

    end subroutine ssxa
    !-Type: double precision
    subroutine dsxa(s, a, m, n) bind(c, name = '_dsxa')
        integer(c_int), intent(in) :: m[value], n[value]
        real(c_double), intent(in) :: s[value]
        real(c_double), intent(inout) :: a(n, m)
        
        !-- body:
        a = s * a;
    
    end subroutine dsxa
    
    !-Purp: sum of binomial scalar times array: x ¡û p ¡Á x(m, n) + q ¡Á y(m, n)
    !-Type: double precision
    subroutine dpxqy(p, x, q, y, m, n) bind(c, name = '_dpxqy')
        integer(c_int), intent(in) :: m[value], n[value]
        real(c_double), intent(in) :: p[value], q[value]
        real(c_double), intent(in) :: y(n, m)
        real(c_double), intent(inout) :: x(n, m)
        
        !-- body:
        x = p * x + q * y;
    
    end subroutine dpxqy
    
    !- Purp: matrix-vector multiplication: x := A(m,n)x(n)
    !- Type: single and double precision
    subroutine smv(A, x, m, n) bind(c, name = '_smv')
        integer(c_int), intent(in):: m[value], n[value]
        real(c_float),  intent(in):: A(n,m)
        real(c_float),  intent(inout):: x(1,n)
      
        x = Matmul(x, A);
        
    end subroutine smv
    subroutine dmv(A, x, m, n) bind(c, name = '_dmv')
        integer(c_int), intent(in):: m[value], n[value]
        real(c_double), intent(in):: A(n, m)
        real(c_double), intent(inout):: x(1,n)
      
        x = Matmul(x, A);
        
    end subroutine dmv
    
    real(c_double) FUNCTION dvmv(y,A,x,m,n) bind(c, name = '_dvmv')
        integer(c_int), intent(in):: m[value], n[value]
        real(c_double), intent(in):: A(n, m)
        real(c_double), intent(in):: y(m), x(1,n)
        real(c_double) :: tx(1,n)
        
        tx = (matmul(x, A));
        dvmv = dot_product(y, tx(1,:));
        
        return
    END FUNCTION dvmv
    
    !- Purp: matrix multiplication: C(m,n) := A(m,k)B(k,n)
    !- Type: single and double precision
    subroutine smm(Ct, At, Bt, m, k, n) bind(c, name = '_smm')
        integer(c_int), intent(in):: m[value], k[value], n[value]
        real(c_float), intent(in):: At(k,m), Bt(n,k) !A->At, B->Bt
        real(c_float), intent(inout):: Ct(n,m)
        
        !C^T = (AB)^T = (B^T)(A^T)
        Ct = MATMUL(Bt, At);
    end subroutine smm
    subroutine dmm(Ct, At, Bt, m, k, n) bind(c, name = '_dmm')
        integer(c_int), intent(in):: m[value], k[value], n[value]
        real(c_double), intent(in):: At(k,m), Bt(n,k) !A->At, B->Bt
        real(c_double), intent(inout):: Ct(n,m)
        
        !C^T = (AB)^T = (B^T)(A^T)
        Ct = MATMUL(Bt, At);
    end subroutine dmm
    
    !- Purp: matrix multiplication: S(n,n) := \matmul{A(m,n)^T, A(m,n)}
    subroutine dmtm(S, A, m, n) bind(c, name = '_dmtm')
        integer(c_int), intent(in):: m[value], n[value]
        real(c_double), intent(in):: A(n, m)
        real(c_double), intent(out):: S(n, n)
        integer i, j
        
        !$OMP DO
        do i = 1, n
            do j = i, n
                S(i,j) = dot_product(A(i,:), A(j,:));
            enddo
        enddo
        !$OMP END DO
        
        !$OMP DO
        do j = 1, n
            S(j+1:n, j) = S(j, j+1:n);
        enddo
        !$OMP END DO
    end subroutine dmtm
    !- Purp: matrix multiplication: S(m,m) := \matmul{A(m,n), A(m,n)^T}
    subroutine dmmt(S, A, m, n) bind(c, name = '_dmmt')
        integer(c_int), intent(in):: m[value], n[value]
        real(c_double), intent(in):: A(n, m)
        real(c_double), intent(out):: S(m, m)
        integer i, j
        
        !$OMP DO
        do i = 1, m
            do j = i, m
                S(i,j) = dot_product(A(:,i), A(:,j));
            enddo
        enddo
        !$OMP END DO
        
        !$OMP DO
        do j = 1, m
            S(j+1:m, j) = S(j, j+1:m);
        enddo
        !$OMP END DO
    end subroutine dmmt
    
    !- Purp: matrix multiplication: S(n,n) := \matmul{A(m,n)^T, A(m,n)}
    subroutine smtm(S, A, m, n) bind(c, name = '_smtm')
        integer(c_int), intent(in):: m[value], n[value]
        real(c_float),  intent(in):: A(n, m)
        real(c_float),  intent(out):: S(n, n)
        integer i, j
        
        !$OMP DO
        do i = 1, n
            do j = i, n
                S(i,j) = dot_product(A(i,:), A(j,:));
            enddo
        enddo
        !$OMP END DO
        
        !$OMP DO
        do j = 1, n
            S(j+1:n, j) = S(j, j+1:n);
        enddo
        !$OMP END DO
    end subroutine smtm
    !- Purp: matrix multiplication: S(m,m) := \matmul{A(m,n), A(m,n)^T}
    subroutine smmt(S, A, m, n) bind(c, name = '_smmt')
        integer(c_int), intent(in):: m[value], n[value]
        real(c_float),  intent(in):: A(n, m)
        real(c_float),  intent(out):: S(m, m)
        integer i, j
        
        !$OMP DO
        do i = 1, m
            do j = i, m
                S(i,j) = dot_product(A(:,i), A(:,j));
            enddo
        enddo
        !$OMP END DO
        
        !$OMP DO
        do j = 1, m
            S(j+1:m, j) = S(j, j+1:m);
        enddo
        !$OMP END DO
    end subroutine smmt
    
    !- Purp: triple square matrices product
    subroutine dtplsm(Rt, At, Bt, Ct, n) bind(c, name = '_dtplsm')
        integer(c_int), intent(in) :: n[value]
        real(c_double), intent(in):: At(n,n), Bt(n,n), Ct(n,n)
        real(c_double), intent(inout):: Rt(n,n)
        
        !R^t = (ABC)^T = C^T B^T A^T
        Rt = MATMUL(Ct, MATMUL(Bt, At));
    end subroutine dtplsm
    
    !- Purp: reduce vector by addition: dredv ¡û sum(v_i^p) 
    !- Type: double precision
    real(c_double) function dredv(v, p, n) bind(c, &
      name = '_dredv')
        integer(c_int), intent(in) :: p[value], n[value]
        real(c_double), intent(in) :: v(n)
        real(c_double) :: vp(n)
        integer :: k
        
        !-- body
        !--1. calculate: v_i^p
        vp = 1;
        do k = 1, p, 1
            vp = vp * v;
        enddo
        !--2. reduce v by addition
        dredv = SUM(vp);
            
        return;
    end function dredv

    !-Purp: generic linear combinations: y ¡û SUM[x(i)*v(:,i)]
    !-Type: double precision
    subroutine dglcv(y, x, v, m, n) bind(c, &
      name = '_dglcv')
        integer(c_int), intent(in) :: m[value], n[value]
        real(c_double), intent(in) :: x(n), v(m, n)
        real(c_double), intent(inout) :: y(n)
        real(c_double) :: sv(m, n)
        !DIR$ ATTRIBUTES ALIGN : 32 :: SV
        integer i
        
        !-- body: smaller dimension first
        y = 0.0;
        If (m < n) Then
            Do i = 1, m, 1
                sv(i, :) = x * v(i, :);
            enddo
            y = SUM(sv, DIM = 2);
        Else
            Do i = 1, n, 1
                y = y + x(i) * v(:, i);
            End do
        End if
    end subroutine dglcv
    
    !-Purp: transpose of matrix
    !-Type: single precision
    subroutine strpsm(a, n) bind(c, name = '_strpsm') !square
        integer(c_int),intent(in):: n[value]
        real(c_float),intent(inout):: a(n, n);
        
        a = transpose(a);
        
    end subroutine strpsm
    subroutine strpm(a, b, m, n) bind(c, & !generic
      name = '_strpm') 
        integer(c_int),intent(in):: m[value], n[value]
        real(c_float),intent(in):: a(n,m);
        real(c_float),intent(inout):: b(m,n);
        
        b = transpose(a);
        
    end subroutine strpm
    
    !-Purp: transpose of matrix
    !-Type: double precision
    subroutine dtrpsm(a, n) bind(c, name = '_dtrpsm') !square
        integer(c_int),intent(in):: n[value]
        real(c_double),intent(inout):: a(n, n);
        
        a = transpose(a);
        
    end subroutine dtrpsm
    subroutine dtrpm(a, b, m, n) bind(c, & !generic
      name = '_dtrpm') 
        integer(c_int),intent(in):: m[value], n[value]
        real(c_double),intent(in):: a(n,m);
        real(c_double),intent(inout):: b(m,n);
        
        b = transpose(a);
        
    end subroutine dtrpm
    
    !-Purp: contraction res = AijBij
    !-- Type: double precision
    real(c_double) function dctrn(A, m1, B, m2, m, n, k) &
      bind(c, name = '_dctr')
        integer(c_int),intent(in):: m[value], n[value]
        real(c_double),intent(in):: m1[value], m2[value]
        real(c_double),intent(in):: A(m*n), B(m*n)
        integer(c_int),intent(in):: k[value]
        
        select case (k)
            case (1)
        	    dctrn = SUM(A - m1);
                return;
            case (2)
                dctrn = SUM((A - m1) * (B - m2));
                return;
        end select
        
    end function dctrn
    !-- Type: single precision
    real(c_double) function sctrn(A, m1, B, m2, m, n, k) &
      bind(c, name = '_sctr')
        integer(c_int),intent(in):: m[value], n[value]
        real(c_double),intent(in):: m1[value], m2[value]
        real(c_float),intent(in):: A(m*n), B(m*n)
        integer(c_int),intent(in):: k[value]
        
        select case (k)
            case (1)
        	    sctrn = SUM(A - m1);
                return;
            case (2)
                sctrn = SUM((A - m1) * (B - m2));
                return;
        end select
        
    end function sctrn
    !-- Type: double precision vector
    real(c_double) function dctrnv(Av, v, m, n, k) & 
        bind(c, name = '_dctrv')
        integer(c_int),intent(in):: m[value], n[value], k[value]
        real(c_double),intent(in):: v(m)
        real(c_double),intent(in):: Av(m, n)
        real(c_double):: A(m, n)
        !DIR$ ATTRIBUTES ALIGN : 32 :: A
        integer i
        
        do i = 1, n, 1
            A(:,i) = Av(:,i) - v;
        enddo
        select case (k)
            case (1)
        	    dctrnv = SUM(A);
                return;
            case (2)
                dctrnv = SUM(A * A);
                return;
        end select
        
    end function dctrnv
    
    
    !-Purp: 2-norm(Frobenius norm) of the input matrix
    !-Type: single and double precision
    real(c_float) function snorm_2(A, m, n) bind(c, name = '_snorm_2')
        integer(c_int),intent(in):: m[value], n[value]
        real(c_float), intent(in):: A(m, n)
        
        snorm_2 = NORM2(A);
        
        return;
    end function snorm_2
    !-Purp: 2-norm(Frobenius norm) of the input matrix
    !-Type: double precision
    real(c_double) function dnorm_2(A, m, n) bind(c, name = '_dnorm_2')
        integer(c_int),intent(in):: m[value], n[value]
        real(c_double), intent(in):: A(m, n)
        
        dnorm_2 = NORM2(A);
        
        return;
    end function dnorm_2
    
end module blak

!module cvg
!    use iso_c_binding
!    implicit none
!    contains
!    
!    !convergence number with normalized square summation
!    real(c_double) function cvgnss(X0, dX, n) bind(c, &
!      name = 'cvg_dnssv')
!        integer(c_int), intent(in):: n[value]
!        real(c_double), intent(in):: X0(n), dX(n)
!        real(c_double) TX0(n), TdX(n)
!        
!    ! Body of dcvgnss
!        TX0 = X0 * X0;
!        TdX = dX * dX;
!        cvgnss = SUM(TdX) / (1 + SUM(TX0));
!        
!        return
!    end function cvgnss 
!end module cvg