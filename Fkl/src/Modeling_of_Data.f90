!------------------------------------------------
!--Description: a module of functions or routines
!  to model input data
!--Author: Dgelom Su
!------------------------------------------------
module dms
    use iso_c_binding
    use f95_precision
    implicit none
    contains
    
    !-- linear fitting: y = a1x + a2
    real(c_double) function dlftab(xy, a, sigma, n) bind(c, &
      name = '_dlftab')
        integer(c_int),intent(in):: n[value]
        real(c_double),intent(in):: xy(2,n)
        real(c_double),intent(inout):: a(2), sigma(2)
        real(c_double):: sx, sy, ax, ay
        real(c_double):: t, st2, sqr(n), sig
        integer:: i
        
        sx = sum(xy(1, :)); sy = sum(xy(2, :));
        ax = sx / n; st2 = 0.0; a = 0.0;
        do i = 1, n, 1
            t = xy(1, i) - ax;
            st2 = st2 + t * t;
            a(1) = a(1) + t*xy(2, i);
        enddo
        
        if (st2 < 1.0D-10 ) then
            a(1) = 1.0; a(2) = -ax;
            dlftab = -1.0;
            sigma = dlftab;
        else
            a(1) = a(1)/st2;
            a(2) = (sy - sx*a(1))/n;
            sigma(1) = SQRT(1.0/st2)
            sigma(2) = SQRT((1.0 + sx*sx/(n*st2))/n);
            sqr = xy(2,:) - a(1)*xy(1,:) - a(2);
            dlftab = SUM(sqr * sqr);
            if (dlftab < 1.0D-20) then
                dlftab = 0.0;
            endif
            sig = 0.0;
            if (n > 2) then
                sig = SQRT(dlftab / (n - 2));
            endif
            sigma = sig * sigma;
        endif
        
        return;
    end function dlftab
    
    !-- Quadratic surface fitting:
    !-- f = [a1,a2,a3,a4,a5,a6,a7,a8,a9,a10][x^2,y^2,z^2,xy,yz,zx,x,y,z,c]
    !-- x(1,:,:), x(2,:,:) and x(3,:,:) are x, y, and z components in order
    subroutine dquadsurft3(x, y, f, m, n) bind(c,name = '_dqsft3')
        use utl, only : dgaussj
        integer(c_int), intent(in):: m[value], n[value]
        real(c_double), intent(in):: x(3,n*m), y(3,n*m)
        real(c_double), intent(out):: f(3,10)
        real(DP):: A(m*n,10), b(m*n,3), xt(3), At(10,m*n), As(10,10)
        integer:: rx, ry, cr, cc, idx, i
        
        !construct Af = b
        A = 1.0_dp;
!$OMP DO        
        do i = 1, m*n
            xt = x(:,i);
            A(i,1:3) = xt*xt;
            A(i,7:9) = xt;
            A(i,4) = xt(1)*xt(2);
            A(i,5) = xt(2)*xt(3);
            A(i,6) = xt(3)*xt(1);
            b(i,:) = y(:,i);
        enddo
!$OMP END DO
        !sove Af = b
        At = transpose(A); As = MATMUL(At,A);
        f(1,:) = dgaussj(As, MATMUL(At,b(:,1)));
        f(2,:) = dgaussj(As, MATMUL(At,b(:,2)));
        f(3,:) = dgaussj(As, MATMUL(At,b(:,3)));
    end subroutine dquadsurft3
    
    !-- Cubic 8-neighbourhood surface fitting
    !-- f = {a1,a2,a3,a4,a5,a6,a6,a8}{1,x,y,x^2,xy,y^2,x^2y,xy^2}
    subroutine dcubsurft2(x, y, f, n) bind(c,name='_dcsft2')
        use utl, only : dgaussj
        integer(c_int),intent(in):: n[value]
        real(c_double),intent(in):: x(2,n), y(n)
        real(c_double),intent(out)::f(8)
        real(dp):: A(8,8), Ao(n,8), At(8,n), b(8)
        integer:: i
        
        Ao = 1.0_dp;
        DO i = 1, n
            Ao(i,2:3) = x(:,i);
            Ao(i,4) = x(1,i)*x(1,i);
            Ao(i,5) = x(1,i)*x(2,i);
            Ao(i,6) = x(2,i)*x(2,i);
            Ao(i,7) = Ao(i,4)*x(2,i);
            Ao(i,8) = x(1,i)*Ao(i,6);
        ENDDO
        At = TRANSPOSE(Ao);
        A = MATMUL(At,Ao); b = MATMUL(At,y);
        f = dgaussj(A,b);
    end subroutine dcubsurft2
    
end module dms