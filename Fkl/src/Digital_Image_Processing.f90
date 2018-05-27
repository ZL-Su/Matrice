!------------------------------------------------
!--Description: a module of functions or routines
!  to impliment basic algorithms for digital image
!  processing.
!--Author: Dgelom Su
!------------------------------------------------   
module dipk
    use iso_c_binding
    use blak
    implicit none
    contains
    
    !-Purp: the mean of input image matrix
    !-Type: double precision
    real(c_double) function mean_f64c1(I, m, n) Bind(c, &
      name = '_mean_f64c1')
        integer(c_int),intent(in):: m[value], n[value]
        real(c_double),intent(in):: I(m*n)
        real(c_double) :: size
        
        size = m * n;
        mean_f64c1 = SUM(I) / size;
        
        return
    end function mean_f64c1
    real(c_double) function mean_f32c1(I, m, n) Bind(c, &
      name = '_mean_f32c1')
        integer(c_int),intent(in):: m[value], n[value]
        real(c_float),intent(in):: I(m*n)
        real(c_double) :: size
            
        size = m * n;
        mean_f32c1 = SUM(I) / size;
        
        return
    end function mean_f32c1
    
    !Purp: convelution of Gaussian:L(x,y;米,考) = G(x,y;米,考)*I(x,y)
    !Kernel: G(x,y;米,考) = 1/sqrt[2羽*sum(考^2)]exp{-[(x-米(1))^2, (y-米(2))^2][2考(1)^2, 2考(2)^2]}
    subroutine convg_u8(I, O, m, n, r, s, u) BIND(c, name = '_dconvgu')
        integer(c_int), intent(in):: m[value], n[value], r[value]
        real(c_double), intent(in), optional:: s(2), u(2)
        integer(kind=1),intent(in):: I(-r+1:n+r,-r+1:m+r)
        real(c_double),intent(out):: O(n,m)
        real(DP):: mu(2), s2(2), a, pi, x(2), kernel(-r:r,-r:r)
        integer(kind=2):: Irr(-r:r,-r:r)
        integer :: j, k, j1, j2
        
        !-- calc. gauss kernel
        mu = 0.0; if (present(u)) mu = u;
        s2 = 2.0; if (present(s)) s2 = 1.0_dp/(2*s*s);
!$OMP DO        
        DO j = -r, r
            DO k = -r, r
                x(1) = k; x(2) = j; x = x-mu; x = x*x;
                a = -dot_product(x, s2);
                kernel(k,j) = exp(a);
            ENDDO
        ENDDO
!$OMP END DO
        
        !convolution
!$OMP DO 
        DO j = 1, m !original image rows
            j1 = j-r; j2 = j+r;
            DO k = 1, n !original image cols
                Irr = I(k-r:k+r,j1:j2);
                where (Irr < 0) Irr = Irr + 256;
                O(k,j) = SUM(Irr*kernel);
            ENDDO
        ENDDO
!$OMP END DO
        pi = acos(-1.0_dp);
        a = (1.0_dp/pi)*sqrt(sum(s2));
        O = a*O;
    end subroutine convg_u8
    
    !-Purp: normalized least square cross correlation criteria, 0 for good correlation
    !-Type: uchar[in], dp[out]
    real(c_double) function dxcorru(fi,gi,n) Bind(c,name = '_dcorru')
        integer(c_int),intent(in):: n[value]
        integer(kind=1),intent(in):: fi(n,n), gi(n,n)
        real(dp):: f(n,n), g(n,n)
        real(dp):: mf, mg, sf, sg
        integer:: size
        
        f = fi; g = gi;
        where(f < 0) f = f + 256;
        where(g < 0) g = g + 256;
        !-- mean
        size = n*n;
        mf = SUM(f)/size; mg = SUM(g)/size;
        !-- variance
        f = f - mf; g = g - mg;
        sf = SUM(f*f); sg = SUM(g*g);
        sf = SQRT(sf); sg = SQRT(sg);
        !-- correlation
        f = f/sf; g = g/sg;
        f = f-g;
        dxcorru = SUM(f*f);
        return;
    end function dxcorru
    
    !-Purp: integral image O(x,y) = SUM(I(i,j)) (i =0,...,x; j=0,...,y)
    !-Type: uchar[in], int[out]
    subroutine iitgu(I, O, m, n) bind(c, name = '_iitguc')
        integer(c_int),intent(in):: m[value], n[value]
        integer(kind=1),intent(in):: I(n,m)
        integer(c_int),intent(out):: O(n,m)
        integer(c_int):: It(n,m), r, c
        
        It = I; where(It < 0) It = It + 256;
        
        O(1,1) = It(1,1);
        !$OMP DO
        DO r = 1, m
            DO c = 1, n
                O(c,r) = SUM(It(1:c,1:r));
            END DO
        END DO
        !$OMP END DO
    
    end subroutine iitgu
    
    !-Purp: reduce sampling according to input steps
    !-Type: uchar[inout]
    subroutine ureduce(I, m, n, O, h, w, inc) bind(c, name='_ucpy')
        integer(c_int),intent(in):: m[value], n[value]
        integer(c_int),intent(in):: h[value], w[value]
        integer(c_int),intent(in), optional:: inc(2)
        integer(kind=1),intent(in):: I(n,m)
        integer(kind=1),intent(out):: O(w,h)
        integer:: x_inc, y_inc
        
        O = 0;
        if (present(inc)) then
            x_inc = inc(1); y_inc = inc(2);
            if (x_inc == 0) x_inc = 1;
            if (inc(2) == 0)y_inc = x_inc;
        else
            x_inc = 1; y_inc = 1;
        endif

        O = I(1:n:x_inc,1:m:y_inc);

    end subroutine ureduce
    
    !-Purp: uniform assignment
    !-Type: uchar[inout]
    subroutine uassign(IO, m, n, val, r) bind(c, name = '_uassign')
        integer(c_int),intent(in):: m[value], n[value]
        integer(kind=1),intent(in), optional::val[value]
        integer(c_int),intent(in), optional:: r(2,2)
        integer(kind=1),intent(inout):: IO(n,m)
        
        if(present(r)) then
            if (present(val)) then
                IO(r(1,1):r(2,1), r(1,2):r(2,2)) = val;
            else 
                IO(r(1,1):r(2,1), r(1,2):r(2,2)) = 0;
            endif
        else
            if (present(val)) IO = val; 
            if (.NOT.(present(val)))IO = 0;
        endif  
    end subroutine uassign
    
    !-Purp : uchar convert to double
    subroutine dcvtuc(O, I, m, n, normlized) bind(c, name = '_dcvtuc')
        integer(c_int),intent(in):: m[value], n[value]
        integer(kind=1),intent(in):: I(n,m)
        real(c_double),intent(out):: O(n,m)
        integer(kind=2) :: temp(n,m)
        integer(c_int), optional:: normlized
        
        where (I < 0) 
            temp = I + 256;
        else where 
            temp = I;
        end where
        
        O = I;
        
        if(present(normlized)) then
            if(normlized == 0) return;
            if(normlized == 1) O = O / 255.0_dp;
        else
            O = O / 255.0_dp;
        endif
        
    end subroutine dcvtuc
    
    function dssdmu(I, m, n, avg) bind(c, name = '_dssdmu')
        integer(c_int),intent(in):: m[value], n[value]
        real(c_double),intent(in):: avg[value]
        integer(kind=1),intent(in):: I(n,m)
        integer(kind=2) :: temp(n,m)
        real(c_double) :: dssdmu
        
        where (I < 0) 
            temp = I + 256;
        else where 
            temp = I;
        end where
        
        dssdmu = SUM((temp - avg)*(temp - avg));
        
        return
    end function dssdmu
end module dipk