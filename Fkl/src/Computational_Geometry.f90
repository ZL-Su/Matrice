!------------------------------------------------
!--Description: a module of functions or routines
!  for computation geometry primative
!--Author: Dgelom Su
!------------------------------------------------
module cgp
use iso_c_binding
    implicit none
    contains
    
    !find the pedal of line through point "p" and perpendicular to line "l"
    integer(c_int) function dpedal(l, p, m) bind(c, &
      name = '_dpedal')
        integer(c_int),intent(in):: m[value]
        real(c_double),intent(in):: l(m, 2)
        real(c_double),intent(inout):: p(m)
        real(c_double):: lambda, d1(m), d2(m), sq
        
        d1 = l(:, 2) - l(:, 1);
        d2 = (l(:, 2) - p) * d1;
        d1 = d1 * d1;
        
        sq = SUM(d1);
        if (sq < 1.0D-6) then
            dpedal = -1;
            return;
        endif
        
        lambda = SUM(d2) / sq;
        d1 = lambda * l(:,1) + (1-lambda) * l(:,2);
        d2 = p - d1;
        sq = SUM(d2 * d2);
        if (sq < 1.0D-6) then
            dpedal = 0;
            return;
        else
            p = d1;
            if (lambda < 0) then
                dpedal = 2;
            else if (lambda > 0) then
                dpedal = 3;
            else
                dpedal = 1;
            endif
        endif
        
        return;
    end function dpedal
    
    !Rodrigues transformation used to convert a rotation matrix to a rotation vector or vice versa
    subroutine srodrgvm(rv, rm) bind(c, name = '_srodrgvm')
        real(c_float), intent(in):: rv(3)
        real(c_float), intent(inout):: rm(3,3)
        real(c_float):: r(3,1), rt(1,3), Im(3, 3), theta
        integer i
        
        theta = norm2(rv);
        r(:,1) = rv / theta;
        
        Im = 0.0;
        forall(i = 1:3)
            IM(i,i) = 1.0;
        end forall
        
        rm = 0.0;
        rm(2,1) = r(3,1);
        rm(1,2) = -rm(2,1);
        rm(1,3) = r(2,1);
        rm(3,1) = -rm(1,3);
        rm(3,2) = r(1,1);
        rm(2,3) = -rm(3,2);
        
        rt(1,:) = r(:,1);
        rm = cos(theta)*IM + (1-cos(theta))*Matmul(r, rt) + sin(theta)*rm;
        rm = transpose(rm);
        
    end subroutine srodrgvm
    real(c_float) function srodrgmv(rm, rv) bind(c, name = '_srodrgmv')
        real(c_float), intent(in):: rm(3,3)
        real(c_float), intent(inout):: rv(3)
        real(c_float):: rmt(3,3), theta, ath
        
        rmt = transpose(rm);
        rmt = (rmt - rm) * 0.5;
        
        rv(1) = rmt(3,2);
        rv(2) = rmt(1,3);
        rv(3) = rmt(2,1);
        
        theta = acos((rm(1,1)+rm(2,2)+rm(3,3)-1.0)*0.5);
        ath = abs(theta);
        if(ath > 1.0D-6 .or. ath /= 3.1415926) then
            rv = rv / sin(theta);
        endif
        srodrgmv = theta;
        return
    end function srodrgmv
    subroutine drodrgvm(rv, rm) bind(c, name = '_drodrgvm')
        real(c_double), intent(in):: rv(3)
        real(c_double), intent(inout):: rm(3,3)
        real(c_double):: r(3,1), rt(1,3), Im(3, 3), theta
        integer i
        
        theta = norm2(rv);
        r(:,1) = rv / theta;
        
        Im = 0.0;
        forall(i = 1:3)
            IM(i,i) = 1.0;
        end forall
        
        rm = 0.0;
        rm(2,1) = r(3,1);
        rm(1,2) = -rm(2,1);
        rm(1,3) = r(2,1);
        rm(3,1) = -rm(1,3);
        rm(3,2) = r(1,1);
        rm(2,3) = -rm(3,2);
        
        rt(1,:) = r(:,1);
        rm = cos(theta)*IM + (1-cos(theta))*Matmul(r, rt) + sin(theta)*rm;
        rm = transpose(rm);
        
    end subroutine drodrgvm
    real(c_double) function drodrgmv(rm, rv) bind(c, name = '_drodrgmv')
        real(c_double), intent(in):: rm(3,3)
        real(c_double), intent(inout):: rv(3)
        real(c_double):: rmt(3,3), theta, ath, trace, s
        
        rmt = transpose(rm);
        rmt = (rmt - rm) * 0.5;
        trace = rm(1,1)+rm(2,2)+rm(3,3);
        if (trace < -1) trace = -1;
        if (trace > 1)  trace = 1;
        theta = acos((trace - 1)*0.5);
        
        rv(1) = rmt(3,2);
        rv(2) = rmt(1,3);
        rv(3) = rmt(2,1);
        ath = abs(theta);
        if(ath > 1.0D-6 .or. ath /= 3.1415926) then
            rv = rv / sin(theta);
        endif
        drodrgmv = theta;
        return
    end function drodrgmv
    
end module cgp