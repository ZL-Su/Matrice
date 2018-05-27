module utl
    include "mkl_direct_call.fi"
    use f95_precision
    use iso_c_binding
    use lapack95
    implicit none
    
    INTERFACE outer_prod
        module procedure douterpd
        module procedure souterpd
    END INTERFACE outer_prod
    INTERFACE pythag
        module procedure dpythag
        module procedure spythag
    END INTERFACE pythag
    INTERFACE deswap    !descending swap
        module procedure dswapdesc
        module procedure sswapdesc
    END INTERFACE deswap
    INTERFACE swapcol   !swap columns
        module procedure dswapcol
        module procedure sswapcol
    END INTERFACE swapcol
    INTERFACE swaprow   !swap rows
        module procedure dswaprow
        module procedure sswaprow
    END INTERFACE swaprow
    INTERFACE imaxloc   !location of array maximum as an integer
        module procedure imaxloc_dp
        module procedure imaxloc_sp
    END INTERFACE imaxloc
    contains
    
    FUNCTION oand(a, na, b, nb)
        integer, intent(in):: na, nb
        logical, intent(in):: a(na), b(nb)
        logical, DIMENSION(na,nb):: oand
       
        oand = spread(a, 2, nb) .AND. spread(b, 1, na);
        return
    END FUNCTION oand
    
    !location of array maximum as an integer
    FUNCTION imaxloc_dp(A)
        real(DP), intent(in):: A(:)
        integer:: imax(1)
        integer:: imaxloc_dp
        
        imax = maxloc(A(:));
        imaxloc_dp = imax(1);
    END FUNCTION imaxloc_dp
    FUNCTION imaxloc_sp(A)
        real(SP), intent(in):: A(:)
        integer:: imax(1)
        integer:: imaxloc_sp
        
        imax = maxloc(A(:));
        imaxloc_sp = imax(1);
    END FUNCTION imaxloc_sp
    
    !outer product
    FUNCTION douterpd(a, b)
        real(c_double), DIMENSION(:), INTENT(IN) :: a, b
        real(c_double), DIMENSION(size(a),size(b)) :: douterpd
        
        douterpd = spread(a, dim=2, ncopies=size(b)) * &
        spread(b, dim=1, ncopies=size(a));
        return
    END FUNCTION douterpd
    FUNCTION souterpd(a, b)
        real(c_float), DIMENSION(:), INTENT(IN) :: a, b
        real(c_float), DIMENSION(size(a),size(b)) :: souterpd
        
        souterpd = spread(a, dim=2, ncopies=size(b)) * &
        spread(b, dim=1, ncopies=size(a));
        return
    END FUNCTION souterpd
    
    !-- solve A(n,n)x(n) = b(n) with partial pivoting
    FUNCTION dgaussj(Ai, bi)
        real(DP), intent(in):: Ai(:,:), bi(:)
        integer:: ipiv(size(bi)), idxr(size(bi)), idxc(size(bi))
        logical:: lpiv(size(bi))
        real(DP):: pivinv, dumc(size(bi)), temp(1,size(bi)), tmp, tp(size(bi))
        real(DP):: A(size(Ai, dim=1),size(Ai, dim=2)), b(size(bi)) 
        real(DP):: dgaussj(size(bi))
        integer, target:: irc(2)
        integer, pointer:: pRow, pCol
        integer:: n, i, j, inc, sts = 1
        
        n = size(bi)
        ipiv = 0; inc = 1;
        pRow => irc(1); pCol => irc(2);
        A = Ai; b = bi;
    !$OMP DO
col_loop: DO i = 1, n, inc
            lpiv = (ipiv == 0);
            irc = MAXLOC(abs(A), spread(lpiv, 2, n) .AND. spread(lpiv, 1, n));
            ipiv(pCol) = ipiv(pCol) + 1;
            
            if(ipiv(pCol) > 1) then
                sts = -1; 
                goto 100;
            endif
            
            !interchange rows
            if(pRow /= pCol) then
                temp(1,:) = A(pRow, :);
                A(pRow,:) = A(pCol,:);
                A(pCol,:) = temp(1,:);
                tmp = b(pRow);
                b(pRow) = b(pCol);
                b(pCol) = tmp;
            endif
            
            idxr(i) = pRow; idxc(i) = pCol;
            if(A(pCol,pCol) == 0.0) then
                sts = -pCol; 
                goto 100; 
            endif
                
            pivinv = 1.0/A(pCol,pCol);
            A(pCol,pCol) = 1.0;
            A(pCol,:) = A(pCol,:)*pivinv;
            b(pCol) = b(pCol)*pivinv;
            dumc = A(:,pCol);
            A(:,pCol) = 0.0;
            A(pCol,pCol) = pivinv;
            
            A(1:pCol-1,:) = A(1:pCol-1,:) - spread(dumc(1:pcol-1), 2, n)*spread(A(pcol,:), 1, size(dumc(1:pcol-1)));
            b(1:pCol-1) = b(1:pCol-1) - dumc(1:pcol-1) * b(pcol);
            A(pCol+1:,:) = A(pCol+1:,:) - spread(dumc(pcol+1:), 2, n)*spread(A(pcol,:), 1, size(dumc(pcol+1:)));
            b(pCol+1:) = b(pCol+1:) - dumc(pCol+1:) * b(pcol);
            
100         if(ipiv(pCol) > 1 .OR. A(pCol,pCol) == 0.0) inc = n - i;

        ENDDO col_loop
    !$OMP END DO
        
    !$OMP DO
        DO j = n, 1, -1
            tp = A(:,idxr(j));
            A(:,idxr(j)) = A(:,idxc(j));
            A(:,idxc(j)) = tp;
        ENDDO
    !$OMP END DO
        dgaussj = b; return;
    END FUNCTION dgaussj
    
    !computes (x^2 + y^2)^(1/2) without destructive under-/overflow
    REAL(DP) FUNCTION dpythag(x, y)
        real(DP), intent(in):: x, y
        real(DP):: absx, absy
        
        absx = abs(x); absy = abs(y);
        if(absx > absy) then 
            dpythag = absx * SQRT(1.0_DP + (absy/absx)**2);
        else
            if(absy == 0) then 
                dpythag = 0.0;
            else
                dpythag = absy * SQRT(1.0_DP + (absx/absy)**2);
            endif
        endif
        return
    END FUNCTION dpythag
    REAL(SP) FUNCTION spythag(x, y)
        real(SP), intent(in):: x, y
        real(SP):: absx, absy
        
        absx = abs(x); absy = abs(y);
        if(absx > absy) then 
            spythag = absx * SQRT(1.0_SP + (absy/absx)**2);
        else
            if(absy == 0) then 
                spythag = 0.0;
            else
                spythag = absy * SQRT(1.0_SP + (absx/absy)**2);
            endif
        endif
        return
    END FUNCTION spythag
    
    !v{...,min,...,max,...} <-- v{max,..., min}
    INTEGER FUNCTION dswapdesc(v, maxmin)
        INTEGER, INTENT(INOUT), OPTIONAL :: maxmin(2)
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: v
        REAL(DP) :: dtemp
        INTEGER :: loc, n
        
        n = size(v);
        if (n <= 0) n = 1;
        dswapdesc = n;
        if (n == 1) return;
        loc = MAXLOC(v, 1); maxmin(1) = loc;
        dtemp = v(1); v(1) = v(loc); v(loc) = dtemp;
        loc = MINLOC(v, 1); maxmin(2) = loc;
        dtemp = v(n); v(n) = v(loc); v(loc) = dtemp;
        
    END FUNCTION dswapdesc
    INTEGER FUNCTION sswapdesc(v, maxmin)
        INTEGER, INTENT(INOUT), OPTIONAL :: maxmin(2)
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: v
        REAL(SP) :: dtemp
        INTEGER :: loc, n
        
        n = size(v); 
        if (n <= 0) n = 1;
        sswapdesc = n;
        if (sswapdesc == 1) return;
        loc = MAXLOC(v, 1); maxmin(1) = loc;
        dtemp = v(1); v(1) = v(loc); v(loc) = dtemp;
        loc = MINLOC(v, 1); maxmin(2) = loc;
        dtemp = v(n); v(n) = v(loc); v(loc) = dtemp;
        
    END FUNCTION sswapdesc
    
    !swap the two columns indicated by input indices 
    SUBROUTINE dswapcol(a, c1, c2)
        INTEGER, INTENT(IN):: c1, c2
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
        REAL(DP) :: temp(size(a, 1))
        
        temp = a(:,c1);
        a(:,c1) = a(:,c2);
        a(:,c2) = temp;
        
    END SUBROUTINE dswapcol
    SUBROUTINE sswapcol(a, c1, c2)
        INTEGER, INTENT(IN):: c1, c2
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
        REAL(SP) :: temp(size(a, 1))
        
        temp = a(:,c1);
        a(:,c1) = a(:,c2);
        a(:,c2) = temp;
        
    END SUBROUTINE sswapcol
    
    !swap the two rows indicated by input indices 
    SUBROUTINE dswaprow(a, r1, r2)
        INTEGER, INTENT(IN):: r1, r2
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
        REAL(DP) :: temp(size(a,2))
        temp = a(r1,:); a(r1,:) = a(r2,:); a(r2,:) = temp;
    END SUBROUTINE dswaprow
    SUBROUTINE sswaprow(a, r1, r2)
        INTEGER, INTENT(IN):: r1, r2
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
        REAL(SP) :: temp(size(a,2))
        temp = a(r1,:); a(r1,:) = a(r2,:); a(r2,:) = temp;
    END SUBROUTINE sswaprow

end module utl