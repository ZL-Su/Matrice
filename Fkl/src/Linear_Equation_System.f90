!------------------------------------------------
!--Description: Solve linear equation system
!--Author: Dgelom Su
!------------------------------------------------
module les
    include "mkl_direct_call.fi"
    use f95_precision
    use iso_c_binding
    use lapack95
    use utl
    implicit none
    contains
    
    !-- solve A(n,n)x(n) = b(n) with partial pivoting
    integer(c_int) function dgesv_gauss(A, b, n) bind(c, name = '_dgesv')
        integer(c_int), intent(in):: n[value]
        real(c_double), intent(inout):: A(n,n), b(n)
        integer(c_int):: ipiv(n), idxr(n), idxc(n)
        logical:: lpiv(n)
        real(c_double):: pivinv, dumc(n), temp(1,n), tmp, tp(n)
        integer(c_int), target:: irc(2)
        integer(c_int), pointer:: pRow, pCol
        integer:: i, j, inc
        
        ipiv = 0; inc = 1;
        pRow => irc(1); pCol => irc(2);
        dgesv_gauss = 1;
        A = transpose(A);
    !$OMP DO
col_loop: DO i = 1, n, inc
            lpiv = (ipiv == 0);
            irc = MAXLOC(abs(A), spread(lpiv, 2, n) .AND. spread(lpiv, 1, n));
            ipiv(pCol) = ipiv(pCol) + 1;
            
            if(ipiv(pCol) > 1) then
                dgesv_gauss = -1; 
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
                dgesv_gauss = -pCol; 
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
        if(dgesv_gauss /= 1) return;
        dgesv_gauss = 1; return;
    end function dgesv_gauss
    integer(c_int) function sgesv_gauss(A, b, n) bind(c, name = '_sgesv')
        integer(c_int), intent(in):: n[value]
        real(c_float), intent(inout):: A(n,n), b(n)
        integer(c_int):: ipiv(n), idxr(n), idxc(n)
        logical:: lpiv(n)
        real(c_float):: pivinv, dumc(n), temp(1,n), tmp, tp(n)
        integer(c_int), target:: irc(2)
        integer(c_int), pointer:: pRow, pCol
        integer:: i, j, inc
        
        ipiv = 0; inc = 1; sgesv_gauss = 1;
        pRow => irc(1); pCol => irc(2);
        A = transpose(A);
    !$OMP DO
col_loop: DO i = 1, n, inc
            lpiv = (ipiv == 0);
            irc = MAXLOC(abs(A), spread(lpiv, 2, n) .AND. spread(lpiv, 1, n));
            ipiv(pCol) = ipiv(pCol) + 1;
            
            if(ipiv(pCol) > 1) then
                sgesv_gauss = -1; 
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
                sgesv_gauss = -pCol; 
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
        if(sgesv_gauss /= 1) return;
        sgesv_gauss = 1; return;
    end function sgesv_gauss
    
    !-- solve A(m,n)x(n) = b(m) with partial pivoting
    integer(c_int) function dgesv2_gauss(AN, bN, m, n) bind(c, name = '_dgeosv')
        integer(c_int), intent(in):: m[value], n[value]
        real(c_double), intent(inout):: AN(n,m), bN(m)
        integer(c_int):: ipiv(n), idxr(n), idxc(n)
        logical:: lpiv(n)
        real(c_double):: pivinv, dumc(n), temp(1,n), tmp, tp(n), A(n,n), b(n)
        integer(c_int), target:: irc(2)
        integer(c_int), pointer:: pRow, pCol
        integer:: i, j, inc
        
        if (m /= n) then
            A = MATMUL(AN, transpose(AN));
            b = MATMUL(AN, bN);
        else 
            A = transpose(AN);
            b = bN;
        endif
        
        ipiv = 0; inc = 1;
        pRow => irc(1); pCol => irc(2);
        dgesv2_gauss = 1;
        
    !$OMP DO
col_loop: DO i = 1, n, inc
            lpiv = (ipiv == 0);
            irc = MAXLOC(abs(A), spread(lpiv, 2, n) .AND. spread(lpiv, 1, n));
            ipiv(pCol) = ipiv(pCol) + 1;
            
            if(ipiv(pCol) > 1) then
                dgesv2_gauss = -1; 
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
                dgesv2_gauss = -pCol; 
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
            bN(j) = b(j);
        ENDDO
    !$OMP END DO
        if(dgesv2_gauss /= 1) return;
        dgesv2_gauss = 1; return;
    end function dgesv2_gauss
    integer(c_int) function sgesv2_gauss(AN, bN, m, n) bind(c, name = '_sgeosv')
        integer(c_int), intent(in):: m[value], n[value]
        real(c_float), intent(inout):: AN(n,m), bN(m)
        integer(c_int):: ipiv(n), idxr(n), idxc(n)
        logical:: lpiv(n)
        real(c_float):: pivinv, dumc(n), temp(1,n), tmp, tp(n), A(n,n), b(n)
        integer(c_int), target:: irc(2)
        integer(c_int), pointer:: pRow, pCol
        integer:: i, j, inc
        
        if (m /= n) then
            A = MATMUL(AN, transpose(AN));
            b = MATMUL(AN, bN);
        else 
            A = transpose(AN);
            b = bN;
        endif
        
        ipiv = 0; inc = 1;
        pRow => irc(1); pCol => irc(2);
        sgesv2_gauss = 1;
        
    !$OMP DO
col_loop: DO i = 1, n, inc
            lpiv = (ipiv == 0);
            irc = MAXLOC(abs(A), spread(lpiv, 2, n) .AND. spread(lpiv, 1, n));
            ipiv(pCol) = ipiv(pCol) + 1;
            
            if(ipiv(pCol) > 1) then
                sgesv2_gauss = -1; 
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
                sgesv2_gauss = -pCol; 
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
            bN(j) = b(j);
        ENDDO
    !$OMP END DO
        if(sgesv2_gauss /= 1) return;
        sgesv2_gauss = 1; return;
    end function sgesv2_gauss
    
    !-- solve inverse of A(n,n)
    integer(c_int) function dginv_gauss(A, n) bind(c, name = '_dginv')
        integer(c_int), intent(in):: n[value]
        real(c_double), intent(inout):: A(n,n)
        integer(c_int):: ipiv(n), idxr(n), idxc(n)
        logical:: lpiv(n)
        real(c_double):: pivinv, dumc(n), temp(1,n), tp(n)
        integer(c_int), target:: irc(2)
        integer(c_int), pointer:: pRow, pCol
        integer:: i, j, inc
        
        ipiv = 0; inc = 1; dginv_gauss = 1;
        pRow => irc(1); pCol => irc(2);
    !$OMP DO
col_loop: DO i = 1, n, inc
            lpiv = (ipiv == 0);
            irc = MAXLOC(abs(A), spread(lpiv, 2, n) .AND. spread(lpiv, 1, n));
            ipiv(pCol) = ipiv(pCol) + 1;
            
            if(ipiv(pCol) > 1) then
                dginv_gauss = -1; 
                goto 100;
            endif
            
            !interchange rows
            if(pRow /= pCol) then
                temp(1,:) = A(pRow, :);
                A(pRow,:) = A(pCol,:);
                A(pCol,:) = temp(1,:);
            endif
            
            idxr(i) = pRow; idxc(i) = pCol;
            if(A(pCol,pCol) == 0.0) then
                dginv_gauss = -pCol; 
                goto 100; 
            endif
                
            pivinv = 1.0/A(pCol,pCol);
            A(pCol,pCol) = 1.0;
            A(pCol,:) = A(pCol,:)*pivinv;
            dumc = A(:,pCol);
            A(:,pCol) = 0.0;
            A(pCol,pCol) = pivinv;
            
            A(1:pCol-1,:) = A(1:pCol-1,:) - spread(dumc(1:pcol-1), 2, n)*spread(A(pcol,:), 1, size(dumc(1:pcol-1)));
            A(pCol+1:,:) = A(pCol+1:,:) - spread(dumc(pcol+1:), 2, n)*spread(A(pcol,:), 1, size(dumc(pcol+1:)));
            
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
        if(dginv_gauss /= 1) return;
        dginv_gauss = 1; return;
    end function dginv_gauss
    integer(c_int) function sginv_gauss(A, n) bind(c, name = '_sginv')
        integer(c_int), intent(in):: n[value]
        real(c_float), intent(inout):: A(n,n)
        integer(c_int):: ipiv(n), idxr(n), idxc(n)
        logical:: lpiv(n)
        real(c_float):: pivinv, dumc(n), temp(1,n), tp(n)
        integer(c_int), target:: irc(2)
        integer(c_int), pointer:: pRow, pCol
        integer:: i, j, inc
        
        ipiv = 0; inc = 1; sginv_gauss = 1;
        pRow => irc(1); pCol => irc(2);
    !$OMP DO
col_loop: DO i = 1, n, inc
            lpiv = (ipiv == 0);
            irc = MAXLOC(abs(A), spread(lpiv, 2, n) .AND. spread(lpiv, 1, n));
            ipiv(pCol) = ipiv(pCol) + 1;
            
            if(ipiv(pCol) > 1) then
                sginv_gauss = -1; 
                goto 100;
            endif
            
            !interchange rows
            if(pRow /= pCol) then
                temp(1,:) = A(pRow, :);
                A(pRow,:) = A(pCol,:);
                A(pCol,:) = temp(1,:);
            endif
            
            idxr(i) = pRow; idxc(i) = pCol;
            if(A(pCol,pCol) == 0.0) then
                sginv_gauss = -pCol; 
                goto 100; 
            endif
                
            pivinv = 1.0/A(pCol,pCol);
            A(pCol,pCol) = 1.0;
            A(pCol,:) = A(pCol,:)*pivinv;
            dumc = A(:,pCol);
            A(:,pCol) = 0.0;
            A(pCol,pCol) = pivinv;
            
            A(1:pCol-1,:) = A(1:pCol-1,:) - spread(dumc(1:pcol-1), 2, n)*spread(A(pcol,:), 1, size(dumc(1:pcol-1)));
            A(pCol+1:,:) = A(pCol+1:,:) - spread(dumc(pcol+1:), 2, n)*spread(A(pcol,:), 1, size(dumc(pcol+1:)));
            
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
        if(sginv_gauss /= 1) return;
        sginv_gauss = 1; return;
    end function sginv_gauss
    
    !-- svd decomposition A = U¦²V, where A should be passed as row-major.
    integer(c_int) function gsvd_sp(A, W, V, row, col) bind(c, name = '_sgesvd')
        USE UTL, ONLY : outer_prod, pythag, deswap, swaprow
        integer(c_int), intent(in):: row[value], col[value]
        real(c_float), intent(inout):: A(col,row), V(col,col), W(col)
        integer :: m, n, i, j, k, l, its, nm, maxits, mmloc(2)
        real(SP) :: anorm, c, f, g, h, s, scale, x, y, z
        real(SP), dimension(size(A,2)) :: tempm
        real(SP), dimension(size(A,1)) :: tempn, rv
        
        maxits = 30; n = size(A,1); m = size(A,2);
        if (n/=size(W) .OR. n /= size(V,1) .OR. n/=size(V,2)) then
            gsvd_sp = -4; goto 100;
        endif
        
        g = 0.0; scale = 0.0;
        !$OMP DO
        DO i = 1, n  ! reduction to bidiagonal form
            l = i + 1; rv(i) = scale * g; g = 0.0; scale = 0.0;
            if (i <= m) then
                scale = SUM(abs(A(i, i:m)));
                if (scale /= 0.0) then
                    A(i, i:m) = A(i, i:m)/scale;
                    s = dot_product(A(i, i:m),A(i, i:m));
                    f = A(i,i); g = -SIGN(sqrt(s),f);
                    h = f*g-s; A(i,i) = f-g;
                    tempn(l:n) = MATMUL(A(l:n, i:m), A(i, i:m))/h;
                    A(l:n, i:m) = A(l:n, i:m) + outer_prod(tempn(l:n), A(i, i:m));
                    A(i, i:m) = scale*A(i, i:m)
                endif
            endif
            W(i) = scale*g; g = 0.0; scale = 0.0;
            if ((i <= m) .AND. (i /= n)) then
                scale = SUM(abs(A(l:n, i)));
                if (scale /= 0.0) then
                    A(l:n, i) = A(l:n, i)/scale
                    s = dot_product(A(l:n, i), A(l:n, i))
                    f = A(l,i); g = -SIGN(sqrt(s),f)
                    h = f*g-s; A(l,i) = f-g
                    rv(l:n) = A(l:n, i)/h
                    tempm(l:m) = MATMUL(A(l:n, i), A(l:n,l:m))
                    A(l:n,l:m) = A(l:n,l:m)+outer_prod(rv(l:n), tempm(l:m))
                    A(l:n, i) = scale*A(l:n, i)
                endif
            endif
        ENDDO
        !$OMP END DO
        anorm = MAXVAL(abs(w) + abs(rv));
        ! accumulation of right-hand transformations
        !$OMP DO
        DO i = n, 1, -1
            if (i < n) then
                if (g /= 0.0) then
                    V(i,l:n) = (A(l:n,i)/A(l,i))/g;
                    tempn(l:n) = MATMUL(V(l:n,l:n), A(l:n,i));
                    V(l:n,l:n) = V(l:n,l:n) + outer_prod(tempn(l:n), V(i, l:n));
                endif
                V(i,l:n) = 0.0; V(l:n,i) = 0.0;
            endif
            V(i,i) = 1.0; g = rv(i); l = i;
        ENDDO
        !$OMP END DO
        ! accumulation of left-hand transformations
        !$OMP DO
        DO i = min(m,n), 1, -1
            l = i + 1; g = W(i); A(l:n, i) = 0.0;
            if (g /= 0.0) then
                g = 1.0_DP/g;
                tempn(l:n) = MATMUL(A(l:n, l:m), A(i, l:m));
                tempn(l:n) = (tempn(l:n)/A(i,i))*g;
                A(l:n,i:m) = A(l:n,i:m) + outer_prod(tempn(l:n), A(i,i:m));
                A(i,i:m) = A(i,i:m) * g;
            else
                A(i,i:m) = 0.0;
            endif
            A(i,i) = A(i,i) + 1.0_SP;
        ENDDO
        !$OMP END DO
        
        !-- diagonalization of the bidiagonal form: loop over singular values and over allowed iterations
        DO k = n, 1, -1
            DO its = 1, maxits
                DO l = k, 1, -1 ! Test for splitting
                    nm = l - 1;
                    if ((ABS(rv(l)) + anorm) == anorm) exit;
                    if ((ABS(W(nm)) + anorm) == anorm) then
                        c = 0.0; s = 1.0;
                        DO i = l, k
                            f = s * rv(i); rv(i) = c * rv(i);
                            if ((ABS(f)+anorm) == anorm) exit;
                            g = W(i); h = pythag(f, g);
                            W(i) = h; h = 1.0_DP/h;
                            c = (g * h); s = -(f*h);
                            tempm(1:m) = A(nm, 1:m);
                            A(nm, 1:m) = A(nm, 1:m)*c + A(i, 1:m)*s;
                            A(i, 1:m) = -tempm(1:m)*s + A(i, 1:m)*c;
                        ENDDO
                        EXIT
                    endif
                ENDDO
                z = W(k);
                if (l == k) then !convergence
                    if (z < 0.0) then
                        W(k) = -z; V(k, 1:n) = -V(k, 1:n);
                    endif
                    EXIT
                endif
                if (its == maxits) then
                    gsvd_sp = -2; goto 100;
                endif
                x = W(l); nm = k - 1; y = W(nm);
                g = rv(nm); h = rv(k);
                f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_SP*h*y);
                g = pythag(f, 1.0_SP);
                f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
                !Next QR transformation
                c = 1.0; s = 1.0;
                !$OMP DO
                DO J = l, nm
                    i = j + 1; 
                    g = rv(i); y = W(i);
                    h = s*g; g = c*g;
                    z = pythag(f, h);
                    rv(j) = z; c = f/z; s = h/z;
                    f = (x*c)+(g*s); g = -(x*s)+(g*c);
                    h = y*s; y = y*c;
                    tempn(1:n) = V(j, 1:n);
                    V(j, 1:n) = V(j, 1:n)*c + V(i, 1:n)*s;
                    V(i, 1:n) = -tempn(1:n)*s + V(i, 1:n)*c;
                    z = pythag(f,h); W(j)=z;
                    if (z /= 0.0) then 
                        z = 1.0_SP/z; c = f*z; s = h*z;
                    endif
                    f = (c*g)+(s*y); x = -(s*g)+(c*y);
                    tempm(1:m) = A(j, 1:m);
                    A(j, 1:m) = A(j, 1:m)*c + A(i, 1:m)*s;
                    A(i, 1:m) = -tempm(1:m)*s + A(i, 1:m)*c
                ENDDO
                !$OMP END DO
                rv(l) = 0.0; rv(k) = f; W(k) = x;
            ENDDO
        ENDDO
        gsvd_sp = 0;
        !reorder
        DO i = 1, n
            l = (n-i+1);
            k = deswap(W(i:l), mmloc);
            if (k == 1) exit;
            call swaprow(A, i, mmloc(1));
            call swaprow(A, l, mmloc(2));
            call swaprow(V, i, mmloc(1));
            call swaprow(V, l, mmloc(2));
        ENDDO
        V = transpose(V)
100     return
    end function gsvd_sp
    integer(c_int) function gsvd_dp(A, W, V, row, col) bind(c, name = '_dgesvd')
        USE UTL, ONLY : outer_prod, pythag, deswap, swaprow
        integer(c_int), intent(in):: row[value], col[value]
        real(c_double), intent(inout):: A(col,row), V(col,col), W(col)
        integer :: m, n, i, j, k, l, its, nm, maxits, mmloc(2)
        real(DP) :: anorm, c, f, g, h, s, scale, x, y, z
        real(DP), dimension(size(A,2)) :: tempm
        real(DP), dimension(size(A,1)) :: tempn, rv
        
        maxits = 30; n = size(A,1); m = size(A,2);
        if (n/=size(W) .OR. n /= size(V,1) .OR. n/=size(V,2)) then
            gsvd_dp = -4; goto 100;
        endif
        
        g = 0.0; scale = 0.0;
        !$OMP DO
        DO i = 1, n  ! reduction to bidiagonal form
            l = i + 1; rv(i) = scale * g; g = 0.0; scale = 0.0;
            if (i <= m) then
                scale = SUM(abs(A(i, i:m)));
                if (scale /= 0.0) then
                    A(i, i:m) = A(i, i:m)/scale;
                    s = dot_product(A(i, i:m),A(i, i:m));
                    f = A(i,i); g = -SIGN(sqrt(s),f);
                    h = f*g-s; A(i,i) = f-g;
                    tempn(l:n) = MATMUL(A(l:n, i:m), A(i, i:m))/h;
                    A(l:n, i:m) = A(l:n, i:m) + outer_prod(tempn(l:n), A(i, i:m));
                    A(i, i:m) = scale*A(i, i:m)
                endif
            endif
            W(i) = scale*g; g = 0.0; scale = 0.0;
            if ((i <= m) .AND. (i /= n)) then
                scale = SUM(abs(A(l:n, i)));
                if (scale /= 0.0) then
                    A(l:n, i) = A(l:n, i)/scale
                    s = dot_product(A(l:n, i), A(l:n, i))
                    f = A(l,i); g = -SIGN(sqrt(s),f)
                    h = f*g-s; A(l,i) = f-g
                    rv(l:n) = A(l:n, i)/h
                    tempm(l:m) = MATMUL(A(l:n, i), A(l:n,l:m))
                    A(l:n,l:m) = A(l:n,l:m)+outer_prod(rv(l:n), tempm(l:m))
                    A(l:n, i) = scale*A(l:n, i)
                endif
            endif
        ENDDO
        !$OMP END DO
        anorm = MAXVAL(abs(w) + abs(rv));
        ! accumulation of right-hand transformations
        !$OMP DO
        DO i = n, 1, -1
            if (i < n) then
                if (g /= 0.0) then
                    V(i,l:n) = (A(l:n,i)/A(l,i))/g;
                    tempn(l:n) = MATMUL(V(l:n,l:n), A(l:n,i));
                    V(l:n,l:n) = V(l:n,l:n) + outer_prod(tempn(l:n), V(i, l:n));
                endif
                V(i,l:n) = 0.0; V(l:n,i) = 0.0;
            endif
            V(i,i) = 1.0; g = rv(i); l = i;
        ENDDO
        !$OMP END DO
        ! accumulation of left-hand transformations
        !$OMP DO
        DO i = min(m,n), 1, -1
            l = i + 1; g = W(i); A(l:n, i) = 0.0;
            if (g /= 0.0) then
                g = 1.0_DP/g;
                tempn(l:n) = MATMUL(A(l:n, l:m), A(i, l:m));
                tempn(l:n) = (tempn(l:n)/A(i,i))*g;
                A(l:n,i:m) = A(l:n,i:m) + outer_prod(tempn(l:n), A(i,i:m));
                A(i,i:m) = A(i,i:m) * g;
            else
                A(i,i:m) = 0.0;
            endif
            A(i,i) = A(i,i) + 1.0_DP;
        ENDDO
        !$OMP END DO
        
        !-- diagonalization of the bidiagonal form: loop over singular values and over allowed iterations
        DO k = n, 1, -1
            DO its = 1, maxits
                DO l = k, 1, -1 ! Test for splitting
                    nm = l - 1;
                    if ((ABS(rv(l)) + anorm) == anorm) exit;
                    if ((ABS(W(nm)) + anorm) == anorm) then
                        c = 0.0; s = 1.0;
                        DO i = l, k
                            f = s * rv(i); rv(i) = c * rv(i);
                            if ((ABS(f)+anorm) == anorm) exit;
                            g = W(i); h = pythag(f, g);
                            W(i) = h; h = 1.0_DP/h;
                            c = (g * h); s = -(f*h);
                            tempm(1:m) = A(nm, 1:m);
                            A(nm, 1:m) = A(nm, 1:m)*c + A(i, 1:m)*s;
                            A(i, 1:m) = -tempm(1:m)*s + A(i, 1:m)*c;
                        ENDDO
                        EXIT
                    endif
                ENDDO
                z = W(k);
                if (l == k) then !convergence
                    if (z < 0.0) then
                        W(k) = -z; V(k, 1:n) = -V(k, 1:n);
                    endif
                    EXIT
                endif
                if (its == maxits) then
                    gsvd_dp = -2; goto 100;
                endif
                x = W(l); nm = k - 1; y = W(nm);
                g = rv(nm); h = rv(k);
                f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_DP*h*y);
                g = pythag(f, 1.0_DP);
                f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
                !Next QR transformation
                c = 1.0; s = 1.0;
                !$OMP DO
                DO J = l, nm
                    i = j + 1; 
                    g = rv(i); y = W(i);
                    h = s*g; g = c*g;
                    z = pythag(f, h);
                    rv(j) = z; c = f/z; s = h/z;
                    f = (x*c)+(g*s); g = -(x*s)+(g*c);
                    h = y*s; y = y*c;
                    tempn(1:n) = V(j, 1:n);
                    V(j, 1:n) = V(j, 1:n)*c + V(i, 1:n)*s;
                    V(i, 1:n) = -tempn(1:n)*s + V(i, 1:n)*c;
                    z = pythag(f,h); W(j)=z;
                    if (z /= 0.0) then 
                        z = 1.0_DP/z; c = f*z; s = h*z;
                    endif
                    f = (c*g)+(s*y); x = -(s*g)+(c*y);
                    tempm(1:m) = A(j, 1:m);
                    A(j, 1:m) = A(j, 1:m)*c + A(i, 1:m)*s;
                    A(i, 1:m) = -tempm(1:m)*s + A(i, 1:m)*c
                ENDDO
                !$OMP END DO
                rv(l) = 0.0; rv(k) = f; W(k) = x;
            ENDDO
        ENDDO
        gsvd_dp = 0;
        !reorder
        DO i = 1, n
            l = (n-i+1);
            k = deswap(W(i:l), mmloc);
            if (k == 1) exit;
            call swaprow(A, i, mmloc(1));
            call swaprow(A, l, mmloc(2));
            call swaprow(V, i, mmloc(1));
            call swaprow(V, l, mmloc(2));
        ENDDO
        V = transpose(V)
100     return
    end function gsvd_dp
    
    !-- svd decomposition A = U¦²V, where A should be passed as column-major if A is not a square matrix
    integer(c_int) function svd1_dp(A, W, V, row, col) bind(c, name = '_dgsvd')
        USE UTL, ONLY : outer_prod, pythag, deswap, swapcol
        integer(c_int), intent(in):: row[value], col[value]
        real(c_double), intent(inout):: A(row,col), V(col,col), W(col)
        integer :: m, n, i, j, k, l, its, nm, maxits, mmloc(2)
        real(DP) :: anorm, c, f, g, h, s, scale, x, y, z
        real(DP), dimension(size(A,1)) :: tempm
        real(DP), dimension(size(A,2)) :: tempn, rv
        
        maxits = 30; m = size(A,1); n = size(A,2);
        if (m == n) A = TRANSPOSE(A);
        if (n/=size(w).OR.n/=size(v,1).OR.n/=size(v,2)) then
            svd1_dp = -4; goto 100;
        endif
        
        g = 0.0; scale = 0.0;
        !$OMP DO
        DO i = 1, n  ! reduction to bidiagonal form
            l = i + 1; rv(i) = scale * g; g = 0.0; scale = 0.0;
            if (i <= m) then
                scale = SUM(abs(A(i:m,i)));
                if (scale /= 0.0) then
                    A(i:m,i) = A(i:m,i)/scale;
                    s = dot_product(A(i:m,i),A(i:m,i));
                    f = A(i,i); g = -SIGN(sqrt(s),f);
                    h = f*g-s; A(i,i) = f-g;
                    tempn(l:n) = MATMUL(A(i:m,i),A(i:m,l:n))/h;
                    A(i:m,l:n) = A(i:m,l:n) + outer_prod(A(i:m,i),tempn(l:n));
                    A(i:m,i) = scale*A(i:m,i)
                endif
            endif
            W(i) = scale*g; g = 0.0; scale = 0.0;
            if ((i <= m) .AND. (i /= n)) then
                scale = SUM(abs(A(i,l:n)));
                if (scale /= 0.0) then
                    A(i,l:n) = A(i,l:n)/scale
                    s = dot_product(A(i,l:n), A(i,l:n))
                    f = A(i,l); g = -SIGN(sqrt(s),f)
                    h = f*g-s; A(i,l) = f-g
                    rv(l:n) = A(i,l:n)/h
                    tempm(l:m) = MATMUL(A(l:m,l:n), A(i,l:n))
                    A(l:m,l:n) = A(l:m,l:n)+outer_prod(tempm(l:m),rv(l:n))
                    A(i,l:n) = scale*A(i,l:n)
                endif
            endif
        ENDDO
        !$OMP END DO
        anorm = MAXVAL(abs(w) + abs(rv));
        ! accumulation of right-hand transformations
        !$OMP DO
        DO i = n, 1, -1
            if (i < n) then
                if (g /= 0.0) then
                    V(l:n,i) = (A(i,l:n)/A(i,l))/g;
                    tempn(l:n) = MATMUL(A(i,l:n), V(l:n,l:n));
                    V(l:n,l:n) = V(l:n,l:n) + outer_prod(V(l:n,i), tempn(l:n));
                endif
                V(i,l:n) = 0.0; V(l:n,i) = 0.0;
            endif
            V(i,i) = 1.0; g = rv(i); l = i;
        ENDDO
        !$OMP END DO
        ! accumulation of left-hand transformations
        !$OMP DO
        DO i = min(m,n), 1, -1
            l = i + 1; g = W(i); A(i, l:n) = 0.0;
            if (g /= 0.0) then
                g = 1.0_DP/g;
                tempn(l:n) = MATMUL(A(l:m,i), A(l:m,l:n));
                tempn(l:n) = (tempn(l:n)/A(i,i))*g;
                A(i:m,l:n) = A(i:m,l:n) + outer_prod(A(i:m,i),tempn(l:n));
                A(i:m,i) = A(i:m,i) * g;
            else
                A(i:m,i) = 0.0;
            endif
            A(i,i) = A(i,i) + 1.0_DP;
        ENDDO
        !$OMP END DO
        
        !-- diagonalization of the bidiagonal form: loop over singular values and over allowed iterations
        DO k = n, 1, -1
            DO its = 1, maxits
                DO l = k, 1, -1 ! Test for splitting
                    nm = l - 1;
                    if ((ABS(rv(l)) + anorm) == anorm) exit;
                    if ((ABS(W(nm)) + anorm) == anorm) then
                        c = 0.0; s = 1.0;
                        DO i = l, k
                            f = s * rv(i); rv(i) = c * rv(i);
                            if ((ABS(f)+anorm) == anorm) exit;
                            g = W(i); h = pythag(f, g);
                            W(i) = h; h = 1.0_DP/h;
                            c = (g * h); s = -(f*h);
                            tempm(1:m) = A(1:m,nm);
                            A(1:m,nm) = A(1:m,nm)*c + A(1:m,i)*s;
                            A(1:m,i) = -tempm(1:m)*s + A(1:m,i)*c;
                        ENDDO
                        EXIT
                    endif
                ENDDO
                z = W(k);
                if (l == k) then !convergence
                    if (z < 0.0) then
                        W(k) = -z; v(1:n,k) = -v(1:n,k);
                    endif
                    EXIT
                endif
                if (its == maxits) then
                    svd1_dp = 0; goto 100;
                endif
                x = W(l); nm = k - 1; y = W(nm);
                g = rv(nm); h = rv(k);
                f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_DP*h*y);
                g = pythag(f, 1.0_DP);
                f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
                !Next QR transformation
                c = 1.0; s = 1.0;
                !$OMP DO
                DO J = l, nm
                    i = j + 1; 
                    g = rv(i); y = W(i);
                    h = s*g; g = c*g;
                    z = pythag(f, h);
                    rv(j) = z; c = f/z; s = h/z;
                    f = (x*c)+(g*s); g = -(x*s)+(g*c);
                    h = y*s; y = y*c;
                    tempn(1:n) = V(1:n,j);
                    V(1:n,j) = V(1:n,j)*c + V(1:n,i)*s;
                    V(1:n,i) = -tempn(1:n)*s + V(1:n,i)*c;
                    z = pythag(f,h); W(j)=z;
                    if (z /= 0.0) then 
                        z = 1.0_DP/z; c = f*z; s = h*z;
                    endif
                    f = (c*g)+(s*y); x = -(s*g)+(c*y);
                    tempm(1:m) = A(1:m,j);
                    A(1:m,j) = A(1:m,j)*c + A(1:m,i)*s;
                    A(1:m,i) = -tempm(1:m)*s + A(1:m,i)*c
                ENDDO
                !$OMP END DO
                rv(l) = 0.0; rv(k) = f; W(k) = x;
            ENDDO
        ENDDO
        svd1_dp = 1;
        !reorder
        DO i = 1, n
            l = (n-i+1);
            k = deswap(W(i:l), mmloc);
            if (k == 1) exit;
            call swapcol(A, i, mmloc(1));
            call swapcol(A, l, mmloc(2));
            call swapcol(V, i, mmloc(1));
            call swapcol(V, l, mmloc(2));
        ENDDO
        
100     return
    end function svd1_dp
    
    !-- fast svd decomposition A = U¦²V, where A should be passed as row-major
    integer(c_int) function svd_dp(A, W, V, Mr, Nc) bind(c, name = '_dfgsvd')
        USE UTL, ONLY : outer_prod, pythag, deswap, swapcol, swaprow
        integer(c_int), intent(in):: Mr[value], Nc[value]
        real(c_double), intent(inout):: A(Nc, Mr)
        real(c_double), intent(out):: V(Nc, Nc), W(Nc)
        real(DP), dimension(size(A,2)) :: tempm
        real(DP), dimension(size(A,1)) :: tempn, rv
        integer :: m, n, i, j, k, l, its, nm, maxits, mmloc(2)
        real(DP) :: anorm, c, f, g, h, s, scale, x, y, z
        
        maxits = 30; n = size(A,1); m = size(A,2);
        !if (m == n) A = TRANSPOSE(A);
        if (n/=size(w).or.n/=size(v,1).or.n/=size(v,2)) then
            svd_dp = -4; goto 100;
        endif
        
        g = 0.0; scale = 0.0;
        !$OMP DO
        DO i = 1, n  ! reduction to bidiagonal form
            l = i + 1; rv(i) = scale * g; g = 0.0; scale = 0.0;
            if (i <= m) then
                scale = SUM(abs(A(i, i:m)));
                if (scale /= 0.0) then
                    A(i, i:m) = A(i, i:m)/scale;
                    s = dot_product(A(i, i:m),A(i, i:m));
                    f = A(i,i); g = -SIGN(sqrt(s),f);
                    h = f*g-s; A(i,i) = f-g;
                    tempn(l:n) = MATMUL(A(l:n, i:m), A(i, i:m))/h;
                    A(l:n, i:m) = A(l:n, i:m) + outer_prod(tempn(l:n), A(i, i:m));
                    A(i, i:m) = scale*A(i, i:m)
                endif
            endif
            W(i) = scale*g; g = 0.0; scale = 0.0;
            if ((i <= m) .AND. (i /= n)) then
                scale = SUM(abs(A(l:n, i)));
                if (scale /= 0.0) then
                    A(l:n, i) = A(l:n, i)/scale
                    s = dot_product(A(l:n, i), A(l:n, i))
                    f = A(l, i); g = -SIGN(sqrt(s),f)
                    h = f*g-s; A(l, i) = f-g
                    rv(l:n) = A(l:n, i)/h
                    tempm(l:m) = MATMUL(A(l:n, i), A(l:n,l:m))
                    A(l:n, l:m) = A(l:n, l:m)+outer_prod(rv(l:n), tempm(l:m))
                    A(l:n, i) = scale*A(l:n, i)
                endif
            endif
        ENDDO
        !$OMP END DO
        anorm = MAXVAL(abs(w) + abs(rv));
        ! accumulation of right-hand transformations
        !$OMP DO
        DO i = n, 1, -1
            if (i < n) then
                if (g /= 0.0) then
                    V(i, l:n) = (A(l:n, i)/A(l, i))/g;
                    tempn(l:n) = MATMUL(V(l:n,l:n), A(l:n, i));
                    V(l:n,l:n) = V(l:n,l:n) + outer_prod(tempn(l:n), V(i, l:n));
                endif
                V(l:n, i) = 0.0; V(i, l:n) = 0.0;
            endif
            V(i,i) = 1.0; g = rv(i); l = i;
        ENDDO
        ! diagonalization of the bidiagonal form: loop over singular values and over allowed iterations
        !$OMP END DO 
        DO k = n, 1, -1
            DO its = 1, maxits
                DO l = k, 1, -1 ! Test for splitting
                    nm = l - 1;
                    if ((ABS(rv(l)) + anorm) == anorm) exit;
                    if ((ABS(W(nm)) + anorm) == anorm) then
                        c = 0.0; s = 1.0;
                        DO i = l, k
                            f = s * rv(i); rv(i) = c * rv(i);
                            if ((ABS(f) + anorm) == anorm) exit;
                            g = W(i); h = pythag(f, g);
                            W(i) = h; h = 1.0_DP/h;
                            c = (g * h); s = -(f * h);
                        ENDDO
                        EXIT
                    endif
                ENDDO
                z = W(k);
                if (l == k) then !convergence
                    if (z < 0.0) then
                        W(k) = -z; V(k, 1:n) = -V(k, 1:n);
                    endif
                    EXIT
                endif
                if (its == maxits) then  
                    svd_dp = -2; goto 100; 
                endif
                    
                x = W(l); nm = k - 1; y = W(nm);
                g = rv(nm); h = rv(k);
                f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_DP*h*y);
                g = pythag(f, 1.0_DP);
                f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
                !Next QR transformation
                c = 1.0; s = 1.0;
                !$OMP DO
                DO J = l, nm
                    i = j + 1; 
                    g = rv(i); y = W(i);
                    h = s*g; g = c*g;
                    z = pythag(f, h);
                    rv(j) = z; c = f/z; s = h/z;
                    f = (x*c)+(g*s); g = -(x*s)+(g*c);
                    h = y*s; y = y*c;
                    tempn(1:n) = V(j, 1:n);
                    V(j, 1:n) = V(j, 1:n)*c + V(i, 1:n)*s;
                    V(i, 1:n) = -tempn(1:n)*s + V(i, 1:n)*c;
                    z = pythag(f,h); W(j)=z;
                    if (z /= 0.0) then 
                        z = 1.0_DP/z; c = f*z; s = h*z;
                    endif
                    f = (c*g)+(s*y); x = -(s*g)+(c*y);
                ENDDO
                !$OMP END DO
                rv(l) = 0.0; rv(k) = f; W(k) = x;
            ENDDO
        ENDDO
        svd_dp = 0;
        !reorder
        DO i = 1, n
            l = (n-i+1); k = deswap(W(i:l), mmloc);
            if (k == 1) exit;
            call swaprow(V, i, mmloc(1));
            call swaprow(V, l, mmloc(2));
        ENDDO
        
100     return
    end function svd_dp
    !-- svd decomposition A = U¦²V, where A should be passed as column-major
    integer(c_int) function svd2_dp(A, W, V, row, col) bind(c, name = '_dfsvd')
        USE UTL, ONLY : outer_prod, pythag, deswap, swapcol
        integer(c_int), intent(in):: row[value], col[value]
        real(c_double), intent(inout):: A(row,col)
        real(c_double), intent(out):: V(col,col), W(col)
        real(DP), dimension(size(A,1)) :: tempm
        real(DP), dimension(size(A,2)) :: tempn, rv
        integer :: m, n, i, j, k, l, its, nm, maxits, mmloc(2)
        real(DP) :: anorm, c, f, g, h, s, scale, x, y, z
        
        maxits = 30; m = size(A,1); n = size(A,2);
        if (m == n) A = TRANSPOSE(A);
        if (n/=size(w).or.n/=size(v,1).or.n/=size(v,2)) then
            svd2_dp = -4; goto 100;
        endif
        
        g = 0.0; scale = 0.0;
        !$OMP DO
        DO i = 1, n  ! reduction to bidiagonal form
            l = i + 1; rv(i) = scale * g; g = 0.0; scale = 0.0;
            if (i <= m) then
                scale = SUM(abs(A(i:m,i)));
                if (scale /= 0.0) then
                    A(i:m,i) = A(i:m,i)/scale;
                    s = dot_product(A(i:m,i),A(i:m,i));
                    f = A(i,i); g = -SIGN(sqrt(s),f);
                    h = f*g-s; A(i,i) = f-g;
                    tempn(l:n) = MATMUL(A(i:m,i),A(i:m,l:n))/h;
                    A(i:m,l:n) = A(i:m,l:n) + outer_prod(A(i:m,i),tempn(l:n));
                    A(i:m,i) = scale*A(i:m,i)
                endif
            endif
            W(i) = scale*g; g = 0.0; scale = 0.0;
            if ((i <= m) .AND. (i /= n)) then
                scale = SUM(abs(A(i,l:n)));
                if (scale /= 0.0) then
                    A(i,l:n) = A(i,l:n)/scale
                    s = dot_product(A(i,l:n), A(i,l:n))
                    f = A(i,l); g = -SIGN(sqrt(s),f)
                    h = f*g-s; A(i,l) = f-g
                    rv(l:n) = A(i,l:n)/h
                    tempm(l:m) = MATMUL(A(l:m,l:n), A(i,l:n))
                    A(l:m,l:n) = A(l:m,l:n)+outer_prod(tempm(l:m),rv(l:n))
                    A(i,l:n) = scale*A(i,l:n)
                endif
            endif
        ENDDO
        !$OMP END DO
        anorm = MAXVAL(abs(w) + abs(rv));
        ! accumulation of right-hand transformations
        !$OMP DO
        DO i = n, 1, -1
            if (i < n) then
                if (g /= 0.0) then
                    V(l:n,i) = (A(i,l:n)/A(i,l))/g;
                    tempn(l:n) = MATMUL(A(i,l:n), V(l:n,l:n));
                    V(l:n,l:n) = V(l:n,l:n) + outer_prod(V(l:n,i), tempn(l:n));
                endif
                V(i,l:n) = 0.0; V(l:n,i) = 0.0;
            endif
            V(i,i) = 1.0; g = rv(i); l = i;
        ENDDO
        !$OMP END DO
        !-- diagonalization of the bidiagonal form: loop over singular values and over allowed iterations
        DO k = n, 1, -1
            DO its = 1, maxits
                DO l = k, 1, -1 ! Test for splitting
                    nm = l - 1;
                    if ((ABS(rv(l)) + anorm) == anorm) exit;
                    if ((ABS(W(nm)) + anorm) == anorm) then
                        c = 0.0; s = 1.0;
                        DO i = l, k
                            f = s * rv(i); rv(i) = c * rv(i);
                            if ((ABS(f)+anorm) == anorm) exit;
                            g = W(i); h = pythag(f, g);
                            W(i) = h; h = 1.0_DP/h;
                            c = (g * h); s = -(f*h);
                        ENDDO
                        EXIT
                    endif
                ENDDO
                z = W(k);
                if (l == k) then !convergence
                    if (z < 0.0) then
                        W(k) = -z; v(1:n,k) = -v(1:n,k);
                    endif
                    EXIT
                endif
                if (its == maxits) then
                    svd2_dp = 0; goto 100;
                endif
                x = W(l); nm = k - 1; y = W(nm);
                g = rv(nm); h = rv(k);
                f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_DP*h*y);
                g = pythag(f, 1.0_DP);
                f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
                !Next QR transformation
                c = 1.0; s = 1.0;
                !$OMP DO
                DO J = l, nm
                    i = j + 1; 
                    g = rv(i); y = W(i);
                    h = s*g; g = c*g;
                    z = pythag(f, h);
                    rv(j) = z; c = f/z; s = h/z;
                    f = (x*c)+(g*s); g = -(x*s)+(g*c);
                    h = y*s; y = y*c;
                    tempn(1:n) = V(1:n,j);
                    V(1:n,j) = V(1:n,j)*c + V(1:n,i)*s;
                    V(1:n,i) = -tempn(1:n)*s + V(1:n,i)*c;
                    z = pythag(f,h); W(j)=z;
                    if (z /= 0.0) then 
                        z = 1.0_DP/z; c = f*z; s = h*z;
                    endif
                    f = (c*g)+(s*y); x = -(s*g)+(c*y);
                ENDDO
                !$OMP END DO
                rv(l) = 0.0; rv(k) = f; W(k) = x;
            ENDDO
        ENDDO
        svd2_dp = 1;
        !reorder
        DO i = 1, n
            l = (n-i+1); k = deswap(W(i:l), mmloc);
            if (k == 1) exit;
            call swapcol(V, i, mmloc(1));
            call swapcol(V, l, mmloc(2));
        ENDDO
        
100     return
    end function svd2_dp
    integer(c_int) function svd2_sp(A, W, V, row, col) bind(c, name = '_sfsvd')
        USE UTL, ONLY : outer_prod, pythag, deswap, swapcol
        integer(c_int), intent(in):: row[value], col[value]
        real(c_float), intent(inout):: A(row,col)
        real(c_float), intent(inout):: V(col,col), W(col)
        real(SP), dimension(size(A,1)) :: tempm
        real(SP), dimension(size(A,2)) :: tempn, rv
        integer :: m, n, i, j, k, l, its, nm, maxits, mmloc(2)
        real(SP) :: anorm, c, f, g, h, s, scale, x, y, z
        
        maxits = 30; m = size(A,1); n = size(A,2);
        if (m == n) A = TRANSPOSE(A);
        if (n/=size(w).or.n/=size(v,1).or.n/=size(v,2)) then
            svd2_sp = -4; goto 100;
        endif
        
        g = 0.0; scale = 0.0;
        !$OMP DO
        DO i = 1, n  ! reduction to bidiagonal form
            l = i + 1; rv(i) = scale * g; g = 0.0; scale = 0.0;
            if (i <= m) then
                scale = SUM(abs(A(i:m,i)));
                if (scale /= 0.0) then
                    A(i:m,i) = A(i:m,i)/scale;
                    s = dot_product(A(i:m,i),A(i:m,i));
                    f = A(i,i); g = -SIGN(sqrt(s),f);
                    h = f*g-s; A(i,i) = f-g;
                    tempn(l:n) = MATMUL(A(i:m,i),A(i:m,l:n))/h;
                    A(i:m,l:n) = A(i:m,l:n) + outer_prod(A(i:m,i),tempn(l:n));
                    A(i:m,i) = scale*A(i:m,i)
                endif
            endif
            W(i) = scale*g; g = 0.0; scale = 0.0;
            if ((i <= m) .AND. (i /= n)) then
                scale = SUM(abs(A(i,l:n)));
                if (scale /= 0.0) then
                    A(i,l:n) = A(i,l:n)/scale
                    s = dot_product(A(i,l:n), A(i,l:n))
                    f = A(i,l); g = -SIGN(sqrt(s),f)
                    h = f*g-s; A(i,l) = f-g
                    rv(l:n) = A(i,l:n)/h
                    tempm(l:m) = MATMUL(A(l:m,l:n), A(i,l:n))
                    A(l:m,l:n) = A(l:m,l:n)+outer_prod(tempm(l:m),rv(l:n))
                    A(i,l:n) = scale*A(i,l:n)
                endif
            endif
        ENDDO
        !$OMP END DO
        anorm = MAXVAL(abs(w) + abs(rv));
        ! accumulation of right-hand transformations
        !$OMP DO
        DO i = n, 1, -1
            if (i < n) then
                if (g /= 0.0) then
                    V(l:n,i) = (A(i,l:n)/A(i,l))/g;
                    tempn(l:n) = MATMUL(A(i,l:n), V(l:n,l:n));
                    V(l:n,l:n) = V(l:n,l:n) + outer_prod(V(l:n,i), tempn(l:n));
                endif
                V(i,l:n) = 0.0; V(l:n,i) = 0.0;
            endif
            V(i,i) = 1.0; g = rv(i); l = i;
        ENDDO
        !$OMP END DO
        !-- diagonalization of the bidiagonal form: loop over singular values and over allowed iterations
        DO k = n, 1, -1
            DO its = 1, maxits
                DO l = k, 1, -1 ! Test for splitting
                    nm = l - 1;
                    if ((ABS(rv(l)) + anorm) == anorm) exit;
                    if ((ABS(W(nm)) + anorm) == anorm) then
                        c = 0.0; s = 1.0;
                        DO i = l, k
                            f = s * rv(i); rv(i) = c * rv(i);
                            if ((ABS(f)+anorm) == anorm) exit;
                            g = W(i); h = pythag(f, g);
                            W(i) = h; h = 1.0_SP/h;
                            c = (g * h); s = -(f*h);
                        ENDDO
                        EXIT
                    endif
                ENDDO
                z = W(k);
                if (l == k) then !convergence
                    if (z < 0.0) then
                        W(k) = -z; v(1:n,k) = -v(1:n,k);
                    endif
                    EXIT
                endif
                if (its == maxits) then
                    svd2_sp = 0; goto 100;
                endif
                x = W(l); nm = k - 1; y = W(nm);
                g = rv(nm); h = rv(k);
                f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_SP*h*y);
                g = pythag(f, 1.0_SP);
                f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
                !Next QR transformation
                c = 1.0; s = 1.0;
                !$OMP DO
                DO J = l, nm
                    i = j + 1; 
                    g = rv(i); y = W(i);
                    h = s*g; g = c*g;
                    z = pythag(f, h);
                    rv(j) = z; c = f/z; s = h/z;
                    f = (x*c)+(g*s); g = -(x*s)+(g*c);
                    h = y*s; y = y*c;
                    tempn(1:n) = V(1:n,j);
                    V(1:n,j) = V(1:n,j)*c + V(1:n,i)*s;
                    V(1:n,i) = -tempn(1:n)*s + V(1:n,i)*c;
                    z = pythag(f,h); W(j)=z;
                    if (z /= 0.0) then 
                        z = 1.0_SP/z; c = f*z; s = h*z;
                    endif
                    f = (c*g)+(s*y); x = -(s*g)+(c*y);
                ENDDO
                !$OMP END DO
                rv(l) = 0.0; rv(k) = f; W(k) = x;
            ENDDO
        ENDDO
        svd2_sp = 1;
        !reorder
        DO i = 1, n
            l = (n-i+1); k = deswap(W(i:l), mmloc);
            if (k == 1) exit;
            call swapcol(V, i, mmloc(1));
            call swapcol(V, l, mmloc(2));
        ENDDO
        
100     return
    end function svd2_sp
    
    !-- solve Ax = 0 using fast svd decomposition
    integer(c_int) function svdsolver_dp(A, xi, row, col, id) bind(c, name = '_dgsvdsv')
        USE UTL, ONLY : outer_prod, pythag, deswap, swapcol
        integer(c_int), intent(in):: row[value], col[value], id[value]
        real(c_double), intent(inout):: A(row,col)
        real(c_double), intent(out):: xi(col)
        real(DP):: V(col, col), W(col)
        real(DP), dimension(size(A,1)) :: tempm
        real(DP), dimension(size(A,2)) :: tempn, rv
        integer :: m, n, i, j, k, l, its, nm, maxits, mmloc(2)
        real(DP) :: anorm, c, f, g, h, s, scale, x, y, z
        
        maxits = 30; m = size(A,1); n = size(A,2);
        if (m == n) A = TRANSPOSE(A);
        if (n/=size(w).or.n/=size(v,1).or.n/=size(v,2)) then
            svdsolver_dp = -4; goto 100;
        endif
        
        g = 0.0; scale = 0.0;
        !$OMP DO
        DO i = 1, n  ! reduction to bidiagonal form
            l = i + 1; rv(i) = scale * g; g = 0.0; scale = 0.0;
            if (i <= m) then
                scale = SUM(abs(A(i:m,i)));
                if (scale /= 0.0) then
                    A(i:m,i) = A(i:m,i)/scale;
                    s = dot_product(A(i:m,i),A(i:m,i));
                    f = A(i,i); g = -SIGN(sqrt(s),f);
                    h = f*g-s; A(i,i) = f-g;
                    tempn(l:n) = MATMUL(A(i:m,i),A(i:m,l:n))/h;
                    A(i:m,l:n) = A(i:m,l:n) + outer_prod(A(i:m,i),tempn(l:n));
                    A(i:m,i) = scale*A(i:m,i)
                endif
            endif
            W(i) = scale*g; g = 0.0; scale = 0.0;
            if ((i <= m) .AND. (i /= n)) then
                scale = SUM(abs(A(i,l:n)));
                if (scale /= 0.0) then
                    A(i,l:n) = A(i,l:n)/scale
                    s = dot_product(A(i,l:n), A(i,l:n))
                    f = A(i,l); g = -SIGN(sqrt(s),f)
                    h = f*g-s; A(i,l) = f-g
                    rv(l:n) = A(i,l:n)/h
                    tempm(l:m) = MATMUL(A(l:m,l:n), A(i,l:n))
                    A(l:m,l:n) = A(l:m,l:n)+outer_prod(tempm(l:m),rv(l:n))
                    A(i,l:n) = scale*A(i,l:n)
                endif
            endif
        ENDDO
        !$OMP END DO
        
        anorm = MAXVAL(abs(w) + abs(rv));
        ! accumulation of right-hand transformations
        !$OMP DO
        DO i = n, 1, -1
            if (i < n) then
                if (g /= 0.0) then
                    V(l:n,i) = (A(i,l:n)/A(i,l))/g;
                    tempn(l:n) = MATMUL(A(i,l:n), V(l:n,l:n));
                    V(l:n,l:n) = V(l:n,l:n) + outer_prod(V(l:n,i), tempn(l:n));
                endif
                V(i,l:n) = 0.0; V(l:n,i) = 0.0;
            endif
            V(i,i) = 1.0; g = rv(i); l = i;
        ENDDO
        !$OMP END DO
        
        !-- diagonalization of the bidiagonal form: loop over singular values and over allowed iterations
        DO k = n, 1, -1
            DO its = 1, maxits
                DO l = k, 1, -1 ! Test for splitting
                    nm = l - 1;
                    if ((ABS(rv(l)) + anorm) == anorm) exit;
                    if ((ABS(W(nm)) + anorm) == anorm) then
                        c = 0.0; s = 1.0;
                        DO i = l, k
                            f = s * rv(i); rv(i) = c * rv(i);
                            if ((ABS(f)+anorm) == anorm) exit;
                            g = W(i); h = pythag(f, g);
                            W(i) = h; h = 1.0_DP/h;
                            c = (g * h); s = -(f*h);
                        ENDDO
                        EXIT
                    endif
                ENDDO
                z = W(k);
                if (l == k) then !convergence
                    if (z < 0.0) then
                        W(k) = -z; v(1:n,k) = -v(1:n,k);
                    endif
                    EXIT
                endif
                if (its == maxits) then
                    svdsolver_dp = 0; goto 100;
                endif
                x = W(l); nm = k - 1; y = W(nm);
                g = rv(nm); h = rv(k);
                f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_DP*h*y);
                g = pythag(f, 1.0_DP);
                f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
                !Next QR transformation
                c = 1.0; s = 1.0;
                !$OMP DO
                DO J = l, nm
                    i = j + 1; 
                    g = rv(i); y = W(i);
                    h = s*g; g = c*g;
                    z = pythag(f, h);
                    rv(j) = z; c = f/z; s = h/z;
                    f = (x*c)+(g*s); g = -(x*s)+(g*c);
                    h = y*s; y = y*c;
                    tempn(1:n) = V(1:n,j);
                    V(1:n,j) = V(1:n,j)*c + V(1:n,i)*s;
                    V(1:n,i) = -tempn(1:n)*s + V(1:n,i)*c;
                    z = pythag(f,h); W(j)=z;
                    if (z /= 0.0) then 
                        z = 1.0_DP/z; c = f*z; s = h*z;
                    endif
                    f = (c*g)+(s*y); x = -(s*g)+(c*y);
                ENDDO
                !$OMP END DO
                rv(l) = 0.0; rv(k) = f; W(k) = x;
            ENDDO
        ENDDO
        svdsolver_dp = 1;
        !reorder
        DO i = 1, n
            l = (n-i+1); k = deswap(W(i:l), mmloc);
            if (k == 1) exit;
            call swapcol(V, i, mmloc(1));
            call swapcol(V, l, mmloc(2));
        ENDDO
        if (id /= -1) then
            xi(:) = V(id, :)
        else 
            xi(:) = V(n, :)
        endif

100     return
    end function svdsolver_dp
    integer(c_int) function svdsolver_sp(A, xi, row, col, id) bind(c, name = '_sgsvdsv')
        USE UTL, ONLY : outer_prod, pythag, deswap, swapcol
        integer(c_int), intent(in):: row[value], col[value], id[value]
        real(c_float), intent(inout):: A(row,col)
        real(c_float), intent(out):: xi(col)
        real(SP):: V(col,col), W(col)
        real(SP), dimension(size(A,1)) :: tempm
        real(SP), dimension(size(A,2)) :: tempn, rv
        integer :: m, n, i, j, k, l, its, nm, maxits, mmloc(2)
        real(SP) :: anorm, c, f, g, h, s, scale, x, y, z
        
        maxits = 30; m = size(A,1); n = size(A,2);
        if (m == n) A = TRANSPOSE(A);
        if (n/=size(w).or.n/=size(v,1).or.n/=size(v,2)) then
            svdsolver_sp = -4; goto 100;
        endif
        
        g = 0.0; scale = 0.0;
        !$OMP DO
        DO i = 1, n  ! reduction to bidiagonal form
            l = i + 1; rv(i) = scale * g; g = 0.0; scale = 0.0;
            if (i <= m) then
                scale = SUM(abs(A(i:m,i)));
                if (scale /= 0.0) then
                    A(i:m,i) = A(i:m,i)/scale;
                    s = dot_product(A(i:m,i),A(i:m,i));
                    f = A(i,i); g = -SIGN(sqrt(s),f);
                    h = f*g-s; A(i,i) = f-g;
                    tempn(l:n) = MATMUL(A(i:m,i),A(i:m,l:n))/h;
                    A(i:m,l:n) = A(i:m,l:n) + outer_prod(A(i:m,i),tempn(l:n));
                    A(i:m,i) = scale*A(i:m,i)
                endif
            endif
            W(i) = scale*g; g = 0.0; scale = 0.0;
            if ((i <= m) .AND. (i /= n)) then
                scale = SUM(abs(A(i,l:n)));
                if (scale /= 0.0) then
                    A(i,l:n) = A(i,l:n)/scale
                    s = dot_product(A(i,l:n), A(i,l:n))
                    f = A(i,l); g = -SIGN(sqrt(s),f)
                    h = f*g-s; A(i,l) = f-g
                    rv(l:n) = A(i,l:n)/h
                    tempm(l:m) = MATMUL(A(l:m,l:n), A(i,l:n))
                    A(l:m,l:n) = A(l:m,l:n)+outer_prod(tempm(l:m),rv(l:n))
                    A(i,l:n) = scale*A(i,l:n)
                endif
            endif
        ENDDO
        !$OMP END DO
        anorm = MAXVAL(abs(w) + abs(rv));
        ! accumulation of right-hand transformations
        !$OMP DO
        DO i = n, 1, -1
            if (i < n) then
                if (g /= 0.0) then
                    V(l:n,i) = (A(i,l:n)/A(i,l))/g;
                    tempn(l:n) = MATMUL(A(i,l:n), V(l:n,l:n));
                    V(l:n,l:n) = V(l:n,l:n) + outer_prod(V(l:n,i), tempn(l:n));
                endif
                V(i,l:n) = 0.0; V(l:n,i) = 0.0;
            endif
            V(i,i) = 1.0; g = rv(i); l = i;
        ENDDO
        !$OMP END DO
        !-- diagonalization of the bidiagonal form: loop over singular values and over allowed iterations
        DO k = n, 1, -1
            DO its = 1, maxits
                DO l = k, 1, -1 ! Test for splitting
                    nm = l - 1;
                    if ((ABS(rv(l)) + anorm) == anorm) exit;
                    if ((ABS(W(nm)) + anorm) == anorm) then
                        c = 0.0; s = 1.0;
                        DO i = l, k
                            f = s * rv(i); rv(i) = c * rv(i);
                            if ((ABS(f)+anorm) == anorm) exit;
                            g = W(i); h = pythag(f, g);
                            W(i) = h; h = 1.0_SP/h;
                            c = (g * h); s = -(f*h);
                        ENDDO
                        EXIT
                    endif
                ENDDO
                z = W(k);
                if (l == k) then !convergence
                    if (z < 0.0) then
                        W(k) = -z; v(1:n,k) = -v(1:n,k);
                    endif
                    EXIT
                endif
                if (its == maxits) then
                    svdsolver_sp = 0; goto 100;
                endif
                x = W(l); nm = k - 1; y = W(nm);
                g = rv(nm); h = rv(k);
                f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_SP*h*y);
                g = pythag(f, 1.0_SP);
                f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
                !Next QR transformation
                c = 1.0; s = 1.0;
                !$OMP DO
                DO J = l, nm
                    i = j + 1; 
                    g = rv(i); y = W(i);
                    h = s*g; g = c*g;
                    z = pythag(f, h);
                    rv(j) = z; c = f/z; s = h/z;
                    f = (x*c)+(g*s); g = -(x*s)+(g*c);
                    h = y*s; y = y*c;
                    tempn(1:n) = V(1:n,j);
                    V(1:n,j) = V(1:n,j)*c + V(1:n,i)*s;
                    V(1:n,i) = -tempn(1:n)*s + V(1:n,i)*c;
                    z = pythag(f,h); W(j)=z;
                    if (z /= 0.0) then 
                        z = 1.0_SP/z; c = f*z; s = h*z;
                    endif
                    f = (c*g)+(s*y); x = -(s*g)+(c*y);
                ENDDO
                !$OMP END DO
                rv(l) = 0.0; rv(k) = f; W(k) = x;
            ENDDO
        ENDDO
        svdsolver_sp = 1;
        !reorder
        DO i = 1, n
            l = (n-i+1);
            k = deswap(W(i:l), mmloc);
            if (k == 1) exit;
            call swapcol(V, i, mmloc(1));
            call swapcol(V, l, mmloc(2));
        ENDDO
        if (id /= -1) then
            xi(:) = V(id, :)
        else 
            xi(:) = V(n, :)
        endif

100     return
    end function svdsolver_sp
    
    !-- LU decomposition with DP type
    integer(c_int) function dLU(A, piv, n) bind(c, name = '_dLU')
        USE utl, ONLY : imaxloc, swaprow, outer_prod
        integer(c_int),intent(in):: n[value]
        real(c_double),intent(inout):: A(n,n)
        integer(c_int),intent(inout):: piv(n)
        real(DP),parameter:: TINY = 1.0E-40_DP
        real(DP):: vv(n)  !stores the implicit scaling of each row
        integer:: j, j1, imax, info
        !DIR$ ATTRIBUTES ALIGN : 32 :: vv
        
        info = 1.0; !no row interchanges yet
        vv = maxval(abs(A), dim = 2); !get the max val of each row
        !loop over rows to get the implicit scaling info
        if (any(vv == 0.0)) then 
            dlu = 0; return; 
        endif
        vv = 1.0_DP/vv; !save the scaling
        !$OMP DO
        DO j = 1, n, 1
            imax = (j-1)+imaxloc(vv(j:n)*abs(A(j:n,j))); !find pivot r
            if (j /= imax) then !interchange rows?
                call swaprow(A, imax, j);
                info = -info;
                vv(imax) = vv(j);
            endif
            piv(j) = imax;
            if(A(j,j) == 0.0) A(j,j) = TINY;
            j1 = j + 1;
            A(j1:n,j) = A(j1:n,j)/A(j,j); !divide by the pivot element
            A(j1:n,j1:n) = A(j1:n,j1:n) - & !reduce remaining submatrix
                outer_prod(A(j1:n,j), A(j,j1:n));
        ENDDO
        !$OMP END DO
        dlu = info;
    end function dLU
    integer(c_int) function sLU(A, piv, n) bind(c, name = '_sLU')
        USE utl, ONLY : imaxloc, swaprow, outer_prod
        integer(c_int),intent(in):: n[value]
        real(c_float),intent(inout):: A(n,n)
        integer(c_int),intent(inout):: piv(n)
        real(SP),parameter:: TINY = 1.0E-40_SP
        real(SP):: vv(n)  !stores the implicit scaling of each row
        integer:: j, j1, imax, info
        !DIR$ ATTRIBUTES ALIGN : 32 :: vv
        
        info = 1.0; !no row interchanges yet
        vv = maxval(abs(A), dim = 2); !get the max val of each row
        !loop over rows to get the implicit scaling info
        if (any(vv == 0.0)) then 
            slu = 0; return; 
        endif
        vv = 1.0_DP/vv; !save the scaling
        !$OMP DO
        DO j = 1, n, 1
            imax = (j-1)+imaxloc(vv(j:n)*abs(A(j:n,j))); !find pivot r
            if (j /= imax) then !interchange rows?
                call swaprow(A, imax, j);
                info = -info;
                vv(imax) = vv(j);
            endif
            piv(j) = imax;
            if(A(j,j) == 0.0) A(j,j) = TINY;
            j1 = j + 1;
            A(j1:n,j) = A(j1:n,j)/A(j,j); !divide by the pivot element
            A(j1:n,j1:n) = A(j1:n,j1:n) - & !reduce remaining submatrix
                outer_prod(A(j1:n,j), A(j,j1:n));
        ENDDO
        !$OMP END DO
        slu = info;
    end function sLU
    
    INTEGER(c_int) FUNCTION dcholesky(A, p, n) bind(c,name='_dcholy2')
        integer(c_int),intent(in):: n[value]
        real(c_double),intent(inout):: A(n,n)
        real(c_double),intent(out)::p(n)
        real(c_double) :: sum
        integer:: i, j, inc
        
        inc = 1;
        j = 1;
        !$OMP DO
        do i = 1, n, inc
            sum = A(i,i) - dot_product(A(i,1:i-1),A(i,1:i-1));
            if (sum <= 0.0) then
                inc = n; j = -i;
            end if
            p(i) = sqrt(sum);
            A(i:n,i) = (A(i,i:n)-matmul(A(i:n,1:i-1),A(i,1:i-1)))/p(i);
        enddo
        !$OMP END DO
        dcholesky = j;
        return;
    END FUNCTION dcholesky
    INTEGER(c_int) FUNCTION dcholesky2(A, n) bind(c,name='_dcholy')
        integer(c_int),intent(in):: n[value]
        real(c_double),intent(inout):: A(n,n)
        real(c_double) :: sum
        integer:: i, j, inc
        
        inc = 1;
        j = 1;
        !$OMP DO
        do i = 1, n, inc
            sum = A(i,i) - dot_product(A(i,1:i-1),A(i,1:i-1));
            if (sum <= 0.0) then
                inc = n; j = -i;
            end if
            A(i:n,i) = (A(i,i:n)-matmul(A(i:n,1:i-1),A(i,1:i-1)))/sqrt(sum);
        enddo
        !$OMP END DO
        dcholesky2 = j;
        return;
    END FUNCTION dcholesky2
    INTEGER(c_int) FUNCTION scholesky2(A, n) bind(c,name='_scholy')
        integer(c_int),intent(in):: n[value]
        real(c_float),intent(inout):: A(n,n)
        real(c_float) :: sum
        integer:: i, j, inc
        
        inc = 1;
        j = 1;
        !$OMP DO
        do i = 1, n, inc
            sum = A(i,i) - dot_product(A(i,1:i-1),A(i,1:i-1));
            if (sum <= 0.0) then
                inc = n; j = -i;
            end if
            A(i:n,i) = (A(i,i:n)-matmul(A(i:n,1:i-1),A(i,1:i-1)))/sqrt(sum);
        enddo
        !$OMP END DO
        scholesky2 = j;
        return;
    END FUNCTION scholesky2
    
    ! solve Ax = b, where A is a positive-definite symmetrix matrix
    ! only the lower triangle of A is accessed.
    SUBROUTINE dcholeskyslv(A, bx, n) bind(c, name = '_dcholsv')
        integer(c_int),intent(in):: n[value]
        real(c_double),intent(in):: A(n,n)
        real(c_double),intent(inout):: bx(n)
        integer:: i
        
        !$OMP DO
        do i = 1, n ! solve L*y = b,storing y in b in place
            bx(i) = bx(i) - dot_product(A(i,1:i-1), bx(1:i-1));
            bx(i) = bx(i) / A(i,i);
        end do
        !$OMP END DO
        !$OMP DO
        do i = n, 1, -1 ! solve L*y = b,storing y in b in place
            bx(i) = bx(i) - dot_product(A(i+1:n,i), bx(i+1:n));
            bx(i) = bx(i) / A(i,i);
        end do
        !$OMP END DO
        
    END SUBROUTINE dcholeskyslv
    
    INTEGER(c_int) FUNCTION dfwsubt(L, b, n) bind(c,name='_dfwsubt')
        integer(c_int),intent(in):: n[value]
        real(c_double),intent(in):: L(n,n)
        real(c_double),intent(inout)::b(n)
        real(c_double) :: temp
        integer:: i, j
        
        if (L(1,1) == 0) then
            dfwsubt = -1; return;
        else
            b(1) = b(1)/L(1,1);
        endif
        temp = 0;
        do i = 2, n
            if (L(i,i) == 0) then
                dfwsubt = -i; return;
            end if
            !$OMP DO
            do j = 1, i
                temp = temp + L(i,j)*b(j);
            enddo
            !$OMP END DO
            b(i) = (b(i)-temp)/L(i,i);
        enddo
        dfwsubt = 1;
    END function dfwsubt
    INTEGER(c_int) FUNCTION dbwsubt(L, b, n) bind(c,name='_dbwsubt')
        integer(c_int),intent(in):: n[value]
        real(c_double),intent(in):: L(n,n)
        real(c_double),intent(inout)::b(n)
        real(c_double) :: temp
        integer:: i, j
        
        if (L(n,n) == 0) then
            dbwsubt = -1; return;
        else
            b(n) = b(n)/L(n,n);
        endif
        temp = 0;
        do i = n, 2, -1
            if (L(i,i) == 0) then
                dbwsubt = -i; return;
            endif
            !$OMP DO
            do j = i + 1, n
                temp = temp + L(i,j)*b(j);
            enddo
            !$OMP END DO
            b(i) = (b(i)-temp)/L(i,i);
        enddo
        dbwsubt = 1;
    END function dbwsubt
    
end module les