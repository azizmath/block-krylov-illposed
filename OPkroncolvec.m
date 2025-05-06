function w = OPkroncolvec(u,A1,A2,tflag)

[m1,n1] = size(A1);
[m2,n2] = size(A2);
if strcmpi(tflag,'size')
    w(1) = m1*m2;
    w(2) = n1*n2;
elseif strcmpi(tflag,'notransp')
    w = zeros(n1*n2, 3);
    u = reshape(u, n1*n2, 3);
    for i = 1:3
        utemp = reshape(u(:,i), n1, n2);
        wtemp = A1*utemp*A2';
        w(:,i) = wtemp(:);
    end
    w = w(:);
else
    w = zeros(m1*m2, 3);
    u = reshape(u, m1*m2, 3);
    for i = 1:3
        utemp = reshape(u(:,i), m1, m2);
        wtemp = A1'*utemp*A2;
        w(:,i) = wtemp(:);
    end
    w = w(:);
end
