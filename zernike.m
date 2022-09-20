% ----------------------------------------------------------------------
% Jichong Zhou  9.20.2022
% 
% Reference: 
% 1. Efficient and robust recurrence relations for the Zernike circle 
%    polynomials and their derivatives in Cartesian coordinates
% ----------------------------------------------------------------------

x = 0.4;
y = 0.8;
MaxOrd = 10;
[Zern, PZPX, PZPY] = ZernikePatial(x, y, MaxOrd);
Norm = NormalizeCoefficient(MaxOrd);


function Norm = NormalizeCoefficient(MaxOrd)

NumItems = (2 + MaxOrd) * (MaxOrd + 1) / 2;
Norm = zeros(NumItems, 1);
% order 0
Norm(1) = 1;
% record index of Norm
CurItem = 1;
for n = 1 : MaxOrd
    for m = 0 : n
        CurItem = CurItem + 1;
        if (n == 2 * m)
            Norm(CurItem) = sqrt(n + 1);
        else
            Norm(CurItem) = sqrt(2 * (n + 1));
        end
    end
end
end


function Zern = ZernikeItem(x, y, MaxOrd)

NumItems = (2 + MaxOrd) * (MaxOrd + 1) / 2;
Zern = zeros(NumItems, 1);

% order 0
Zern(1) = 1;
% order 1
Zern(2) = y;
Zern(3) = x;

CurItem = 3;

% start from order 2
for n = 2 : MaxOrd

    cond1 = (n - 1) / 2;
    cond2 = n / 2;
    for m = 0 : n
        CurItem = CurItem + 1;  % update the current index of Zern
        if (m == 0)
            Zern(CurItem) = x * Zern(CurItem - n) + y * Zern(CurItem - 1);

        elseif (m == n)
            Zern(CurItem) = x * Zern(CurItem - n - 1) - y * Zern(CurItem - 2 * n);

        elseif (mod(n,2) == 0 && m == cond2)
            Zern(CurItem) = 2 * x * Zern(CurItem - n) + 2 * y * Zern(CurItem - n -1) - Zern(CurItem - 2 * n);

        elseif (mod(n,2) == 1 && m == cond1)
            Zern(CurItem) = y * Zern(CurItem - 2 * m - 1) + x * Zern(CurItem - n - 1) - y * Zern(CurItem - 2 * m) - Zern(CurItem - 2 * n);

        elseif (mod(n,2) == 1 && m == (cond1 + 1))
            Zern(CurItem) = x * Zern(CurItem - n) + y * Zern(CurItem - 2 * m - 1) + x * Zern(CurItem - n - 1) - Zern(CurItem - 2 * n);
        
        else 
            Zern(CurItem) = x * Zern(CurItem - n) + y * Zern(CurItem - 2 * m - 1) + x * Zern(CurItem - n - 1) - y * Zern(CurItem - 2 * m) - Zern(CurItem - 2 * n);

        end

    end
end
end


function [Zern, PZPX, PZPY] = ZernikePatial(x, y, MaxOrd)
NumItems = (2 + MaxOrd) * (MaxOrd + 1) / 2;
Zern = zeros(NumItems, 1);
PZPX = zeros(NumItems, 1);
PZPY = zeros(NumItems, 1);

% order 0
Zern(1) = 1;
PZPX(1) = 0;
PZPY(1) = 0;
% order 1
Zern(2) = y;
Zern(3) = x;
PZPX(2) = 0;
PZPX(3) = 1;
PZPY(2) = 1;
PZPY(3) = 0;

CurItem = 3;

% start from order 2
for n = 2 : MaxOrd
    
    cond1 = (n - 1) / 2;
    cond2 = n / 2;
    for m = 0 : n
        CurItem = CurItem + 1;  % update the current index of Zern
        if (m == 0)
            Zern(CurItem) = x * Zern(CurItem - n) + y * Zern(CurItem - 1);
            PZPX(CurItem) = n * Zern(CurItem - n);
            PZPY(CurItem) = n * Zern(CurItem - 1);
            
        elseif (m == n)
            Zern(CurItem) = x * Zern(CurItem - n - 1) - y * Zern(CurItem - 2 * n);
            PZPX(CurItem) = n * Zern(CurItem - n - 1);
            PZPY(CurItem) = - n * Zern(CurItem - 2 * n);
            
        elseif (mod(n,2) == 0 && m == cond2)
            Zern(CurItem) = 2 * x * Zern(CurItem - n) + 2 * y * Zern(CurItem - n -1) - Zern(CurItem - 2 * n);
            PZPX(CurItem) = 2 * n * Zern(CurItem - n) + PZPX(CurItem - 2 * n);
            PZPY(CurItem) = 2 * n * Zern(CurItem - n -1) + PZPY(CurItem - 2 * n);
            
        elseif (mod(n,2) == 1 && m == cond1)
            Zern(CurItem) = y * Zern(CurItem - 2 * m - 1) + x * Zern(CurItem - n - 1) - y * Zern(CurItem - 2 * m) - Zern(CurItem - 2 * n);
            PZPX(CurItem) = n * Zern(CurItem - n - 1) - PZPX(CurItem - 2 * n);
            PZPY(CurItem) = n * Zern(CurItem - 2 * m - 1) - n * Zern(CurItem - 2 * m) + PZPY(CurItem - 2 * n);
            
        elseif (mod(n,2) == 1 && m == (cond1 + 1))
            Zern(CurItem) = x * Zern(CurItem - n) + y * Zern(CurItem - 2 * m - 1) + x * Zern(CurItem - n - 1) - Zern(CurItem - 2 * n);
            PZPX(CurItem) = n * Zern(CurItem - n) + n * Zern(CurItem - n - 1) + PZPX(CurItem - 2 * n);
            PZPY(CurItem) = n * Zern(CurItem - 2 * m - 1) + PZPY(CurItem - 2 * n);
            
        else 
            Zern(CurItem) = x * Zern(CurItem - n) + y * Zern(CurItem - 2 * m - 1) + x * Zern(CurItem - n - 1) - y * Zern(CurItem - 2 * m) - Zern(CurItem - 2 * n);
            PZPX(CurItem) = n * Zern(CurItem - n) + n * Zern(CurItem - n - 1) + PZPX(CurItem - 2 * n);
            PZPY(CurItem) = n * Zern(CurItem - 2 * m - 1) - n * Zern(CurItem - 2 * m) + PZPY(CurItem - 2 * n);
            
        end

    end
end
end