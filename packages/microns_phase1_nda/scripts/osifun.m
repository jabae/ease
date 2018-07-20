function osi = osifun(s) 
numerator = (1-[s.r]).*([s.w2] - [s.w3].*[s.r]);
denominator = 2*[s.w1] + (1+[s.r]).*([s.w2] + [s.w3].*[s.r]);
osi = numerator./denominator;
end