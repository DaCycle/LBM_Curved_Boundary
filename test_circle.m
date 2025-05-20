function s = test_circle(x,y,R,center_x,center_y)

mag = sqrt((x-center_x)^2+(y-center_y)^2);
if mag > R
    s = 0;
else
    s = 1;
end