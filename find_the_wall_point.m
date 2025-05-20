function C = find_the_wall_point(x1,y1,x2,y2,R,center_x,center_y)
%%% This function finds the coordinates [x; y] of the point that is the intercept point
%%% Given any given circle (radius is R, coordinates of the center is [center_x; center_y]),
%%% a point located inside or on the circle ([x1; y1]), and
%%% another point located outside of the circle ([x2; y2]), this function
%%% finds the coordinates of the intercept point ([x; y]) between the
%%% circle and the straight line connecting point [x1; y1] and point [x2; y2].

epsilon=1e-6;
iter = 0;

x_low=x1;
y_low=y1;
x_high=x2;  %Differ by directions
y_high=y2;      %Differ by directions

x_w=(x_low+x_high)/2;
y_w=(y_low+y_high)/2;
R_w=sqrt((x_w-center_x)^2+(y_w-center_y)^2);
i=0;
while abs(R_w-R)>=epsilon || iter<1e3
    iter = iter + 1;
    if R_w>R
        x_high=x_w;
        y_high=y_w;
    else
        x_low=x_w;
        y_low=y_w;
    end
    x_w=(x_low+x_high)/2;
    y_w=(y_low+y_high)/2;
    R_w=sqrt((x_w-center_x)^2+(y_w-center_y)^2);
    i=i+1;
    R_w;    
end

C=[x_w;y_w];