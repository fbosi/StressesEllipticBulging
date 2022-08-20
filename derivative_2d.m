function dz = derivative_2d(x,y,z,m,n)

% Function for calculating mixed numerical derivatives d^(m+n)z/dx^m dy^n upto 4th order

dx = x(2,1)-x(1,1);
dy = y(1,2)-y(1,1);

[length_x,length_y]=size(z);

if n==0
    if m==1
        for j=3:length_y-2
            for i=3:length_x-2
                dz(i,j) = (1/2).*dx.^(-1).*((-1).*z((-1)+i,j)+z(1+i,j));
            end
        end
    end
    
    if m==2
        for j=3:length_y-2
            for i=3:length_x-2
                dz(i,j) = dx.^(-2).*(z(-1+i,j) + (-2).*z(i,j) + z(1+i,j));
            end
        end
    end
    
    if m==3
        for j=3:length_y-2
            for i=3:length_x-2
                dz(i,j) = (1/2).*dx.^(-3).*((-1).*z((-2)+i,j)+2.*z((-1)+i,j)+(-2).*z(1+i,j)+z(2+i,j));
            end
        end
    end

    if m==4
        for j=3:length_y-2
            for i=3:length_x-2
                dz(i,j) = dx.^(-4).*(z((-2)+i,j)+(-4).*z((-1)+i,j)+6.*z(i,j)+(-4).*z(1+i,j)+z(2+i,j));
            end
        end
    end
    
end


if n==1
    if m==0
        for j=3:length_y-2
            for i=3:length_x-2
                dz(i,j) = (1/2).*dy.^(-1).*((-1).*z(i,(-1)+j)+z(i,1+j));
            end
        end
    end
    
    if m==1
        for j=3:length_y-2
            for i=3:length_x-2
                dz(i,j) = (1/4).*dx.^(-1).*dy.^(-1).*(z((-1)+i,(-1)+j)+(-1).*z((-1)+i,1+j)+(-1).*z(1+i,(-1)+j)+z(1+i,1+j));
            end
        end
    end
    
    if m==2
        for j=3:length_y-2
            for i=3:length_x-2
                dz(i,j) = (1/2).*dx.^(-2).*dy.^(-1).*((-1).*z((-1)+i,(-1)+j)+z((-1)+i,1+j)+2.*z(i,(-1)+j)+(-2).*z(i,1+j)+(-1).*z(1+i,(-1)+j)+z(1+i,1+j));
            end
        end
    end

    if m==3
        for j=3:length_y-2
            for i=3:length_x-2
                dz(i,j) = (1/4).*dx.^(-3).*dy.^(-1).*(z((-2)+i,(-1)+j)+(-1).*z((-2)+i,1+j)+(-2).*z((-1)+i,(-1)+j)+2.*z((-1)+i,1+j)+2.*z(1+i,(-1)+j)+(-2).*z(1+i,1+j)+(-1).*z(2+i,(-1)+j)+z(2+i,1+j));
            end
        end
    end
    
    
end

if n==2
    if m==0
        for j=3:length_y-2
            for i=3:length_x-2
                dz(i,j) = dy.^(-2).*(z(i,(-1)+j)+(-2).*z(i,j)+z(i,1+j));
            end
        end
    end
    
    if m==1
        for j=3:length_y-2
            for i=3:length_x-2
                dz(i,j) = (1/2).*dx.^(-1).*dy.^(-2).*((-1).*z((-1)+i,(-1)+j)+2.*z((-1)+i,j)+(-1).*z((-1)+i,1+j)+z(1+i,(-1)+j)+(-2).*z(1+i,j)+z(1+i,1+j));
            end
        end
    end
    
    if m==2
        for j=3:length_y-2
            for i=3:length_x-2
                dz(i,j) = dx.^(-2).*dy.^(-2).*(z((-1)+i,(-1)+j)+(-2).*z((-1)+i,j)+z((-1)+i,1+j)+(-2).*z(i,(-1)+j)+4.*z(i,j)+(-2).*z(i,1+j)+z(1+i,(-1)+j)+(-2).*z(1+i,j)+z(1+i,1+j));
            end
        end
    end
end

if n==3
    if m==0
        for j=3:length_y-2
            for i=3:length_x-2
                dz(i,j) = (1/2).*dy.^(-3).*((-1).*z(i,(-2)+j)+2.*z(i,(-1)+j)+(-2).*z(i,1+j)+z(i,2+j));
            end
        end
    end
    
    if m==1
        for j=3:length_y-2
            for i=3:length_x-2
                dz(i,j) = (1/4).*dx.^(-1).*dy.^(-3).*(z((-1)+i,(-2)+j)+(-2).*z((-1)+i,(-1)+j)+2.*z((-1)+i,1+j)+(-1).*z((-1)+i,2+j)+(-1).*z(1+i,(-2)+j)+2.*z(1+i,(-1)+j)+(-2).*z(1+i,1+j)+z(1+i,2+j));
            end
        end
    end
end

if n==4
    if m==0
        for j=3:length_y-2
            for i=3:length_x-2
                dz(i,j) = dy.^(-4).*(z(i,(-2)+j)+(-4).*z(i,(-1)+j)+6.*z(i,j)+(-4).*z(i,1+j)+z(i,2+j));
            end
        end
    end
    
end

dz = dz(3:length_x-2,3:length_y-2);

end

