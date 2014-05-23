close all; clear all; clc;
image = imread('Penguins.jpg');
image2= image;
[Irows, Icols, chan]=size(image);
if(chan == 3)
    image = rgb2gray(image); %convert the input image to gray level (if it not gray already)
end
I = double(image);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scale space peak selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigma = 1.6;
k=sqrt(2);
DOG_Times = 4;
Octaves = 3;
Octave_outputs = ScaleSpace( I,DOG_Times,Octaves,Sigma,k );
%Peaks Detection
counter = 1;
r_tr = 10;
OriBins = 36;
winfactor  = 1.5 ;
for i=1:Octaves
    for j=1:DOG_Times
        if(j~=DOG_Times && j~=1) %avoid min and max scale output because they dont have adjacent scales
            current = Octave_outputs(i).DOG{j}.DOG;
            below = Octave_outputs(i).DOG{j-1}.DOG;
            above = Octave_outputs(i).DOG{j+1}.DOG;
            L = Octave_outputs(i).DOG{j}.L1;
            S = Octave_outputs(i).DOG{j}.S;
            [rows, cols]=size(Octave_outputs(i).DOG{j}.DOG);
            for r=1:rows
                for c=1:cols
                   if (c~=1 && c~=cols && r~=1 && r~=rows) %avoid edges becasue they dont have 26 pixels around them in the current and adjecnt scales
                    %compare the current pixle with the other adjacent 26 pixles
                    X = [current(r,c) current(r+1,c) current(r-1,c) current(r,c+1) current(r,c-1) current(r+1,c+1) current(r-1,c-1) current(r+1,c-1) current(r-1,c+1) below(r,c) below(r+1,c) below(r-1,c) below(r,c+1) below(r,c-1) below(r+1,c+1) below(r-1,c-1) below(r+1,c-1) below(r-1,c+1) above(r,c) above(r+1,c) above(r-1,c) above(r,c+1) above(r,c-1) above(r+1,c+1) above(r-1,c-1) above(r+1,c-1) above(r-1,c+1)];
                    [max_X,I_max_X]= max(X);
                    [min_X,I_min_X]= min(X);
                    if(I_max_X == 1 || I_min_X == 1) %Pick the point if it is the maximum or the minimum
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Key point localization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        Dx= (current(r,c+1) - current(r,c-1))/2;
                        Dy= (current(r+1,c) - current(r-1,c))/2;
                        Dz= (above(r,c) - below(r,c))/2;
                        D_firstDerivative = [Dx; Dy; Dz];
                        Dxx = current(r,c-1) - 2*current(r,c) + current(r,c+1);
                        Dyy = current(r-1,c) - 2*current(r,c) + current(r+1,c);
                        Dzz = below(r,c) - 2*current(r,c) + above(r,c);
                        Dxy = ((current(r+1,c+1) - current(r+1,c-1)) - (current(r-1,c+1) - current(r-1,c-1))) / 4;
                        Dxz = ((above(r,c+1) - above(r,c-1)) - (below(r,c+1) - below(r,c-1))) / 4;
                        Dyz = ((above(r+1,c) - above(r-1,c)) - (below(r+1,c) - below(r-1,c))) / 4;
                        D_secondDerivative = [Dxx Dxy Dxz;Dxy Dyy Dyz;Dxz Dyz Dzz];
                        extremum = -inv(D_secondDerivative)*D_firstDerivative;
                        if(I_max_X == 1)
                            D=max_X;
                        else
                            D=min_X;
                        end
                        %Maxima and Minima with offest
                        %Interpolated estimate for the location of the extremum.
                        x=c;
                        y=r;
                        z=j;
                        if(extremum(1) > 0.5 && ceil(c+extremum(1))<cols)
                            extremum(1) = ceil(x+extremum(1));
                            x= extremum(1);
                        end
                        if(extremum(2) > 0.5 && ceil(r+extremum(2))<rows)
                            extremum(2) = ceil(y+extremum(2));
                            y= extremum(2);
                        end
                        if(extremum(3) > 0.5 && ceil(j+extremum(3))<=DOG_Times)
                            extremum(3) = ceil(z+extremum(3));
                            z= extremum(3);
                            L = Octave_outputs(i).DOG{z}.L1;
                            S = Octave_outputs(i).DOG{z}.S;
                        end 
                        D_extremum = D + 0.5 * transp(D_firstDerivative) * extremum;
                        Dxx = current(y,x-1) - 2*current(y,x) + current(y,x+1);
                        Dyy = current(y-1,x) - 2*current(y,x) + current(y+1,x);
                        Dxy = ((current(y+1,x+1) - current(y+1,x-1)) - (current(y-1,x+1) - current(y-1,x-1))) / 4;
                        if (abs(D_extremum) >= (0.03*255)) %all extrema with a value of |D(ˆx)| less than 0.03 were discarded
                            H = [Dxx Dxy;Dxy Dyy]; %Eliminating edge responses
                            if ( (trace(H).^2)/det(H) < ((r_tr+1).^2)/r_tr ) %Eliminate key points if r_tr >10
                               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Orientation Assignment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                               sigma = winfactor * Sigma * (2.^(z/k));
                               radius = round(3*sigma);
                               hist = zeros(1,OriBins);
                               angle = atan((L(y,x+1)-L(y,x-1))/(L(y+1,x)-L(y-1,x)));
                               %Create a window around the key point and create histogram around this window
                               for new_r = y-radius:y+radius
                                for new_c = x-radius:x+radius
                                    if (new_r > 1 && new_c > 1 && new_r < rows - 2 && new_c < cols - 2)
                                        gradVal = sqrt(power(L(new_r+1,new_c)-L(new_r-1,new_c),2)+power(L(new_r,new_c+1)-L(new_r,new_c-1),2));
                                        distsq = (y - new_r).^2 + (x - new_c).^2;
                                        if (gradVal > 0.0  &&  distsq < radius * radius + 0.5)
                                            weight = exp(- distsq / (2.0 * sigma * sigma));
                                            % Ori is in range of -pi to pi
                                            angle = atan((L(new_r,new_c+1)-L(new_r,new_c-1))/(L(new_r+1,new_c)-L(new_r-1,new_c)));
                                            bin =  ceil(OriBins * (angle + pi + 0.001) / (2.0 * pi));
                                            bin = min(bin, OriBins);
                                            hist(bin) = hist(bin) + weight * gradVal;
                                        end
                                    end
                                end
                               end
                               %extract dominant orientation
                               [~, in] = max(hist);
                               %conversion from bins to angles
                               angle = 2*pi*(in)/OriBins;
                               temp_Peaks(counter) = [struct('r_peak',y,'c_peak',x,'i_octave',i,'ori',angle,'S',S)];
                               counter = counter +1;
                            end   
                        end            
                    end
                   end
                end
            end
        end
    end
end

%Creat the final matrix to draw the output
for i=1:length(temp_Peaks)
    c = temp_Peaks(i).c_peak * power(2,(temp_Peaks(i).i_octave)-1); %to get the original position in case the octave was not the first
    r = temp_Peaks(i).r_peak * power(2,(temp_Peaks(i).i_octave)-1);
    locs(i,:) = [r c temp_Peaks(i).S temp_Peaks(i).ori];
end

image = image2;
computed_features = counter - 1
showkeys(image, locs);

