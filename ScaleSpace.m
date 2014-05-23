function Octave_outputs = ScaleSpace( I_input,DOG_Times,Octaves,Sigma,k )
%A funcion that computes a scal space using Difference of Gaussian and returns multiple output Images based on Number of Octaves and Number of DOF computed 
%I_input: the orginal input Image
%DOG_Times: the number of times that Difference of Gaussian should be computed
%Octaves: How many octaves should be computed
%Sigma: is the sigma value for Gaussain filter
%k: is the scale value that will be multiplied by sigma to compute the Gaussian and find DOG 


DOG = cell(1,DOG_Times); %Create cell array of empty matrices.
counter = 0; %for mutlible Octaves different k

for j=1:Octaves
       %First Gaussian Filter
        GaussainSize = ceil(7*(Sigma * (k^counter)));
        H = fspecial('gaussian' , [GaussainSize GaussainSize], Sigma * (k^counter));
        if(j>1) %Down scale the image after the first octave
            I_input = I_input(1:2:end,1:2:end);
        end
    for i=1:DOG_Times
      L1 = conv2(I_input,H,'same');
      GaussainSize = ceil(7*(Sigma*(k^(i+counter))));
      H = fspecial('gaussian' , [GaussainSize GaussainSize], Sigma*(k^(i+counter)));
      L2 = conv2(I_input,H,'same');
      %DOG and save it
      Octave_outputs(j).DOG{i} = struct('DOG',L2 - L1,'L2',L2,'L1',L1,'S',Sigma*(k^(i+counter)));
    end
    counter = counter + 1;
end

end

