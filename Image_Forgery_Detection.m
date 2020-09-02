%Clearing All Variables in MATLAB Worksapce.
clear all
 
%Clearing MATLAB Screen.
clc
 
BLOCK_S = 20;
OVERLAP_S = 10;
MSB_TH = 60;
 
%Reading an Image & saving it in 'original_image' variable.
original_image=imread('Original.jpg');
 
%Cropping The Read Image by 200x200 Pixels & saving it in 'cropped_image'
%variable.
cropped_image = imcrop(original_image,[0 0 200 200]);
imwrite(cropped_image,'Cropped.jpg')
 
%Copy a small part(40x40 Pixels)from the Cropped image(200x200 Pixels)
%starting from pixel#100 on x-axis & on y-axis.
image_copy=imcrop(cropped_image,[90 90 39 39]);
 
%Paste the copied image on the cropped image starting from pixel#20
%on x-axis & on y-axis
for i=1:length(image_copy)
    for j=1:length(image_copy)
        for k=1:3
            cropped_image(i+20,j+20,k)=image_copy(i,j,k);
        end
    end
end
 
%Saving the forged image.
imwrite(cropped_image,'forged.jpg')
 
%Display the cropped & forged image.
subplot(1,2,1)
imshow(imread('cropped.jpg'))
xlabel('Original Cropped Image')
subplot(1,2,2)
imshow(cropped_image)
xlabel('Forged Image')
 
%Decomposing the forged image into R , G & B colors.
R  = cropped_image(:,:,1);
G  = cropped_image(:,:,2);
B  = cropped_image(:,:,3);
 
%Delete No more necessary variables.
clear original_image cropped_image image_copy1 image_copy
 
 
%--------------------------------------------------------------------------
% Dividing R,G & B Chanells into blocks , each block has a 
% dimension of BLOCK_S with overlapping of OVERLAP_S pixels .
%--------------------------------------------------------------------------
 
%Using Block Size of 40x40 Pixels.
block_size = [BLOCK_S BLOCK_S];
 
%Using overlapping size of 20x20 Pixels.
overlapping_size = [OVERLAP_S OVERLAP_S];
 
tot_block_row = ((size(R, 1) - BLOCK_S) / OVERLAP_S) + 1;
 
tot_block = tot_block_row * tot_block_row;
tot_block_pair = tot_block * (tot_block - 1) / 2;
retain_block_pair = floor(tot_block_pair / 16);
 
startPos = 0: overlapping_size: size(R, 1) - block_size;
 
%--------------------------------------------------------------------------
% Dividing R into overlpping blocks
%--------------------------------------------------------------------------
Sub_R = cell(numel(startPos), numel(startPos));
for i = 1:numel(startPos)
    for j = 1:numel(startPos)
        Sub_R{i,j} = R(startPos(i)+(1:block_size), startPos(j)+(1:block_size));
    end
end
 
 
%--------------------------------------------------------------------------
% Finding Histogram of each block for "Sub_R" color by finding the historam of
% each block.
%--------------------------------------------------------------------------
 
for i = 1:tot_block_row
    for j = 1:tot_block_row
        temp = single(cell2mat(Sub_R(i,j)));
        Hist_R{i,j}=lbp(temp);
    end
end
%Delete No more necessary variables.
clear Sub_R temp
 
%--------------------------------------------------------------------------
% Dividing G into overlpping blocks
%--------------------------------------------------------------------------
Sub_G = cell(numel(startPos), numel(startPos));
for i = 1:numel(startPos)
    for j = 1:numel(startPos)
        Sub_G{i,j} = G(startPos(i)+(1:block_size), startPos(j)+(1:block_size));
    end
end
 
%--------------------------------------------------------------------------
% Finding Histogram of each block for "Sub_G" color by finding the historam of
% each block.
%--------------------------------------------------------------------------
for i = 1:tot_block_row
    for j = 1:tot_block_row
        temp = single(cell2mat(Sub_G(i,j)));
        Hist_G{i,j}=lbp(temp);
    end
end
%Delete No more necessary variables.
clear Sub_G temp
 
%--------------------------------------------------------------------------
% Dividing B into overlpping blocks
%--------------------------------------------------------------------------
Sub_B = cell(numel(startPos), numel(startPos));
for i = 1:numel(startPos)
    for j = 1:numel(startPos)
        Sub_B{i,j} = B(startPos(i)+(1:block_size), startPos(j)+(1:block_size));
    end
end
 
%--------------------------------------------------------------------------
% Finding Histogram of each block for "Sub_B" color by finding the historam of
% each block.
%--------------------------------------------------------------------------
for i = 1:tot_block_row
    for j = 1:tot_block_row
        temp = single(cell2mat(Sub_B(i,j)));
        Hist_B{i,j}=lbp(temp);
    end
end
%Delete No more necessary variables.
clear temp Sub_B
 
%--------------------------------------------------------------------------
%Finding the Distance('norm') for each pair of blocks in the 'Sub_R' color 
%and saving the result in 'Distance_R' variable.
%--------------------------------------------------------------------------
%Converting the Sub_R from 'Matrix' shape to 'Single Column' shape.
Hist_R = reshape(Hist_R,1,[])';
counter=1;
x_R=1;
y_R=1;
index_R=1;
for i=1:(tot_block)-1           
    for j=(i+1):(tot_block)     
       temp = kolmogorov_smirnov_distance(double(Hist_R{i,1}), double(Hist_R{j,1}));
       x_R(counter)=i;
       y_R(counter)=j;
       Distance_R(counter)=temp;
       index_R(counter)=counter;
       counter=counter+1;
    end
end
R_Final=[index_R' x_R' y_R' Distance_R'];
 
%Delete No more necessary variables.
clear Hist_R temp x_R y_R Distance_R index_R
 
%--------------------------------------------------------------------------
%Finding the Distance('norm') for each pair of blocks in the 'Sub_G' color 
%and saving the result in 'Distance_G' variable.
%--------------------------------------------------------------------------
%Converting the Sub_R from 'Matrix' shape to 'Single Column' shape.
Hist_G = reshape(Hist_G,1,[])';
counter=1;
x_G=1;
y_G=1;
index_G=1;
for i=1:(tot_block)-1           %Counter for Block X-index.
    for j=(i+1):(tot_block)     %Counter for Block Y-index.
       temp= kolmogorov_smirnov_distance(double(Hist_G{i,1}), double(Hist_G{j,1}));
       x_G(counter)=i;
       y_G(counter)=j;
       Distance_G(counter)=temp;
       index_G(counter)=counter;
       counter=counter+1;
    end
end
G_Final=[index_G' x_G' y_G' Distance_G'];
 
%Delete No more necessary variables.
clear Hist_G temp x_G y_G Distance_G index_G
 
%--------------------------------------------------------------------------
%Finding the Distance('norm') for each pair of blocks in the 'Sub_B' color 
%and saving the result in 'Distance_B' variable.
%--------------------------------------------------------------------------
%Converting the Sub_R from 'Matrix' shape to 'Single Column' shape.
Hist_B = reshape(Hist_B,1,[])';
counter=1;
x_B=1;
y_B=1;
index_B=1;
for i=1:(tot_block)-1           %Counter for Block X-index.
    for j=(i+1):(tot_block)     %Counter for Block Y-index.
       temp= kolmogorov_smirnov_distance(double(Hist_B{i,1}), double(Hist_B{j,1}));
       x_B(counter)=i;
       y_B(counter)=j;
       Distance_B(counter)=temp;
       index_B(counter)=counter;
       counter=counter+1;
    end
end
B_Final=[index_B' x_B' y_B' Distance_B'];
 
%Delete No more necessary variables.
clear Hist_B temp x_B y_B Distance_B index_B
 
%--------------------------------------------------------------------------
% Sorting the distances in ascent way for Distance_R .
%--------------------------------------------------------------------------
R_Final=sortrows(R_Final,4);
 
%--------------------------------------------------------------------------
% Sorting the distances in ascent way for Distance_G .
%--------------------------------------------------------------------------
G_Final=sortrows(G_Final,4);
 
%--------------------------------------------------------------------------
% Sorting the distances in ascent way for Distance_B .
%--------------------------------------------------------------------------
B_Final=sortrows(B_Final,4);
 
%--------------------------------------------------------------------------
% Keeping the first 1/5 from the total R_Final with indices.
%--------------------------------------------------------------------------
R_Final=R_Final(1:retain_block_pair,:);
 
%--------------------------------------------------------------------------
% Keeping the first 1/5 from the total G_Final with indices.
%--------------------------------------------------------------------------
G_Final=G_Final(1:retain_block_pair,:);
 
%--------------------------------------------------------------------------
% Keeping the first 1/5 from the total B_Final with indices.
%--------------------------------------------------------------------------
B_Final=B_Final(1:retain_block_pair,:);
 
%-----------------------------------------------------------------------------------
% Checking if any block pair indices exists (Matched) in all the three 
% Fianl_R,Final_G & Final_B sorted lists.
%-----------------------------------------------------------------------------------
x=1;
temp=0;
for i=1:retain_block_pair
    temp=R_Final(i);
    for j=1:retain_block_pair
        if isequal(temp,G_Final(j))
            for k=1:retain_block_pair
                if isequal(temp,B_Final(k))
                    Matching_Indices_R(x) = i; %Saving existing indices.
                    Matching_Indices_G(x) = j;
                    Matching_Indices_B(x) = k;
                    x=x+1;
                end
            end
        end
    end
end
                   
 
if exist('Matching_Indices_R','var')
    disp('Number of Existing Indices :')
    length(Matching_Indices_R)
    %------------------------------------------------------------------
    % Cut the Lists to keep the matching indices only.
    %------------------------------------------------------------------
    for w = 1:length(Matching_Indices_R)
       Final_Result_R(w,1:4) = R_Final(Matching_Indices_R(w),1:4);
       Final_Result_G(w,1:4) = G_Final(Matching_Indices_G(w),1:4);
       Final_Result_B(w,1:4) = B_Final(Matching_Indices_B(w),1:4);
    end
    %Delete No more necessary variables.
    clear R_Final G_Final B_Final Matching_Indices_B Matching_Indices_G
    
    
   %------------------------------------------------------------------
   % Finding MSB of the First Final Selected Blocks 
   % (MSB_R1 , MSB_G1 & MSB_B1).
   %------------------------------------------------------------------
   Final_Image=imread('forged.jpg');
   nBins = [0:7];
   
   for w = 1:length(Matching_Indices_R) % Starting by the 1st Selected Block in R,G&B Color.
       counter=1;
       R1 = (rem(Final_Result_R(w,2)-1 , tot_block_row)) * OVERLAP_S + 1;
       R2 = (ceil(Final_Result_R(w,2) / tot_block_row) - 1) * OVERLAP_S + 1;
 
       for i=R1:(R1+(BLOCK_S-1))
           for j=R2:(R2+(BLOCK_S-1))
               MSB_RR1 = bitget(Final_Image(i,j,1),8);
               MSB_RR2 = bitget(Final_Image(i,j,1),7);
               MSB_RR3 = bitget(Final_Image(i,j,1),6);
               MSB_R1(w,counter)  = MSB_RR1*4 + MSB_RR2*2 + MSB_RR3*1;
               
               MSB_GG1 = bitget(Final_Image(i,j,2),8);
               MSB_GG2 = bitget(Final_Image(i,j,2),7);
               MSB_GG3 = bitget(Final_Image(i,j,2),6);
               MSB_G1(w,counter)  = MSB_GG1*4 + MSB_GG2*2 + MSB_GG3*1;
               
               MSB_BB1 = bitget(Final_Image(i,j,3),8);
               MSB_BB2 = bitget(Final_Image(i,j,3),7);
               MSB_BB3 = bitget(Final_Image(i,j,3),6);
               MSB_B1(w,counter)  = MSB_BB1*4 + MSB_BB2*2 + MSB_BB3*1;
                                             
               counter=counter+1;
           end
       end
   end   
   for w = 1:length(Matching_Indices_R) 
        Hist_R1{w}=histc(MSB_R1(w,:),nBins);
        Hist_G1{w}=histc(MSB_G1(w,:),nBins);
        Hist_B1{w}=histc(MSB_B1(w,:),nBins);
   end
   
   %------------------------------------------------------------------
   % Finding MSB of the Second Final Selected Blocks 
   % (MSB_R2 , MSB_G2 & MSB_B2).
   %------------------------------------------------------------------
   
   for w = 1:length(Matching_Indices_R) % Starting by the 1st Selected Block in R,G&B Color.
       counter=1;
       R1 = (rem(Final_Result_R(w,3)-1 , tot_block_row)) * OVERLAP_S + 1;
       R2 = (ceil(Final_Result_R(w,3) / tot_block_row) - 1) * OVERLAP_S + 1;
       for i=R1:(R1+(BLOCK_S-1))
           for j=R2:(R2+(BLOCK_S-1))
               MSB_RR1 = bitget(Final_Image(i,j,1),8);
               MSB_RR2 = bitget(Final_Image(i,j,1),7);
               MSB_RR3 = bitget(Final_Image(i,j,1),6);
               MSB_R2(w,counter)  = MSB_RR1*4 + MSB_RR2*2 + MSB_RR3*1;
               
               MSB_GG1 = bitget(Final_Image(i,j,2),8);
               MSB_GG2 = bitget(Final_Image(i,j,2),7);
               MSB_GG3 = bitget(Final_Image(i,j,2),6);
               MSB_G2(w,counter)  = MSB_GG1*4 + MSB_GG2*2 + MSB_GG3*1;
               
               MSB_BB1 = bitget(Final_Image(i,j,3),8);
               MSB_BB2 = bitget(Final_Image(i,j,3),7);
               MSB_BB3 = bitget(Final_Image(i,j,3),6);
               MSB_B2(w,counter)  = MSB_BB1*4 + MSB_BB2*2 + MSB_BB3*1;
                            
               counter=counter+1;
           end
       end
   end  
   for w = 1:length(Matching_Indices_R)
       Hist_R2{w}=histc(MSB_R2(w,:),nBins);
       Hist_G2{w}=histc(MSB_G2(w,:),nBins);
       Hist_B2{w}=histc(MSB_B2(w,:),nBins);
   end
   %------------------------------------------------------------------
   % Finding the Distances Between the Selected Blocks Using the
   % Results from Histrogram of MSB's.
   %------------------------------------------------------------------
   for w = 1:length(Matching_Indices_R)
       %Finding distances between the Blocks in MSB_R1 & MSB_R2.
       dR1_2(w) = kolmogorov_smirnov_distance(double(cell2mat(Hist_R1(w))), double(cell2mat(Hist_R2(w))));
       %Finding distances between the Blocks in MSB_G1 & MSB_G2.
       dG1_2(w) = kolmogorov_smirnov_distance(double(cell2mat(Hist_G1(w))), double(cell2mat(Hist_G2(w))));
       %Finding distances between the Blocks in MSB_B1 & MSB_B2.
       dB1_2(w) = kolmogorov_smirnov_distance(double(cell2mat(Hist_B1(w))), double(cell2mat(Hist_B2(w))));
   end
   
   %------------------------------------------------------------------
   % Checking the filterred blocks for deleting the unwanted blocks.
   %------------------------------------------------------------------ 
  counter2 = 0; 
  for w = 1:length(Matching_Indices_R)
       if ((dR1_2(w) < MSB_TH) && (dG1_2(w) < MSB_TH) && (dB1_2(w) < MSB_TH))
            counter2 = counter2 + 1;
            Final_Result_R1(counter2) = Final_Result_R(w,2);
            Final_Result_R2(counter2) = Final_Result_R(w,3);
       end
   end 
   
   
   %------------------------------------------------------------------
   % Filling common blocks with white color.
   %------------------------------------------------------------------
   Final_Image=imread('forged.jpg');
   
  for w = 1:counter2 % Filling 1st Block in R,G&B with White Color.
%      White_Blocks_R(1,w) = Final_Result_R(w,2);
       R1 = (rem(Final_Result_R1(w)-1 , tot_block_row)) * OVERLAP_S + 1;
       R2 = (ceil(Final_Result_R1(w) / tot_block_row) - 1) * OVERLAP_S + 1;
       R11(w) = R1;
       R12(w) = R2;
       
       R1 = (rem(Final_Result_R2(w)-1 , tot_block_row)) * OVERLAP_S + 1;
       R2 = (ceil(Final_Result_R2(w) / tot_block_row) - 1) * OVERLAP_S + 1;
       R21(w) = R1;
       R22(w) = R2;
       
%      Neighborhood blocks in the block pairs are removed 
       N_dist(w) = power(((R11(w)-R21(w))^2 + (R12(w)-R22(w))^2), 0.5); 
       
       if (N_dist(w) > MSB_TH)
           for i=R11(w):(R11(w)+(BLOCK_S-1))
                for j=R12(w):(R12(w)+(BLOCK_S-1))
                    Final_Image(i,j,1)=255;
                    Final_Image(i,j,2)=255;
                    Final_Image(i,j,3)=255;
                end
           end
           for i=R21(w):(R21(w)+(BLOCK_S-1))
                for j=R22(w):(R22(w)+(BLOCK_S-1))
                    Final_Image(i,j,1)=255;
                    Final_Image(i,j,2)=255;
                    Final_Image(i,j,3)=255;
                end
           end
       end
   end    
   
      figure
      imshow(Final_Image)
      xlabel('Forged Image')
      
     else
       disp(' No matching Blocks !')
end
