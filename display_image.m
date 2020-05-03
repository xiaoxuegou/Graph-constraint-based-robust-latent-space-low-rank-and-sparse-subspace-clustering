function [] = display_image(X, database_name, showmode)
% display_image(data_train, 'PIE', 'imshow');
% [D N] = size(X);
% randArray = randperm(N);
% X = X(:, randArray);
is_digit = false;
%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp( database_name, 'ORL' ))
    faceW = 26;
    faceH = 32;
elseif ( strcmp( database_name, 'UMIST' ) )
    faceW = 32;
    faceH = 32;
elseif ( strcmp( database_name, 'JAFFE' ) )
    faceW = 32;
    faceH = 32;
elseif ( strcmp( database_name, 'CMU' ) )
    faceW = 32;
    faceH = 30;
else
    disp('inexact databasename');
    return;
end
numPerLine = 20; 
ShowLine = 2; 

Y = zeros(faceH*ShowLine,faceW*numPerLine); 
for i=0:ShowLine-1 
  	for j=0:numPerLine-1 
        if (is_digit == true)
            Y(i*faceH+1:(i+1)*faceH,j*faceW+1:(j+1)*faceW) = reshape(X(:, i*numPerLine+j+1),[faceH,faceW])';
        else
            Y(i*faceH+1:(i+1)*faceH,j*faceW+1:(j+1)*faceW) = reshape(X(:, i*numPerLine+j+1),[faceH,faceW]);
        end
  	end 
end 

if ( strcmp( showmode, 'imshow' ) )
    imshow(Y);
else
    imagesc(Y);colormap(gray);
end

return;