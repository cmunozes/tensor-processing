function [eigenvalues, eigenvectors] = compute_eig_of_tensor3d(t11, t12, t13, t22, t23, t33, mask)

[n_lines, n_cols, n_slices] = size(t11);

% Initialise eigenvalue matrix
eigenvalues  = zeros([size(t11) 3]);
eigenvectors = zeros([size(t11) 3 3]);


for slice_index=1:n_slices
    for line_index=1:n_lines
        for col_index=1:n_cols

            if mask(line_index,col_index,slice_index)
                % Calculate the eigenvectors and eigenvalues of the structural
                % tensor
                st = [t11(line_index,col_index,slice_index) t12(line_index,col_index,slice_index) t13(line_index,col_index,slice_index); ...
                    t12(line_index,col_index,slice_index) t22(line_index,col_index,slice_index) t23(line_index,col_index,slice_index); ...
                    t13(line_index,col_index,slice_index) t23(line_index,col_index,slice_index) t33(line_index,col_index,slice_index)];
                [eigenvectors_temp,eigenvalues_temp] = eig(st);

                %each column is an eigenvector
                eigenvectors(line_index,col_index,slice_index,:,:)=eigenvectors_temp;
                %the values are in the diagonal of the matrix
                eigenvalues(line_index,col_index,slice_index,:)=diag(eigenvalues_temp);

                %sort dti_maps.eigenvalues by descending order
                [t,i]=sort(squeeze(eigenvalues(line_index,col_index,slice_index,:)),'descend');
                eigenvalues(line_index,col_index,slice_index,:)=t;

                %sort the dti_maps.eigenvectors by the same order
                eigenvectors_old(line_index,col_index,slice_index,:,:)=eigenvectors(line_index,col_index,slice_index,:,:);
                eigenvectors(line_index,col_index,slice_index,:,1)=eigenvectors_old(line_index,col_index,slice_index,:,i(1));
                eigenvectors(line_index,col_index,slice_index,:,2)=eigenvectors_old(line_index,col_index,slice_index,:,i(2));
                eigenvectors(line_index,col_index,slice_index,:,3)=eigenvectors_old(line_index,col_index,slice_index,:,i(3));

                %%%number of negative eigenvalues from the tensor fit
                %negative_eigenvalues(line_index,col_index,slice_index)=size(find(t<0),1);

            end

        end
    end
end

