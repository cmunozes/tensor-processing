function [HA,long_axis_info]=...
    helix_angle_approx(mask,eigenvectors,centroid,region_option)

[n_lines, n_cols, n_slices] = size(mask);
prot.n_lines  = n_lines;
prot.n_cols   = n_cols;
prot.n_slices = n_slices; 

if isempty(centroid)
    centroid = [n_cols/2 n_lines/2];
end

first_eigenvectors=squeeze(eigenvectors(:,:,:,3));

if region_option == "LV_SA"
    
    %these vectors are not used for short-axis slices but they need to exist
    circ_vector=zeros(prot.n_lines,prot.n_cols,3);
    long_axis_info.vector=zeros(1,3);
    long_axis_info.points=zeros(3,2);
    long_axis_info.curve=zeros(2,2);
    long_axis_info.divider=zeros(prot.n_lines,prot.n_cols);
        
    %store theta values for all the pixels in the image considering the centre
    %to be the centroid
    theta_matrix=zeros(prot.n_lines,prot.n_cols);
    for iindex=1:prot.n_lines
        for jindex=1:prot.n_cols
            [THETA,RHO] = cart2pol(jindex-centroid(1),iindex-centroid(2));
            theta_matrix(iindex,jindex)=-THETA;
        end
    end
    
    % calculate the three vectors that will define the epicardial-tangent azimuthal and radial planes.
    % One of the planes is always [0 0 1], i.strain.E. through-plane the other is
    % in-plane and will depend on the theta value. Also important, we are
    % looking to the image from the feet.
    tp_component=zeros(prot.n_lines,prot.n_cols,3);
    circ_component=zeros(prot.n_lines,prot.n_cols,3);
    radial_component=zeros(prot.n_lines,prot.n_cols,3);
    for iindex=1:prot.n_lines
        for jindex=1:prot.n_cols
            tp_component(iindex,jindex,:)=[0 0 1];
            circ_component(iindex,jindex,:)=[sin(theta_matrix(iindex,jindex)),cos(theta_matrix(iindex,jindex)),0];
            radial_component(iindex,jindex,:)=[cos(theta_matrix(iindex,jindex)),-sin(theta_matrix(iindex,jindex)),0];
        end
    end
    
    %calculating the helical angle. The azimuthal component of the vector
    %points towards the patient's head/ hearts base. Assign signs to the
    %helical angle accordingly.
    first_eigenvector_proj=zeros(prot.n_lines,prot.n_cols,3); %projection in the epicardial plane
    HA=zeros(prot.n_lines,prot.n_cols);
    for iindex=1:prot.n_lines
        for jindex=1:prot.n_cols
            
            if mask(iindex,jindex) ~= 0
                
                current_plane_vector_1=squeeze(circ_component(iindex,jindex,:));
                current_plane_vector_2=squeeze(tp_component(iindex,jindex,:));
                A=squeeze(first_eigenvectors(iindex,jindex,:));
                B=[current_plane_vector_1 current_plane_vector_2];
                
                first_eigenvector_proj(iindex,jindex,:)=B*inv(B'*B)*B'*A;
                
                HA(iindex,jindex)=acosd(dot(squeeze(first_eigenvector_proj(iindex,jindex,:))...
                    ,squeeze(circ_component(iindex,jindex,:)))...
                    /(norm(squeeze(first_eigenvector_proj(iindex,jindex,:)))...
                    *norm(squeeze(circ_component(iindex,jindex,:)))));
                
                if first_eigenvector_proj(iindex,jindex,3)>0
                    
                    if HA(iindex,jindex)>90
                        HA(iindex,jindex)=(-1)*(180-HA(iindex,jindex));
                    end
                else
                    if HA(iindex,jindex)>90
                        HA(iindex,jindex)=(180-HA(iindex,jindex));
                    else
                        HA(iindex,jindex)=(-1)*HA(iindex,jindex);
                    end
                end
                
            end
            
        end
        
    end
    
end