function Und_subjects
 %subjects=['CK'; 'GM'; 'GS'; 'PP'; 'TC'];
 subjects='xx';
 th=5;
 n_sub=length(subjects);
 
 for subj=1:n_sub
    %names= subjects(subj);
    file_mat=['xx_TOTAL_connmatrix.txt'];
    C=load(file_mat);
    file_vec=['xx_ROIs_vector.txt'];
    voxvect=load(file_vec);
    voxmat=repmat(voxvect,1,90);
    C=C./voxmat;
       
    % Set diagonal to zero
    for n=1:90
        C(n,n)=0;
    end
    
     % Create undirectional matrix
     for i=1:90, 
         for j=1:90, 
             new(i,j)=mean([C(i,j) C(j,i)],2); 
         end;
     end;
     
     C=new;
     
    % Threshold to th=5
    I=find(C<th);
    C(I)=0;
    
 
    
    file_save=['Und_PRE_xx'];
    save(file_save,'C')
    
    
 end